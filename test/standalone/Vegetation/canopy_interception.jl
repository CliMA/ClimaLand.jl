using Test
import ClimaParams
import ClimaComms
ClimaComms.@import_required_backends
using ClimaCore
using Thermodynamics
using Dates
using ClimaLand
using ClimaLand: PrescribedAtmosphere, PrescribedRadiativeFluxes
using ClimaUtilities.TimeVaryingInputs: TimeVaryingInput
using ClimaLand.Canopy
using ClimaLand.Canopy:
    PrescribedBiomassModel,
    PlantHydraulicsModel,
    BeerLambertModel,
    ConstantGFunction,
    FarquharModel,
    MedlynConductanceModel,
    AutotrophicRespirationModel,
    BigLeafEnergyModel,
    NoMoistureStressModel,
    CanopyInterceptionModel,
    CanopyInterceptionParameters,
    NoCanopyInterception
using ClimaLand.Domains: Point
import ClimaLand.Parameters as LP

@testset "Canopy Interception - unit formulas" begin
    FT = Float64
    # Test free throughfall fraction (Beer-Lambert analogy)
    k_int = FT(0.5)
    LAI = FT(4.0)
    p_free = exp(-k_int * LAI)
    @test p_free ≈ exp(-2.0)
    @test 0 < p_free < 1

    # Test maximum storage capacity
    S_leaf = FT(0.0002)
    S_stem = FT(0.00025)
    SAI = FT(1.0)
    S_max = S_leaf * LAI + S_stem * SAI
    @test S_max ≈ 0.0002 * 4.0 + 0.00025 * 1.0

    # Test wet fraction (Deardorff 1978)
    W_int = FT(0.5) * S_max
    f_wet = min(FT(1), (W_int / S_max)^FT(2 / 3))
    @test f_wet ≈ 0.5^(2 / 3)
    @test 0 < f_wet < 1

    # f_wet = 1 when W_int >= S_max
    f_wet_full = min(FT(1), (S_max / S_max)^FT(2 / 3))
    @test f_wet_full ≈ 1.0

    # f_wet = 0 when W_int = 0
    f_wet_empty = min(FT(1), (FT(0) / S_max)^FT(2 / 3))
    @test f_wet_empty ≈ 0.0

    # Test drainage (Rutter exponential)
    D_s = FT(2.4e-5)
    b = FT(3700.0)
    # Below capacity: D ≈ D_s
    W_below = FT(0.9) * S_max
    D_below = D_s * exp(b * max(W_below - S_max, FT(0)))
    @test D_below ≈ D_s  # since W_below < S_max, exponent is 0

    # Above capacity: D > D_s
    W_above = S_max + FT(0.001)
    D_above = D_s * exp(b * max(W_above - S_max, FT(0)))
    @test D_above > D_s
end

@testset "Canopy Interception - NoCanopyInterception default" begin
    for FT in (Float32, Float64)
        toml_dict = LP.create_toml_dict(FT)
        domain = Point(; z_sfc = FT(0), longlat = (FT(-180), FT(1)))

        # Build a canopy model with default (NoCanopyInterception)
        radiation_parameters = (;
            α_PAR_leaf = FT(0.1),
            α_NIR_leaf = FT(0.4),
            Ω = FT(1),
            G_Function = ConstantGFunction(FT(0.5)),
        )
        rt_model = BeerLambertModel{FT}(domain, toml_dict; radiation_parameters)
        photosynthesis_model = FarquharModel{FT}(
            domain,
            toml_dict;
            photosynthesis_parameters = (;
                fractional_c3 = FT(1),
                Vcmax25 = FT(9e-5),
            ),
        )
        energy_model = BigLeafEnergyModel{FT}(toml_dict)
        AR_model = AutotrophicRespirationModel{FT}(toml_dict)
        stomatal_model =
            MedlynConductanceModel{FT}(domain, toml_dict; g1 = FT(790))
        LAI = TimeVaryingInput(t -> FT(4))
        hydraulics = PlantHydraulicsModel{FT}(domain, toml_dict)
        biomass = PrescribedBiomassModel{FT}(domain, LAI, toml_dict)

        earth_param_set = LP.LandParameters(toml_dict)
        start_date = DateTime(2005)

        atmos = PrescribedAtmosphere(
            TimeVaryingInput(t -> FT(1e-5)),  # liquid_precip
            TimeVaryingInput(t -> FT(0)),
            TimeVaryingInput(t -> FT(290)),
            TimeVaryingInput(t -> FT(10)),
            TimeVaryingInput(t -> FT(0.001)),
            TimeVaryingInput(t -> FT(1e5)),
            start_date,
            FT(30),
            toml_dict;
            c_co2 = TimeVaryingInput(t -> FT(4.11e-4)),
        )
        radiation = PrescribedRadiativeFluxes(
            FT,
            TimeVaryingInput(t -> FT(500)),
            TimeVaryingInput(t -> FT(200)),
            start_date;
            toml_dict,
        )

        ground = PrescribedGroundConditions{FT}()
        forcing = (; atmos, radiation, ground)

        canopy = CanopyModel{FT}(
            domain,
            forcing,
            LAI,
            toml_dict;
            autotrophic_respiration = AR_model,
            radiative_transfer = rt_model,
            photosynthesis = photosynthesis_model,
            conductance = stomatal_model,
            soil_moisture_stress = NoMoistureStressModel{FT}(),
            hydraulics,
            energy = energy_model,
            biomass,
        )

        # Default interception is NoCanopyInterception
        @test canopy.interception isa NoCanopyInterception{FT}

        # Initialize and check state
        Y, p, coords = ClimaLand.initialize(canopy)
        @test :interception in propertynames(p.canopy)
        @test :W_int in propertynames(p.canopy.interception)
        @test :f_wet in propertynames(p.canopy.interception)
        @test :throughfall in propertynames(p.canopy.interception)

        # NoCanopyInterception should have no prognostic vars
        @test !hasproperty(Y.canopy, :interception)
    end
end

@testset "Canopy Interception - CanopyInterceptionModel component" begin
    for FT in (Float32, Float64)
        toml_dict = LP.create_toml_dict(FT)
        domain = Point(; z_sfc = FT(0), longlat = (FT(-180), FT(1)))

        # Build parameters
        params = CanopyInterceptionParameters{FT}()
        @test params.k_int ≈ FT(0.5)
        @test params.S_leaf ≈ FT(0.0002)

        # Build from TOML
        params_toml = CanopyInterceptionParameters(toml_dict)
        @test params_toml.k_int ≈ FT(0.5)

        # Build model
        int_model = CanopyInterceptionModel{FT}(params)
        @test int_model isa CanopyInterceptionModel{FT}
        @test ClimaLand.Canopy.name(int_model) == :interception
        @test ClimaLand.prognostic_vars(int_model) == (:W_int,)
        @test ClimaLand.auxiliary_vars(int_model) == (:f_wet, :throughfall)

        # Build a canopy model with active interception
        radiation_parameters = (;
            α_PAR_leaf = FT(0.1),
            α_NIR_leaf = FT(0.4),
            Ω = FT(1),
            G_Function = ConstantGFunction(FT(0.5)),
        )
        rt_model = BeerLambertModel{FT}(domain, toml_dict; radiation_parameters)
        photosynthesis_model = FarquharModel{FT}(
            domain,
            toml_dict;
            photosynthesis_parameters = (;
                fractional_c3 = FT(1),
                Vcmax25 = FT(9e-5),
            ),
        )
        energy_model = BigLeafEnergyModel{FT}(toml_dict)
        AR_model = AutotrophicRespirationModel{FT}(toml_dict)
        stomatal_model =
            MedlynConductanceModel{FT}(domain, toml_dict; g1 = FT(790))
        LAI = TimeVaryingInput(t -> FT(4))
        hydraulics = PlantHydraulicsModel{FT}(domain, toml_dict)
        biomass = PrescribedBiomassModel{FT}(domain, LAI, toml_dict)

        earth_param_set = LP.LandParameters(toml_dict)
        start_date = DateTime(2005)

        atmos = PrescribedAtmosphere(
            TimeVaryingInput(t -> FT(1e-5)),
            TimeVaryingInput(t -> FT(0)),
            TimeVaryingInput(t -> FT(290)),
            TimeVaryingInput(t -> FT(10)),
            TimeVaryingInput(t -> FT(0.001)),
            TimeVaryingInput(t -> FT(1e5)),
            start_date,
            FT(30),
            toml_dict;
            c_co2 = TimeVaryingInput(t -> FT(4.11e-4)),
        )
        radiation = PrescribedRadiativeFluxes(
            FT,
            TimeVaryingInput(t -> FT(500)),
            TimeVaryingInput(t -> FT(200)),
            start_date;
            toml_dict,
        )

        ground = PrescribedGroundConditions{FT}()
        forcing = (; atmos, radiation, ground)

        canopy = CanopyModel{FT}(
            domain,
            forcing,
            LAI,
            toml_dict;
            autotrophic_respiration = AR_model,
            radiative_transfer = rt_model,
            photosynthesis = photosynthesis_model,
            conductance = stomatal_model,
            soil_moisture_stress = NoMoistureStressModel{FT}(),
            hydraulics,
            energy = energy_model,
            biomass,
            interception = int_model,
        )

        @test canopy.interception isa CanopyInterceptionModel{FT}

        Y, p, coords = ClimaLand.initialize(canopy)

        # Check prognostic variable exists
        @test hasproperty(Y.canopy, :interception)
        @test :W_int in propertynames(Y.canopy.interception)

        # Check auxiliary variables
        @test :f_wet in propertynames(p.canopy.interception)
        @test :throughfall in propertynames(p.canopy.interception)

        # Check structure of Y is valid
        @test !isnothing(zero(Y))
    end
end
