using Test
using Dates
using StaticArrays
using Insolation
import ClimaComms
ClimaComms.@import_required_backends
using ClimaCore
import ClimaParams as CP
using ClimaLand
using ClimaLand.Canopy
using ClimaLand.Canopy.PlantHydraulics
using ClimaLand.Domains: Column, HybridBox

import ClimaLand
import ClimaLand.Parameters as LP

for FT in (Float32, Float64)
    cmax = FT(0)
    cmin = FT(-2)
    nelems = 10
    col = Column(; zlim = (cmin, cmax), nelements = nelems)
    box = ClimaLand.Domains.HybridBox(;
        xlim = (cmin, cmax),
        ylim = (cmin, cmax),
        zlim = (cmin, cmax),
        nelements = (nelems, nelems, nelems),
    )

    domains = [col, box]
    t0 = 0.0
    earth_param_set = LP.LandParameters(FT)
    start_date = DateTime(2005)
    Δt = FT(180.0)
    "Radiation"
    SW_d = TimeVaryingInput((t) -> eltype(t)(20.0))
    LW_d = TimeVaryingInput((t) -> eltype(t)(20.0))
    zenith_angle =
        (t, s) -> default_zenith_angle(
            t,
            s;
            insol_params = earth_param_set.insol_params,
            latitude = FT(40.0),
            longitude = FT(-120.0),
        )

    rad = PrescribedRadiativeFluxes(
        FT,
        SW_d,
        LW_d,
        start_date;
        θs = zenith_angle,
        earth_param_set = earth_param_set,
    )
    "Atmos"
    precip = TimeVaryingInput((t) -> eltype(t)(0)) # no precipitation
    T_atmos = TimeVaryingInput((t) -> eltype(t)(290.0))
    u_atmos = TimeVaryingInput((t) -> eltype(t)(2.0))
    q_atmos = TimeVaryingInput((t) -> eltype(t)(0.011))
    h_atmos = FT(3)
    P_atmos = TimeVaryingInput((t) -> eltype(t)(101325))
    atmos = ClimaLand.PrescribedAtmosphere(
        precip,
        precip,
        T_atmos,
        u_atmos,
        q_atmos,
        P_atmos,
        start_date,
        h_atmos,
        earth_param_set,
    )

    ld = FT(0.5)
    α_PAR_leaf = FT(0.2)
    α_NIR_leaf = α_PAR_leaf
    is_c3 = FT(1.0)
    Vcmax25 = FT(5e-6)
    g1 = FT(120)
    G_Function = ConstantGFunction(ld)
    RTparams = BeerLambertParameters(FT; α_PAR_leaf, α_NIR_leaf, G_Function)
    photosynthesis_params = FarquharParameters(FT, is_c3; Vcmax25)
    stomatal_g_params = MedlynConductanceParameters(FT; g1)
    stomatal_model = MedlynConductanceModel{FT}(stomatal_g_params)
    photosynthesis_model = FarquharModel{FT}(photosynthesis_params)
    rt_model = BeerLambertModel{FT}(RTparams)
    energy_model = BigLeafEnergyModel{FT}(BigLeafEnergyParameters{FT}())
    thermo_params = LP.thermodynamic_parameters(earth_param_set)
    LAI = FT(8.0) # m2 [leaf] m-2 [ground]
    z_0m = FT(0.1) # m, Roughness length for momentum - value from tall forest ChatGPT
    z_0b = FT(0.1) # m, Roughness length for scalars - value from tall forest ChatGPT
    h_int = FT(30.0) # m, "where measurements would be taken at a typical flux tower of a 20m canopy"
    shared_params = SharedCanopyParameters{FT, typeof(earth_param_set)}(
        z_0m,
        z_0b,
        earth_param_set,
    )
    autotrophic_parameters = AutotrophicRespirationParameters(FT)
    autotrophic_respiration_model =
        AutotrophicRespirationModel{FT}(autotrophic_parameters)
    RAI = FT(1)
    SAI = FT(1)
    lai_fun = t -> LAI
    ai_parameterization = PlantHydraulics.PrescribedSiteAreaIndex{FT}(
        TimeVaryingInput(lai_fun),
        SAI,
        RAI,
    )
    K_sat_plant = FT(1.8e-8) # m/s
    ψ63 = FT(-4 / 0.0098) # / MPa to m, Holtzman's original parameter value
    Weibull_param = FT(4) # unitless, Holtzman's original c param value
    a = FT(0.05 * 0.0098) # Holtzman's original parameter for the bulk modulus of elasticity
    plant_ν = FT(0.7) # m3/m3
    plant_S_s = FT(1e-2 * 0.0098) # m3/m3/MPa to m3/m3/m
    conductivity_model =
        PlantHydraulics.Weibull{FT}(K_sat_plant, ψ63, Weibull_param)
    retention_model = PlantHydraulics.LinearRetentionCurve{FT}(a)
    rooting_depth = FT(0.5)
    param_set = PlantHydraulics.PlantHydraulicsParameters(;
        ai_parameterization = ai_parameterization,
        ν = plant_ν,
        S_s = plant_S_s,
        rooting_depth = rooting_depth,
        conductivity_model = conductivity_model,
        retention_model = retention_model,
    )
    Δz = FT(1.0) # height of compartments
    n_stem = Int64(2) # number of stem elements
    n_leaf = Int64(1) # number of leaf elements
    compartment_centers =
        FT.(
            Vector(
                range(
                    start = Δz / 2,
                    step = Δz,
                    stop = Δz * (n_stem + n_leaf) - (Δz / 2),
                ),
            ),
        )
    compartment_faces =
        FT.(
            Vector(
                range(start = 0.0, step = Δz, stop = Δz * (n_stem + n_leaf)),
            )
        )

    ψ_soil0 = FT(0.0)
    T_soil0 = FT(290)
    root_depths = SVector{10, FT}(-(10:-1:1.0) ./ 10.0 * 2.0 .+ 0.2 / 2.0)

    soil_driver = PrescribedGroundConditions(
        root_depths,
        (t) -> ψ_soil0,
        (t) -> T_soil0,
        FT(0.2),
        FT(0.4),
        FT(0.98),
    )
    plant_hydraulics = PlantHydraulics.PlantHydraulicsModel{FT}(;
        parameters = param_set,
        n_stem = n_stem,
        n_leaf = n_leaf,
        compartment_surfaces = compartment_faces,
        compartment_midpoints = compartment_centers,
    )
    @testset "Canopy model total energy and water, FT = $FT" begin
        for domain in domains
            canopy = ClimaLand.Canopy.CanopyModel{FT}(;
                parameters = shared_params,
                domain = ClimaLand.Domains.obtain_surface_domain(domain),
                radiative_transfer = rt_model,
                photosynthesis = photosynthesis_model,
                conductance = stomatal_model,
                autotrophic_respiration = autotrophic_respiration_model,
                energy = energy_model,
                hydraulics = plant_hydraulics,
                boundary_conditions = Canopy.AtmosDrivenCanopyBC(
                    atmos,
                    rad,
                    soil_driver,
                ),
            )

            Y, p, cds = initialize(canopy)
            ϑ0 = plant_ν / 2
            Temp0 = FT(290.5)

            Y.canopy.hydraulics.ϑ_l.:1 .= ϑ0
            Y.canopy.hydraulics.ϑ_l.:2 .= ϑ0
            Y.canopy.hydraulics.ϑ_l.:3 .= ϑ0
            Y.canopy.energy.T .= Temp0

            set_initial_cache! = make_set_initial_cache(canopy)
            set_initial_cache!(p, Y, t0)
            total_energy = ClimaCore.Fields.zeros(domain.space.surface)
            ClimaLand.total_energy_per_area!(total_energy, canopy, Y, p, t0)
            @test all(
                parent(total_energy) .≈
                (LAI + SAI) * canopy.energy.parameters.ac_canopy * Temp0,
            )

            total_water = ClimaCore.Fields.zeros(domain.space.surface)
            ClimaLand.total_liq_water_vol_per_area!(
                total_water,
                canopy,
                Y,
                p,
                t0,
            )
            @test all(
                parent(total_water) .≈
                (LAI * n_leaf * Δz * ϑ0 + SAI * n_stem * Δz * ϑ0),
            )
        end
    end
end
