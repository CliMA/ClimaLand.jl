
using Test
using ClimaLand
import ClimaComms
ClimaComms.@import_required_backends
using ClimaLand.Canopy
using DelimitedFiles
using ClimaLand.Domains: Point
import ClimaLand.Parameters as LP
import ClimaParams

function setup_standalone_canopy(FT::Type{<:AbstractFloat}) 
    domain = Point(; z_sfc = FT(0.0))

    # Autotrophic respiration 
    AR_params = AutotrophicRespirationParameters(FT)
    AR_model = AutotrophicRespirationModel{FT}(AR_params)

    # Radiative transfer 
    RTparams = BeerLambertParameters(FT)
    rt_model = BeerLambertModel{FT}(RTparams)

    # Photosynthesis 
    is_c3 = FT(1) 
    photosynthesis_params = FarquharParameters(FT, is_c3)
    photosynthesis_model = FarquharModel{FT}(photosynthesis_params)

    # Stomatal conductance
    stomatal_g_params = MedlynConductanceParameters(FT)
    stomatal_model = MedlynConductanceModel{FT}(stomatal_g_params)
    
    # Plant hydraulics 
    n_stem = Int64(0) # number of stem elements
    n_leaf = Int64(1) # number of leaf elements
    SAI = FT(0) # m2/m2
    RAI = FT(0) # m2/m2
    lai_fun = t -> 0
    ai_parameterization = PlantHydraulics.PrescribedSiteAreaIndex{FT}(
        TimeVaryingInput(lai_fun),
        SAI,
        RAI,
    )
    K_sat_plant = 0 # m/s.
    ψ63 = FT(-4 / 0.0098) # / MPa to m
    Weibull_param = FT(4) # unitless
    a = FT(0.05 * 0.0098) # 1/m
    plant_ν = FT(0.1) # m3/m3
    plant_S_s = FT(1e-2)
    conductivity_model =
        PlantHydraulics.Weibull{FT}(K_sat_plant, ψ63, Weibull_param)
    retention_model = PlantHydraulics.LinearRetentionCurve{FT}(a)
    root_depths = [FT(0.0)]# 1st element is the deepest root depth
    compartment_midpoints = [h_canopy]
    compartment_surfaces = [FT(0.0), h_canopy]
    # set rooting_depth param to largest possible value to test no roots
    param_set = PlantHydraulics.PlantHydraulicsParameters(;
        ai_parameterization = ai_parameterization,
        ν = plant_ν,
        S_s = plant_S_s,
        rooting_depth = maxintfloat(FT),
        conductivity_model = conductivity_model,
        retention_model = retention_model,
    )

    transpiration = DiagnosticTranspiration{FT}()
    soil_driver = PrescribedGroundConditions(FT)
    plant_hydraulics = PlantHydraulics.PlantHydraulicsModel{FT}(;
        parameters = param_set,
        transpiration = transpiration,
        n_stem = n_stem,
        n_leaf = n_leaf,
        compartment_surfaces = compartment_surfaces,
        compartment_midpoints = compartment_midpoints,
    )

    autotrophic_parameters = AutotrophicRespirationParameters(FT)
    autotrophic_respiration_model =
        AutotrophicRespirationModel(autotrophic_parameters)


    earth_param_set = LP.LandParameters(FT)
    thermo_params = LP.thermodynamic_parameters(earth_param_set)
    LAI = FT(0.0) # m2 [leaf] m-2 [ground]
    z_0m = FT(2.0) # m, Roughness length for momentum
    z_0b = FT(0.1) # m, Roughness length for scalars
    h_canopy = FT(0.0) # m, canopy height
    h_sfc = FT(0.0) # m
    h_int = FT(30.0) # m, "where measurements would be taken at a typical flux tower of a 20m canopy"
    shared_params = SharedCanopyParameters{FT, typeof(earth_param_set)}(
        z_0m,
        z_0b,
        earth_param_set,
    )
    lat = FT(0.0) # degree
    long = FT(-180) # degree
    start_date = DateTime(2005)

    # Dummy radiative forcing
    zenith_angle =
        (t, s) -> default_zenith_angle(
            t,
            s;
            insol_params = earth_param_set.insol_params,
            latitude = lat,
            longitude = long,
        )
    shortwave_radiation = (t; latitude=lat, longitude=long, insol_params=earth_param_set.insol_params) -> 1000
    longwave_radiation = t -> 200
    radiation = PrescribedRadiativeFluxes(
        FT,
        TimeVaryingInput(shortwave_radiation),
        TimeVaryingInput(longwave_radiation),
        start_date;
        θs = zenith_angle,
        earth_param_set = earth_param_set,
    )

    # Dummy atmospheric forcing 
    u_atmos = t -> 10 #m.s-1
    liquid_precip = (t) -> 0 # m
    snow_precip = (t) -> 0 # m
    T_atmos = t -> 290 # Kelvin
    q_atmos = t -> 0.001 # kg/kg
    P_atmos = t -> 1e5 # Pa
    h_atmos = h_int # m
    c_atmos = (t) -> 4.11e-4 # mol/mol
    atmos = PrescribedAtmosphere(
        TimeVaryingInput(liquid_precip),
        TimeVaryingInput(snow_precip),
        TimeVaryingInput(T_atmos),
        TimeVaryingInput(u_atmos),
        TimeVaryingInput(q_atmos),
        TimeVaryingInput(P_atmos),
        start_date,
        h_atmos,
        earth_param_set;
        c_co2 = TimeVaryingInput(c_atmos),
    )


    model = ClimaLand.Canopy.CanopyModel{FT}(;
        parameters = shared_params,
        domain = domain,
        autotrophic_respiration = autotrophic_respiration_model,
        radiative_transfer = rt_model,
        photosynthesis = photosynthesis_model,
        conductance = stomatal_model,
        hydraulics = plant_hydraulics,
        boundary_conditions = Canopy.AtmosDrivenCanopyBC(
            atmos,
            radiation,
            soil_driver,
        ),
    )

    return model
end 


@testset "Initialize canopy with soil moisture stress" begin
    for FT in (Float32, Float64)
        model = setup_canopy(FT)
        @testset "Initialize canopy with soil moisture stress for float type $FT" begin
            Y, p, coords = initialize(model)
            dY = similar(Y)
            for i in 1:(n_stem + n_leaf)
                Y.canopy.hydraulics.ϑ_l.:($i) .= FT(0.1)
                p.canopy.hydraulics.ψ.:($i) .= NaN
                p.canopy.hydraulics.fa.:($i) .= NaN
                dY.canopy.hydraulics.ϑ_l.:($i) .= NaN
            end
            set_initial_cache! = make_set_initial_cache(model)
            set_initial_cache!(p, Y, FT(0.0))
            @test all(parent(p.canopy.hydraulics.fa) .≈ FT(0.0))
            @test all(parent(p.canopy.hydraulics.fa_roots) .≈ FT(0.0))
            @test all(parent(p.canopy.turbulent_fluxes.transpiration) .≈ FT(0.0))
            @test all(parent(p.canopy.radiative_transfer.par.abs) .≈ FT(0.0))
            exp_tend! = make_exp_tendency(model)
            exp_tend!(dY, Y, p, FT(0))
            @test all(parent(dY.canopy.hydraulics.ϑ_l.:1) .≈ FT(0.0))

        end
    end
end