
using Test
using ClimaLand
import ClimaComms
ClimaComms.@import_required_backends
using ClimaLand.Canopy
using ClimaLand.Soil
using ClimaLand.Canopy.PlantHydraulics
using DelimitedFiles
using ClimaLand.Domains: Point
import ClimaLand.Parameters as LP
import ClimaParams
import Insolation
using Dates

function setup_standalone_canopy(FT::Type{<:AbstractFloat}, moisture_stress_model) 
    domain = Point(; z_sfc = FT(0.0))

    # Autotrophic respiration 
    ar_params = AutotrophicRespirationParameters(FT)
    ar_model = AutotrophicRespirationModel{FT}(ar_params)

    # Radiative transfer 
    rt_params = BeerLambertParameters(FT)
    rt_model = BeerLambertModel{FT}(rt_params)

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
    h_canopy = FT(0) # m, height of the canopy
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

    earth_param_set = LP.LandParameters(FT)
    z_0m = FT(2.0) # m, Roughness length for momentum
    z_0b = FT(0.1) # m, Roughness length for scalars
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
        autotrophic_respiration = ar_model,
        radiative_transfer = rt_model,
        photosynthesis = photosynthesis_model,
        conductance = stomatal_model,
        soil_moisture_stress = moisture_stress_model, 
        hydraulics = plant_hydraulics,
        boundary_conditions = Canopy.AtmosDrivenCanopyBC(
            atmos,
            radiation,
            soil_driver,
        ),
    )

    return model
end 


for FT in (Float32, Float64)
    soil_moisture_stress_models = (
        TuzetMoistureStressModel{FT}(TuzetMoistureStressParameters{FT}(
            sc = FT(0.01), # Pa^-1
            pc = FT(-0.5e6), # Pa
        )),
        NoMoistureStressModel{FT}()
    )

    soil_moisture_stress_names = (
        "Tuzet",
        "No Stress",
    )

    for (sms_model, name) in zip(soil_moisture_stress_models, soil_moisture_stress_names)
        canopy = setup_standalone_canopy(FT, sms_model)
        @testset "Initialize canopy with type $name for float type $FT" begin
            Y, p, coords = initialize(canopy)
            dY = similar(Y)

            @test all(parent(p.canopy.soil_moisture_stress.βm) .≈ FT(0.0))

            if name == "Piecewise"
                @test all(parent(p.canopy.soil_moisture_stress.ϑ_root) .≈ FT(0.0))
            end

            
            # set initial conditions 
            n_stem = 0
            n_leaf = 1
            for i in 1:(n_stem + n_leaf)
                Y.canopy.hydraulics.ϑ_l.:($i) .= FT(0.1)
                p.canopy.hydraulics.ψ.:($i) .= NaN
                p.canopy.hydraulics.fa.:($i) .= NaN
                dY.canopy.hydraulics.ϑ_l.:($i) .= NaN
            end

            set_initial_cache! = make_set_initial_cache(canopy)
            set_initial_cache!(p, Y, FT(0.0))
            
            @test all(parent(p.canopy.soil_moisture_stress.βm) .≈ FT(1.0))
        end
    end

    @testset "Test correctness for piecewise soil moisture stress for float type $FT" begin
        # construct piecewise moisture stress parameters from hydrology
        hydrology_cm = vanGenuchten{FT}(; α = FT(0.04), n = FT(2))
        ν = FT(0.5)
        θ_r = FT(0.1)
        c = FT(0.5)
        β0 = FT(0.87)
        sms_params = PiecewiseMoistureStressParametersFromHydrology(
            hydrology_cm,
            ν,
            θ_r,
            c = c,
            β0 = β0
        )
        θ_w = sms_params.θ_w
        θ_c = sms_params.θ_c
        println("Wilting point = $θ_w, field capacity = $θ_c, c = $c, β0 = $β0")
        
        # check that \theta_w and \theta_c are physical
        @test θ_w < θ_c
        @test θ_w > FT(0.0)
        @test θ_c > FT(0.0)
        @test θ_w < FT(1.0)
        @test θ_c < FT(1.0)

        # check correctness
        θ_root = (θ_w + θ_c) / 2 # midpoint
        βm_expected = β0 * ((θ_root - θ_w) / (θ_c - θ_w))^c
        βm_computed = compute_piecewise_moisture_stress(
            sms_params,
            θ_root
        )
        @test βm_computed ≈ βm_expected

        θ_root = (θ_w) / 2 # drier than wilting point
        βm_expected = FT(0)
        βm_computed = compute_piecewise_moisture_stress(
            sms_params,
            θ_root
        )
        @test βm_computed ≈ βm_expected

        θ_root = (1 + θ_c) / 2 # wetter than field capacity
        βm_expected = β0
        βm_computed = compute_piecewise_moisture_stress(
            sms_params,
            θ_root
        )
        @test βm_computed ≈ βm_expected
    end 
    
    @testset "Test parameter domain validation for piecewise soil moisture stress for float type $FT" begin
        # Test direct constructor with invalid parameters
        θ_c = FT(0.3)
        θ_w = FT(0.1)
        
        # Test negative c parameter
        @test_throws ArgumentError PiecewiseMoistureStressParameters(
            θ_c = θ_c,
            θ_w = θ_w,
            c = FT(-1.0),
            β0 = FT(1.0)
        )
        
        # Test zero c parameter
        @test_throws ArgumentError PiecewiseMoistureStressParameters(
            θ_c = θ_c,
            θ_w = θ_w,
            c = FT(0.0),
            β0 = FT(1.0)
        )
        
        # Test negative β0 parameter
        @test_throws ArgumentError PiecewiseMoistureStressParameters(
            θ_c = θ_c,
            θ_w = θ_w,
            c = FT(1.0),
            β0 = FT(-0.5)
        )
        
        # Test zero β0 parameter
        @test_throws ArgumentError PiecewiseMoistureStressParameters(
            θ_c = θ_c,
            θ_w = θ_w,
            c = FT(1.0),
            β0 = FT(0.0)
        )
        
        # Test validation in helper function
        hydrology_cm = vanGenuchten{FT}(; α = FT(0.04), n = FT(2))
        ν = FT(0.5)
        θ_r = FT(0.1)
        
        # Test negative c parameter in helper function
        @test_throws ArgumentError PiecewiseMoistureStressParametersFromHydrology(
            hydrology_cm,
            ν,
            θ_r;
            c = FT(-1.0)
        )
        
        # Test negative β0 parameter in helper function
        @test_throws ArgumentError PiecewiseMoistureStressParametersFromHydrology(
            hydrology_cm,
            ν,
            θ_r;
            β0 = FT(-0.5)
        )
        
        # Test that valid parameters work correctly
        @test isa(PiecewiseMoistureStressParameters(
            θ_c = θ_c,
            θ_w = θ_w,
            c = FT(1.0),
            β0 = FT(1.0)
        ), PiecewiseMoistureStressParameters)
        
        @test isa(PiecewiseMoistureStressParametersFromHydrology(
            hydrology_cm,
            ν,
            θ_r;
            c = FT(1.0),
            β0 = FT(1.0)
        ), PiecewiseMoistureStressParameters)
    end
end

