
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

for FT in (Float32, Float64)
    soil_moisture_stress_models = (
        TuzetMoistureStressModel{FT}(
            TuzetMoistureStressParameters{FT}(
                sc = FT(0.01), # Pa^-1
                pc = FT(-0.5e6), # Pa
            ),
        ),
        NoMoistureStressModel{FT}(),
    )

    soil_moisture_stress_names = ("Tuzet", "No Stress")

    lat = FT(38.7441)
    long = FT(-92.2000)
    domain = Point(; z_sfc = FT(0.0), longlat = (long, lat))
    earth_param_set = LP.LandParameters(FT)

    # dummy forcing
    atmos, radiation = prescribed_analytic_forcing(FT; h_atmos = FT(1))
    ground = PrescribedGroundConditions{FT}()
    LAI = TimeVaryingInput(t -> FT(1))
    forcing = (; atmos, radiation, ground)

    for (sms_model, name) in
        zip(soil_moisture_stress_models, soil_moisture_stress_names)
        canopy = Canopy.CanopyModel{FT}(
            domain,
            forcing,
            LAI,
            earth_param_set;
            soil_moisture_stress = sms_model,
        )
        @testset "Initialize canopy with type $name for float type $FT" begin
            Y, p, coords = initialize(canopy)

            @test all(parent(p.canopy.soil_moisture_stress.βm) .≈ FT(0.0))

            if name == "Piecewise"
                @test all(
                    parent(p.canopy.soil_moisture_stress.ϑ_root) .≈ FT(0.0),
                )
            end

            # # set initial conditions 
            n_stem = 0
            n_leaf = 1
            for i in 1:(n_stem + n_leaf)
                Y.canopy.hydraulics.ϑ_l.:($i) .= FT(0.1)
                p.canopy.hydraulics.ψ.:($i) .= NaN
                p.canopy.hydraulics.fa.:($i) .= NaN
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
        soil_zmin = FT(-10)
        sms_params = PiecewiseMoistureStressParametersFromHydrology(
            FT,
            hydrology_cm,
            ν,
            θ_r,
            soil_zmin,
            c = c,
            β0 = β0,
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
        βm_computed = compute_piecewise_moisture_stress(θ_c, θ_w, c, β0, θ_root)
        @test βm_computed ≈ βm_expected

        θ_root = (θ_w) / 2 # drier than wilting point
        βm_expected = FT(0)
        βm_computed = compute_piecewise_moisture_stress(θ_c, θ_w, c, β0, θ_root)
        @test βm_computed ≈ βm_expected

        θ_root = (θ_w) / 2 # drier than wilting point
        βm_expected = FT(0)
        βm_computed = compute_piecewise_moisture_stress(θ_c, θ_w, c, β0, θ_root)
        @test βm_computed ≈ βm_expected

        θ_root = (1 + θ_c) / 2 # wetter than field capacity
        βm_expected = β0
        βm_computed = compute_piecewise_moisture_stress(θ_c, θ_w, c, β0, θ_root)
        @test βm_computed ≈ βm_expected
    end

end
