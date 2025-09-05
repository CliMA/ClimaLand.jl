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
    default_params_filepath =
        joinpath(pkgdir(ClimaLand), "toml", "default_parameters.toml")
    toml_dict = LP.create_toml_dict(FT, default_params_filepath)
    soil_moisture_stress_models = (
        TuzetMoistureStressModel{FT}(; sc = FT(0.01), pc = FT(-0.5e6)),
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
            toml_dict;
            soil_moisture_stress = sms_model,
        )
        @testset "Initialize canopy with type $name for float type $FT" begin
            Y, p, coords = initialize(canopy)

            @test all(parent(p.canopy.soil_moisture_stress.βm) .≈ FT(0.0))


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
        α = FT(0.04)
        n = FT(2)
        ν = FT(0.5)
        θ_r = FT(0.1)
        sms = PiecewiseMoistureStressModel{FT}(
            toml_dict;
            porosity_residual = true,
            soil_params = (; ν, θ_r),
        )
        @test θ_r == sms.θ_low
        @test ν == sms.θ_high
        βm_expected = FT(0.5)
        βm_computed = compute_piecewise_moisture_stress(
            sms.θ_high,
            sms.θ_low,
            sms.c,
            θ_r + (ν - θ_r) / 2,
        )
        @test βm_computed ≈ βm_expected
        sms = PiecewiseMoistureStressModel{FT}(
            toml_dict;
            porosity_residual = false,
            soil_params = (; ν, θ_r, α, n),
        )
        βm_expected =
            (θ_r + (ν - θ_r) / 2 - sms.θ_low) / (sms.θ_high - sms.θ_low)
        βm_computed = compute_piecewise_moisture_stress(
            sms.θ_high,
            sms.θ_low,
            sms.c,
            θ_r + (ν - θ_r) / 2,
        )
        @test βm_computed ≈ βm_expected
    end

end
