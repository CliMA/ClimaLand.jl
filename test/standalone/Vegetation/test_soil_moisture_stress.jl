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
    toml_dict = LP.create_toml_dict(FT)
    soil_moisture_stress_models = (
        TuzetMoistureStressModel{FT}(; sc = FT(0.01), pc = FT(-0.5e6)),
        NoMoistureStressModel{FT}(),
    )
    soil_moisture_stress_names = ("Tuzet", "No Stress")

    lat = FT(38.7441)
    long = FT(-92.2000)
    domain = Point(; z_sfc = FT(0.0), longlat = (long, lat))
    earth_param_set = LP.LandParameters(toml_dict)

    # dummy forcing
    atmos, radiation =
        prescribed_analytic_forcing(FT; toml_dict, h_atmos = FT(1))
    ground = PrescribedGroundConditions{FT}()
    LAI = TimeVaryingInput(t -> FT(1))
    forcing = (; atmos, radiation, ground)

    # Set up optimal LAI model
    lai_model =
        Canopy.OptimalLAIModel{FT}(Canopy.OptimalLAIParameters{FT}(toml_dict))

    for (sms_model, name) in
        zip(soil_moisture_stress_models, soil_moisture_stress_names)
        canopy = Canopy.CanopyModel{FT}(
            domain,
            forcing,
            LAI,
            toml_dict;
            soil_moisture_stress = sms_model,
            lai_model,
        )
        @testset "Initialize canopy with type $name for float type $FT" begin
            Y, p, coords = initialize(canopy)

            @test all(parent(p.canopy.soil_moisture_stress.βm) .≈ FT(0.0))
            # # set initial conditions
            Y.canopy.hydraulics.ϑ_l.:1 .= canopy.hydraulics.parameters.ν
            p.canopy.hydraulics.ψ.:1 .= NaN
            set_initial_cache! = make_set_initial_cache(canopy)
            set_initial_cache!(p, Y, FT(0.0))

            @test all(parent(p.canopy.soil_moisture_stress.βm) .≈ FT(1.0))

            # Repeat with different IC
            p.canopy.soil_moisture_stress.βm .= 0
            Y.canopy.hydraulics.ϑ_l.:1 .= canopy.hydraulics.parameters.ν / 2
            p.canopy.hydraulics.ψ.:1 .= NaN
            set_initial_cache! = make_set_initial_cache(canopy)
            set_initial_cache!(p, Y, FT(0.0))
            grav = LP.grav(earth_param_set)
            ρ_water = LP.ρ_cloud_liq(earth_param_set)
            if canopy.soil_moisture_stress isa TuzetMoistureStressModel
                p_leaf = @. p.canopy.hydraulics.ψ.:1 * ρ_water * grav
                @test all(
                    parent(p.canopy.soil_moisture_stress.βm) .≈ parent(
                        compute_tuzet_moisture_stress.(
                            p_leaf,
                            canopy.soil_moisture_stress.pc,
                            canopy.soil_moisture_stress.sc,
                        ),
                    ),
                )
            else
                @test all(parent(p.canopy.soil_moisture_stress.βm) .≈ FT(1.0))
            end
        end
    end

    @testset "Test correctness for piecewise soil moisture stress for float type $FT" begin
        # construct piecewise moisture stress parameters from hydrology

        ν = FT(0.5)
        θ_r = FT(0.1)
        sms = PiecewiseMoistureStressModel{FT}(
            domain,
            toml_dict;
            soil_params = (; ν, θ_r),
        )
        @test θ_r == sms.θ_low
        @test ν == sms.θ_high
        βm_expected = FT(0.5)^sms.c
        βm_computed = compute_piecewise_moisture_stress(
            sms.θ_high,
            sms.θ_low,
            sms.c,
            θ_r + (ν - θ_r) / 2,
        )
        @test βm_computed ≈ βm_expected
    end

end
