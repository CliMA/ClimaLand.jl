using Test
import ForwardDiff
import SurfaceFluxes
import ClimaLand
import ClimaLand: Soil, Snow
import ClimaLand.Parameters as LP

# Adapted from Thermodynamics.jl's AD tests
function check_derivative(
    f,
    x;
    rtol = x isa Float32 ? 2e-2 : 1e-4,
    atol = 1e-8,
    ε = 1e-3,
)
    ad = ForwardDiff.derivative(f, x)
    fd = (f(x + ε) - f(x - ε)) / (2ε)
    @test isapprox(ad, fd; rtol, atol)
end

function check_inferred(f, x)
    result = @inferred f(x)
    for name in propertynames(result)
        @test getproperty(result, name) isa typeof(x)
    end
    @test @inferred(f(ForwardDiff.Dual(x, one(x)))) isa NamedTuple
end

# Check for type stability and whether the central finite difference is close to
# the derivative by ForwardDiff
# This assumes f(x) returns a NamedTuple
function check_ad(f, x)
    check_inferred(f, x)
    for name in keys(f(x))
        @testset "$name" begin
            check_derivative(y -> getproperty(f(y), name), x)
        end
    end
end

# Selects the default get_∂T_sfc∂T_function and get_∂q_sfc∂T_function, which
# the soil and snow models use
struct DummyModel{FT} <: ClimaLand.AbstractModel{FT} end

@testset "Soil skin-temperature solve, FT = $FT" for FT in (Float32, Float64)
    earth_param_set = LP.LandParameters(LP.create_toml_dict(FT))
    z_0m = FT(0.001)
    z_0b = FT(0.001)
    roughness_model = SurfaceFluxes.ConstantRoughnessParams{FT}(z_0m, z_0b)
    gustiness = SurfaceFluxes.FlooredDeardorffGustinessSpec(FT(1))
    T_top = FT(280)
    r = FT(0.0154)
    ϵ = FT(0.99)
    albedo = Soil.ConstantTwoBandSoilAlbedo{FT}()
    α_sfc = (albedo.PAR_albedo + albedo.NIR_albedo) / 2
    SW_n = -(1 - α_sfc) * FT(500)
    LW_d = FT(5.67e-8 * 280.0^4)
    g_liq = FT(2.98e-3)
    β_ice = FT(0)
    ψ_sfc = FT(-1.47)
    Tf_depressed = FT(273.138)
    h_sfc = FT(0)
    displ = FT(0)
    P_atmos = FT(101325)
    q_atmos = FT(0.005)
    u_atmos = FT(3)
    atmos_h = FT(3)
    update_∂T_sfc∂T =
        ClimaLand.get_∂T_sfc∂T_function(DummyModel{FT}(), nothing, nothing)
    update_∂q_sfc∂T =
        ClimaLand.get_∂q_sfc∂T_function(DummyModel{FT}(), nothing, nothing)

    soil_fluxes =
        T_atmos -> Soil.solve_soil_surface_temperature_at_a_point(
            Val(false),
            T_top,
            r,
            ϵ,
            SW_n,
            LW_d,
            g_liq,
            β_ice,
            ψ_sfc,
            Tf_depressed,
            h_sfc,
            displ,
            P_atmos,
            T_atmos,
            q_atmos,
            u_atmos,
            roughness_model,
            atmos_h,
            gustiness,
            update_∂T_sfc∂T,
            update_∂q_sfc∂T,
            earth_param_set,
        )
    check_ad(soil_fluxes, FT(285))

    @testset "Screen-level reconstruction" begin
        z_screen = ClimaLand.screen_height(FT)
        z_anemometer = ClimaLand.anemometer_height(FT)
        screen_values =
            T_atmos -> begin
                fluxes = soil_fluxes(T_atmos)
                ClimaLand.screen_level_values(
                    fluxes.T_sfc,
                    fluxes.q_sfc,
                    fluxes.ustar,
                    fluxes.ζ,
                    fluxes.Δz_eff,
                    z_0m,
                    z_0b,
                    T_atmos,
                    q_atmos,
                    z_screen,
                    z_anemometer,
                    earth_param_set,
                )
            end
        check_ad(screen_values, FT(285))
    end
end

@testset "Snow skin-temperature solve, FT = $FT" for FT in (Float32, Float64)
    earth_param_set = LP.LandParameters(LP.create_toml_dict(FT))
    z_0m = FT(0.0024)
    z_0b = FT(8e-2)
    roughness_model = SurfaceFluxes.ConstantRoughnessParams{FT}(z_0m, z_0b)
    gustiness = SurfaceFluxes.ConstantGustinessSpec(FT(1))
    T_initial_guess = FT(272)
    T_bulk = FT(270)
    κ_snow = FT(0.08)
    ρ_snow = FT(500)
    ϵ_snow = FT(0.97)
    α_snow = FT(0.8)
    SW_net = (α_snow - FT(1)) * FT(20)
    LW_d = FT(20)
    q_l = FT(0)
    h_sfc = FT(0)
    displ = FT(0)
    P_atmos = FT(101325)
    q_atmos = FT(0.003)
    u_atmos = FT(3)
    atmos_h = FT(1)
    update_∂T_sfc∂T =
        ClimaLand.get_∂T_sfc∂T_function(DummyModel{FT}(), nothing, nothing)
    update_∂q_sfc∂T =
        ClimaLand.get_∂q_sfc∂T_function(DummyModel{FT}(), nothing, nothing)
    surf_temp = Snow.EquilibriumGradientTemperatureModel{FT}()

    @testset "z_snow = $z_snow" for z_snow in (FT(0), FT(1))
        snow_fluxes =
            T_atmos -> Snow.solve_for_surface_temp_at_a_point(
                Val(false),
                T_initial_guess,
                T_bulk,
                z_snow,
                κ_snow,
                ρ_snow,
                ϵ_snow,
                SW_net,
                LW_d,
                q_l,
                h_sfc,
                displ,
                P_atmos,
                T_atmos,
                q_atmos,
                u_atmos,
                roughness_model,
                atmos_h,
                gustiness,
                update_∂T_sfc∂T,
                update_∂q_sfc∂T,
                earth_param_set,
                surf_temp,
            )
        check_ad(snow_fluxes, FT(245))
    end
end
