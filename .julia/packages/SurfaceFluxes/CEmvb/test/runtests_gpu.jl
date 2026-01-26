# Testing and Device
using Test
#using ClimaComms
using CUDA
using KernelAbstractions

# Physics
import Thermodynamics as TD
import SurfaceFluxes as SF
import SurfaceFluxes.Parameters as SFP
import SurfaceFluxes.UniversalFunctions.BusingerParams
import ClimaParams as CP

ArrayType = CUDA.CuArray

@info "GPU Compatibility Tests"
@info ArrayType

@kernel function test_surface_conditions_kernel!(
    param_set,
    output::AbstractArray{FT},
    z,
    θ,
    θ_sfc,
    z0,
    speed,
    u_star,
    b_star,
    κ,
) where {FT}
    ii = @index(Global)
    thermo_params = SFP.thermodynamics_params(param_set)
    # Assumption on surface density
    ρ_sfc = FT(1.15)
    ρ_in = FT(1.13)
    qt_sfc = FT(0.01)
    qt_in = FT(0.009)
    @inbounds begin
        L_MO = u_star[ii]^2 / κ / b_star[ii]
        ts_sfc = TD.PhaseEquil_ρθq(thermo_params, ρ_sfc, θ_sfc[ii], qt_sfc)
        ts_in = TD.PhaseEquil_ρθq(thermo_params, ρ_in, θ[ii], qt_in)
        state_sfc = SF.StateValues(FT(0), (FT(0), FT(0)), ts_sfc)
        state_in = SF.StateValues(z[ii], (FT(speed[ii]), FT(0)), ts_in)
        z0m = z0[ii]
        z0b = FT(0.001)
        sc = (SF.ValuesOnly(state_in, state_sfc, z0m, z0b),)
        output[1, ii] = SF.surface_conditions(param_set, sc[ii]).L_MO
    end
end

function test_surfacefluxes_gpu(FT)
    param_set = SFP.SurfaceFluxesParameters(FT, BusingerParams)
    data_length = 4
    output = ArrayType{FT}(undef, 1, data_length)
    fill!(output, FT(-99999.99))

    ndrange = (data_length,)
    z = ArrayType([
        FT(29.432779269303),
        FT(30.0497139076724),
        FT(31.6880000418153),
        FT(34.1873479240475),
    ])
    θ = ArrayType([
        FT(268.559120403867),
        FT(269.799228886728),
        FT(277.443023238556),
        FT(295.79192777341),
    ])
    θ_sfc = ArrayType([
        FT(273.42369841804),
        FT(272.551410044203),
        FT(278.638168565727),
        FT(298.133068766049),
    ])
    z0 = ArrayType([
        FT(5.86144925739178e-05),
        FT(0.0001),
        FT(0.000641655193293549),
        FT(3.23383768877187e-05),
    ])
    speed = ArrayType([
        FT(2.9693638452068),
        FT(2.43308757772094),
        FT(5.69418282305367),
        FT(9.5608693754561),
    ])
    u_star = ArrayType([
        FT(0.109462510724615),
        FT(0.0932942802513508),
        FT(0.223232887323184),
        FT(0.290918439028557),
    ])
    b_star = ArrayType([
        FT(0.00690834676781433),
        FT(0.00428178089592372),
        FT(0.00121229800895103),
        FT(0.00262353784027441),
    ])
    κ = SFP.von_karman_const(param_set)

    backend = KernelAbstractions.get_backend(output)
    kernel! = test_surface_conditions_kernel!(backend)
    event = kernel!(
        param_set,
        output,
        z,
        θ,
        θ_sfc,
        z0,
        speed,
        u_star,
        b_star,
        κ,
        ndrange = ndrange,
    )

    KernelAbstractions.synchronize(backend)
    return output
end

out_float32 = test_surfacefluxes_gpu(Float32)
out_float64 = test_surfacefluxes_gpu(Float64)
@test !any(x -> x == -99999.99, out_float32)
@test !any(x -> x == -99999.99, out_float64)
@test out_float32 ≈ out_float32
