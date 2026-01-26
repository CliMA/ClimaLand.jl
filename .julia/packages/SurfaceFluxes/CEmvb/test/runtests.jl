using Test

import Random
Random.seed!(1234)
import KernelAbstractions: CPU

import Thermodynamics as TD
using SurfaceFluxes
import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF
import SurfaceFluxes.Parameters as SFP
import SurfaceFluxes.UniversalFunctions.BusingerParams
import ClimaParams as CP

device(::T) where {T <: Array} = CPU()

ArrayType = Array
@info "CPU Tests"
@info ArrayType
FloatType = Float32

@testset "SurfaceFluxes - Recovery Profiles" begin
    param_set = SFP.SurfaceFluxesParameters(FloatType, BusingerParams)
    thermo_params = param_set.thermo_params
    uft = UF.universal_func_type(typeof(param_set.ufp))
    ρ_sfc = FloatType(1.15)
    ρ_in = FloatType(1.13)
    qt_sfc = FloatType(0.01)
    qt_in = FloatType(0.009)
    ## Discretisation altitude z
    z = ArrayType(
        FloatType[
            29.432779269303,
            30.0497139076724,
            31.6880000418153,
            34.1873479240475,
        ],
    )
    ## Virtual Potential Temperature at height z
    θ = ArrayType(
        FloatType[
            268.559120403867,
            269.799228886728,
            277.443023238556,
            295.79192777341,
        ],
    )
    ## Surface Pottemp
    θ_sfc = ArrayType(
        FloatType[
            273.42369841804,
            272.551410044203,
            278.638168565727,
            298.133068766049,
        ],
    )
    ## Roughness lengths
    z0 = ArrayType(
        FloatType[
            5.86144925739178e-05,
            0.0001,
            0.000641655193293549,
            3.23383768877187e-05,
        ],
    )
    ## Speed
    speed = ArrayType(
        FloatType[
            2.9693638452068,
            2.43308757772094,
            5.69418282305367,
            9.5608693754561,
        ],
    )
    ## Scale velocity and moisture
    u_star = ArrayType(
        FloatType[
            0.109462510724615,
            0.0932942802513508,
            0.223232887323184,
            0.290918439028557,
        ],
    )
    # No explicit buoyancy terms in ClimateMachine
    b_star = ArrayType([
        0.00690834676781433,
        0.00428178089592372,
        0.00121229800895103,
        0.00262353784027441,
    ])

    κ = SFP.von_karman_const(param_set)
    for ii in 1:length(b_star)
        # Compute L_MO given u_star and b_star
        L_MO = u_star[ii]^2 / κ / b_star[ii]

        ts_sfc = TD.PhaseEquil_ρθq(thermo_params, ρ_sfc, θ_sfc[ii], qt_sfc)
        ts_in = TD.PhaseEquil_ρθq(thermo_params, ρ_in, θ[ii], qt_in)

        state_sfc =
            SF.StateValues(FloatType(0), (FloatType(0), FloatType(0)), ts_sfc)
        state_in =
            SF.StateValues(z[ii], (FloatType(speed[ii]), FloatType(0)), ts_in)

        # State containers
        z0m = z0[ii]
        z0b = FloatType(0.001)
        sc = (
            SF.Fluxes(
                state_in,
                state_sfc,
                FloatType(0),
                FloatType(0),
                z0m,
                z0b,
            ),
            SF.FluxesAndFrictionVelocity(
                state_in,
                state_sfc,
                FloatType(0),
                FloatType(0),
                u_star[ii],
                z0m,
                z0b,
            ),
            SF.ValuesOnly(state_in, state_sfc, z0m, z0b),
        )
        for jj in 1:length(sc)
            u_scale_fd = SF.compute_physical_scale_coeff(
                param_set,
                sc[jj],
                L_MO,
                UF.MomentumTransport(),
                uft,
                SF.PointValueScheme(),
            )
            Δu_fd = u_star[ii] / u_scale_fd
            u_scale_fv = SF.compute_physical_scale_coeff(
                param_set,
                sc[jj],
                L_MO,
                UF.MomentumTransport(),
                uft,
                SF.LayerAverageScheme(),
            )
            Δu_fv = u_star[ii] / u_scale_fv
            @test (Δu_fd - Δu_fv) ./ Δu_fd * 100 <= FloatType(50)
        end
    end
end


@testset "Near-zero Obukhov length (Floating Point Consistency)" begin
    FloatTypes = (Float32, Float64)
    z_levels = [1, 5, 10, 20, 40, 80, 160, 320, 640] # [m] level of first interior grid point
    for (i, FloatType) in enumerate(FloatTypes)
        sf_params = SFP.SurfaceFluxesParameters(FloatType, BusingerParams)
        for (jj, z_int) in enumerate(z_levels)
            ts_int_test = TD.PhaseEquil{FloatType}(
                1.1751807f0,
                97086.64f0,
                10541.609f0,
                0.0f0,
                287.85202f0,
            )
            ts_sfc_test = TD.PhaseEquil{FloatType}(
                1.2176297f0,
                102852.51f0,
                45087.812f0,
                0.013232904f0,
                291.96683f0,
            )
            sc = SF.ValuesOnly(
                SF.StateValues(
                    FloatType(z_int),
                    (FloatType(0), FloatType(0)),
                    ts_int_test,
                ),
                SF.StateValues(
                    FloatType(0),
                    (FloatType(0), FloatType(0)),
                    ts_sfc_test,
                ),
                FloatType(1e-5),
                FloatType(1e-5),
            )

            sfc_output = SF.surface_conditions(sf_params, sc; maxiter = 20)
            @test abs(SF.obukhov_length(sfc_output)) > FloatType(0)
            @test sign(SF.non_zero(1.0)) == 1
            @test sign(SF.non_zero(-1.0)) == -1
            @test sign(SF.non_zero(-0.0)) == 1
            @test sign(SF.non_zero(0.0)) == 1
        end
    end
end

@testset "Identical thermodynamic states (Floating Point Consistency)" begin
    FloatTypes = (Float32, Float64)
    z_levels = [1, 5, 10, 20, 40, 80, 160] # [m] level of first interior grid point
    z0_m = [1e-5, 1e-4, 1e-3] # roughness length [momentum]
    z0_b = [1e-5, 1e-4, 1e-3] # roughness length [heat] 
    sol_mat =
        Array{Any, 4}(undef, 2, length(z_levels), length(z0_m), length(z0_b))
    for (ii, FloatType) in enumerate(FloatTypes)
        sf_params = SFP.SurfaceFluxesParameters(FloatType, BusingerParams)
        for (jj, z_int) in enumerate(z_levels)
            for (kk, z0m) in enumerate(z0_m)
                for (ll, z0b) in enumerate(z0_b)
                    # Test case with identical interior and surface states
                    ts_int_test = TD.PhaseEquil{FloatType}(
                        1.1751807f0,
                        97086.64f0,
                        10541.609f0,
                        0.0f0,
                        287.85202f0,
                    )
                    ts_sfc_test = TD.PhaseEquil{FloatType}(
                        1.1751807f0,
                        97086.64f0,
                        10541.609f0,
                        0.0f0,
                        287.85202f0,
                    )
                    sc = SF.ValuesOnly(
                        SF.StateValues(
                            FloatType(z_int),
                            (FloatType(0), FloatType(0)),
                            ts_int_test,
                        ),
                        SF.StateValues(
                            FloatType(0),
                            (FloatType(0), FloatType(0)),
                            ts_sfc_test,
                        ),
                        FloatType(z0m),
                        FloatType(z0b),
                    )
                    sfc_output = SF.surface_conditions(
                        sf_params,
                        sc;
                        maxiter = 20,
                        noniterative_stable_sol = true,
                    )
                    sol_mat[ii, jj, kk, ll] =
                        isinf(sfc_output.L_MO) ? FloatType(1e6) :
                        sfc_output.L_MO
                end
            end
        end
    end
    rdiff_sol =
        (sol_mat[1, :, :, :] .- sol_mat[2, :, :, :]) ./ sol_mat[2, :, :, :]
    @test all(x -> x <= FloatType(0.005), abs.(rdiff_sol))
end

@testset "Test profiles" begin
    include("test_profiles.jl")
end
@testset "Test Float32 solver convergence" begin
    include("test_cases.jl")
end
@testset "Test universal functions" begin
    include("test_universal_functions.jl")
end
@testset "Test generated thermodynamic states" begin
    include("test_convergence.jl")
end

@testset "Quality assurance" begin
    include("aqua.jl")
end
