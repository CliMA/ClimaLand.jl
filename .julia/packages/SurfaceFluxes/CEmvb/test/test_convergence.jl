using Test

import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF
import SurfaceFluxes.Parameters.SurfaceFluxesParameters
import ClimaParams as CP

using Statistics
using Thermodynamics
using Thermodynamics.TemperatureProfiles
using Thermodynamics.TestedProfiles
Thermodynamics.print_warning() = false

using RootSolvers
const RS = RootSolvers

import Thermodynamics.TestedProfiles: input_config, PhaseEquilProfiles

abstract type TestProfiles end
struct DryProfiles <: TestProfiles end
struct MoistEquilProfiles <: TestProfiles end

function input_config(
    ArrayType;
    n = 5,
    n_RS1 = 5,
    n_RS2 = 5,
    T_surface = 290,
    T_min = 150,
)
    n_RS = n_RS1 + n_RS2
    z_range = ArrayType(range(0, stop = 80, length = n))
    relative_sat1 = ArrayType(range(0, stop = 1, length = n_RS1))
    relative_sat2 = ArrayType(range(1, stop = 1.02, length = n_RS2))
    relative_sat = vcat(relative_sat1, relative_sat2)
    return z_range, relative_sat, T_surface, T_min
end

function generate_profiles(FT, ::DryProfiles; uf_params = UF.BusingerParams)
    param_set = SurfaceFluxesParameters(FT, uf_params)
    thermo_params = param_set.thermo_params
    profiles = collect(
        Thermodynamics.TestedProfiles.PhaseDryProfiles(
            thermo_params,
            Array{FT},
        ),
    )
    profiles_sfc = filter(p -> iszero(p.z), profiles)
    profiles_int = filter(p -> !iszero(p.z), profiles)
    ## Properties contained in `profiles_<sfc, int>`
    ## :z, :T, :p, :RS, :e_int, :h, :ρ, 
    ## :θ_liq_ice, :q_tot, :q_liq, :q_ice, :q_pt, :RH, 
    ## :e_pot, :u, :v, :w, :e_kin, :phase_type
    return profiles_sfc, profiles_int, param_set
end

function generate_profiles(
    FT,
    ::MoistEquilProfiles;
    uf_params = UF.BusingerParams,
)
    param_set = SurfaceFluxesParameters(FT, uf_params)
    thermo_params = param_set.thermo_params
    profiles = collect(
        Thermodynamics.TestedProfiles.PhaseEquilProfiles(
            thermo_params,
            Array{FT},
        ),
    )
    profiles_sfc = filter(p -> iszero(p.z), profiles)
    profiles_int = filter(p -> !iszero(p.z), profiles)
    ## Properties contained in `profiles_<sfc, int>`
    ## :z, :T, :p, :RS, :e_int, :h, :ρ, 
    ## :θ_liq_ice, :q_tot, :q_liq, :q_ice, :q_pt, :RH, 
    ## :e_pot, :u, :v, :w, :e_kin, :phase_type
    return profiles_sfc, profiles_int, param_set
end

function compute_dse_diff(thermo_params, ts_int, ts_sfc, prof_int, prof_sfc)
    DSEᵥ_int =
        TD.virtual_dry_static_energy(thermo_params, ts_int, prof_int.e_pot)
    DSEᵥ_sfc =
        TD.virtual_dry_static_energy(thermo_params, ts_sfc, prof_sfc.e_pot)
    return DSEᵥ_int - DSEᵥ_sfc
end

function assemble_surface_conditions(
    prof_int,
    prof_sfc,
    ts_int,
    ts_sfc,
    z0m,
    z0b,
)
    FT = eltype(z0m)
    state_in =
        SF.StateValues(prof_int.z, (prof_int.u / 10, prof_int.v / 10), ts_int)
    state_sfc = SF.StateValues(FT(0), (FT(0), FT(0)), ts_sfc)
    return SF.ValuesOnly(state_in, state_sfc, z0m, z0b)
end

function check_over_dry_states(
    ::Type{FT},
    profiles_int,
    profiles_sfc,
    scheme,
    z0_momentum,
    z0_thermal,
    maxiter,
    tol_neutral,
    gryanik_noniterative::Bool,
) where {FT}
    param_set = SFP.SurfaceFluxesParameters(FT, UF.BusingerParams)
    counter = [0, 0, 0] # Stable, Unstable, Neutral
    @inbounds for (ii, prof_int) in enumerate(profiles_int)
        @inbounds for (jj, prof_sfc) in enumerate(profiles_sfc)
            @inbounds for (kk, z0m) in enumerate(z0_momentum)
                @inbounds for (ll, z0b) in enumerate(z0_thermal)
                    @inbounds for sch in scheme
                        if z0m / z0b >= FT(1e2) || z0b / z0m >= FT(1e2)
                            nothing
                        else
                            thermo_params =
                                SF.Parameters.thermodynamics_params(param_set)
                            ts_sfc = Thermodynamics.PhaseDry{FT}(
                                prof_sfc.e_int,
                                prof_sfc.ρ,
                            )
                            ts_int = Thermodynamics.PhaseDry{FT}(
                                prof_int.e_int,
                                prof_int.ρ,
                            )
                            sc = assemble_surface_conditions(
                                prof_int,
                                prof_sfc,
                                ts_int,
                                ts_sfc,
                                z0m,
                                z0b,
                            )
                            ΔDSEᵥ = compute_dse_diff(
                                thermo_params,
                                ts_int,
                                ts_sfc,
                                prof_int,
                                prof_sfc,
                            )
                            if abs(ΔDSEᵥ) <= tol_neutral
                                counter[3] += 1
                            else
                                sign(ΔDSEᵥ) == 1 ? counter[1] += 1 :
                                counter[2] += 1
                            end
                            sfcc = SF.surface_conditions(
                                param_set,
                                sc,
                                sch;
                                maxiter,
                                tol_neutral,
                                soltype = RS.VerboseSolution(),
                                noniterative_stable_sol = gryanik_noniterative,
                            )
                            if abs(ΔDSEᵥ) <= tol_neutral &&
                               gryanik_noniterative == false
                                @test isinf(sfcc.L_MO)
                            else
                                @test sign.(sfcc.ρτxz) == -sign(prof_int.u)
                                @test sign.(sfcc.ρτyz) == -sign(prof_int.v)
                                @test sign(sfcc.L_MO) == sign(ΔDSEᵥ)
                                @test sign(sfcc.shf) == -sign(ΔDSEᵥ)
                                @test sign(
                                    SF.compute_bstar(
                                        param_set,
                                        sfcc.L_MO,
                                        sc,
                                        SF.Parameters.universal_func_type(
                                            param_set,
                                        ),
                                        sch,
                                    ),
                                ) == sign(ΔDSEᵥ)
                                @test sign(sfcc.buoy_flux) == -sign(ΔDSEᵥ)
                            end

                        end
                    end
                end
            end
        end
    end
    return counter
end

function check_over_moist_states(
    ::Type{FT},
    profiles_int,
    profiles_sfc,
    scheme,
    z0_momentum,
    z0_thermal,
    maxiter,
    tol_neutral,
    gryanik_noniterative::Bool,
) where {FT}
    param_set = SFP.SurfaceFluxesParameters(FT, BusingerParams)
    counter = [0, 0, 0] # Stable, Unstable, Neutral
    @inbounds for (ii, prof_int) in enumerate(profiles_int)
        @inbounds for (jj, prof_sfc) in enumerate(profiles_sfc)
            @inbounds for (kk, z0m) in enumerate(z0_momentum)
                @inbounds for (ll, z0b) in enumerate(z0_thermal)
                    @inbounds for sch in scheme
                        if z0m / z0b >= FT(1e2) || z0b / z0m >= FT(1e2)
                            nothing
                        else
                            thermo_params =
                                SF.Parameters.thermodynamics_params(param_set)
                            ts_sfc = Thermodynamics.PhaseEquil{FT}(
                                prof_sfc.ρ,
                                prof_sfc.p,
                                prof_sfc.e_int,
                                prof_sfc.q_tot,
                                prof_sfc.T,
                            )
                            ts_int = Thermodynamics.PhaseEquil{FT}(
                                prof_int.ρ,
                                prof_int.p,
                                prof_int.e_int,
                                prof_int.q_tot,
                                prof_int.T,
                            )
                            sc = assemble_surface_conditions(
                                prof_int,
                                prof_sfc,
                                ts_int,
                                ts_sfc,
                                z0m,
                                z0b,
                            )
                            ΔDSEᵥ = compute_dse_diff(
                                thermo_params,
                                ts_int,
                                ts_sfc,
                                prof_int,
                                prof_sfc,
                            )
                            if abs(ΔDSEᵥ) <= tol_neutral
                                counter[3] += 1
                            else
                                sign(ΔDSEᵥ) == 1 ? counter[1] += 1 :
                                counter[2] += 1
                            end
                            sfcc = SF.surface_conditions(
                                param_set,
                                sc,
                                sch;
                                maxiter,
                                tol_neutral,
                                soltype = RS.VerboseSolution(),
                                noniterative_stable_sol = gryanik_noniterative,
                            )
                            if abs(ΔDSEᵥ) <= tol_neutral &&
                               gryanik_noniterative == false
                                @test isinf(sfcc.L_MO)
                            else
                                @test sign.(sfcc.evaporation) ==
                                      -sign(prof_int.q_tot - prof_sfc.q_tot)
                                @test sign.(sfcc.ρτxz) == -sign(prof_int.u)
                                @test sign.(sfcc.ρτyz) == -sign(prof_int.v)
                                @test sign(sfcc.L_MO) == sign(ΔDSEᵥ)
                                @test sign(
                                    SF.compute_bstar(
                                        param_set,
                                        sfcc.L_MO,
                                        sc,
                                        SF.Parameters.universal_func_type(
                                            param_set,
                                        ),
                                        sch,
                                    ),
                                ) == sign(ΔDSEᵥ)
                            end
                        end
                    end
                end
            end
        end
    end
    return counter
end

@testset "Check convergence (dry thermodynamic states): Stable/Unstable" begin
    for FT in [Float32, Float64]
        for uf_params in [UF.BusingerParams, UF.GryanikParams, UF.GrachevParams]
            for profile_type in [DryProfiles(), MoistEquilProfiles()]
                profiles_sfc, profiles_int =
                    generate_profiles(FT, profile_type; uf_params)
                scheme = [SF.LayerAverageScheme(), SF.PointValueScheme()]
                z0_momentum = Array{FT}(range(1e-6, stop = 1e-1, length = 2))
                z0_thermal = Array{FT}(range(1e-6, stop = 1e-1, length = 2))
                maxiter = 10
                tol_neutral = FT(SF.Parameters.cp_d(param_set) / 10)
                for iteration_option in [true]
                    counter = check_over_dry_states(
                        FT,
                        profiles_int,
                        profiles_sfc,
                        scheme,
                        z0_momentum,
                        z0_thermal,
                        maxiter,
                        tol_neutral,
                        iteration_option,
                    )
                    counter = check_over_moist_states(
                        FT,
                        profiles_int,
                        profiles_sfc,
                        scheme,
                        z0_momentum,
                        z0_thermal,
                        maxiter,
                        tol_neutral,
                        iteration_option,
                    )
                end
            end
        end
    end
    @info "Tested Businger/Gryanik/Grachev u.f., Float32/Float64, Dry/EquilMoist profiles, non-iterative solver on/off"

end
