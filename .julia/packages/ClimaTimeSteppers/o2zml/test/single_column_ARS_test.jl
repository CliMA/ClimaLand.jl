using DiffEqBase, ClimaTimeSteppers, LinearAlgebra, StaticArrays, Test
import ClimaTimeSteppers as CTS
# Unit test for ARS schemes on the single column hydrostatic problem
# Here we run 1 time step and 1 Newton iteration to compare the result
mutable struct Single_Stack
    L::Float64
    N::Int64
    t::Float64
    Δz::Array{Float64, 1}

    ρ::Array{Float64, 1}
    w::Array{Float64, 1}
    ρθ::Array{Float64, 1}

    ρ_res::Array{Float64, 1}
    w_res::Array{Float64, 1}
    ρθ_res::Array{Float64, 1}

    _grav::Float64
    γ::Float64
    _R_m::Float64

    P::Function
    dP::Function
end

function combine!(ρ, w, ρθ, u)
    N = length(ρ)
    u[1:N] .= ρ
    u[(N + 1):(2N)] .= ρθ
    u[(2N + 1):(3N + 1)] .= w
end

function split!(u, ρ, w, ρθ)

    N = length(ρ)
    ρ .= u[1:N]
    ρθ .= u[(N + 1):(2N)]
    w .= u[(2N + 1):(3N + 1)]
end

function DecayingTemperatureProfile(z::Float64; T_virt_surf, T_min_ref, _R_m, _grav, _MSLP)
    # Scale height for surface temperature
    H_sfc = _R_m * T_virt_surf / _grav
    H_t = H_sfc

    z′ = z / H_t
    tanh_z′ = tanh(z′)

    ΔTv = T_virt_surf - T_min_ref
    Tv = T_virt_surf - ΔTv * tanh_z′

    ΔTv′ = ΔTv / T_virt_surf
    p = -H_t * (z′ + ΔTv′ * (log(1 - ΔTv′ * tanh_z′) - log(1 + tanh_z′) + z′))
    p /= H_sfc * (1 - ΔTv′^2)
    p = _MSLP * exp(p)
    ρ = p / (_R_m * Tv)
    return (Tv, p, ρ)
end

function spatial_residual!(du, u, ss::Single_Stack, t)
    # u is a 3N+1 vector
    N, L = ss.N, ss.L
    Δz = L / N
    _grav = ss._grav
    P = ss.P

    ρ, w, ρθ = ss.ρ, ss.w, ss.ρθ
    ρ_res, w_res, ρθ_res = ss.ρ_res, ss.w_res, ss.ρθ_res

    split!(u, ρ, w, ρθ)


    for i in 1:N
        ρ₊ = (i == N ? 0.0 : (ρ[i] + ρ[i + 1]) / 2.0)
        ρ₋ = (i == 1 ? 0.0 : (ρ[i] + ρ[i - 1]) / 2.0)

        ρ_res[i] = -(w[i + 1]ρ₊ - w[i]ρ₋) / Δz

        ρθ₊ = (i == N ? 0.0 : (ρθ[i] + ρθ[i + 1]) / 2.0)
        ρθ₋ = (i == 1 ? 0.0 : (ρθ[i] + ρθ[i - 1]) / 2.0)

        ρθ_res[i] = -(w[i + 1]ρθ₊ - w[i]ρθ₋) / Δz
    end

    w_res[1] = w_res[N + 1] = 0
    for i in 2:N
        w_res[i] = -(P(ρθ[i]) - P(ρθ[i - 1])) / (1 / 2 * (ρ[i - 1] + ρ[i])) / Δz - _grav
    end

    combine!(ρ_res, w_res, ρθ_res, du)
end

function jacobian!(J, u, ss, t)
    P, dP = ss.P, ss.dP
    Δz = ss.Δz

    # N cells
    N = div(length(u) - 1, 3)

    ρ, ρθ, w = u[1:N], u[(N + 1):(2N)], u[(2N + 1):(3N + 1)]

    # interpolate ρ and ρθ to cell face values
    ρh = [ρ[1]; (ρ[1:(N - 1)] + ρ[2:N]) / 2.0; ρ[N]]
    ρθh = [ρθ[1]; (ρθ[1:(N - 1)] + ρθ[2:N]) / 2.0; ρθ[N]]

    # pressure and its derivative with respect to ρθ at cell center
    pc = P.(ρθ)
    dpc_dρθ = dP.(ρθ)

    # half cell size
    Δzh = [NaN64; (Δz[1:(N - 1)] + Δz[2:N]) / 2.0; NaN64]

    # P = ([zeros(N,N)     zeros(N,N)      ∂ρ/∂w ;
    #       zeros(N,N)     zeros(N,N)      ∂ρθ/∂w
    #       ∂w/∂ρ          ∂w/∂ρθ      zeros(N+1,N+1)])

    for i in 1:N
        J[i, i + 2N] = ρh[i] / Δz[i]
        J[i, i + 2N + 1] = -ρh[i + 1] / Δz[i]
    end

    for i in 1:N
        J[i + N, i + 2N] = ρθh[i] / Δz[i]
        J[i + N, i + 2N + 1] = -ρθh[i + 1] / Δz[i]
    end

    for i in 2:N #(dw_i+1/2)
        J[i + 2N, (i - 1)] = (pc[i] - pc[i - 1]) / Δzh[i] / (2 * ρh[i]^2)
        J[i + 2N, (i - 1) + 1] = (pc[i] - pc[i - 1]) / Δzh[i] / (2 * ρh[i]^2)

        J[i + 2N, (i - 1) + N] = dpc_dρθ[i - 1] ./ (ρh[i] * Δzh[i])
        J[i + 2N, (i - 1) + 1 + N] = -dpc_dρθ[i] ./ (ρh[i] * Δzh[i])
    end
end

function single_column_Wfact!(J, u, ss, γ, t)
    J .= 0.0
    jacobian!(J, u, ss, t)
    J .*= γ
    N = ss.N
    for i in 1:(3N + 1)
        J[i, i] -= 1.0
    end
    return J
end

function linsolve_direct(::Type{Val{:init}}, f, u0; kwargs...)
    function _linsolve!(x, A, b, update_matrix = false; kwargs...)
        x .= A \ b
    end
end

@testset "IMEX ARS schemes on a single column" begin
    N = 30
    L = 30e3
    N_iter = 1
    dt = 100.0

    _MSLP = 1e5
    _grav = 9.80616
    _R_m = 287.0024093890231
    γ = 1.4
    _C_p = _R_m * γ / (γ - 1)
    _C_v = _R_m / (γ - 1)


    P(ρθ) = _MSLP^(-_R_m / _C_v) * _R_m^(_C_p / _C_v) * (ρθ)^(_C_p / _C_v)
    dP(ρθ) = _C_p / _C_v * _MSLP^(-_R_m / _C_v) * _R_m^(_C_p / _C_v) * (ρθ)^(_C_p / _C_v - 1)


    # Initialization
    ρ = zeros(Float64, N)
    ρθ = zeros(Float64, N)
    w = zeros(Float64, N + 1)
    Δz = zeros(N + 1) .+ L / N
    for i in 1:N
        z = (i - 0.5) * L / N
        Tvᵢ, pᵢ, ρᵢ = DecayingTemperatureProfile(
            z;
            T_virt_surf = 300.0,
            T_min_ref = 300.0,
            _R_m = _R_m,
            _grav = _grav,
            _MSLP = _MSLP,
        )
        ρ[i] = ρᵢ
        ρθ[i] = ρᵢ * Tvᵢ * (_MSLP / pᵢ)^(_R_m / _C_p)
    end
    ss = Single_Stack(L, N, 0.0, Δz, copy(ρ), copy(w), copy(ρθ), copy(ρ), copy(w), copy(ρθ), _grav, γ, _R_m, P, dP)
    u0 = zeros(3 * N + 1)
    combine!(ρ, w, ρθ, u0)






    # Generate reference values
    begin
        N_iter = 1
        function f_exp(du, u)
            du .= 0
        end
        function f_imp(du, u)
            spatial_residual!(du, u, ss, 0.0)

        end
        function jac_imp(J, u)
            jacobian!(J, u, ss, 0.0)
        end

        function f_imp_solve(f_imp, jac_imp, u0, γ)
            N = 30
            k0 = zeros(3 * N + 1)
            J0 = zeros(3 * N + 1, 3 * N + 1)
            f_imp(k0, u0)
            jac_imp(J0, u0)

            k = (I - γ * J0) \ k0
            return k
        end

        # TIME_INTEGRATOR == "ARS343"
        γ = 0.4358665215084590
        a42 = 0.5529291480359398
        a43 = 0.5529291480359398
        b1 = -3 / 2 * γ^2 + 4 * γ - 1 / 4
        b2 = 3 / 2 * γ^2 - 5 * γ + 5 / 4
        a31 =
            (1 - 9 / 2 * γ + 3 / 2 * γ^2) * a42 + (11 / 4 - 21 / 2 * γ + 15 / 4 * γ^2) * a43 - 7 / 2 + 13 * γ -
            9 / 2 * γ^2
        a32 =
            (-1 + 9 / 2 * γ - 3 / 2 * γ^2) * a42 + (-11 / 4 + 21 / 2 * γ - 15 / 4 * γ^2) * a43 + 4 - 25 / 2 * γ +
            9 / 2 * γ^2
        a41 = 1 - a42 - a43
        a_exp = [
            0 0 0 0
            γ 0 0 0
            a31 a32 0 0
            a41 a42 a43 0
        ]
        b_exp = [0, b1, b2, γ]
        a_imp = [
            0 0 0 0
            0 γ 0 0
            0 (1 - γ)/2 γ 0
            0 b1 b2 γ
        ]
        b_imp = [0, b1, b2, γ]


        n_stage = size(a_imp)[1] - 1
        K_exp = zeros(3 * N + 1, n_stage + 1)
        K_imp = zeros(3 * N + 1, n_stage + 1)
        u = copy(u0)
        for i in 1:N_iter
            u_stage = copy(u)
            for i_stage in 1:n_stage
                f_exp(@view(K_exp[:, i_stage]), u_stage)
                u_stage_temp =
                    u +
                    dt * (
                        K_imp[:, 1:(i_stage - 1)] * a_imp[i_stage + 1, 2:i_stage] +
                        K_exp[:, 1:i_stage] * a_exp[i_stage + 1, 1:i_stage]
                    )

                K_imp[:, i_stage] = f_imp_solve(f_imp, jac_imp, u_stage_temp, dt * a_imp[i_stage + 1, i_stage + 1])
                u_stage = u_stage_temp + dt * K_imp[:, i_stage] * a_imp[i_stage + 1, i_stage + 1]
            end
            f_exp(@view(K_exp[:, n_stage + 1]), u_stage)
            u .+=
                dt * (K_exp[:, 1:(n_stage + 1)] * b_exp[1:(n_stage + 1)] + K_imp[:, 1:n_stage] * b_imp[2:(n_stage + 1)])
        end
        ref_ARS343 = norm(u)
    end






    algorithms = (
        CTS.IMEXAlgorithm(ARS111(), NewtonsMethod(; max_iters = 1)),
        CTS.IMEXAlgorithm(ARS121(), NewtonsMethod(; max_iters = 1)),
        CTS.IMEXAlgorithm(ARS122(), NewtonsMethod(; max_iters = 1)),
        CTS.IMEXAlgorithm(ARS232(), NewtonsMethod(; max_iters = 1)),
        CTS.IMEXAlgorithm(ARS222(), NewtonsMethod(; max_iters = 1)),
        CTS.IMEXAlgorithm(ARS343(), NewtonsMethod(; max_iters = 1)),
    )
    reference_sol_norm = [
        860.2745315698107
        860.2745315698107
        860.4393569534262
        860.452530117785
        860.452530117785
        ref_ARS343
    ]



    # ARS solve
    for (i, algo) in (enumerate(algorithms))

        func_kwargs = (; jac_prototype = zeros(Float64, 3N + 1, 3N + 1), Wfact = single_column_Wfact!)
        tendency_func = ClimaODEFunction(;
            T_exp! = ODEFunction((du, u, p, t) -> (du .= 0.0)),
            T_imp! = ODEFunction(spatial_residual!; func_kwargs...),
        )

        single_column_prob_wfact_split = ODEProblem(tendency_func, copy(u0), (0.0, N_iter * dt), ss)

        u = solve(single_column_prob_wfact_split, algo; dt = dt)
        # @info norm(u.u[end])
        @test norm(u.u[end]) ≈ reference_sol_norm[i] atol = 2e3eps()
    end

end
