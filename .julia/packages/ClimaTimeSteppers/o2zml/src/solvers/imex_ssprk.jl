struct IMEXSSPRKCache{U, SCI, B, Γ, NMC}
    U::U
    U_exp::U
    U_lim::U
    T_lim::U
    T_exp::U
    T_imp::SCI # sparse container of length s
    temp::U
    β::B
    γ::Γ
    newtons_method_cache::NMC
end

function init_cache(prob, alg::IMEXAlgorithm{SSP}; kwargs...)
    (; u0, f) = prob
    (; T_imp!) = f
    (; tableau, newtons_method) = alg
    (; a_exp, b_exp, a_imp, b_imp) = tableau
    s = length(b_exp)
    inds = ntuple(i -> i, s)
    inds_T_imp = filter(i -> !all(iszero, a_imp[:, i]) || !iszero(b_imp[i]), inds)
    U = zero(u0)
    U_exp = zero(u0)
    T_lim = zero(u0)
    T_exp = zero(u0)
    U_lim = zero(u0)
    T_imp = SparseContainer(map(i -> zero(u0), collect(1:length(inds_T_imp))), inds_T_imp)
    temp = zero(u0)
    â_exp = SparseCoeffs(vcat(a_exp.coeffs, b_exp.coeffs'))
    β = SparseCoeffs(diag(â_exp, -1))
    for i in 1:length(β)
        if â_exp.coeffs[(i + 1):end, i] != cumprod(β.coeffs[i:end])
            error("The SSP IMEXAlgorithm currently only supports an \
                   IMEXTableau that specifies a \"low-storage\" IMEX SSPRK \
                   algorithm, where the canonical Shu-Osher representation of \
                   the i-th explicit stage for i > 1 must have the form U[i] = \
                   (1 - β[i-1]) * u + β[i-1] * (U[i-1] + dt * T_exp(U[i-1])). \
                   So, it must be possible to express vcat(a_exp, b_exp') as\n \
                   0                  0           0    …\n \
                   β[1]               0           0    …\n \
                   β[1] * β[2]        β[2]        0    …\n \
                   β[1] * β[2] * β[3] β[2] * β[3] β[3] …\n \
                   ⋮                  ⋮            ⋮    ⋱\n \
                   The given IMEXTableau does not satisfy this property.")
        end
    end
    γs = unique(filter(!iszero, diag(a_imp)))
    γ = length(γs) == 1 ? γs[1] : nothing # TODO: This could just be a constant.
    jac_prototype = has_jac(T_imp!) ? T_imp!.jac_prototype : nothing
    newtons_method_cache =
        isnothing(T_imp!) || isnothing(newtons_method) ? nothing : allocate_cache(newtons_method, u0, jac_prototype)
    return IMEXSSPRKCache(U, U_exp, U_lim, T_lim, T_exp, T_imp, temp, β, γ, newtons_method_cache)
end

function step_u!(integrator, cache::IMEXSSPRKCache)
    (; u, p, t, dt, alg) = integrator
    (; f) = integrator.sol.prob
    (; cache!, cache_imp!) = f
    (; T_exp_T_lim!, T_lim!, T_exp!, T_imp!, lim!, dss!) = f
    (; tableau, newtons_method) = alg
    (; a_imp, b_imp, c_exp, c_imp) = tableau
    (; U, U_lim, U_exp, T_lim, T_exp, T_imp, temp, β, γ, newtons_method_cache) = cache
    s = length(b_imp)

    if !isnothing(T_imp!) && !isnothing(newtons_method)
        (; update_j) = newtons_method
        jacobian = newtons_method_cache.j
        if (!isnothing(jacobian)) && needs_update!(update_j, NewTimeStep(t))
            if γ isa Nothing
                sdirk_error(name)
            else
                T_imp!.Wfact(jacobian, u, p, dt * γ, t)
            end
        end
    end

    @. U = u

    for i in 1:s
        t_exp = t + dt * c_exp[i]
        t_imp = t + dt * c_imp[i]
        dtγ = dt * a_imp[i, i]

        if i == 1
            @. U_exp = u
        elseif !iszero(β[i - 1])
            if has_T_lim(f)
                @. U_lim = U_exp + dt * T_lim
                lim!(U_lim, p, t_exp, U_exp)
                @. U_exp = U_lim
            end
            if has_T_exp(f)
                @. U_exp += dt * T_exp
            end
            @. U_exp = (1 - β[i - 1]) * u + β[i - 1] * U_exp
        end

        @. U = U_exp
        if !isnothing(T_imp!) # Update based on implicit tendencies from previous stages
            for j in 1:(i - 1)
                iszero(a_imp[i, j]) && continue
                @. U += dt * a_imp[i, j] * T_imp[j]
            end
        end

        # Run the implicit solver, apply DSS, and update the cache. When γ == 0,
        # the implicit solver does not need to be run. On stage i == 1, we do
        # not need to apply DSS and update the cache because we did that at the
        # end of the previous timestep.
        i ≠ 1 && dss!(U, p, t_exp)
        if isnothing(T_imp!) || iszero(a_imp[i, i])
            i ≠ 1 && cache!(U, p, t_exp)
        else
            @assert !isnothing(newtons_method)
            i ≠ 1 && cache_imp!(U, p, t_imp)
            @. temp = U
            implicit_equation_residual! = (residual, U′) -> begin
                T_imp!(residual, U′, p, t_imp)
                @. residual = temp + dtγ * residual - U′
            end
            solve_newton!(
                newtons_method,
                newtons_method_cache,
                U,
                implicit_equation_residual!,
                (jacobian, U′) -> T_imp!.Wfact(jacobian, U′, p, dtγ, t_imp),
                U′ -> cache_imp!(U′, p, t_imp),
            )
            dss!(U, p, t_imp)
            cache!(U, p, t_imp)
        end

        if !iszero(β[i])
            isnothing(T_exp_T_lim!) || T_exp_T_lim!(T_lim, T_exp, U, p, t_exp)
            isnothing(T_lim!) || T_lim!(T_lim, U, p, t_exp)
            isnothing(T_exp!) || T_exp!(T_exp, U, p, t_exp)
        end
        if !all(iszero, a_imp[:, i]) || !iszero(b_imp[i])
            if iszero(a_imp[i, i])
                # When γ == 0, T_imp[i] is treated explicitly.
                isnothing(T_imp!) || T_imp!(T_imp[i], U, p, t_imp)
            else
                # When γ != 0, T_imp[i] is treated implicitly, so it must satisfy
                # the implicit equation. To ensure that T_imp[i] only includes the
                # effect of applying DSS to T_imp!, U and temp must be DSSed.
                isnothing(T_imp!) || @. T_imp[i] = (U - temp) / dtγ
            end
        end
    end

    t_final = t + dt

    if !iszero(β[s])
        if has_T_lim(f)
            @. U_lim = U_exp + dt * T_lim
            lim!(U_lim, p, t_final, U_exp)
            @. U_exp = U_lim
        end
        if has_T_exp(f)
            @. U_exp += dt * T_exp
        end
        @. u = (1 - β[s]) * u + β[s] * U_exp
    end

    if !isnothing(T_imp!) # Update based on implicit tendencies from previous stages
        for j in 1:s
            iszero(b_imp[j]) && continue
            @. u += dt * b_imp[j] * T_imp[j]
        end
    end

    dss!(u, p, t_final)
    cache!(u, p, t_final)

    return u
end
