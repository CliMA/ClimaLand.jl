import NVTX

has_jac(T_imp!) =
    hasfield(typeof(T_imp!), :Wfact) &&
    hasfield(typeof(T_imp!), :jac_prototype) &&
    !isnothing(T_imp!.Wfact) &&
    !isnothing(T_imp!.jac_prototype)

sdirk_error(name) = error("$(isnothing(name) ? "The given IMEXTableau" : name) \
                           has implicit stages with distinct coefficients (it \
                           is not SDIRK), and an update is required whenever a \
                           stage has a different coefficient from the previous \
                           stage. Do not update on the NewTimeStep signal when \
                           using $(isnothing(name) ? "this tableau" : name).")

struct IMEXARKCache{SCU, SCE, SCI, T, Γ, NMC}
    U::SCU     # sparse container of length s
    T_lim::SCE # sparse container of length s
    T_exp::SCE # sparse container of length s
    T_imp::SCI # sparse container of length s
    temp::T
    γ::Γ
    newtons_method_cache::NMC
end

function init_cache(prob, alg::IMEXAlgorithm{Unconstrained}; kwargs...)
    (; u0, f) = prob
    (; T_imp!) = f
    (; tableau, newtons_method) = alg
    (; a_exp, b_exp, a_imp, b_imp) = tableau
    s = length(b_exp)
    inds = ntuple(i -> i, s)
    inds_T_exp = filter(i -> !all(iszero, a_exp[:, i]) || !iszero(b_exp[i]), inds)
    inds_T_imp = filter(i -> !all(iszero, a_imp[:, i]) || !iszero(b_imp[i]), inds)
    U = zero(u0)
    T_lim = SparseContainer(map(i -> zero(u0), collect(1:length(inds_T_exp))), inds_T_exp)
    T_exp = SparseContainer(map(i -> zero(u0), collect(1:length(inds_T_exp))), inds_T_exp)
    T_imp = SparseContainer(map(i -> zero(u0), collect(1:length(inds_T_imp))), inds_T_imp)
    temp = zero(u0)
    γs = unique(filter(!iszero, diag(a_imp)))
    γ = length(γs) == 1 ? γs[1] : nothing # TODO: This could just be a constant.
    jac_prototype = has_jac(T_imp!) ? T_imp!.jac_prototype : nothing
    newtons_method_cache =
        isnothing(T_imp!) || isnothing(newtons_method) ? nothing : allocate_cache(newtons_method, u0, jac_prototype)
    return IMEXARKCache(U, T_lim, T_exp, T_imp, temp, γ, newtons_method_cache)
end

# generic fallback
function step_u!(integrator, cache::IMEXARKCache)
    (; u, p, t, dt, alg) = integrator
    (; f) = integrator.sol.prob
    (; cache!, cache_imp!) = f
    (; T_lim!, T_exp!, T_imp!, lim!, dss!) = f
    (; tableau, newtons_method) = alg
    (; a_exp, b_exp, a_imp, b_imp, c_exp, c_imp) = tableau
    (; U, T_lim, T_exp, T_imp, temp, γ, newtons_method_cache) = cache
    s = length(b_exp)

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

    update_stage!(integrator, cache, ntuple(i -> i, Val(s)))

    t_final = t + dt

    if has_T_lim(f) # Update based on limited tendencies from previous stages
        assign_fused_increment!(temp, u, dt, b_exp, T_lim, Val(s))
        lim!(temp, p, t_final, u)
        @. u = temp
    end

    # Update based on tendencies from previous stages
    has_T_exp(f) && fused_increment!(u, dt, b_exp, T_exp, Val(s))
    isnothing(T_imp!) || fused_increment!(u, dt, b_imp, T_imp, Val(s))

    dss!(u, p, t_final)
    cache!(u, p, t_final)

    return u
end


@inline update_stage!(integrator, cache, ::Tuple{}) = nothing
@inline update_stage!(integrator, cache, is::Tuple{Int}) = update_stage!(integrator, cache, first(is))
@inline function update_stage!(integrator, cache, is::Tuple)
    update_stage!(integrator, cache, first(is))
    update_stage!(integrator, cache, Base.tail(is))
end
@inline function update_stage!(integrator, cache::IMEXARKCache, i::Int)
    (; u, p, t, dt, alg) = integrator
    (; f) = integrator.sol.prob
    (; cache!, cache_imp!) = f
    (; T_exp_T_lim!, T_lim!, T_exp!, T_imp!, lim!, dss!) = f
    (; tableau, newtons_method) = alg
    (; a_exp, b_exp, a_imp, b_imp, c_exp, c_imp) = tableau
    (; U, T_lim, T_exp, T_imp, temp, newtons_method_cache) = cache
    s = length(b_exp)

    t_exp = t + dt * c_exp[i]
    t_imp = t + dt * c_imp[i]
    dtγ = float(dt) * a_imp[i, i]

    if has_T_lim(f) # Update based on limited tendencies from previous stages
        assign_fused_increment!(U, u, dt, a_exp, T_lim, Val(i))
        i ≠ 1 && lim!(U, p, t_exp, u)
    else
        @. U = u
    end

    # Update based on tendencies from previous stages
    has_T_exp(f) && fused_increment!(U, dt, a_exp, T_exp, Val(i))
    isnothing(T_imp!) || fused_increment!(U, dt, a_imp, T_imp, Val(i))

    # Run the implicit solver, apply DSS, and update the cache. When γ == 0,
    # the implicit solver does not need to be run. On stage i == 1, we do not
    # need to apply DSS and update the cache because we did that at the end of
    # the previous timestep.
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

    if !all(iszero, a_exp[:, i]) || !iszero(b_exp[i])
        isnothing(T_exp_T_lim!) || T_exp_T_lim!(T_exp[i], T_lim[i], U, p, t_exp)
        isnothing(T_lim!) || T_lim!(T_lim[i], U, p, t_exp)
        isnothing(T_exp!) || T_exp!(T_exp[i], U, p, t_exp)
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

    return nothing
end
