# This code has gotten a bit stale, but is still useful for
# understanding the order of operations.
function step_u!(integrator, cache::IMEXARKCache, ::ARS343)
    (; u, p, t, dt, sol, alg) = integrator
    (; f) = sol.prob
    (; T_imp!, lim!, dss!) = f
    (; cache!, cache_imp!) = f
    (; tableau, newtons_method) = alg
    (; a_exp, b_exp, a_imp, b_imp, c_exp, c_imp) = tableau
    (; U, T_lim, T_exp, T_imp, temp, γ, newtons_method_cache) = cache
    T_lim! = !isnothing(f.T_lim!) ? f.T_lim! : (args...) -> nothing
    T_exp! = !isnothing(f.T_exp!) ? f.T_exp! : (args...) -> nothing
    dtγ = dt * γ

    if !isnothing(newtons_method_cache)
        jacobian = newtons_method_cache.j
        if (!isnothing(jacobian)) && needs_update!(newtons_method.update_j, NewTimeStep(t))
            T_imp!.Wfact(jacobian, u, p, dtγ, t)
        end
    end

    s = 4

    i::Int = 1
    t_exp = t
    @. U = u # TODO: This is unnecessary; we can just pass u to T_exp and T_lim
    T_lim!(T_lim[i], U, p, t_exp)
    T_exp!(T_exp[i], U, p, t_exp)

    i = 2
    t_exp = t + dt * c_exp[i]
    @. U = u + dt * a_exp[i, 1] * T_lim[1]
    lim!(U, p, t_exp, u)
    @. U += dt * a_exp[i, 1] * T_exp[1]
    dss!(U, p, t_exp)
    @. temp = U # used in closures
    let i = i
        t_imp = t + dt * c_imp[i]
        cache_imp!(U, p, t_imp)
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
    T_lim!(T_lim[i], U, p, t_exp)
    T_exp!(T_exp[i], U, p, t_exp)
    @. T_imp[i] = (U - temp) / dtγ

    i = 3
    t_exp = t + dt * c_exp[i]
    @. U = u + dt * a_exp[i, 1] * T_lim[1] + dt * a_exp[i, 2] * T_lim[2]
    lim!(U, p, t_exp, u)
    @. U += dt * a_exp[i, 1] * T_exp[1] + dt * a_exp[i, 2] * T_exp[2] + dt * a_imp[i, 2] * T_imp[2]
    dss!(U, p, t_exp)
    @. temp = U # used in closures
    let i = i
        t_imp = t + dt * c_imp[i]
        cache_imp!(U, p, t_imp)
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
    T_lim!(T_lim[i], U, p, t_exp)
    T_exp!(T_exp[i], U, p, t_exp)
    @. T_imp[i] = (U - temp) / dtγ

    i = 4
    t_exp = t + dt
    @. U = u + dt * a_exp[i, 1] * T_lim[1] + dt * a_exp[i, 2] * T_lim[2] + dt * a_exp[i, 3] * T_lim[3]
    lim!(U, p, t_exp, u)
    @. U +=
        dt * a_exp[i, 1] * T_exp[1] +
        dt * a_exp[i, 2] * T_exp[2] +
        dt * a_exp[i, 3] * T_exp[3] +
        dt * a_imp[i, 2] * T_imp[2] +
        dt * a_imp[i, 3] * T_imp[3]
    dss!(U, p, t_exp)
    @. temp = U # used in closures
    let i = i
        t_imp = t + dt * c_imp[i]
        cache_imp!(U, p, t_imp)
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
    T_lim!(T_lim[i], U, p, t_exp)
    T_exp!(T_exp[i], U, p, t_exp)
    @. T_imp[i] = (U - temp) / dtγ

    t_final = t + dt
    @. temp = u + dt * b_exp[2] * T_lim[2] + dt * b_exp[3] * T_lim[3] + dt * b_exp[4] * T_lim[4]
    lim!(temp, p, t_final, u)
    @. u =
        temp +
        dt * b_exp[2] * T_exp[2] +
        dt * b_exp[3] * T_exp[3] +
        dt * b_exp[4] * T_exp[4] +
        dt * b_imp[2] * T_imp[2] +
        dt * b_imp[3] * T_imp[3] +
        dt * b_imp[4] * T_imp[4]
    dss!(u, p, t_final)
    cache!(u, p, t_final)
    return u
end
