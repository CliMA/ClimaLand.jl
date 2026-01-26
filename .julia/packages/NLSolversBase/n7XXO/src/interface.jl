"""
Force (re-)evaluation of the objective value at `x`.

Returns `f(x)` and stores the value in `obj.F`
"""
function value!!(obj::AbstractObjective, x)
    obj.f_calls .+= 1
    copyto!(obj.x_f, x)
    obj.F = obj.f(x)
    value(obj)
end
"""
Evaluates the objective value at `x`.

Returns `f(x)`, but does *not* store the value in `obj.F`
"""
function value(obj::AbstractObjective, x)
    obj.f_calls .+= 1
    return obj.f(x)
end
"""
Evaluates the objective value at `x`.

Returns `f(x)` and stores the value in `obj.F`
"""
function value!(obj::AbstractObjective, x)
    if x != obj.x_f
        value!!(obj, x)
    end
    value(obj)
end

"""
Evaluates the gradient value at `x`

This does *not* update `obj.DF` or `obj.x_df`.
"""
function gradient(obj::AbstractObjective, x)
    newdf = copy(obj.DF)
    obj.df(newdf, x)
    obj.df_calls .+= 1
    return newdf
end
"""
Evaluates the gradient value at `x`.

Stores the value in `obj.DF`.
"""
function gradient!(obj::AbstractObjective, x)
    if x != obj.x_df
        gradient!!(obj, x)
    end
    gradient(obj)
end
"""
Force (re-)evaluation of the gradient value at `x`.

Stores the value in `obj.DF`.
"""
function gradient!!(obj::AbstractObjective, x)
    obj.df_calls .+= 1
    copyto!(obj.x_df, x)
    obj.df(obj.DF, x)
    gradient(obj)
end

function value_gradient!(obj::AbstractObjective, x)
    if x != obj.x_f && x != obj.x_df
        value_gradient!!(obj, x)
    elseif x != obj.x_f
        value!!(obj, x)
    elseif x != obj.x_df
        gradient!!(obj, x)
    end
    value(obj), gradient(obj)
end
function value_gradient!!(obj::AbstractObjective, x)
    obj.f_calls .+= 1
    obj.df_calls .+= 1
    copyto!(obj.x_f, x)
    copyto!(obj.x_df, x)
    obj.F = obj.fdf(gradient(obj), x)
    value(obj), gradient(obj)
end

function hessian!(obj::AbstractObjective, x)
    if x != obj.x_h
        hessian!!(obj, x)
    end
    hessian(obj)
end
function hessian!!(obj::AbstractObjective, x)
    obj.h_calls .+= 1
    copyto!(obj.x_h, x)
    obj.h(obj.H, x)
    hessian(obj)
end

# Getters are without ! and accept only an objective and index or just an objective
"Get the most recently evaluated objective value of `obj`."
value(obj::AbstractObjective) = obj.F
"Get the most recently evaluated gradient of `obj`."
gradient(obj::AbstractObjective) = obj.DF
"Get the most recently evaluated Jacobian of `obj`."
jacobian(obj::AbstractObjective) = obj.DF
"Get the `i`th element of the most recently evaluated gradient of `obj`."
gradient(obj::AbstractObjective, i::Integer) = obj.DF[i]
"Get the most recently evaluated Hessian of `obj`"
hessian(obj::AbstractObjective) = obj.H

value_jacobian!(obj, x) = value_jacobian!(obj, obj.F, obj.DF, x)
function value_jacobian!(obj, F, J, x)
    if x != obj.x_f && x != obj.x_df
        value_jacobian!!(obj, F, J, x)
    elseif x != obj.x_f
        value!!(obj, F, x)
    elseif x != obj.x_df
        jacobian!!(obj, J, x)
    end
    F, J
end
value_jacobian!!(obj, x) = value_jacobian!!(obj, obj.F, obj.DF, x)
function value_jacobian!!(obj, F, J, x)
    obj.fdf(F, J, x)
    copyto!(obj.x_f, x)
    copyto!(obj.x_df, x)
    obj.f_calls .+= 1
    obj.df_calls .+= 1
    obj.df_calls
    F, J
end

function jacobian!(obj, x)
    if x != obj.x_df
        jacobian!!(obj, x)
    end
    jacobian(obj)
end

jacobian!!(obj, x) = jacobian!!(obj, obj.DF, x)
function jacobian!!(obj, J, x)
    obj.df(J, x)
    copyto!(obj.x_df, x)
    obj.df_calls .+= 1
    obj.df_calls
    J
end
function jacobian(obj::AbstractObjective, x)
    tmp = copy(obj.DF)
    jacobian!!(obj, x)
    newdf = copy(obj.DF)
    copyto!(obj.DF, tmp)
    return newdf
end

value(obj::NonDifferentiable{TF, TX}, x) where {TF<:AbstractArray, TX} = value(obj, copy(obj.F), x)
value(obj::OnceDifferentiable{TF, TDF, TX}, x) where {TF<:AbstractArray, TDF, TX} = value(obj, copy(obj.F), x)
function value(obj::AbstractObjective, F, x)
    obj.f_calls .+= 1
    return obj.f(F, x)
end

value!!(obj::NonDifferentiable{TF, TX}, x) where {TF<:AbstractArray, TX} = value!!(obj, obj.F, x)
value!!(obj::OnceDifferentiable{TF, TDF, TX}, x) where {TF<:AbstractArray, TDF, TX} = value!!(obj, obj.F, x)
function value!!(obj::AbstractObjective, F, x)
    obj.f(F, x)
    copyto!(obj.x_f, x)
    obj.f_calls .+= 1
    obj.f_calls
    F
end

function _clear_f!(d::NLSolversBase.AbstractObjective)
    d.f_calls .= 0
    if typeof(d.F) <: AbstractArray
        d.F .= eltype(d.F)(NaN)
    else
        d.F = typeof(d.F)(NaN)
    end
    d.x_f .= eltype(d.x_f)(NaN)
    nothing
end

function _clear_df!(d::NLSolversBase.AbstractObjective)
    d.df_calls .= 0
    d.DF .= eltype(d.DF)(NaN)
    d.x_df .= eltype(d.x_df)(NaN)
    nothing
end

function _clear_h!(d::NLSolversBase.AbstractObjective)
    d.h_calls .= 0
    d.H .= eltype(d.H)(NaN)
    d.x_h .= eltype(d.x_h)(NaN)
    nothing
end

function _clear_hv!(d::NLSolversBase.AbstractObjective)
    d.hv_calls .= 0
    d.Hv .= eltype(d.Hv)(NaN)
    d.x_hv .= eltype(d.x_hv)(NaN)
    d.v_hv .= eltype(d.v_hv)(NaN)
    nothing
end

clear!(d::NonDifferentiable)  = _clear_f!(d)

function clear!(d::OnceDifferentiable)
    _clear_f!(d)
    _clear_df!(d)
    nothing
end

function clear!(d::TwiceDifferentiable)
    _clear_f!(d)
    _clear_df!(d)
    _clear_h!(d)
    nothing
end

function clear!(d::TwiceDifferentiableHV)
    _clear_f!(d)
    _clear_df!(d)
    _clear_hv!(d)
    nothing
end

g_calls(d::NonDifferentiable) = 0
h_calls(d::Union{NonDifferentiable, OnceDifferentiable}) = 0
f_calls(d) = first(d.f_calls)
g_calls(d) = first(d.df_calls)
h_calls(d) = first(d.h_calls)
hv_calls(d) = 0
h_calls(d::TwiceDifferentiableHV) = 0
hv_calls(d::TwiceDifferentiableHV) = first(d.hv_calls)
