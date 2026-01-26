# Vector interface implementation: general fallbacks?
#######################################################
# Are these desirable; should they exist?
# They might lead to plenty of invalidations?
# They do not enforce any compatibility between the different vector types

@noinline function _warn_message(fun, type)
    return string("The function `", fun,
                  "` is not implemented for (values of) type `", type, "`;\n",
                  "this fallback will disappear in future versions of VectorInterface.jl")
end
@noinline function _error_message(fun, type)
    return string("No fallback for applying `", fun, "` to (values of) type `", type,
                  "` could be determined")
end

# scalartype
#------------
function scalartype(T::Type)
    @warn _warn_message(scalartype, T) maxlog = 1
    # the following does not work because `eltype` is always applicable.
    # applicable(eltype, T) && return scalartype(eltype(T))
    elT = eltype(T)
    elT === T && throw(ArgumentError(_error_message(scalartype, T)))
    return scalartype(elT)
end
# should this try to use `eltype` instead? e.g. scalartype(T) = scalartype(eltype(T))

# zerovector & zerovector!!
#---------------------------
function zerovector(x, ::Type{S}) where {S<:Number}
    T = Tuple{typeof(x),S}
    @warn _warn_message(zerovector, T) maxlog = 1
    if applicable(similar, x, S)
        return zerovector!(similar(x, S))
    else
        throw(ArgumentError(_error_message(zerovector, T)))
    end
end
function zerovector!(x)
    T = typeof(x)
    @warn _warn_message(zerovector!, T) maxlog = 1
    if applicable(LinearAlgebra.rmul!, x, false)
        return LinearAlgebra.rmul!(x, false)
    else
        throw(ArgumentError(_error_message(zerovector!, T)))
    end
end
function zerovector!!(x)
    T = typeof(x)
    @warn _warn_message(zerovector!!, T) maxlog = 1
    if applicable(LinearAlgebra.rmul!, x, false)
        return LinearAlgebra.rmul!(x, false)
    elseif applicable(zero, x)
        return zero(x)
    elseif applicable(similar, x) && applicable(LinearAlgebra.rmul!, x, false)
        return LinearAlgebra.rmul!(similar(x), false)
    else
        throw(ArgumentError(_error_message(zerovector!!, T)))
    end
end

# scale, scale! & scale!!
#-------------------------
function scale(x, α::Number)
    (α === One()) && return x
    T = Tuple{typeof(x),typeof(α)}
    @warn _warn_message(scale, T) maxlog = 1
    if applicable(*, x, α)
        return x * α
    elseif applicable(*, α, x)
        return α * x
    else
        throw(ArgumentError(_error_message(scale, T)))
    end
end

function scale!(x, α::Number)
    (α === One()) && return x
    T = Tuple{typeof(x),typeof(α)}
    @warn _warn_message(scale!, T) maxlog = 1
    if applicable(LinearAlgebra.rmul!, x, α)
        return LinearAlgebra.rmul!(x, α)
    else
        throw(ArgumentError(_error_message(scale!, T)))
    end
end

function scale!!(x, α::Number)
    (α === One()) && return x
    T = Tuple{typeof(x),typeof(α)}
    @warn _warn_message(scale!!, T) maxlog = 1
    if applicable(LinearAlgebra.rmul!, x, α) && promote_scale(x, α) <: scalartype(x)
        return LinearAlgebra.rmul!(x, α)
    elseif applicable(*, x, α)
        return x * α
    else
        throw(ArgumentError(_error_message(scale!!, T)))
    end
end

function scale!(y, x, α::Number)
    T = Tuple{typeof(y),typeof(x),typeof(α)}
    @warn _warn_message(scale!, T) maxlog = 1
    if applicable(LinearAlgebra.mul!, y, x, α)
        return LinearAlgebra.mul!(y, x, α)
    else
        throw(ArgumentError(_error_message(scale!, T)))
    end
end

function scale!!(y, x, α::Number)
    T = Tuple{typeof(y),typeof(x),typeof(α)}
    @warn _warn_message(scale!!, T) maxlog = 1
    if applicable(LinearAlgebra.mul!, y, x, α) && promote_scale(y, x, α) <: scalartype(y)
        return LinearAlgebra.mul!(y, x, α)
    else
        α_Ty = α * one(scalartype(y))
        if applicable(*, x, α_Ty)
            return x * α_Ty
        else
            throw(ArgumentError(_error_message(scale!!, T)))
        end
    end
end

# add, add! & add!!
#-------------------
function add(y, x, α::Number, β::Number)
    T = Tuple{typeof(y),typeof(x),typeof(α),typeof(β)}
    @warn _warn_message(add, T) maxlog = 1
    yb = scale(y, β)
    xa = scale(x, α)
    if applicable(+, yb, xa)
        return yb + xa
    else
        throw(ArgumentError(_error_message(add, T)))
    end
end

function add!(y, x, α::Number, β::Number)
    T = Tuple{typeof(y),typeof(x),typeof(α),typeof(β)}
    @warn _warn_message(add!, T) maxlog = 1

    α′ = (α === One()) ? true : α
    if β === One()
        if applicable(LinearAlgebra.axpy!, α′, x, y)
            return LinearAlgebra.axpy!(α′, x, y)
        else
            throw(ArgumentError(_error_message(add!, T)))
        end
    else
        if applicable(LinearAlgebra.axpby!, α′, x, β, y)
            return LinearAlgebra.axpby!(α′, x, β, y)
        else
            throw(ArgumentError(_error_message(add!, T)))
        end
    end
end

function add!!(y, x, α::Number, β::Number)
    T = Tuple{typeof(y),typeof(x),typeof(α),typeof(β)}
    @warn _warn_message(add!!, T) maxlog = 1

    if β === One() && α === One()
        if applicable(LinearAlgebra.axpy!, true, x, y) && promote_add(y, x) <: scalartype(y)
            return LinearAlgebra.axpy!(true, x, y)
        elseif applicable(+, y, x)
            return y + x
        else
            throw(ArgumentError(_error_message(add!!, T)))
        end
    elseif β === One()
        α′ = (α === One()) ? true : α
        if applicable(LinearAlgebra.axpy!, α′, x, y) &&
           promote_add(y, x, α) <: scalartype(y)
            return LinearAlgebra.axpy!(α′, x, y)
        else
            return add!!(y, scale(x, α))
        end
    else
        α′ = (α === One()) ? true : α
        if applicable(LinearAlgebra.axpby!, α′, x, β, y) &&
           promote_add(y, x, α, β) <: scalartype(y)
            return LinearAlgebra.axpby!(α′, x, β, y)
        else
            return add!!(scale!!(y, β), scale(x, α))
        end
    end
end

# inner
#-------
function inner(x, y)
    T = Tuple{typeof(x),typeof(y)}
    @warn _warn_message(inner, T) maxlog = 1

    if applicable(LinearAlgebra.dot, x, y)
        return LinearAlgebra.dot(x, y)
    else
        throw(ArgumentError(_error_message(inner, T)))
    end
end
