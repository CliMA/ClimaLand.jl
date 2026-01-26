# Vector interface implementation for subtypes of `Number`
###########################################################

# scalartype
#------------
scalartype(::Type{T}) where {T<:Number} = T

# zerovector & zerovector!!
#---------------------------
@inline zerovector(::Number, ::Type{S}) where {S<:Number} = zero(S)
@inline zerovector!!(x::Number) = zero(x)

# scale & scale!!
#-----------------
# note: required to make scale(NaN, 0) = 0
@inline scale(x::Number, α::Number) = (iszero(α) ? zero(x) : x) * α
@inline function scale(y::Number, x::Number, α::Number)
    T = Base.promote_op(*, typeof(y), typeof(x), typeof(α))
    return convert(T, scale(x, α))
end
@inline scale!!(x::Number, α::Number) = scale(x, α)
@inline scale!!(y::Number, x::Number, α::Number) = scale(y, x, α)

# add & add!!
#-------------
@inline add(y::Number, x::Number, ::One, ::One) = y + x
# @inline add(y::Number, x::Number, ::One, β::Number) = add(x, y, β)
@inline function add(y::Number, x::Number, α::Number, ::One)
    return iszero(α) ? muladd(zero(x), α, y) : muladd(x, α, y) # required to make zero coefficients kill NaN values in x or y
end
@inline add(y::Number, x::Number, α::Number, β::Number) = add(scale(y, β), x, α)

@inline add!!(y::Number, x::Number, α::Number, β::Number) = add(y, x, α, β)

# inner
#-------
@inline inner(x::Number, y::Number) = conj(x) * y
