# Vector interface implementation for subtypes of `Tuple`
###########################################################
# (as long as they have a homogeneous scalar type)
# general comment: `map` on tuples is better at preserving type information than broadcast

# scalartype
#------------
scalartype(::Type{NamedTuple{names,T}}) where {names,T<:Tuple} = scalartype(T)

# zerovector & zerovector!!
#---------------------------
function zerovector(x::NamedTuple, ::Type{S}) where {S<:Number}
    return NamedTuple{keys(x)}(zerovector(values(x), S))
end
function zerovector!!(x::NamedTuple)
    return NamedTuple{keys(x)}(zerovector!!(values(x)))
end

# scale & scale!!
#-----------------
scale(x::NamedTuple, α::Number) = NamedTuple{keys(x)}(map(xᵢ -> scale(xᵢ, α), values(x)))
function scale!!(x::NamedTuple, α::Number)
    return NamedTuple{keys(x)}(map(xᵢ -> scale!!(xᵢ, α), values(x)))
end
function scale!!(y::NamedTuple{names}, x::NamedTuple{names}, α::Number) where {names}
    xvals = values(x)
    yvals = values(y)
    yxvals = ntuple(i -> (yvals[i], xvals[i]), length(x))
    return NamedTuple{names}(map(yxᵢ -> scale!!(yxᵢ[1], yxᵢ[2], α), yxvals))
end

# add & add!!
#-------------
function add(y::NamedTuple{names}, x::NamedTuple{names},
             α::Number, β::Number) where {names}
    xvals = values(x)
    yvals = values(y)
    yxvals = ntuple(i -> (yvals[i], xvals[i]), length(x))
    return NamedTuple{names}(map(yxᵢ -> add(yxᵢ[1], yxᵢ[2], α, β), yxvals))
end
function add!!(y::NamedTuple{names}, x::NamedTuple{names},
               α::Number, β::Number) where {names}
    xvals = values(x)
    yvals = values(y)
    yxvals = ntuple(i -> (yvals[i], xvals[i]), length(x))
    return NamedTuple{names}(map(yxᵢ -> add!!(yxᵢ[1], yxᵢ[2], α, β), yxvals))
end

# inner
#-------
function inner(x::NamedTuple{names}, y::NamedTuple{names}) where {names}
    xvals = values(x)
    yvals = values(y)
    xyvals = ntuple(i -> (xvals[i], yvals[i]), length(x))
    return sum(map(xyᵢ -> inner(xyᵢ[1], xyᵢ[2]), xyvals))
end
