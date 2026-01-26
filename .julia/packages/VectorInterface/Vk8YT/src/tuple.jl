# Vector interface implementation for subtypes of `Tuple`
###########################################################
# (as long as they have a homogeneous scalar type)
# general comment: `map` on tuples is better at preserving type information than broadcast

# scalartype
#------------
function scalartype(::Type{Tuple{}})
    throw(ArgumentError("no scalar type is defined for empty tuple"))
end
scalartype(::Type{Tuple{T}}) where {T} = scalartype(T)
scalartype(::Type{<:Tuple{T,Vararg{T}}}) where {T} = scalartype(T)
function scalartype(::Type{TT}) where {TT<:Tuple}
    S = scalartype(Base.tuple_type_head(TT))
    S2 = scalartype(Base.tuple_type_tail(TT))
    S == S2 || throw(ArgumentError("no (unique) scalar type is defined for type $TT"))
    return S
end

# zerovector & zerovector!!
#---------------------------
function zerovector(x::Tuple, ::Type{S}) where {S<:Number}
    y = map(xᵢ -> zerovector(xᵢ, S), x)
    return y
end
function zerovector!!(x::Tuple)
    y = map(xᵢ -> zerovector!!(xᵢ), x)
    return y
end

# scale & scale!!
#-----------------
scale(x::Tuple, α::Number) = map(xᵢ -> scale(xᵢ, α), x)
scale!!(x::Tuple, α::Number) = map(xᵢ -> scale!!(xᵢ, α), x)
function scale!!(y::Tuple, x::Tuple, α::Number)
    lx = length(x)
    ly = length(y)
    lx == ly || throw(DimensionMismatch("non-matching tuple lengths $lx and $ly"))
    yx = ntuple(i -> (y[i], x[i]), lx)
    return map(yxᵢ -> scale!!(yxᵢ[1], yxᵢ[2], α), yx)
end

# add & add!!
#-------------
function add(y::Tuple, x::Tuple, α::Number=One(), β::Number=One())
    lx = length(x)
    ly = length(y)
    lx == ly || throw(DimensionMismatch("non-matching tuple lengths $lx and $ly"))
    yx = ntuple(i -> (y[i], x[i]), lx)
    return map(yxᵢ -> add(yxᵢ[1], yxᵢ[2], α, β), yx)
end
function add!!(y::Tuple, x::Tuple, α::Number=One(), β::Number=One())
    lx = length(x)
    ly = length(y)
    lx == ly || throw(DimensionMismatch("non-matching tuple lengths $lx and $ly"))
    yx = ntuple(i -> (y[i], x[i]), lx)
    return map(yxᵢ -> add!!(yxᵢ[1], yxᵢ[2], α, β), yx)
end

# inner
#-------
function inner(x::Tuple, y::Tuple)
    lx = length(x)
    ly = length(y)
    lx == ly || throw(DimensionMismatch("non-matching tuple lengths $lx and $ly"))
    xy = ntuple(i -> (x[i], y[i]), lx)
    return sum(map(xyᵢ -> inner(xyᵢ[1], xyᵢ[2]), xy))
end
