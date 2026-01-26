module VectorInterface

using LinearAlgebra: LinearAlgebra, norm, BlasFloat

# General interface: export
export scalartype
export zerovector, zerovector!, zerovector!!
export scale, scale!, scale!!
export add, add!, add!!
export inner, norm
export One, Zero

include("interface.jl")

# Auxiliary type for representing constant 1 or 0 in vector addition / linear combinations
include("onezero.jl")

# Auxiliary methods for determining types
promote_scale(x, α::Number) = promote_scale(scalartype(x), typeof(α))
promote_scale(::Type{Tx}, ::Type{Tα}) where {Tx,Tα<:Number} = Base.promote_op(scale, Tx, Tα)
promote_scale(x, y, α::Number) = promote_scale(scalartype(x), scalartype(y), typeof(α))
function promote_scale(::Type{Tx}, ::Type{Ty}, ::Type{Tα}) where {Tx,Ty,Tα<:Number}
    return Base.promote_op(scale, Tx, Ty, Tα)
end

function promote_add(x, y, α::Number=One(), β::Number=One())
    return promote_add(scalartype(x), scalartype(y), typeof(α), typeof(β))
end
function promote_add(::Type{Ty}, ::Type{Tx},
                     ::Type{Tα}=One, ::Type{Tβ}=One) where {Ty,Tx,Tα<:Number,Tβ<:Number}
    return Base.promote_op(add, Tx, Ty, Tα, Tβ)
end

promote_inner(x, y) = promote_inner(scalartype(x), scalartype(y))
promote_inner(::Type{Tx}, ::Type{Ty}) where {Tx,Ty} = Base.promote_op(inner, Tx, Ty)

# Specific implementations for Base types / type hierarchies
include("number.jl")
include("abstractarray.jl")
include("tuple.jl")
include("namedtuple.jl")

# General fallback implementation: comes with warning and some overhead
include("fallbacks.jl")

end
