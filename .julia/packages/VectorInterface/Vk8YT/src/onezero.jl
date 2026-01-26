"""
    struct One end

Singleton type for representing a hard-coded constant 1 in vector addition / linear
combinations.
"""
struct One <: Number end

"""
    struct Zero end
    
Singleton type for representing a hard-coded constant 0 in vector addition / linear
combinations.

```quote
    Now I am become Zero, the destroyer of NaN. - Vishnu
```
"""
struct Zero <: Number end

# Base Arithmetic
# ---------------

Base.:(-)(::One) = -1
Base.:(-)(::Zero) = Zero()

Base.:(+)(::Zero, x::Number) = x
Base.:(+)(x::Number, ::Zero) = x
Base.:(+)(::Zero, ::Zero) = Zero()
Base.:(+)(::One, ::One) = 2

Base.:(-)(x::Number, ::Zero) = x
Base.:(-)(::Zero, x::Number) = -x
Base.:(-)(::Zero, ::Zero) = Zero()
Base.:(-)(::One, ::One) = Zero()

Base.:(*)(::One, x::Number) = x
Base.:(*)(::Zero, x::Number) = zero(x)
Base.:(*)(x::Number, ::One) = x
Base.:(*)(x::Number, ::Zero) = zero(x)
Base.:(*)(::Zero, ::Zero) = Zero()
Base.:(*)(::One, ::One) = One()
Base.:(*)(::Zero, ::One) = Zero()
Base.:(*)(::One, ::Zero) = Zero()

Base.:(/)(::Zero, ::Zero) = throw(DivideError())
Base.:(/)(::One, ::One) = One()
Base.:(/)(::Zero, ::One) = Zero()
Base.:(/)(::One, ::Zero) = throw(DivideError())
Base.:(/)(::Number, ::Zero) = throw(DivideError())
Base.:(/)(::Zero, x::Number) = iszero(x) ? throw(DivideError()) : zero(x)
Base.:(/)(x::Number, ::One) = x
Base.:(/)(::One, x::Number) = inv(x)

Base.inv(::One) = One()

Base.conj(::One) = One()
Base.conj(::Zero) = Zero()

Base.one(::Type{One}) = One()
Base.one(::Type{Zero}) = One()
Base.one(::Union{Zero,One}) = One()
Base.zero(::Type{One}) = Zero()
Base.zero(::Type{Zero}) = Zero()
Base.zero(::Union{Zero,One}) = Zero()

Base.:(==)(::One, ::One) = true
Base.:(==)(::Zero, ::Zero) = true
Base.:(==)(::One, ::Zero) = false
Base.:(==)(::Zero, ::One) = false

# Promotion
# ---------
Base.promote_rule(::Type{Zero}, ::Type{One}) = Bool
Base.promote_rule(::Type{One}, ::Type{Zero}) = Bool
Base.promote_rule(::Type{One}, ::Type{T}) where {T<:Number} = T
Base.promote_rule(::Type{Zero}, ::Type{T}) where {T<:Number} = T
# disambiguate:
Base.promote_rule(::Type{Bool}, ::Type{One}) = Bool
Base.promote_rule(::Type{Bool}, ::Type{Zero}) = Bool
Base.convert(::Type{T}, ::One) where {T<:Number} = one(T)
(T::Type{<:Number})(::One) = one(T)
Base.convert(::Type{T}, ::Zero) where {T<:Number} = zero(T)
(T::Type{<:Number})(::Zero) = zero(T)
