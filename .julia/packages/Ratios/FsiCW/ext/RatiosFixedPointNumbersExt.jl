module RatiosFixedPointNumbersExt

using Ratios
isdefined(Base, :get_extension) ? (using FixedPointNumbers) : (using ..FixedPointNumbers)

using .FixedPointNumbers: FixedPoint, Fixed, Normed, rawone

rawone_noerr(::Type{Fixed{T,f}}) where {T,f} = widen(oneunit(T)) << f
rawone_noerr(::Type{N}) where N<:Normed = rawone(N)
rawone_noerr(x::FixedPoint) = rawone_noerr(typeof(x))
Base.promote_rule(::Type{SimpleRatio{S}}, ::Type{<:FixedPoint{T}}) where {S<:Integer,T<:Integer} = SimpleRatio{promote_type(S, T)}
Ratios.SimpleRatio{S}(x::FixedPoint) where S<:Integer = SimpleRatio{S}(reinterpret(x), rawone_noerr(x))
Ratios.SimpleRatio(x::FixedPoint) = SimpleRatio(reinterpret(x), rawone_noerr(x))
Base.convert(::Type{S}, x::FixedPoint) where S<:SimpleRatio = S(x)

end
