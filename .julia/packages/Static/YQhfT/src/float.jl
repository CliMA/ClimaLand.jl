
"""
    StaticFloat64(F::Float64)::StaticFloat64{F}

A statically sized `Float64`.
Use `StaticFloat64(N)` instead of `Val(N)` when you want it to behave like a number.
"""
struct StaticFloat64{N} <: Real
    StaticFloat64{N}() where {N} = new{N::Float64}()
    StaticFloat64(x::Float64) = new{x}()
    StaticFloat64(x::Int) = new{Base.sitofp(Float64, x)::Float64}()
    StaticFloat64(x::StaticInt{N}) where {N} = StaticFloat64(convert(Float64, N))
    StaticFloat64(x::Complex) = StaticFloat64(convert(Float64, x))
    StaticFloat64(@nospecialize x::StaticFloat64) = x
end

Base.zero(@nospecialize T::Type{<:StaticFloat64}) = Float64(0.0)
Base.one(@nospecialize T::Type{<:StaticFloat64}) = Float64(1.0)

Base.round(::StaticFloat64{M}) where {M} = StaticFloat64(round(M))
roundtostaticint(::StaticFloat64{M}) where {M} = StaticInt(round(Int, M))
roundtostaticint(x::AbstractFloat) = round(Int, x)
floortostaticint(::StaticFloat64{M}) where {M} = StaticInt(Base.fptosi(Int, M))
floortostaticint(x::AbstractFloat) = Base.fptosi(Int, x)

Base.rad2deg(::StaticFloat64{M}) where {M} = StaticFloat64(rad2deg(M))
Base.deg2rad(::StaticFloat64{M}) where {M} = StaticFloat64(deg2rad(M))
@generated Base.cbrt(::StaticFloat64{M}) where {M} = StaticFloat64(cbrt(M))
Base.mod2pi(::StaticFloat64{M}) where {M} = StaticFloat64(mod2pi(M))
@generated Base.sinpi(::StaticFloat64{M}) where {M} = StaticFloat64(sinpi(M))
@generated Base.cospi(::StaticFloat64{M}) where {M} = StaticFloat64(cospi(M))
Base.exp(::StaticFloat64{M}) where {M} = StaticFloat64(exp(M))
Base.exp2(::StaticFloat64{M}) where {M} = StaticFloat64(exp2(M))
Base.exp10(::StaticFloat64{M}) where {M} = StaticFloat64(exp10(M))
@generated Base.expm1(::StaticFloat64{M}) where {M} = StaticFloat64(expm1(M))
@generated Base.log(::StaticFloat64{M}) where {M} = StaticFloat64(log(M))
@generated Base.log2(::StaticFloat64{M}) where {M} = StaticFloat64(log2(M))
@generated Base.log10(::StaticFloat64{M}) where {M} = StaticFloat64(log10(M))
@generated Base.log1p(::StaticFloat64{M}) where {M} = StaticFloat64(log1p(M))
@generated Base.sin(::StaticFloat64{M}) where {M} = StaticFloat64(sin(M))
@generated Base.cos(::StaticFloat64{M}) where {M} = StaticFloat64(cos(M))
@generated Base.tan(::StaticFloat64{M}) where {M} = StaticFloat64(tan(M))
Base.sec(x::StaticFloat64{M}) where {M} = inv(cos(x))
Base.csc(x::StaticFloat64{M}) where {M} = inv(sin(x))
Base.cot(x::StaticFloat64{M}) where {M} = inv(tan(x))
@generated Base.asin(::StaticFloat64{M}) where {M} = StaticFloat64(asin(M))
@generated Base.acos(::StaticFloat64{M}) where {M} = StaticFloat64(acos(M))
@generated Base.atan(::StaticFloat64{M}) where {M} = StaticFloat64(atan(M))
@generated Base.sind(::StaticFloat64{M}) where {M} = StaticFloat64(sind(M))
@generated Base.cosd(::StaticFloat64{M}) where {M} = StaticFloat64(cosd(M))
Base.tand(x::StaticFloat64{M}) where {M} = sind(x) / cosd(x)
Base.secd(x::StaticFloat64{M}) where {M} = inv(cosd(x))
Base.cscd(x::StaticFloat64{M}) where {M} = inv(sind(x))
Base.cotd(x::StaticFloat64{M}) where {M} = inv(tand(x))
Base.asind(x::StaticFloat64{M}) where {M} = rad2deg(asin(x))
Base.acosd(x::StaticFloat64{M}) where {M} = rad2deg(acos(x))
Base.asecd(x::StaticFloat64{M}) where {M} = rad2deg(asec(x))
Base.acscd(x::StaticFloat64{M}) where {M} = rad2deg(acsc(x))
Base.acotd(x::StaticFloat64{M}) where {M} = rad2deg(acot(x))
Base.atand(x::StaticFloat64{M}) where {M} = rad2deg(atan(x))
@generated Base.sinh(::StaticFloat64{M}) where {M} = StaticFloat64(sinh(M))
Base.cosh(::StaticFloat64{M}) where {M} = StaticFloat64(cosh(M))
Base.tanh(x::StaticFloat64{M}) where {M} = StaticFloat64(tanh(M))
Base.sech(x::StaticFloat64{M}) where {M} = inv(cosh(x))
Base.csch(x::StaticFloat64{M}) where {M} = inv(sinh(x))
Base.coth(x::StaticFloat64{M}) where {M} = inv(tanh(x))
@generated Base.asinh(::StaticFloat64{M}) where {M} = StaticFloat64(asinh(M))
@generated Base.acosh(::StaticFloat64{M}) where {M} = StaticFloat64(acosh(M))
@generated Base.atanh(::StaticFloat64{M}) where {M} = StaticFloat64(atanh(M))
Base.asech(x::StaticFloat64{M}) where {M} = acosh(inv(x))
Base.acsch(x::StaticFloat64{M}) where {M} = asinh(inv(x))
Base.acoth(x::StaticFloat64{M}) where {M} = atanh(inv(x))
Base.asec(x::StaticFloat64{M}) where {M} = acos(inv(x))
Base.acsc(x::StaticFloat64{M}) where {M} = asin(inv(x))
Base.acot(x::StaticFloat64{M}) where {M} = atan(inv(x))

Base.rem(x::Real, ::StaticFloat64{Y}) where {Y} = rem(x, Y)
Base.rem(::StaticFloat64{X}, y::Real) where {X} = rem(X, y)
Base.rem(::StaticFloat64{X}, ::StaticFloat64{Y}) where {X, Y} = StaticFloat64(rem(X, Y))

Base.min(x::StaticFloat64{X}, y::StaticFloat64{Y}) where {X, Y} = X > Y ? y : x
Base.min(x::Real, ::StaticFloat64{Y}) where {Y} = min(x, Y)
Base.min(::StaticFloat64{X}, y::Real) where {X} = min(X, y)

Base.max(x::StaticFloat64{X}, y::StaticFloat64{Y}) where {X, Y} = X > Y ? x : y
Base.max(x::Real, ::StaticFloat64{Y}) where {Y} = max(x, Y)
Base.max(::StaticFloat64{X}, y::Real) where {X} = max(X, y)

Base.isless(::StaticFloat64{X}, ::StaticFloat64{Y}) where {X, Y} = isless(X, Y)
Base.isless(x::AbstractFloat, ::StaticFloat64{Y}) where {Y} = isless(x, Y)
Base.isless(::StaticFloat64{X}, y::AbstractFloat) where {X} = isless(X, y)

Base.inv(::StaticFloat64{N}) where {N} = StaticFloat64(1.0 / N)

Base.:(*)(::StaticFloat64{X}, ::StaticFloat64{Y}) where {X, Y} = static(X * Y)
Base.:(/)(::StaticFloat64{X}, ::StaticFloat64{Y}) where {X, Y} = static(X / Y)
Base.:(-)(::StaticFloat64{X}, ::StaticFloat64{Y}) where {X, Y} = static(X - Y)
Base.:(+)(::StaticFloat64{X}, ::StaticFloat64{Y}) where {X, Y} = static(X + Y)

Base.isinteger(x::StaticFloat64{X}) where {X} = isinteger(X)
