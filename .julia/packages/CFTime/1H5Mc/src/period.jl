# Period constructors or methods with Period as first (or main) argument

# sadly Dates.CompoundPeriod allocates a vector otherwise
# Dates.CompoundPeriod could have been used instead of CFTime.Period

Period(duration::Number, factor, exponent = -3) = Period{typeof(duration), Val(factor), Val(exponent)}(duration)

@inline function Period(duration::T, units::Val) where {T <: Number}
    Tfactor = _Tfactor(units)
    Texponent = _Texponent(units)
    return Period{T, Tfactor, Texponent}(duration)
end

@inline function Period(duration::Number, units::Union{Symbol, AbstractString})
    return Period(duration, Val(Symbol(units)))
end

_type(p::Period{T, Tfactor, Texponent}) where {T, Tfactor, Texponent} = T
_factor(p::Period{T, Tfactor, Texponent}) where {T, Tfactor, Texponent} = unwrap(Tfactor)
_exponent(p::Period{T, Tfactor, Texponent}) where {T, Tfactor, Texponent} = unwrap(Texponent)


_type(::Type{Period{T, Tfactor, Texponent}}) where {T, Tfactor, Texponent} = T
_factor(::Type{Period{T, Tfactor, Texponent}}) where {T, Tfactor, Texponent} = unwrap(Tfactor)
_exponent(::Type{Period{T, Tfactor, Texponent}}) where {T, Tfactor, Texponent} = unwrap(Texponent)

Dates.value(p::Period) = p.duration

function Base.zero(p::Period{T, Tfactor, Texponent}) where {T, Tfactor, Texponent}
    return Period{T, Tfactor, Texponent}(0)
end

function Base.one(p::Period{T, Tfactor, Texponent}) where {T, Tfactor, Texponent}
    return Period{T, Tfactor, Texponent}(1)
end

function Base.abs(p::Period{T, Tfactor, Texponent}) where {T, Tfactor, Texponent}
    return Period{T, Tfactor, Texponent}(abs(value(p)))
end

Base.sign(x::Period) = sign(value(x))
Base.signbit(x::Period) = signbit(value(x))

# helper functions for _timetuple
@inline __tf(result, time) = result
@inline function __tf(result, time::T, d1, dn...) where {T}
    Ti = promote_type(Int64, T)

    return if d1 == 0
        __tf((result..., Ti(0)), 0, dn...)
    else
        p, time2 = divrem(time, d1, RoundDown)
        __tf((result..., Ti(p)), time2, dn...)
    end
end

"""
    _timetuple(time,divi)

Recursively divides `time` into the tuple `divi`. For example

```julia
divi = (24*60*60,60*60,60,1)
_timetuple(1234567,divi)
# output
# (14, 6, 56, 7)
# 14 days, 6 hours, 56 minutes, 7 seconds
sum(_timetuple(1234567,divi) .* divi) == 1234567
# output
# true
```
"""
@inline function _timetuple(time, divi)
    return __tf((), time, divi...)
end


# T(10)^(-n) for n <= 0 avoiding problems with type-inference on the GPU (CUDA)
# but slow on CPU
#@inline ten_tom(T,n::Integer) = (n >= 0 ? T(1) : T(10) * ten_tom(T,n+1))

# rescale the time units for the ratio factor/exponent
@inline function division(::Type{T}, factor, exponent) where {T}
    return ntuple(length(TIME_DIVISION)) do i
        (T(10)^(-exponent) * TIME_DIVISION[i][2]) ÷
            (T(10)^(-TIME_DIVISION[i][3]) * factor)
        #        (ten_tom(T,exponent) * TIME_DIVISION[i][2]) ÷
        #            (ten_tom(T,TIME_DIVISION[i][3]) * factor)
    end
end

@inline function division(::Type{T}, factor, exponent) where {T <: AbstractFloat}
    return ntuple(length(TIME_DIVISION)) do i
        (T(10)^(-exponent) * TIME_DIVISION[i][2]) /
            (T(10)^(-TIME_DIVISION[i][3]) * factor)
    end
end

@inline function _datenum(tuf::Tuple, factor, exponent)
    T = promote_type(typeof.(tuf)...)
    divi = division(T, factor, exponent)

    for (_tuf, _divi) in zip(tuf, divi)
        if (_tuf != 0) && (_divi == 0)
            error("Cannot express the time tuple $tuf (representing days,[hours,minutes,seconds,milliseconds...]) in $factor × 10^($exponent) seconds without loss of precision")
        end
    end
    return sum(divi[1:length(tuf)] .* tuf)
end

"""
    days,h,mi,s,ms,... = timetuplefrac(t::Period)

Return a tuple with the number of whole days, hours (`h`), minutes (`mi`),
seconds (`s`) and millisecods (`ms`),... from the time period `t`.
"""
function timetuplefrac(t::Period{T, Tfactor}) where {T, Tfactor}
    # for integers
    factor = _factor(t)
    exponent = _exponent(t)
    divi = division(T, factor, exponent)
    time = t.duration
    return _timetuple(time, divi)
end


@inline function Period(tuf::Tuple, factor, exponent = -3)
    duration = _datenum(tuf, factor, exponent)
    return Period{typeof(duration), Val(factor), Val(exponent)}(duration)
end

# must be inlined for type-stability using the CUDA compiler
@inline function Period(T::DataType, tuf::Tuple, factor, exponent = -3)
    duration = T(_datenum(tuf, factor, exponent))
    return Period{typeof(duration), Val(factor), Val(exponent)}(duration)
end

function promote_rule(
        ::Type{Period{T1, Tfactor1, Texponent1}},
        ::Type{Period{T2, Tfactor2, Texponent2}}
    ) where
    {T1, Tfactor1, Texponent1, T2, Tfactor2, Texponent2}

    factor1 = unwrap(Tfactor1)
    factor2 = unwrap(Tfactor2)
    exponent1 = unwrap(Texponent1)
    exponent2 = unwrap(Texponent2)
    T = promote_type(T1, T2)

    # the 10^exp is problematic for the CUDA compiler and type-stability
    # all these clauses expect the last one redundant but necessary to
    # workaround limitation of the CUDA compiler

    if factor1 == factor2
        if exponent1 < exponent2
            return Period{T, Tfactor1, Texponent1}
        else
            return Period{T, Tfactor2, Texponent2}
        end
    elseif exponent1 == exponent2
        if factor1 < factor2
            return Period{T, Tfactor1, Texponent1}
        else
            return Period{T, Tfactor2, Texponent2}
        end
    elseif (factor2 > 1) && (factor1 == 1) && (exponent1 < exponent2)
        return Period{T, Tfactor1, Texponent1}
    elseif (factor1 > 1) && (factor2 == 1) && (exponent2 < exponent1)
        return Period{T, Tfactor2, Texponent2}
    else
        if factor1 / 10^(-exponent1) <= factor2 / 10^(-exponent2)
            return Period{T, Tfactor1, Texponent1}
        else
            return Period{T, Tfactor2, Texponent2}
        end
    end
end

for Tfactor1 in (SOLAR_YEAR, SOLAR_YEAR ÷ 12)
    @eval function promote_rule(
            ::Type{Period{T1, Val($Tfactor1), Val(-3)}},
            ::Type{Period{T2, Tfactor2, Texponent2}}
        ) where
        {T1, T2, Tfactor2, Texponent2}
        return promote_rule(
            Period{T1, Val(1), Val(-3)},
            Period{T2, Tfactor2, Texponent2}
        )
    end
end


@inline function convert(
        ::Type{Period{T1, Tfactor1, Texponent1}},
        p::Period{T2, Tfactor2, Texponent2}
    ) where
    {T1, Tfactor1, Texponent1, T2, Tfactor2, Texponent2}

    factor1 = unwrap(Tfactor1)
    factor2 = unwrap(Tfactor2)
    exponent1 = unwrap(Texponent1)
    exponent2 = unwrap(Texponent2)

    exp = exponent2 - exponent1

    if T1 <: AbstractFloat
        duration = (T1(p.duration) * factor2 * 10^(exp)) / factor1
    else
        (num, denom) =
        if exp > 0
            ((T1(p.duration) * factor2 * 10^(exp)), factor1)
        else
            ((T1(p.duration) * factor2), (factor1 * 10^(-exp)))
        end

        if num % denom != 0
            throw(InexactError(:convert, Type{Period{T1, Tfactor1, Texponent1}}, p))
        end

        duration = num ÷ denom
    end

    return Period{T1, Tfactor1, Texponent1}(duration)
end

function convert(::Type{Period{T, Tfactor, Texponent}}, t::Union{Dates.Day, Dates.TimePeriod}) where {T, Tfactor, Texponent}
    p = convert(Period, t)
    return convert(Period{T, Tfactor, Texponent}, p)
end

# operations between two CFTime.Periods (or Dates.Period)

# operators returning a CFTime.Period
for op in (:+, :-, :mod, :lcm, :gcd, :rem)
    @eval begin
        function $op(p1::Period{T, Tfactor, Texponent}, p2::Period{T, Tfactor, Texponent}) where {T, Tfactor, Texponent}
            return Period{T, Tfactor, Texponent}($op(p1.duration, p2.duration))
        end
    end
end

function Base.gcdx(a::T, b::T) where {T <: Period}
    (g, x, y) = gcdx(value(a), value(b))
    return (T(g), x, y)
end

# operators not returning a CFTime.Period
for op in (:/, :div, :(==), :isless)
    @eval begin
        function $op(p1::Period{T, Tfactor, Texponent}, p2::Period{T, Tfactor, Texponent}) where {T, Tfactor, Texponent}
            return $op(p1.duration, p2.duration)
        end
    end
end


for op in (:+, :-, :/, :div, :mod, :(==), :isless, :lcm, :gcd, :gcdx, :rem)
    @eval begin
        $op(p1::Period, p2::Period) = $op(promote(p1, p2)...)
        $op(p1::Period, p2::Dates.Period) = $op(promote(p1, p2)...)
        $op(p1::Dates.Period, p2::Period) = $op(promote(p1, p2)...)
    end
end

# operations between CFTime.Period and a number
for op in (:*, :/, :div)
    @eval begin
        function $op(p::Period{T, Tfactor, Texponent}, v::Number) where {T, Tfactor, Texponent}
            pv = $op(p.duration, v)
            return Period{typeof(pv), Tfactor, Texponent}(pv)
        end
    end
end

*(v::Number, p::Period) = p * v

function -(p::Period{T, Tfactor, Texponent}) where {T, Tfactor, Texponent}
    return Period{T, Tfactor, Texponent}(-p.duration)
end


for T in (:Day, :Hour, :Minute, :Second, :Millisecond, :Microsecond, :Nanosecond)
    unit = Symbol(lowercase(string(T)))
    factor, exponent = filter(td -> td[1] == unit, TIME_DIVISION)[1][2:end]

    @eval begin
        convert(::Type{Period}, t::Dates.$T) =
            Period{Int64, Val($factor), Val($exponent)}(Dates.value(t))

        function promote_rule(
                ::Type{Period{T, Tfactor, Texponent}},
                ::Type{Dates.$T}
            ) where {T, Tfactor, Texponent}
            return promote_type(
                Period{T, Tfactor, Texponent},
                Period{Int64, Val($factor), Val($exponent)}
            )
        end

        # Can throw an InexactError
        @inline function Dates.$T(p::Period{T, Val($factor), Val($exponent)}) where {T}
            return Dates.$T(Int64(p.duration))
        end

        @inline function Dates.$T(p::Period)
            return Dates.$T(convert(Period{Int64, Val($factor), Val($exponent)}, p))
        end
    end
end

# picoseconds and smaller
for (name, factor, exponent) in TIME_DIVISION[8:end]
    name_str = string(name)
    name_titlecase = Symbol(titlecase(name_str))

    @eval begin
        @doc """
             CFTime.$($name_titlecase)(v::Number)

        Construct an object representing the duration of `v` $($name_str)s.
        """ $name_titlecase(d::T) where {T <: Number} = Period{T, Val($factor), Val($exponent)}(d)
    end

    if VERSION >= v"1.11"
        eval(Meta.parse("public $name_titlecase"))
    end
end


function units(p::Period{T, Tfactor, Texponent}) where {T, Tfactor, Texponent}

    for (name, factor, exponent) in TIME_DIVISION
        if (Val(factor) == Tfactor) && (Val(exponent) == Texponent)
            # always append s for plural
            return string(name, "s")
        end
    end

    return string(unwrap(Tfactor), " × 10^", unwrap(Texponent), " s")
end


function Base.show(io::IO, p::Period)
    exp = _exponent(p)
    fact = _factor(p)

    time_divisions = (
        (:solar_year, SOLAR_YEAR, -3),
        (:solar_month, SOLAR_YEAR ÷ 12, -3),
        TIME_DIVISION...,
    )

    for (name, factor, exponent) in time_divisions
        if (fact == factor) && (exp == exponent)
            print(io, "$(p.duration) $name")
            if p.duration != 1
                print(io, "s")
            end
            return
        end
    end
    print(io, "$(p.duration * fact) ")
    if exp != 0
        print(io, "× 10^($(exp)) ")
    end
    return print(io, "s")
end


# Missing support
(==)(x::Period, y::Missing) = missing
(==)(x::Missing, y::Period) = missing


typemax(::Type{Period{T, Tfactor, Texponent}}) where {T, Tfactor, Texponent} = Period{T, Tfactor, Texponent}(typemax(T))

typemin(::Type{Period{T, Tfactor, Texponent}}) where {T, Tfactor, Texponent} = Period{T, Tfactor, Texponent}(typemin(T))

typemax(p::T) where {T <: Period} = typemax(T)
typemin(p::T) where {T <: Period} = typemin(T)
