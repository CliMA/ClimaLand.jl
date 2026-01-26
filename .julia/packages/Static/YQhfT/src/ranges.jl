
"""
    OptionallyStaticUnitRange(start, stop) <: AbstractUnitRange{Int}

Similar to `UnitRange` except each field may be an `Int` or `StaticInt`. An
`OptionallyStaticUnitRange` is intended to be constructed internally from other valid
indices. Therefore, users should not expect the same checks are used to ensure construction
of a valid `OptionallyStaticUnitRange` as a `UnitRange`.
"""
struct OptionallyStaticUnitRange{F <: IntType, L <: IntType} <:
       AbstractUnitRange{Int}
    start::F
    stop::L

    function OptionallyStaticUnitRange(start::IntType,
            stop::IntType)
        new{typeof(start), typeof(stop)}(start, stop)
    end
    function OptionallyStaticUnitRange(start, stop)
        OptionallyStaticUnitRange(IntType(start), IntType(stop))
    end
    OptionallyStaticUnitRange(@nospecialize x::OptionallyStaticUnitRange) = x
    function OptionallyStaticUnitRange(x::AbstractRange)
        step(x) == 1 && return OptionallyStaticUnitRange(static_first(x), static_last(x))

        errmsg(x) = throw(ArgumentError("step must be 1, got $(step(x))")) # avoid GC frame
        errmsg(x)
    end
    function OptionallyStaticUnitRange{F, L}(x::AbstractRange) where {F, L}
        OptionallyStaticUnitRange(x)
    end
    function OptionallyStaticUnitRange{StaticInt{F}, StaticInt{L}}() where {F, L}
        new{StaticInt{F}, StaticInt{L}}()
    end
end

"""
    OptionallyStaticStepRange(start, step, stop) <: OrdinalRange{Int,Int}

Similarly to [`OptionallyStaticUnitRange`](@ref), `OptionallyStaticStepRange` permits
a combination of static and standard primitive `Int`s to construct a range. It
specifically enables the use of ranges without a step size of 1. It may be constructed
through the use of `OptionallyStaticStepRange` directly or using static integers with
the range operator (i.e., `:`).

```julia
julia> using Static

julia> x = static(2);

julia> x:x:10
static(2):static(2):10

julia> Static.OptionallyStaticStepRange(x, x, 10)
static(2):static(2):10

```
"""
struct OptionallyStaticStepRange{F <: IntType, S <: IntType,
    L <: IntType} <: OrdinalRange{Int, Int}
    start::F
    step::S
    stop::L

    global function _OptionallyStaticStepRange(@nospecialize(start::IntType),
            @nospecialize(step::IntType),
            @nospecialize(stop::IntType))
        new{typeof(start), typeof(step), typeof(stop)}(start, step, stop)
    end
end
@noinline function OptionallyStaticStepRange(@nospecialize(start::IntType),
        ::StaticInt{0},
        @nospecialize(stop::IntType))
    throw(ArgumentError("step cannot be zero"))
end
# we don't need to check the `stop` if we know it acts like a unit range
function OptionallyStaticStepRange(@nospecialize(start::IntType),
        step::StaticInt{1},
        @nospecialize(stop::IntType))
    _OptionallyStaticStepRange(start, step, stop)
end
function OptionallyStaticStepRange(@nospecialize(start::IntType),
        @nospecialize(step::StaticInt),
        @nospecialize(stop::IntType))
    _OptionallyStaticStepRange(start, step, _steprange_last(start, step, stop))
end
function OptionallyStaticStepRange(start, step, stop)
    OptionallyStaticStepRange(IntType(start), IntType(step), IntType(stop))
end
function OptionallyStaticStepRange(@nospecialize(start::IntType),
        step::Int,
        @nospecialize(stop::IntType))
    if step === 0
        throw(ArgumentError("step cannot be zero"))
    else
        _OptionallyStaticStepRange(start, step, _steprange_last(start, step, stop))
    end
end
OptionallyStaticStepRange(@nospecialize x::OptionallyStaticStepRange) = x
function OptionallyStaticStepRange(x::AbstractRange)
    _OptionallyStaticStepRange(IntType(static_first(x)), IntType(static_step(x)),
        IntType(static_last(x)))
end

# to make StepRange constructor inlineable, so optimizer can see `step` value
@inline function _steprange_last(start::StaticInt, step::StaticInt, stop::StaticInt)
    StaticInt(_steprange_last(Int(start), Int(step), Int(stop)))
end
@inline function _steprange_last(start::Union{StaticInt, Int},
        step::Union{StaticInt, Int},
        stop::Union{StaticInt, Int})
    _steprange_last(Int(start), Int(step), Int(stop))
end
@inline function _steprange_last(start::Int, step::Int, stop::Int)
    if stop === start
        return stop
    elseif step > 0
        if stop > start
            return stop - rem(stop - start, step)
        else
            return start - 1
        end
    else
        if stop > start
            return start + 1
        else
            return stop + rem(start - stop, -step)
        end
    end
end

"""
    SUnitRange(start::Int, stop::Int)

An alias for `OptionallyStaticUnitRange` where both the start and stop are known statically.
"""
const SUnitRange{F, L} = OptionallyStaticUnitRange{StaticInt{F}, StaticInt{L}}
SUnitRange(start::Int, stop::Int) = SUnitRange{start, stop}()

"""
    SOneTo(n::Int)

An alias for `OptionallyStaticUnitRange` usfeul for statically sized axes.
"""
const SOneTo{L} = SUnitRange{1, L}
SOneTo(n::Int) = SOneTo{n}()
Base.oneto(::StaticInt{N}) where {N} = SOneTo{N}()

const OptionallyStaticRange{
    F, L} = Union{OptionallyStaticUnitRange{F, L},
    OptionallyStaticStepRange{F, <:Any, L}}

"""
    static_first(x::AbstractRange)

Attempt to return `static(first(x))`, if known at compile time. Otherwise, return
`first(x)`.

See also: [`static_step`](@ref), [`static_last`](@ref)

# Examples

```julia
julia> static_first(static(2):10)
static(2)

julia> static_first(1:10)
1

julia> static_first(Base.OneTo(10))
static(1)

```
"""
static_first(x::Base.OneTo) = StaticInt(1)
static_first(x::Union{Base.Slice, Base.IdentityUnitRange}) = static_first(x.indices)
static_first(x::OptionallyStaticRange) = getfield(x, :start)
static_first(x) = first(x)

"""
    static_step(x::AbstractRange)

Attempt to return `static(step(x))`, if known at compile time. Otherwise, return
`step(x)`.

See also: [`static_first`](@ref), [`static_last`](@ref)

# Examples

```julia
julia> static_step(static(1):static(3):9)
static(3)

julia> static_step(1:3:9)
3

julia> static_step(1:9)
static(1)

```
"""
static_step(@nospecialize x::AbstractUnitRange) = StaticInt(1)
static_step(x::OptionallyStaticStepRange) = getfield(x, :step)
static_step(x) = step(x)

"""
    static_last(x::AbstractRange)

Attempt to return `static(last(x))`, if known at compile time. Otherwise, return
`last(x)`.

See also: [`static_first`](@ref), [`static_step`](@ref)

# Examples

```julia
julia> static_last(static(1):static(10))
static(10)

julia> static_last(static(1):10)
10

```
"""
static_last(x::OptionallyStaticRange) = getfield(x, :stop)
static_last(x) = last(x)
static_last(x::Union{Base.Slice, Base.IdentityUnitRange}) = static_last(x.indices)

Base.first(x::OptionallyStaticRange{Int}) = getfield(x, :start)
Base.first(::OptionallyStaticRange{StaticInt{F}}) where {F} = F
Base.step(x::OptionallyStaticStepRange{<:Any, Int}) = getfield(x, :step)
Base.step(::OptionallyStaticStepRange{<:Any, StaticInt{S}}) where {S} = S
Base.last(x::OptionallyStaticRange{<:Any, Int}) = getfield(x, :stop)
Base.last(::OptionallyStaticRange{<:Any, StaticInt{L}}) where {L} = L

# FIXME this line causes invalidations
Base.:(:)(L::Integer, ::StaticInt{U}) where {U} = OptionallyStaticUnitRange(L, StaticInt(U))
Base.:(:)(::StaticInt{L}, U::Integer) where {L} = OptionallyStaticUnitRange(StaticInt(L), U)
function Base.:(:)(::StaticInt{L}, ::StaticInt{U}) where {L, U}
    OptionallyStaticUnitRange(StaticInt(L), StaticInt(U))
end
function Base.:(:)(::StaticInt{F}, ::StaticInt{S}, ::StaticInt{L}) where {F, S, L}
    OptionallyStaticStepRange(StaticInt(F), StaticInt(S), StaticInt(L))
end
function Base.:(:)(start::Integer, ::StaticInt{S}, ::StaticInt{L}) where {S, L}
    OptionallyStaticStepRange(start, StaticInt(S), StaticInt(L))
end
function Base.:(:)(::StaticInt{F}, ::StaticInt{S}, stop::Integer) where {F, S}
    OptionallyStaticStepRange(StaticInt(F), StaticInt(S), stop)
end
function Base.:(:)(::StaticInt{F}, step::Integer, ::StaticInt{L}) where {F, L}
    OptionallyStaticStepRange(StaticInt(F), step, StaticInt(L))
end
function Base.:(:)(start::Integer, step::Integer, ::StaticInt{L}) where {L}
    OptionallyStaticStepRange(start, step, StaticInt(L))
end
function Base.:(:)(start::Integer, ::StaticInt{S}, stop::Integer) where {S}
    OptionallyStaticStepRange(start, StaticInt(S), stop)
end
function Base.:(:)(::StaticInt{F}, step::Integer, stop::Integer) where {F}
    OptionallyStaticStepRange(StaticInt(F), step, stop)
end
Base.:(:)(start::StaticInt{F}, ::StaticInt{1}, stop::StaticInt{L}) where {F, L} = start:stop
Base.:(:)(start::Integer, ::StaticInt{1}, stop::StaticInt{L}) where {L} = start:stop
Base.:(:)(start::StaticInt{F}, ::StaticInt{1}, stop::Integer) where {F} = start:stop
function Base.:(:)(start::Integer, ::StaticInt{1}, stop::Integer)
    OptionallyStaticUnitRange(start, stop)
end

Base.isempty(r::OptionallyStaticUnitRange) = first(r) > last(r)
@inline function Base.isempty(x::OptionallyStaticStepRange)
    start = first(x)
    stop = last(x)
    if start === stop
        return false
    else
        s = step(x)
        s > 0 ? start > stop : start < stop
    end
end

function Base.checkindex(::Type{Bool},
        ::SUnitRange{F1, L1},
        ::SUnitRange{F2, L2}) where {F1, L1, F2, L2}
    (F1::Int <= F2::Int) && (L1::Int >= L2::Int)
end

function Base.getindex(r::OptionallyStaticUnitRange,
        s::AbstractUnitRange{<:Integer})
    @boundscheck checkbounds(r, s)
    f = static_first(r)
    fnew = f - one(f)
    return (fnew + static_first(s)):(fnew + static_last(s))
end

function Base.getindex(x::OptionallyStaticUnitRange{StaticInt{1}}, i::Int)
    @boundscheck checkbounds(x, i)
    i
end
function Base.getindex(x::OptionallyStaticUnitRange, i::Int)
    val = first(x) + (i - 1)
    @boundscheck ((i < 1) || val > last(x)) && throw(BoundsError(x, i))
    val::Int
end

## length
@inline function Base.length(x::OptionallyStaticUnitRange)
    start = first(x)
    stop = last(x)
    start > stop ? 0 : stop - start + 1
end
Base.length(r::OptionallyStaticStepRange) = _range_length(first(r), step(r), last(r))
@inline function _range_length(start::Int, s::Int, stop::Int)
    if s > 0
        stop < start ? 0 : div(stop - start, s) + 1
    else
        stop > start ? 0 : div(start - stop, -s) + 1
    end
end

Base.AbstractUnitRange{Int}(r::OptionallyStaticUnitRange) = r
function Base.AbstractUnitRange{T}(r::OptionallyStaticUnitRange) where {T}
    start = static_first(r)
    if isa(start, StaticInt{1}) && T <: Integer
        return Base.OneTo{T}(T(static_last(r)))
    else
        return UnitRange{T}(T(static_first(r)), T(static_last(r)))
    end
end

Base.isdone(x::OptionallyStaticRange, state::Int) = state === last(x)
function _next(x::OptionallyStaticRange)
    new_state = first(x)
    (new_state, new_state)
end
@inline function _next(@nospecialize(x::OptionallyStaticUnitRange), state::Int)
    new_state = state + 1
    (new_state, new_state)
end
@inline function _next(x::OptionallyStaticStepRange, state::Int)
    new_state = state + step(x)
    (new_state, new_state)
end
@inline Base.iterate(x::OptionallyStaticRange) = isempty(x) ? nothing : _next(x)
@inline function Base.iterate(x::OptionallyStaticRange, s::Int)
    Base.isdone(x, s) ? nothing : _next(x, s)
end

Base.to_shape(x::OptionallyStaticRange) = length(x)
Base.to_shape(x::Base.Slice{T}) where {T <: OptionallyStaticRange} = length(x)
Base.axes(S::Base.Slice{<:OptionallyStaticUnitRange{StaticInt{1}}}) = (S.indices,)
Base.axes(S::Base.Slice{<:OptionallyStaticRange}) = (Base.IdentityUnitRange(S.indices),)

Base.axes(x::OptionallyStaticRange) = (Base.axes1(x),)
function Base.axes1(x::OptionallyStaticUnitRange)
    OptionallyStaticUnitRange(StaticInt(1), length(x))
end
Base.axes1(x::OptionallyStaticUnitRange{StaticInt{1}}) = x
function Base.axes1(x::OptionallyStaticUnitRange{StaticInt{F}, StaticInt{L}}) where {F, L}
    OptionallyStaticUnitRange(StaticInt(1), StaticInt(L - F + 1))
end
function Base.axes1(x::OptionallyStaticStepRange)
    OptionallyStaticUnitRange(StaticInt(1), length(x))
end
function Base.axes1(x::OptionallyStaticStepRange{
        StaticInt{F}, StaticInt{S}, StaticInt{L}}) where {
        F,
        S,
        L
}
    OptionallyStaticUnitRange(StaticInt(1), StaticInt(_range_length(F, S, L)))
end
Base.axes1(x::Base.Slice{<:OptionallyStaticUnitRange{One}}) = x.indices
Base.axes1(x::Base.Slice{<:OptionallyStaticRange}) = Base.IdentityUnitRange(x.indices)

Base.:(-)(r::OptionallyStaticRange) = (-static_first(r)):(-static_step(r)):(-static_last(r))

function Base.reverse(x::OptionallyStaticUnitRange)
    _OptionallyStaticStepRange(getfield(x, :stop), StaticInt(-1), getfield(x, :start))
end
function Base.reverse(x::OptionallyStaticStepRange)
    _OptionallyStaticStepRange(getfield(x, :stop), -getfield(x, :step), getfield(x, :start))
end

Base.show(io::IO, @nospecialize(x::OptionallyStaticRange)) = show(io, MIME"text/plain"(), x)
function Base.show(io::IO, ::MIME"text/plain", @nospecialize(r::OptionallyStaticUnitRange))
    print(io, "$(getfield(r, :start)):$(getfield(r, :stop))")
end
function Base.show(io::IO, ::MIME"text/plain", @nospecialize(r::OptionallyStaticStepRange))
    print(io, "$(getfield(r, :start)):$(getfield(r, :step)):$(getfield(r, :stop))")
end

# we overload properties because occasionally Base assumes that abstract range types have
# the same exact same set up as native types where `x.start === first(x)`
@inline function Base.getproperty(x::OptionallyStaticRange, s::Symbol)
    if s === :start
        return first(x)
    elseif s === :step
        return step(x)
    elseif s === :stop
        return last(x)
    else
        error("$x has no property $s")
    end
end

function Base.Broadcast.axistype(r::OptionallyStaticUnitRange{StaticInt{1}}, _)
    Base.OneTo(last(r))
end
function Base.Broadcast.axistype(_, r::OptionallyStaticUnitRange{StaticInt{1}})
    Base.OneTo(last(r))
end
function Base.Broadcast.axistype(r::OptionallyStaticUnitRange{StaticInt{1}},
        ::OptionallyStaticUnitRange{StaticInt{1}})
    Base.OneTo(last(r))
end
function Base.similar(::Type{<:Array{T}},
        axes::Tuple{OptionallyStaticUnitRange{StaticInt{1}},
            Vararg{
                Union{Base.OneTo,
                OptionallyStaticUnitRange{StaticInt{1}}}}}) where {
        T
}
    Array{T}(undef, map(last, axes))
end
function Base.similar(::Type{<:Array{T}},
        axes::Tuple{Base.OneTo, OptionallyStaticUnitRange{StaticInt{1}},
            Vararg{
                Union{Base.OneTo,
                OptionallyStaticUnitRange{StaticInt{1}}}}}) where {
        T
}
    Array{T}(undef, map(last, axes))
end

function Base.first(x::OptionallyStaticUnitRange, n::IntType)
    n < 0 && throw(ArgumentError("Number of elements must be nonnegative"))
    start = static_first(x)
    OptionallyStaticUnitRange(start, min(start - one(start) + n, static_last(x)))
end
function Base.first(x::OptionallyStaticStepRange, n::IntType)
    n < 0 && throw(ArgumentError("Number of elements must be nonnegative"))
    start = static_first(x)
    s = static_step(x)
    stop = min(((n - one(n)) * s) + static_first(x), static_last(x))
    OptionallyStaticStepRange(start, s, stop)
end
function Base.last(x::OptionallyStaticUnitRange, n::IntType)
    n < 0 && throw(ArgumentError("Number of elements must be nonnegative"))
    stop = static_last(x)
    OptionallyStaticUnitRange(max(stop + one(stop) - n, static_first(x)), stop)
end
function Base.last(x::OptionallyStaticStepRange, n::IntType)
    n < 0 && throw(ArgumentError("Number of elements must be nonnegative"))
    start = static_first(x)
    s = static_step(x)
    stop = static_last(x)
    OptionallyStaticStepRange(max(stop + one(stop) - (n * s), start), s, stop)
end
