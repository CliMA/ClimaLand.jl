module Static

import IfElse: ifelse
using SciMLPublic: @public

export StaticInt, StaticFloat64, StaticSymbol, True, False, StaticBool, NDIndex
export dynamic, is_static, known, static, static_promote, static_first, static_step,
       static_last

@public OptionallyStaticRange,
OptionallyStaticUnitRange, OptionallyStaticStepRange, SUnitRange, SOneTo
@public eachop, eachop_tuple, reduce_tup, eq, ne, gt, ge, le, lt, mul, add

import PrecompileTools: @recompile_invalidations

@recompile_invalidations begin
    import CommonWorldInvalidations
end

"""
    StaticSymbol(S::Symbol)::StaticSymbol{S}

A statically typed `Symbol`.
"""
struct StaticSymbol{s}
    StaticSymbol{s}() where {s} = new{s::Symbol}()
    StaticSymbol(s::Symbol) = new{s}()
    StaticSymbol(@nospecialize x::StaticSymbol) = x
    StaticSymbol(x) = StaticSymbol(Symbol(x))
    StaticSymbol(x, y) = StaticSymbol(Symbol(x, y))
    @generated function StaticSymbol(::StaticSymbol{X}, ::StaticSymbol{Y}) where {X, Y}
        StaticSymbol(Symbol(X, Y))
    end
    StaticSymbol(x, y, z...) = StaticSymbol(StaticSymbol(x, y), z...)
end

Base.Symbol(@nospecialize(s::StaticSymbol)) = known(s)

abstract type StaticInteger{N} <: Number end

"""
    StaticBool(x::Bool)::Union{True, False}

A statically typed `Bool`.
"""
abstract type StaticBool{bool} <: StaticInteger{bool} end

struct True <: StaticBool{true} end

struct False <: StaticBool{false} end

StaticBool{true}() = True()
StaticBool{false}() = False()
StaticBool(x::StaticBool) = x
function StaticBool(x::Bool)
    if x
        return True()
    else
        return False()
    end
end

"""
    StaticInt(N::Int)::StaticInt{N}

A statically sized `Int`.
Use `StaticInt(N)` instead of `Val(N)` when you want it to behave like a number.
"""
struct StaticInt{N} <: StaticInteger{N}
    StaticInt{N}() where {N} = new{N::Int}()
    StaticInt(N::Int) = new{N}()
    StaticInt(@nospecialize N::StaticInt) = N
    StaticInt(::Val{N}) where {N} = StaticInt(N)
end

"""
    IntType(x::Integer)::Union{Int, StaticInt}

`IntType` is a union of `Int` and `StaticInt`. As a function, it ensures that `x` one of the
two.
"""
const IntType = Union{StaticInt, Int}
IntType(x::Integer) = Int(x)
IntType(@nospecialize x::Union{Int, StaticInt}) = x

Base.isinteger(@nospecialize x::StaticInteger) = true

include("float.jl")

const StaticNumber{N} = Union{StaticInt{N}, StaticBool{N}, StaticFloat64{N}}

Base.promote_rule(::Type{Bool}, T::Type{<:StaticNumber}) = promote_rule(Bool, eltype(T))

Base.getindex(x::Tuple, ::StaticInt{N}) where {N} = getfield(x, N)

Base.getindex(x::NamedTuple, ::StaticSymbol{N}) where {N} = getfield(x, N)
function Base.getindex(x::NamedTuple, symbols::Tuple{StaticSymbol, Vararg{StaticSymbol}})
    return x[map(known, symbols)]
end

Base.zero(@nospecialize(::StaticInt)) = StaticInt{0}()

Base.to_index(x::StaticInt) = known(x)
function Base.checkindex(::Type{Bool}, inds::AbstractUnitRange, ::StaticNumber{N}) where {N}
    checkindex(Bool, inds, N)
end
function Base.checkindex(::Type{Bool}, inds::Base.IdentityUnitRange,
        ::StaticFloat64{N}) where {N}
    checkindex(Bool, inds.indices, N)
end

ifelse(::True, @nospecialize(x), @nospecialize(y)) = x
ifelse(::False, @nospecialize(x), @nospecialize(y)) = y

const Zero = StaticInt{0}
const One = StaticInt{1}
const FloatOne = StaticFloat64{one(Float64)}
const FloatZero = StaticFloat64{zero(Float64)}

const StaticType{T} = Union{StaticNumber{T}, StaticSymbol{T}}

StaticInt(x::False) = Zero()
StaticInt(x::True) = One()
Base.Bool(::True) = true
Base.Bool(::False) = false

Base.eltype(@nospecialize(T::Type{<:StaticFloat64})) = Float64
Base.eltype(@nospecialize(T::Type{<:StaticInt})) = Int
Base.eltype(@nospecialize(T::Type{<:StaticBool})) = Bool

"""
    NDIndex(i, j, k...)   -> I
    NDIndex((i, j, k...)) -> I

A multidimensional index that refers to a single element. Each dimension is represented by
a single `Int` or `StaticInt`.

```julia
julia> using Static

julia> i = NDIndex(static(1), 2, static(3))
NDIndex(static(1), 2, static(3))

julia> i[static(1)]
static(1)

julia> i[1]
1

```
"""
struct NDIndex{N, I <: Tuple{Vararg{Union{StaticInt, Int}, N}}} <:
       Base.AbstractCartesianIndex{N}
    index::I

    NDIndex{N}(i::Tuple{Vararg{Union{StaticInt, Int}, N}}) where {N} = new{N, typeof(i)}(i)
    NDIndex{N}(index::Tuple) where {N} = _ndindex(static(N), _flatten(index...))
    NDIndex{N}(index...) where {N} = NDIndex{N}(index)

    NDIndex{0}(::Tuple{}) = new{0, Tuple{}}(())
    NDIndex{0}() = NDIndex{0}(())

    NDIndex(i::Tuple{Vararg{Union{StaticInt, Int}, N}}) where {N} = new{N, typeof(i)}(i)
    NDIndex(i::Vararg{Union{StaticInt, Int}, N}) where {N} = NDIndex(i)

    NDIndex(index::Tuple) = NDIndex(_flatten(index...))
    NDIndex(index...) = NDIndex(index)
end

_ndindex(n::StaticInt{N}, i::Tuple{Vararg{Union{Int, StaticInt}, N}}) where {N} = NDIndex(i)
function _ndindex(n::StaticInt{N}, i::Tuple{Vararg{Any, M}}) where {N, M}
    M > N && throw(ArgumentError("input tuple of length $M, requested $N"))
    return NDIndex(_fill_to_length(i, n))
end
_fill_to_length(x::Tuple{Vararg{Any, N}}, n::StaticInt{N}) where {N} = x
@inline function _fill_to_length(x::Tuple{Vararg{Any, M}}, n::StaticInt{N}) where {M, N}
    return _fill_to_length((x..., static(1)), n)
end

_flatten(i::StaticInt{N}) where {N} = (i,)
_flatten(i::Integer) = (Int(i),)
_flatten(i::Base.AbstractCartesianIndex) = _flatten(Tuple(i)...)
@inline _flatten(i::StaticInt, I...) = (i, _flatten(I...)...)
@inline _flatten(i::Integer, I...) = (Int(i), _flatten(I...)...)
@inline function _flatten(i::Base.AbstractCartesianIndex, I...)
    return (_flatten(Tuple(i)...)..., _flatten(I...)...)
end
Base.Tuple(@nospecialize(x::NDIndex)) = getfield(x, :index)

"""
    known(T::Type)

Returns the known value corresponding to a static type `T`. If `T` is not a static type then
`nothing` is returned.

`known` ensures that the type of the returned value is always inferred, even if the
compiler fails to infer the exact value.

See also: [`static`](@ref), [`is_static`](@ref), [`dynamic`](@ref)

# Examples
```julia
julia> known(StaticInt{1})
1

julia> known(Int)

```
"""
known(@nospecialize(x)) = known(typeof(x))
_get_known(::Type{T}, dim::StaticInt{D}) where {T, D} = known(field_type(T, dim))
known(@nospecialize(T::Type{<:Tuple})) = eachop(_get_known, nstatic(Val(fieldcount(T))), T)
known(T::DataType) = nothing
known(@nospecialize(T::Type{<:NDIndex})) = known(T.parameters[2])
known(::Union{True, Type{True}}) = true
known(::Union{False, Type{False}}) = false
known(::Union{Val{V}, Type{Val{V}}}) where {V} = V
known(::Union{StaticInt{N}, Type{StaticInt{N}}}) where {N} = _return_int(N)
_return_int(x::Int) = x
known(::Union{StaticSymbol{S}, Type{StaticSymbol{S}}}) where {S} = _return_symbol(S)
_return_symbol(x::Symbol) = x
known(::Union{StaticFloat64{F}, Type{StaticFloat64{F}}}) where {F} = _return_float(F)
_return_float(x::Float64) = x

"""
    static(x)

Returns a static form of `x`. If `x` is already in a static form then `x` is returned. If
there is no static alternative for `x` then an error is thrown.

See also: [`is_static`](@ref), [`known`](@ref)

# Examples

```julia
julia> using Static

julia> static(1)
static(1)

julia> static(true)
True()

julia> static(:x)
static(:x)

```
"""
static(@nospecialize(x::Union{StaticSymbol, StaticNumber})) = x
static(x::Integer) = StaticInt(x)
function static(x::Union{AbstractFloat, Complex, Rational, AbstractIrrational})
    StaticFloat64(Float64(x))
end
static(x::Bool) = StaticBool(x)
static(x::Union{Symbol, AbstractChar, AbstractString}) = StaticSymbol(x)
static(x::Tuple{Vararg{Any}}) = map(static, x)
static(::Val{V}) where {V} = static(V)
static(x::CartesianIndex) = NDIndex(static(Tuple(x)))
function static(x::X) where {X}
    Base.issingletontype(X) && return x
    error("There is no static alternative for type $(typeof(x)).")
end

"""
    is_static(::Type{T})::Union{True, False}

If `T` is a static type return `static(true)::True` and otherwise returns
`static(false)::False`

See also: [`static`](@ref), [`known`](@ref)
"""
is_static(@nospecialize(x)) = is_static(typeof(x))
is_static(@nospecialize(x::Type{<:StaticType})) = True()
is_static(@nospecialize(x::Type{<:Val})) = True()
_tuple_static(::Type{T}, i) where {T} = is_static(field_type(T, i))
@inline function is_static(@nospecialize(T::Type{<:Tuple}))
    if all(eachop(_tuple_static, nstatic(Val(fieldcount(T))), T))
        return True()
    else
        return False()
    end
end
is_static(T::DataType) = False()

"""
    dynamic(x)

Returns the "dynamic" or non-static form of `x`. If `x` is not a static type, then it is
returned unchanged.

`dynamic` ensures that the type of the returned value is always inferred, even if the
compiler fails to infer the exact value.

See also: [`known`](@ref)

# Examples

```julia
julia> dynamic(static(1))
1

julia> dynamic(1)
1

```
"""
@inline dynamic(@nospecialize x::StaticInt) = known(x)
@inline dynamic(@nospecialize x::StaticFloat64) = known(x)
@inline dynamic(@nospecialize x::StaticSymbol) = known(x)
@inline dynamic(@nospecialize x::Union{True, False}) = known(x)
@inline dynamic(@nospecialize x::Tuple) = map(dynamic, x)
dynamic(@nospecialize(x::NDIndex)) = CartesianIndex(dynamic(Tuple(x)))
dynamic(@nospecialize x) = x

"""
    static_promote(x, y)

Throws an error if `x` and `y` are not equal, preferentially returning the one that is known
at compile time.
"""
@inline static_promote(x::StaticType{X}, ::StaticType{X}) where {X} = x
@noinline function static_promote(::StaticType{X}, ::StaticType{Y}) where {X, Y}
    error("$X and $Y are not equal")
end
Base.@propagate_inbounds function static_promote(::StaticType{N}, x) where {N}
    static(static_promote(N, x))
end
Base.@propagate_inbounds function static_promote(x, ::StaticType{N}) where {N}
    static(static_promote(N, x))
end
Base.@propagate_inbounds static_promote(x, y) = _static_promote(x, y)
Base.@propagate_inbounds function _static_promote(x, y)
    @boundscheck x === y || error("$x and $y are not equal")
    x
end
_static_promote(::Nothing, ::Nothing) = nothing
_static_promote(x, ::Nothing) = x
_static_promote(::Nothing, y) = y

"""
    static_promote(x::AbstractRange{<:Integer}, y::AbstractRange{<:Integer})

A type stable method for combining two equal ranges into a new range that preserves static
parameters. Throws an error if `x != y`.

# Examples

```julia
julia> static_promote(static(1):10, 1:static(10))
static(1):static(10)

julia> static_promote(1:2:9, static(1):static(2):static(9))
static(1):static(2):static(9)
```
"""
Base.@propagate_inbounds @inline function static_promote(x::AbstractUnitRange{<:Integer},
        y::AbstractUnitRange{<:Integer})
    fst = static_promote(static_first(x), static_first(y))
    lst = static_promote(static_last(x), static_last(y))
    return OptionallyStaticUnitRange(fst, lst)
end
Base.@propagate_inbounds @inline function static_promote(x::AbstractRange{<:Integer},
        y::AbstractRange{<:Integer})
    fst = static_promote(static_first(x), static_first(y))
    stp = static_promote(static_step(x), static_step(y))
    lst = static_promote(static_last(x), static_last(y))
    return _OptionallyStaticStepRange(fst, stp, lst)
end
function static_promote(x::Base.Slice, y::Base.Slice)
    Base.Slice(static_promote(x.indices, y.indices))
end

Base.@propagate_inbounds function _promote_shape(a::Tuple{A, Vararg{Any}},
        b::Tuple{B, Vararg{Any}}) where {A, B}
    (static_promote(getfield(a, 1), getfield(b, 1)),
        _promote_shape(Base.tail(a), Base.tail(b))...)
end
_promote_shape(::Tuple{}, ::Tuple{}) = ()
Base.@propagate_inbounds function _promote_shape(::Tuple{}, b::Tuple{B}) where {B}
    (static_promote(static(1), getfield(b, 1)),)
end
Base.@propagate_inbounds function _promote_shape(a::Tuple{A}, ::Tuple{}) where {A}
    (static_promote(static(1), getfield(a, 1)),)
end
Base.@propagate_inbounds function Base.promote_shape(
        a::Tuple{
            Vararg{Union{Int, StaticInt}},
        },
        b::Tuple{Vararg{Union{Int, StaticInt}}
        })
    _promote_shape(a, b)
end

function Base.promote_rule(@nospecialize(T1::Type{<:StaticNumber}),
        @nospecialize(T2::Type{<:StaticNumber}))
    promote_rule(eltype(T1), eltype(T2))
end
function Base.promote_rule(::Type{<:Base.TwicePrecision{R}},
        @nospecialize(T::Type{<:StaticNumber})) where {R <: Number}
    promote_rule(Base.TwicePrecision{R}, eltype(T))
end
function Base.promote_rule(@nospecialize(T1::Type{<:StaticNumber}),
        T2::Type{<:Union{Rational, AbstractFloat, Signed}})
    promote_rule(T2, eltype(T1))
end

Base.:(~)(::StaticInteger{N}) where {N} = static(~N)

Base.inv(x::StaticInteger{N}) where {N} = one(x) / x

Base.zero(@nospecialize T::Type{<:StaticInt}) = StaticInt(0)
Base.zero(@nospecialize T::Type{<:StaticBool}) = False()

Base.one(@nospecialize T::Type{<:StaticBool}) = True()
Base.one(@nospecialize T::Type{<:StaticInt}) = StaticInt(1)

@inline Base.iszero(::Union{StaticInt{0}, StaticFloat64{0.0}, False}) = true
@inline Base.iszero(@nospecialize x::Union{StaticInt, True, StaticFloat64}) = false

@inline Base.isone(::Union{One, True, StaticFloat64{1.0}}) = true
@inline Base.isone(@nospecialize x::Union{StaticInt, False, StaticFloat64}) = false

@inline Base.iseven(@nospecialize x::StaticNumber) = iseven(known(x))
@inline Base.isodd(@nospecialize x::StaticNumber) = isodd(known(x))

Base.AbstractFloat(x::StaticNumber) = StaticFloat64(x)

Base.abs(::StaticNumber{N}) where {N} = static(abs(N))
Base.abs2(::StaticNumber{N}) where {N} = static(abs2(N))
Base.sign(::StaticNumber{N}) where {N} = static(sign(N))

Base.widen(@nospecialize(x::StaticNumber)) = widen(known(x))

function Base.convert(::Type{T}, @nospecialize(N::StaticNumber)) where {T <: Number}
    convert(T, known(N))
end

#Base.Bool(::StaticInt{N}) where {N} = Bool(N)

Base.Integer(@nospecialize(x::StaticInt)) = x
(::Type{T})(x::StaticInteger) where {T <: Real} = T(known(x))
function (@nospecialize(T::Type{<:StaticNumber}))(x::Union{AbstractFloat,
        AbstractIrrational, Integer,
        Rational})
    static(convert(eltype(T), x))
end

@inline Base.:(-)(::StaticNumber{N}) where {N} = static(-N)
Base.:(*)(::Union{AbstractFloat, AbstractIrrational, Integer, Rational}, y::Zero) = y
Base.:(*)(x::Zero, ::Union{AbstractFloat, AbstractIrrational, Integer, Rational}) = x
Base.:(*)(::StaticInteger{X}, ::StaticInteger{Y}) where {X, Y} = static(X * Y)
Base.:(/)(::StaticInteger{X}, ::StaticInteger{Y}) where {X, Y} = static(X / Y)
Base.:(-)(::StaticInteger{X}, ::StaticInteger{Y}) where {X, Y} = static(X - Y)
Base.:(+)(::StaticInteger{X}, ::StaticInteger{Y}) where {X, Y} = static(X + Y)

@generated Base.sqrt(::StaticNumber{N}) where {N} = :($(static(sqrt(N))))

function Base.div(::StaticNumber{X}, ::StaticNumber{Y}, m::RoundingMode) where {X, Y}
    static(div(X, Y, m))
end
Base.div(x::Real, ::StaticInteger{Y}, m::RoundingMode) where {Y} = div(x, Y, m)
Base.div(::StaticNumber{X}, y::Real, m::RoundingMode) where {X} = div(X, y, m)
Base.div(x::StaticBool, y::False) = throw(DivideError())
Base.div(x::StaticBool, y::True) = x

Base.rem(@nospecialize(x::StaticNumber), T::Type{<:Integer}) = rem(known(x), T)
Base.rem(::StaticNumber{X}, ::StaticNumber{Y}) where {X, Y} = static(rem(X, Y))
Base.rem(x::Real, ::StaticInteger{Y}) where {Y} = rem(x, Y)
Base.rem(::StaticInteger{X}, y::Real) where {X} = rem(X, y)
Base.mod(::StaticNumber{X}, ::StaticNumber{Y}) where {X, Y} = static(mod(X, Y))
Base.mod(x::Real, ::StaticNumber{Y}) where {Y} = mod(x, Y)
Base.mod(::StaticNumber{X}, y::Real) where {X} = mod(X, y)

Base.:(==)(::StaticNumber{X}, ::StaticNumber{Y}) where {X, Y} = ==(X, Y)

Base.:(<)(::StaticNumber{X}, ::StaticNumber{Y}) where {X, Y} = <(X, Y)

Base.isless(::StaticInteger{X}, ::StaticInteger{Y}) where {X, Y} = isless(X, Y)
Base.isless(::StaticInteger{X}, y::Real) where {X} = isless(X, y)
Base.isless(x::Real, ::StaticInteger{Y}) where {Y} = isless(x, Y)

Base.min(::StaticInteger{X}, ::StaticInteger{Y}) where {X, Y} = static(min(X, Y))
Base.min(::StaticInteger{X}, y::Number) where {X} = min(X, y)
Base.min(x::Number, ::StaticInteger{Y}) where {Y} = min(x, Y)

Base.max(::StaticInteger{X}, ::StaticInteger{Y}) where {X, Y} = static(max(X, Y))
Base.max(::StaticInteger{X}, y::Number) where {X} = max(X, y)
Base.max(x::Number, ::StaticInteger{Y}) where {Y} = max(x, Y)

Base.minmax(::StaticNumber{X}, ::StaticNumber{Y}) where {X, Y} = static(minmax(X, Y))

Base.:(<<)(::StaticInteger{X}, ::StaticInteger{Y}) where {X, Y} = static(<<(X, Y))
Base.:(<<)(::StaticInteger{X}, n::Integer) where {X} = <<(X, n)
Base.:(<<)(x::Integer, ::StaticInteger{N}) where {N} = <<(x, N)

Base.:(>>)(::StaticInteger{X}, ::StaticInteger{Y}) where {X, Y} = static(>>(X, Y))
Base.:(>>)(::StaticInteger{X}, n::Integer) where {X} = >>(X, n)
Base.:(>>)(x::Integer, ::StaticInteger{N}) where {N} = >>(x, N)

Base.:(>>>)(::StaticInteger{X}, ::StaticInteger{Y}) where {X, Y} = static(>>>(X, Y))
Base.:(>>>)(::StaticInteger{X}, n::Integer) where {X} = >>>(X, n)
Base.:(>>>)(x::Integer, ::StaticInteger{N}) where {N} = >>>(x, N)

Base.:(&)(::StaticInteger{X}, ::StaticInteger{Y}) where {X, Y} = static(X & Y)
Base.:(&)(::StaticInteger{X}, y::Union{Integer, Missing}) where {X} = X & y
Base.:(&)(x::Union{Integer, Missing}, ::StaticInteger{Y}) where {Y} = x & Y
Base.:(&)(x::Bool, y::True) = x
Base.:(&)(x::Bool, y::False) = y
Base.:(&)(x::True, y::Bool) = y
Base.:(&)(x::False, y::Bool) = x

Base.:(|)(::StaticInteger{X}, ::StaticInteger{Y}) where {X, Y} = static(|(X, Y))
Base.:(|)(::StaticInteger{X}, y::Union{Integer, Missing}) where {X} = X | y
Base.:(|)(x::Union{Integer, Missing}, ::StaticInteger{Y}) where {Y} = x | Y
Base.:(|)(x::Bool, y::True) = y
Base.:(|)(x::Bool, y::False) = x
Base.:(|)(x::True, y::Bool) = x
Base.:(|)(x::False, y::Bool) = y

Base.xor(::StaticInteger{X}, ::StaticInteger{Y}) where {X, Y} = static(xor(X, Y))
Base.xor(::StaticInteger{X}, y::Union{Integer, Missing}) where {X} = xor(X, y)
Base.xor(x::Union{Integer, Missing}, ::StaticInteger{Y}) where {Y} = xor(x, Y)

Base.:(!)(::True) = False()
Base.:(!)(::False) = True()

Base.all(::Tuple{True, Vararg{True}}) = true
Base.all(::Tuple{Union{True, False}, Vararg{Union{True, False}}}) = false

Base.any(::Tuple{False, Vararg{False}}) = false
Base.any(::Tuple{Union{True, False}, Vararg{Union{True, False}}}) = true

Base.real(@nospecialize(x::StaticNumber)) = x
Base.real(@nospecialize(T::Type{<:StaticNumber})) = eltype(T)
Base.imag(@nospecialize(x::StaticNumber)) = zero(x)

"""
    field_type(::Type{T}, f)

Functionally equivalent to `fieldtype(T, f)` except `f` may be a static type.
"""
@inline field_type(T::Type, f::Union{Int, Symbol}) = fieldtype(T, f)
@inline field_type(::Type{T}, ::StaticInt{N}) where {T, N} = fieldtype(T, N)
@inline field_type(::Type{T}, ::StaticSymbol{S}) where {T, S} = fieldtype(T, S)

@inline Base.exponent(::StaticNumber{M}) where {M} = static(exponent(M))

Base.:(^)(::StaticFloat64{x}, y::Float64) where {x} = exp2(log2(x) * y)

Base.:(+)(x::True) = One()
Base.:(+)(x::False) = Zero()

# from `^(x::Bool, y::Bool) = x | !y`
Base.:(^)(x::StaticBool, y::False) = True()
Base.:(^)(x::StaticBool, y::True) = x
Base.:(^)(x::Integer, y::False) = one(x)
Base.:(^)(x::Integer, y::True) = x
Base.:(^)(x::BigInt, y::False) = one(x)
Base.:(^)(x::BigInt, y::True) = x

@inline function Base.ntuple(f::F, ::StaticInt{N}) where {F, N}
    (N >= 0) || throw(ArgumentError(string("tuple length should be ≥ 0, got ", N)))
    if @generated
        quote
            Base.Cartesian.@ntuple $N i->f(i)
        end
    else
        Tuple(f(i) for i in 1:N)
    end
end

@inline function invariant_permutation(@nospecialize(x::Tuple), @nospecialize(y::Tuple))
    if y === x === ntuple(static, StaticInt(nfields(x)))
        return True()
    else
        return False()
    end
end

@inline nstatic(::Val{N}) where {N} = ntuple(StaticInt, Val(N))

permute(@nospecialize(x::Tuple), @nospecialize(perm::Val)) = permute(x, static(perm))
@inline function permute(@nospecialize(x::Tuple), @nospecialize(perm::Tuple))
    if invariant_permutation(nstatic(Val(nfields(x))), perm) === False()
        return eachop(getindex, perm, x)
    else
        return x
    end
end

"""
    Static.eachop(op, args...; iterator::Tuple{Vararg{StaticInt}})::Tuple

Produces a tuple of `(op(args..., iterator[1]), op(args..., iterator[2]),...)`.
"""
@inline function eachop(op::F, itr::Tuple{T, Vararg{Any}}, args::Vararg{Any}) where {F, T}
    (op(args..., first(itr)), eachop(op, Base.tail(itr), args...)...)
end
eachop(::F, ::Tuple{}, args::Vararg{Any}) where {F} = ()

"""
    Static.eachop_tuple(op, arg, args...; iterator::Tuple{Vararg{StaticInt}})::Type{Tuple}

Produces a tuple type of `Tuple{op(arg, args..., iterator[1]), op(arg, args..., iterator[2]),...}`.
Note that if one of the arguments passed to `op` is a `Tuple` type then it should be the first argument
instead of one of the trailing arguments, ensuring type inference of each element of the tuple.
"""
eachop_tuple(op, itr, arg, args...) = _eachop_tuple(op, itr, arg, args)
@generated function _eachop_tuple(op, ::I, arg, args::A) where {A, I}
    t = Expr(:curly, Tuple)
    narg = length(A.parameters)
    for p in I.parameters
        call_expr = Expr(:call, :op, :arg)
        if narg > 0
            for i in 1:narg
                push!(call_expr.args, :(getfield(args, $i)))
            end
        end
        push!(call_expr.args, :(StaticInt{$(p.parameters[1])}()))
        push!(t.args, call_expr)
    end
    Expr(:block, Expr(:meta, :inline), t)
end

#=
    find_first_eq(x, collection::Tuple)

Finds the position in the tuple `collection` that is exactly equal (i.e. `===`) to `x`.
If `x` and `collection` are static (`is_static`) and `x` is in `collection` then the return
value is a `StaticInt`.
=#
@generated function find_first_eq(x::X, itr::I) where {X, N, I <: Tuple{Vararg{Any, N}}}
    # we avoid incidental code gen when evaluated a tuple of known values by iterating
    #  through `I.parameters` instead of `known(I)`.
    index = known(X) === nothing ? nothing : findfirst(==(X), I.parameters)
    if index === nothing
        :(Base.Cartesian.@nif $(N + 1) d->(dynamic(x) == dynamic(getfield(itr, d))) d->(d) d->(nothing))
    else
        :($(static(index)))
    end
end

function Base.invperm(p::Tuple{StaticInt, Vararg{StaticInt, N}}) where {N}
    map(Base.Fix2(find_first_eq, p), ntuple(static, StaticInt(N + 1)))
end

"""
  reduce_tup(f::F, inds::Tuple{Vararg{Any,N}}) where {F,N}

An optimized `reduce` for tuples. `Base.reduce`'s `afoldl` will often not inline.
Additionally, `reduce_tup` attempts to order the reduction in an optimal manner.

```julia
julia> using StaticArrays, Static, BenchmarkTools

julia> rsum(v::SVector) = Static.reduce_tup(+, v.data)
rsum (generic function with 2 methods)

julia> for n ∈ 2:16
           @show n
           v = @SVector rand(n)
           s1 = @btime  sum(\$(Ref(v))[])
           s2 = @btime rsum(\$(Ref(v))[])
       end
n = 2
  0.863 ns (0 allocations: 0 bytes)
  0.863 ns (0 allocations: 0 bytes)
n = 3
  0.862 ns (0 allocations: 0 bytes)
  0.863 ns (0 allocations: 0 bytes)
n = 4
  0.862 ns (0 allocations: 0 bytes)
  0.862 ns (0 allocations: 0 bytes)
n = 5
  1.074 ns (0 allocations: 0 bytes)
  0.864 ns (0 allocations: 0 bytes)
n = 6
  0.864 ns (0 allocations: 0 bytes)
  0.862 ns (0 allocations: 0 bytes)
n = 7
  1.075 ns (0 allocations: 0 bytes)
  0.864 ns (0 allocations: 0 bytes)
n = 8
  1.077 ns (0 allocations: 0 bytes)
  0.865 ns (0 allocations: 0 bytes)
n = 9
  1.081 ns (0 allocations: 0 bytes)
  0.865 ns (0 allocations: 0 bytes)
n = 10
  1.195 ns (0 allocations: 0 bytes)
  0.867 ns (0 allocations: 0 bytes)
n = 11
  1.357 ns (0 allocations: 0 bytes)
  1.400 ns (0 allocations: 0 bytes)
n = 12
  1.543 ns (0 allocations: 0 bytes)
  1.074 ns (0 allocations: 0 bytes)
n = 13
  1.702 ns (0 allocations: 0 bytes)
  1.077 ns (0 allocations: 0 bytes)
n = 14
  1.913 ns (0 allocations: 0 bytes)
  0.867 ns (0 allocations: 0 bytes)
n = 15
  2.076 ns (0 allocations: 0 bytes)
  1.077 ns (0 allocations: 0 bytes)
n = 16
  2.273 ns (0 allocations: 0 bytes)
  1.078 ns (0 allocations: 0 bytes)
```

More importantly, `reduce_tup(_pick_range, inds)` often performs better than `reduce(_pick_range, inds)`.
```julia
julia> using ArrayInterface, BenchmarkTools, Static

julia> inds = (Base.OneTo(100), 1:100, 1:static(100))
(Base.OneTo(100), 1:100, 1:static(100))

julia> @btime reduce(ArrayInterface._pick_range, \$(Ref(inds))[])
  6.405 ns (0 allocations: 0 bytes)
Base.Slice(static(1):static(100))

julia> @btime Static.reduce_tup(ArrayInterface._pick_range, \$(Ref(inds))[])
  2.570 ns (0 allocations: 0 bytes)
Base.Slice(static(1):static(100))

julia> inds = (Base.OneTo(100), 1:100, 1:UInt(100))
(Base.OneTo(100), 1:100, 0x0000000000000001:0x0000000000000064)

julia> @btime reduce(ArrayInterface._pick_range, \$(Ref(inds))[])
  6.411 ns (0 allocations: 0 bytes)
Base.Slice(static(1):100)

julia> @btime Static.reduce_tup(ArrayInterface._pick_range, \$(Ref(inds))[])
  2.592 ns (0 allocations: 0 bytes)
Base.Slice(static(1):100)

julia> inds = (Base.OneTo(100), 1:100, 1:UInt(100), Int32(1):Int32(100))
(Base.OneTo(100), 1:100, 0x0000000000000001:0x0000000000000064, 1:100)

julia> @btime reduce(ArrayInterface._pick_range, \$(Ref(inds))[])
  9.048 ns (0 allocations: 0 bytes)
Base.Slice(static(1):100)

julia> @btime Static.reduce_tup(ArrayInterface._pick_range, \$(Ref(inds))[])
  2.569 ns (0 allocations: 0 bytes)
Base.Slice(static(1):100)
```
"""
@generated function reduce_tup(f::F, inds::Tuple{Vararg{Any, N}}) where {F, N}
    q = Expr(:block, Expr(:meta, :inline, :propagate_inbounds))
    if N == 1
        push!(q.args, :(inds[1]))
        return q
    end
    syms = Vector{Symbol}(undef, N)
    i = 0
    for n in 1:N
        syms[n] = iₙ = Symbol(:i_, (i += 1))
        push!(q.args, Expr(:(=), iₙ, Expr(:ref, :inds, n)))
    end
    W = 1 << (8sizeof(N) - 2 - leading_zeros(N))
    while W > 0
        _N = length(syms)
        for _ in (2W):W:_N
            for w in 1:W
                new_sym = Symbol(:i_, (i += 1))
                push!(q.args, Expr(:(=), new_sym, Expr(:call, :f, syms[w], syms[w + W])))
                syms[w] = new_sym
            end
            deleteat!(syms, (1 + W):(2W))
        end
        W >>>= 1
    end
    q
end

# This method assumes that `f` uetrieves compile time information and `g` is the fall back
# for the corresponding dynamic method. If the `f(x)` doesn't return `nothing` that means
# the value is known and compile time and returns `static(f(x))`.
@inline function maybe_static(f::F, g::G, x) where {F, G}
    L = f(x)
    if L === nothing
        return g(x)
    else
        return static(L)
    end
end

"""
    Static.eq(x, y)::Union{Bool, True, False}

Equivalent to `==` but if `x` and `y` are static the return value is a `StaticBool.
"""
eq(x::X, y::Y) where {X, Y} = ifelse(is_static(X) & is_static(Y), static, identity)(x == y)

"""
    Static.eq(x)::Base.Fix2{typeof(Static.eq}}

Create a function that compares `x` to other values using `Static.eq` (i.e. a
function equivalent to `y -> Static.eq(y, x)`).
"""
eq(x) = Base.Fix2(eq, x)

"""
    Static.ne(x, y)::Union{Bool, True, False}

Equivalent to `!=` but if `x` and `y` are static the return value is a `StaticBool.
"""
ne(x::X, y::Y) where {X, Y} = !eq(x, y)

"""
    Static.ne(x)::Base.Fix2{typeof(Static.ne}}

Create a function that compares `x` to other values using `Static.ne` (i.e. a
function equivalent to `y -> Static.ne(y, x))`.
"""
ne(x) = Base.Fix2(ne, x)

"""
    Static.ne(x, y)::Union{Bool, True, False}

Equivalent to `>` but if `x` and `y` are static the return value is a `StaticBool.
"""
gt(x::X, y::Y) where {X, Y} = ifelse(is_static(X) & is_static(Y), static, identity)(x > y)

"""
    Static.gt(x)::Base.Fix2{typeof(Static.gt}}

Create a function that compares `x` to other values using `Static.gt` (i.e. a
function equivalent to `y -> Static.gt(y, x))`.
"""
gt(x) = Base.Fix2(gt, x)

"""
    Static.ge(x, y)::Union{Bool, True, False}

Equivalent to `>=` but if `x` and `y` are static the return value is a `StaticBool.
"""
ge(x::X, y::Y) where {X, Y} = ifelse(is_static(X) & is_static(Y), static, identity)(x >= y)

"""
    Static.ge(x)::Base.Fix2{typeof(Static.ge}}

Create a function that compares `x` to other values using `Static.ge` (i.e. a
function equivalent to `y -> Static.ge(y, x)`).
"""
ge(x) = Base.Fix2(ge, x)

"""
    Static.le(x, y)::Union{Bool, True, False}

Equivalent to `<=` but if `x` and `y` are static the return value is a `StaticBool.
"""
le(x::X, y::Y) where {X, Y} = ifelse(is_static(X) & is_static(Y), static, identity)(x <= y)

"""
    Static.le(x)::Base.Fix2{typeof(Static.le}}

Create a function that compares `x` to other values using `Static.le` (i.e. a
function equivalent to `y -> Static.le(y, x)`).
"""
le(x) = Base.Fix2(le, x)

"""
    Static.lt(x, y)::Union{Bool, True, False}

Equivalent to `<` but if `x` and `y` are static the return value is a `StaticBool.`
"""
lt(x::X, y::Y) where {X, Y} = ifelse(is_static(X) & is_static(Y), static, identity)(x < y)

"""
    Static.lt(x)::Base.Fix2{typeof(Static.lt}}

Create a function that compares `x` to other values using `Static.lt` (i.e. a
function equivalent to y -> Static.lt(y, x)).
"""
lt(x) = Base.Fix2(lt, x)

"""
    Static.mul(x)::Base.Fix2{typeof(*)}

Create a function that multiplies `x` with other values (i.e. a function
equivalent to `y -> y * x`).
"""
mul(x) = Base.Fix2(*, x)

"""
    Static.add(x) -> Base.Fix2(+, x)
    Static.add(x, y)

Create a function that adds `x` to other values (i.e. a function equivalent to
`y -> y + x`).
"""
add(x) = Base.Fix2(+, x)

const Mul{X} = Base.Fix2{typeof(*), X}
const Add{X} = Base.Fix2{typeof(+), X}

Base.:∘(::Add{StaticInt{X}}, ::Add{StaticInt{Y}}) where {X, Y} = Base.Fix2(+, static(X + Y))
Base.:∘(x::Mul{Int}, ::Add{StaticInt{0}}) = x
Base.:∘(x::Mul{StaticInt{X}}, ::Add{StaticInt{0}}) where {X} = x
Base.:∘(x::Mul{StaticInt{0}}, ::Add{StaticInt{0}}) = x
Base.:∘(x::Mul{StaticInt{1}}, ::Add{StaticInt{0}}) = x
Base.:∘(x::Mul{StaticInt{0}}, y::Add{StaticInt{Y}}) where {Y} = x
Base.:∘(::Mul{StaticInt{1}}, y::Add{StaticInt{Y}}) where {Y} = y
Base.:∘(x::Mul{StaticInt{0}}, y::Add{Int}) = x
Base.:∘(::Mul{StaticInt{1}}, y::Add{Int}) = y
Base.:∘(::Mul{StaticInt{X}}, ::Mul{StaticInt{Y}}) where {X, Y} = Base.Fix2(*, static(X * Y))
Base.:∘(x::Mul{StaticInt{0}}, y::Mul{Int}) = x
Base.:∘(::Mul{Int}, y::Mul{StaticInt{0}}) = y
Base.:∘(::Mul{StaticInt{1}}, y::Mul{Int}) = y
Base.:∘(x::Mul{Int}, ::Mul{StaticInt{1}}) = x

# length
Base.length(@nospecialize(x::NDIndex))::Int = length(Tuple(x))
Base.length(::Type{<:NDIndex{N}}) where {N} = N

# indexing
Base.@propagate_inbounds function Base.getindex(x::NDIndex{N, T}, i::Int)::Int where {N, T}
    return Int(getfield(Tuple(x), i))
end
Base.@propagate_inbounds function Base.getindex(x::NDIndex{N, T},
        i::StaticInt{I}) where {N, T, I}
    return getfield(Tuple(x), I)
end

# Base.get(A::AbstractArray, I::CartesianIndex, default) = get(A, I.I, default)
# eltype(::Type{T}) where {T<:CartesianIndex} = eltype(fieldtype(T, :I))

Base.setindex(x::NDIndex, i, j) = NDIndex(Base.setindex(Tuple(x), i, j))

# equality
Base.:(==)(@nospecialize(x::NDIndex), @nospecialize(y::NDIndex)) = ==(Tuple(x), Tuple(y))

# zeros and ones
Base.zero(@nospecialize(x::NDIndex)) = zero(typeof(x))
function Base.zero(@nospecialize(T::Type{<:NDIndex}))
    NDIndex(ntuple(_ -> static(0), Val(length(T))))
end
Base.oneunit(@nospecialize(x::NDIndex)) = oneunit(typeof(x))
function Base.oneunit(@nospecialize(T::Type{<:NDIndex}))
    NDIndex(ntuple(_ -> static(1), Val(length(T))))
end

@inline function Base.IteratorsMD.split(i::NDIndex, V::Val)
    i, j = Base.IteratorsMD.split(Tuple(i), V)
    return NDIndex(i), NDIndex(j)
end

# arithmetic, min/max
@inline Base.:(-)(@nospecialize(i::NDIndex)) = NDIndex(map(-, Tuple(i)))
@inline function Base.:(+)(@nospecialize(i1::NDIndex), @nospecialize(i2::NDIndex))
    NDIndex(map(+, Tuple(i1), Tuple(i2)))
end
@inline function Base.:(-)(@nospecialize(i1::NDIndex), @nospecialize(i2::NDIndex))
    NDIndex(map(-, Tuple(i1), Tuple(i2)))
end
@inline function Base.min(@nospecialize(i1::NDIndex), @nospecialize(i2::NDIndex))
    NDIndex(map(min, Tuple(i1), Tuple(i2)))
end
@inline function Base.max(@nospecialize(i1::NDIndex), @nospecialize(i2::NDIndex))
    NDIndex(map(max, Tuple(i1), Tuple(i2)))
end
@inline function Base.:(*)(a::Integer, @nospecialize(i::NDIndex))
    NDIndex(map(x -> a * x, Tuple(i)))
end
@inline Base.:(*)(@nospecialize(i::NDIndex), a::Integer) = *(a, i)

Base.CartesianIndex(@nospecialize(x::NDIndex)) = dynamic(x)

# comparison
@inline function Base.isless(@nospecialize(x::NDIndex), @nospecialize(y::NDIndex))
    Bool(_isless(static(0), Tuple(x), Tuple(y)))
end

function lt(@nospecialize(x::NDIndex), @nospecialize(y::NDIndex))
    _isless(static(0), Tuple(x), Tuple(y))
end

_final_isless(c::Int) = c === 1
_final_isless(::StaticInt{N}) where {N} = static(false)
_final_isless(::StaticInt{1}) = static(true)
_isless(c::C, x::Tuple{}, y::Tuple{}) where {C} = _final_isless(c)
function _isless(c::C, x::Tuple, y::Tuple) where {C}
    _isless(icmp(c, x, y), Base.front(x), Base.front(y))
end
icmp(::StaticInt{0}, x::Tuple, y::Tuple) = icmp(last(x), last(y))
icmp(::StaticInt{N}, x::Tuple, y::Tuple) where {N} = static(N)
function icmp(cmp::Int, x::Tuple, y::Tuple)
    if cmp === 0
        return icmp(Int(last(x)), Int(last(y)))
    else
        return cmp
    end
end
icmp(a, b) = _icmp(lt(a, b), a, b)
_icmp(x::StaticBool, a, b) = ifelse(x, static(1), __icmp(eq(a, b)))
_icmp(x::Bool, a, b) = ifelse(x, 1, __icmp(a == b))
__icmp(x::StaticBool) = ifelse(x, static(0), static(-1))
__icmp(x::Bool) = ifelse(x, 0, -1)

#  Necessary for compatibility with Base
# In simple cases, we know that we don't need to use axes(A). Optimize those
# until Julia gets smart enough to elide the call on its own:
@inline function Base.to_indices(A, inds, I::Tuple{NDIndex, Vararg{Any}})
    to_indices(A, inds, (Tuple(I[1])..., Base.tail(I)...))
end
# But for arrays of CartesianIndex, we just skip the appropriate number of inds
@inline function Base.to_indices(A, inds,
        I::Tuple{AbstractArray{NDIndex{N, J}}, Vararg{Any}}) where {
        N,
        J
}
    _, indstail = Base.IteratorsMD.split(inds, Val(N))
    return (Base.to_index(A, I[1]), to_indices(A, indstail, Base.tail(I))...)
end

function Base.show(@nospecialize(io::IO), @nospecialize(x::Union{
        StaticNumber, StaticSymbol, NDIndex}))
    show(io, MIME"text/plain"(), x)
end
function Base.show(
        @nospecialize(io::IO), ::MIME"text/plain", @nospecialize(x::Union{
            StaticNumber, StaticSymbol}))
    print(io, "static(" * repr(known(typeof(x))) * ")")
    nothing
end
function Base.show(@nospecialize(io::IO), m::MIME"text/plain", @nospecialize(x::NDIndex))
    print(io, "NDIndex")
    show(io, m, Tuple(x))
    nothing
end

include("ranges.jl")

end
