module StructUtils

using Dates, UUIDs

export @noarg, @defaults, @tags, @kwarg, @nonstruct, Selectors

"""
    StructUtils.StructStyle

Abstract type that all concrete struct styles must subtype.
Custom struct styles allow fine-grained control over various
StructUtils.jl interface methods like `fieldtags`, `fielddefaults`,
`lift`, `lower`, etc.
"""
abstract type StructStyle end

"""
    StructUtils.DefaultStyle

Default struct style that all StructUtils.jl interface methods
are defined for by default.
"""
struct DefaultStyle <: StructStyle end

include("macros.jl")

"""
  StructUtils.dictlike(x) -> Bool
  StructUtils.dictlike(::StructStyle, x) -> Bool
  StructUtils.dictlike(::StructStyle, ::Type{T}) -> Bool

Returns `true` if `x` or type `T` is dictionary-like, `false` otherwise.
When `StructUtils.make(T, source)` is called, if `dictlike(T)` is `true`,
an instance will be `initialize`d, and then `addkeyval!`ed for each
key-value pair in `source`.
"""
function dictlike end

dictlike(st::StructStyle, x) = dictlike(st, typeof(x))
dictlike(::StructStyle, T::Type) = dictlike(T)
dictlike(::Type{<:AbstractDict}) = true
dictlike(::Type{<:AbstractVector{<:Pair}}) = true
dictlike(@nospecialize(T)) = false

"""
  StructUtils.noarg(x) -> Bool
  StructUtils.noarg(::StructStyle, x) -> Bool
  StructUtils.noarg(::StructStyle, ::Type{T}) -> Bool

Signals that `x` or type `T` is a mutable type that can be constructed by calling an empty
constructor, like `t = T()`. Automatically overloaded when structs use the
`@noarg` macro in their struct definition. The default value is `false` unless
explicitly overloaded.
"""
function noarg end

noarg(st::StructStyle, x) = noarg(st, typeof(x))
noarg(::StructStyle, T::Type) = noarg(T)
noarg(@nospecialize(T)) = false

"""
  StructUtils.kwarg(x) -> Bool
  StructUtils.kwarg(::StructStyle, x) -> Bool
  StructUtils.kwarg(::StructStyle, ::Type{T}) -> Bool

Signals that `x` or type `T` can be constructed by passing struct fields as keyword arguments
to the constructor, like `t = T(field1=a, field2=b, ...)`. Automatically overloaded
when structs use the `StructUtils.@kwarg` macro in their struct definition. The default value
is `false` unless explicitly overloaded.

Note that `StructUtils.@kwarg` is a separate implementation of `Base.@kwdef`, yet should
be a drop-in replacement for it.
"""
function kwarg end

kwarg(st::StructStyle, x) = kwarg(st, typeof(x))
kwarg(::StructStyle, T::Type) = kwarg(T)
kwarg(@nospecialize(T)) = false

"""
    StructUtils.fieldtagkey(::StructStyle) -> Symbol

Field tags defined on struct fields can be grouped by keys that are associated with
a particular struct style. This function returns the key that should be used to
retrieve field tags for a given struct style. By default, this function returns
`nothing`. An example overload might look like:

```julia
struct MySQLStyle <: StructStyle end

StructUtils.fieldtagkey(::MySQLStyle) = :mysql

@tags struct Foo
    a::Int &(mysql=(name="foo_a",),)
    b::String
end
```

In this example, when `StructUtils.make` is called on `Foo` with the `MySQLStyle` style,
only `(name="foo_a",)` will be retrieved from the field tags for `a` because the
`mysql` key is associated with the `MySQLStyle` struct style. In other words, fieldtag keys
allow custom struct styles to "namespace" field tags so structs can overload specific tags
in multiple ways for different namespaces, i.e. `a::Int &(mysql=(name="foo_a",), json=(name="json_a",))`.
"""
function fieldtagkey end

fieldtagkey(::StructStyle) = nothing

"""
    StructUtils.defaultstate(::StructStyle) -> Any

Returns the default state for a given struct style. This is used to initialize
the state of a struct when no state is provided. The default implementation
returns `nothing`.
"""
defaultstate(::StructStyle) = nothing

"""
    StructUtils.fieldtags(::StructStyle, ::Type{T}) -> NamedTuple
    StructUtils.fieldtags(::StructStyle, ::Type{T}, fieldname) -> NamedTuple

Returns a `NamedTuple` of field tags for the struct `T`. Field tags can be
added manually by overloading `fieldtags`, or included via convenient syntax
using the StructUtils.jl macros: `@tags`, `@noarg`, `@defaults`, or `@kwarg`.
Note this function returns the tags of *all* fields as a single NamedTuple.
"""
function fieldtags end

fieldtags(::StructStyle, T::Type)::NamedTuple{(),Tuple{}} = (;)

function fieldtags(st::StructStyle, T::Type, field)
    ft = fieldtags(st, T)
    isempty(ft) && return (;)
    fft = get(() -> (;), ft, field)
    ftk = fieldtagkey(st)
    return ftk === nothing ? fft : get(fft, ftk, fft)
end

"""
    StructUtils.fielddefaults(::StructStyle, ::Type{T}) -> NamedTuple
    StructUtils.fielddefault(::StructStyle, ::Type{T}, fieldname) -> NamedTuple

Returns a `NamedTuple` of field defaults for the struct `T`. Field defaults can be
added manually by overloading `fielddefaults`, or included via convenient syntax
using the StructUtils.jl macros: `@tags`, `@noarg`, `@defaults`, or `@kwarg`.
"""
function fielddefaults end

fielddefaults(::StructStyle, T::Type)::NamedTuple{(),Tuple{}} = (;)
fielddefault(st::StructStyle, T::Type, key) = get(fielddefaults(st, T), key, nothing)
@doc (@doc fielddefaults) fielddefault

"""
    StructUtils.initialize(::StructStyle, T, source) -> T

In `StructUtils.make`, this function is called to initialize a new instance of `T`,
when `T` is `dictlike`, `arraylike`, or `noarg`. The `source` is passed from the call to `make`,
and can be used for initialization if appropriate.
The default implementation of `initialize` is to call `T()` or `T(undef, 0)`
for `<:AbstractArray` types.
"""
function initialize end

initialize(st::StructStyle, T::Type, @nospecialize(source)) =
    arraylike(st, T) ? T(undef, 0) : T()

function initialize(st::StructStyle, ::Type{A}, source) where {A<:AbstractArray}
    if ndims(A) > 1
        dims = discover_dims(st, source)
        return A(undef, dims)
    else
        return A(undef, 0)
    end
end

initialize(::StructStyle, ::Type{T}, source) where {T<:AbstractSet} = T()

"""
    StructUtils.addkeyval!(d, k, v)

Add a key-value pair to a dictionary-like object `d`. This function is called
by `StructUtils.make` when `d` is `dictlike`. The default implementation is to
call `d[k] = v` for `AbstractDict`.
"""
function addkeyval! end

addkeyval!(d::AbstractDict, k, v) = d[k] = v
addkeyval!(d::AbstractVector, k, v) = push!(d, k => v)

_keytype(d) = keytype(d)
_keytype(::AbstractVector{Pair{A,B}}) where {A,B} = A
_valtype(d) = valtype(d)
_valtype(::AbstractVector{Pair{A,B}}) where {A,B} = B

"""
  StructUtils.arraylike(x) -> Bool
  StructUtils.arraylike(::StructStyle, x) -> Bool
  StructUtils.arraylike(::StructStyle, ::Type{T}) -> Bool

Returns `true` if `x` or type `T` is array-like, `false` otherwise. This function is
called by `StructUtils.make` to determine if `T` is array-like. The default
implementation returns `true` for `<:AbstractArray`, `<:AbstractSet`, `<:Tuple`,
`<:Base.Generator`, and `<:Core.SimpleVector` types, and `false` for `<:AbstractArray{T,0}`.

Once `initialize` is called, `StructUtils.make` will call `push!` to add values
to the array-like object.
"""
function arraylike end

arraylike(st::StructStyle, x) = arraylike(st, typeof(x))
arraylike(::StructStyle, T::Type) = arraylike(T)
arraylike(::Type{<:AbstractArray{T,0}}) where {T} = false
arraylike(::Type{<:AbstractArray}) = true
arraylike(::Type{<:AbstractSet}) = true
arraylike(::Type{<:Tuple}) = true
arraylike(::Type{<:Base.Generator}) = true
arraylike(::Type{<:Core.SimpleVector}) = true
arraylike(@nospecialize(T)) = false

"""
  StructUtils.structlike(x) -> Bool
  StructUtils.structlike(::StructStyle, x) -> Bool
  StructUtils.structlike(::StructStyle, ::Type{T}) -> Bool

Returns `true` if `x` or type `T` is struct-like, `false` otherwise. This function is
called by `StructUtils.make` to determine if `T` is struct-like. The default
implementation returns `true` for `isstructtype(T)` and `!Base.issingletontype(T)`.

`structlike` structs are expected to be able to be constructed by the default constructor
like `T(field1, field2, ...)`.

Due to how `StructUtils.make` works, `structlike` is often overloaded to `false` by "unit"/"atom" types
where fields should be considered private to the `make` process and should instead attempt to
`lift` the `source` object into the `unit` type.
"""
function structlike end

structlike(st::StructStyle, x) = structlike(st, typeof(x))
structlike(::StructStyle, T::Type) = structlike(T)
structlike(::Type{<:Function}) = false
structlike(::Type{<:Module}) = false
structlike(::Type{<:AbstractArray{T,0}}) where {T} = false
structlike(::Type{<:AbstractChar}) = false
structlike(::Type{<:AbstractString}) = false
structlike(::Type{Symbol}) = false
structlike(::Type{Regex}) = false
structlike(::Type{<:Dates.TimeType}) = false
structlike(::Type{Number}) = false
structlike(::Type{Nothing}) = false
structlike(::Type{Missing}) = false
structlike(::Type{UUID}) = false
structlike(::Type{VersionNumber}) = false
structlike(@nospecialize(T)) = isstructtype(T) && !Base.issingletontype(T)

"""
  StructUtils.nulllike(x) -> Bool
  StructUtils.nulllike(::StructStyle, x) -> Bool
  StructUtils.nulllike(::StructStyle, ::Type{T}) -> Bool

Returns `true` if `x` or type `T` is null-like, `false` otherwise. This function is
mainly used in the `make!` implementation to determine if a
`Union` type can be narrowed by excluding `nulllike` types like `Nothing` and `Missing`.
"""
function nulllike end

nulllike(st::StructStyle, x) = nulllike(st, typeof(x))
nulllike(::StructStyle, T::Type) = nulllike(T)
nulllike(@nospecialize(T)) = T === Missing || T === Nothing

"""
  StructUtils.lower(x) -> x
  StructUtils.lower(::StructStyle, x) -> x

Domain value transformation function. This function is called by
`StructUtils.applyeach` on each value in the `source` object before
calling the apply function. By default, `lower` is the identity function.
This allows a domain transformation of values according to the
style used.
"""
function lower end

lower(::StructStyle, x) = lower(x)
lower(x) = x

function lower(st::StructStyle, x, tags)
    # there are a few builtin tags supported
    if x isa Dates.TimeType && haskey(tags, :dateformat)
        return Dates.format(x, tags.dateformat)
    elseif haskey(tags, :lower)
        return tags.lower(x)
    else
        return lower(st, x)
    end
end

"""
  StructUtils.lowerkey(x) -> x
  StructUtils.lowerkey(style::StructUtils.StructStyle, x) -> x

Allows customizing how a value is lowered when used specifically as a key.
By default, calls [`StructUtils.lower`](@ref). Called from [`StructUtils.applyeach`](@ref)
on the key or index before passed to the key-value function.

### Example

```julia
struct Point
    x::Int; y::Int
end

# lower a Point as a single string value
StructUtils.lowerkey(::StructUtils.StructStyle, p::Point) = "\$(p.x)_\$(p.y)"

d = Dict(Point(1, 2) => 99)

StructUtils.make(Dict{String, Dict{String, Point}}, Dict(Point(1, 2) => Dict(Point(3, 4) => Point(5, 6))))
# Dict{String, Dict{String, Point}} with 1 entry:
#   "1_2" => Dict("3_4"=>Point(5, 6))
```

For loss-less round-tripping also provide a [`StructUtils.liftkey`](@ref) overload to "lift" the key back.
"""
lowerkey(::StructStyle, x) = lowerkey(x)
lowerkey(x) = x

"""
  StructUtils.lift(::Type{T}, x) -> T
  StructUtils.lift(::StructStyle, ::Type{T}, x) -> Tuple{T, Any}

Lifts a value `x` to a type `T`. This function is called by `StructUtils.make`
to lift unit/atom values to the appropriate type. The default implementation is
the identity function for most types, but it also includes special cases
for `Symbol`, `Char`, `UUID`, `VersionNumber`, `Regex`, and `TimeType` types to be
constructed from strings.
Allows transforming a "domain value" that may be some primitive representation
into a more complex Julia type.

The method with a `StructStyle` argument should return a tuple of the lifted value and any side-effect state
derived from lifting the value.
"""
function lift end

lift(::Type{Symbol}, x) = Symbol(x)
lift(::Type{String}, x::Symbol) = String(x)
lift(::Type{T}, x) where {T} = Base.issingletontype(T) ? T() : convert(T, x)
lift(::Type{>:Missing}, ::Nothing) = missing
lift(::Type{>:Nothing}, ::Nothing) = nothing
lift(::Type{>:Union{Missing,Nothing}}, ::Nothing) = nothing
lift(::Type{>:Union{Missing,Nothing}}, ::Missing) = missing
lift(::Type{Char}, x::AbstractString) = length(x) == 1 ? x[1] : throw(ArgumentError("expected single character, got $x"))
lift(::Type{UUID}, x::AbstractString) = UUID(x)
lift(::Type{VersionNumber}, x::AbstractString) = VersionNumber(x)
lift(::Type{Regex}, x::AbstractString) = Regex(x)
lift(::Type{T}, x::AbstractString) where {T<:Dates.TimeType} = T(x)

function lift(::Type{T}, x::AbstractString) where {T<:Enum}
    sym = Symbol(x)
    for (k, v) in Base.Enums.namemap(T)
        v === sym && return T(k)
    end
    throw(ArgumentError("invalid `$T` string value: \"$sym\""))
end

lift(st::StructStyle, ::Type{T}, x) where {T} = lift(T, x), defaultstate(st)

# bit of an odd case, but support 0-dimensional array lifting from scalar value
function lift(st::StructStyle, ::Type{A}, x) where {A<:AbstractArray{T,0}} where {T}
    m = A(undef)
    m[1] = lift(st, T, x)
    return m, defaultstate(st)
end

function lift(st::StructStyle, ::Type{T}, x, tags) where {T}
    if haskey(tags, :lift)
        return tags.lift(x), defaultstate(st)
    elseif T <: Dates.TimeType && haskey(tags, :dateformat)
        if tags.dateformat isa String
            return parse(T, x, Dates.DateFormat(tags.dateformat)), defaultstate(st)
        else
            return parse(T, x, tags.dateformat), defaultstate(st)
        end
    else
        return lift(st, T, x)
    end
end

"""
  StructUtils.liftkey(::Type{T}, x) -> x
  StructUtils.liftkey(style::StructStyle, ::Type{T}, x) -> x

Allows customizing how a key is lifted before being passed to [`addkeyval!`](@ref)
in `dictlike` construction.

By default, calls [`StructUtils.lift`](@ref).

### Example

```julia
struct Point
    x::Int; y::Int
end

# lift a Point from a string value
StructUtils.liftkey(::StructUtils.StructStyle, x::String) = Point(parse(Int, split(x, "_")[1]), parse(Int, split(x, "_")[2]))

d = Dict("1_2" => 99)
StructUtils.make(Dict{Point, Int}, Dict("1_2" => 99))
# Dict{Point, Int} with 1 entry:
#   Point(1, 2) => 99
```

For loss-less round-tripping also provide a [`StructUtils.lowerkey`](@ref) overload to "lower" the key.
"""
function liftkey end

liftkey(::StructStyle, ::Type{T}, x) where {T} = liftkey(T, x)
liftkey(::Type{T}, x) where {T} = lift(T, x)
liftkey(f, st::StructStyle, ::Type{T}, x) where {T} = f(liftkey(st, T, x))

"""
    StructUtils.applyeach(style, f, x) -> Union{StructUtils.EarlyReturn, Nothing}

A custom `foreach`-like function that operates specifically on `(key, val)` or `(ind, val)` pairs,
and supports short-circuiting (via `StructUtils.EarlyReturn`). It also supports a `StructStyle` argument
to allow for style-specific behavior for non-owned types.

For each key-value or index-value pair in `x`, call `f(k, v)`.
If `f` returns a `StructUtils.EarlyReturn` instance, `applyeach` should
return the `EarlyReturn` immediately and stop iterating (i.e. short-circuit).
Otherwise, the return value of `f` can be ignored and iteration continues.

Key types are generally expected to be Symbols, Strings, or Integers.

An example overload of `applyeach` for a generic iterable would be:

```julia
function StructUtils.applyeach(style::StructUtils.StructStyle, f, x::MyIterable)
    for (i, v) in enumerate(x)
        ret = f(StructUtils.lowerkey(style, i), StructUtils.lower(style, v))
        # if `f` returns EarlyReturn, return immediately
        ret isa StructUtils.EarlyReturn && return ret
    end
    return
end
```

Note that `applyeach` must include the `style` argument when overloading.

Also note that before applying `f`, the key or index is passed through `StructUtils.lowerkey(style, k)`,
and the value `v` is passed through `StructUtils.lower(style, v)`.

If a value is `#undef` or otherwise not defined, the `f` function should generally be called with `nothing` or skipped.
"""
function applyeach end

"""
    StructUtils.EarlyReturn{T}

A wrapper type that can be used in function arguments to `applyeach`
to short-circuit iteration and return a value from `applyeach`.

Example usage:

```julia
function find_needle_in_haystack(haystack, needle)
    ret = applyeach(haystack) do k, v
        k == needle && return StructUtils.EarlyReturn(v)
    end
    ret isa StructUtils.EarlyReturn && return ret.value
    throw(ArgumentError("needle not found in haystack")
end
````
"""
struct EarlyReturn{T}
    value::T
end

applyeach(f, x) = applyeach(DefaultStyle(), f, x)
applyeach(f, st::StructStyle, x) = applyeach(st, f, x)

function applyeach(st::StructStyle, f, x::AbstractArray)
    for i in eachindex(x)
        ret = if @inbounds(isassigned(x, i))
            f(lowerkey(st, i), lower(st, @inbounds(x[i])))
        else
            f(lowerkey(st, i), lower(st, nothing))
        end
        ret isa EarlyReturn && return ret
    end
    return defaultstate(st)
end

# special-case Pair vectors to act like Dicts
function applyeach(st::StructStyle, f, x::AbstractVector{Pair{K,V}}) where {K,V}
    for (k, v) in x
        ret = f(lowerkey(st, k), lower(st, v))
        ret isa EarlyReturn && return ret
    end
    return defaultstate(st)
end

# appropriate definition for iterables that
# can't have #undef values
function applyeach(st::StructStyle, f, x::Union{AbstractSet,Base.Generator,Core.SimpleVector})
    for (i, v) in enumerate(x)
        ret = f(lowerkey(st, i), lower(st, v))
        ret isa EarlyReturn && return ret
    end
    return defaultstate(st)
end

# generic definition for Tuple, NamedTuple, structs
function applyeach(st::StructStyle, f, x::T) where {T}
    if @generated
        N = fieldcount(T)
        ex = quote
            defs = fielddefaults(st, T)
        end
        for i = 1:N
            fname = Meta.quot(fieldname(T, i))
            push!(ex.args, quote
                ftags = fieldtags(st, T, $fname)
                if !haskey(ftags, :ignore) || !ftags.ignore
                    fname = get(ftags, :name, $fname)
                    ret = if isdefined(x, $i)
                        f(lowerkey(st, fname), lower(st, getfield(x, $i), ftags))
                    elseif haskey(defs, $fname)
                        # this branch should be really rare because we should
                        # have applied a field default in the struct constructor
                        f(lowerkey(st, fname), lower(st, defs[$fname], ftags))
                    else
                        f(lowerkey(st, fname), lower(st, nothing, ftags))
                    end
                    ret isa EarlyReturn && return ret
                end
            end)
        end
        push!(ex.args, :(return defaultstate(st)))
        return ex
    else
        defs = fielddefaults(st, T)
        for i = 1:fieldcount(T)
            fname = fieldname(T, i)
            ftags = fieldtags(st, T, fname)
            if !haskey(ftags, :ignore) || !ftags.ignore
                fname = get(ftags, :name, fname)
                ret = if isdefined(x, i)
                    f(lowerkey(st, fname), lower(st, getfield(x, i), ftags))
                elseif haskey(defs, fname)
                    f(lowerkey(st, fname), lower(st, defs[fname], ftags))
                else
                    f(lowerkey(st, fname), lower(st, nothing, ftags))
                end
                ret isa EarlyReturn && return ret
            end
        end
        return defaultstate(st)
    end
end

function applyeach(st::StructStyle, f, x::AbstractDict)
    for (k, v) in x
        ret = f(lowerkey(st, k), lower(st, v))
        ret isa EarlyReturn && return ret
    end
    return defaultstate(st)
end

@static if VERSION < v"1.10"
    function _isfieldatomic(t::Type, s::Int)
        t = Base.unwrap_unionall(t)
        # TODO: what to do for `Union`?
        isa(t, DataType) || return false # uncertain
        ismutabletype(t) || return false # immutable structs are never atomic
        1 <= s <= length(t.name.names) || return false # OOB reads are not atomic (they always throw)
        atomicfields = t.name.atomicfields
        atomicfields === C_NULL && return false
        s -= 1
        return unsafe_load(Ptr{UInt32}(atomicfields), 1 + s ÷ 32) & (1 << (s % 32)) != 0
    end
else
    const _isfieldatomic = Base.isfieldatomic
end

_setfield!(x, i, v) = setfield!(x, i, v, _isfieldatomic(typeof(x), i) ? :sequentially_consistent : :not_atomic)

keyeq(a::Symbol, b::String) = a === Symbol(b)
keyeq(a::String, b::Symbol) = Symbol(a) == b
keyeq(a, b::String) = string(a) == b
keyeq(a::AbstractString, b::String) = String(a) == b
keyeq(a, b) = isequal(a, b)
keyeq(x) = y -> keyeq(x, y)

macro _f(i)
    esc(:(x = f($i); x isa EarlyReturn && return x.value))
end

@generated function _foreach(f, ::Type{T}) where {T}
    # marked inline since this benefits from constant propagation of `n`
    n = fieldcount(T)
    ex = Expr(:block)
    push!(ex.args, :(Base.@_inline_meta))
    for i = 1:n
        push!(ex.args, :(@_f($i)))
    end
    return ex
    return
end

# helper closure that computes the length of an applyeach source
# note that it should be used sparingly/carefully since it consumes
# the source object and we generally want to do a single pass
if VERSION < v"1.10"
    mutable struct LengthClosure
        len::Int
    end
    (f::LengthClosure)(_, _) = f.len += 1
    function applylength(x)
        lc = LengthClosure(0)
        StructUtils.applyeach(lc, x)
        return lc.len
    end
else
    struct LengthClosure
        len::Ptr{Int}
    end

    (f::LengthClosure)(_, _) = unsafe_store!(f.len, unsafe_load(f.len) + 1)

    function applylength(x)
        ref = Ref(0)
        GC.@preserve ref begin
            lc = LengthClosure(Base.unsafe_convert(Ptr{Int}, ref))
            StructUtils.applyeach(lc, x)
            return unsafe_load(lc.len)
        end
    end
end # VERSION < v"1.10"

# recursively build up multidimensional array dimensions
# "[[1.0],[2.0]]" => (1, 2)
# "[[1.0,2.0]]" => (2, 1)
# "[[[1.0]],[[2.0]]]" => (1, 1, 2)
# "[[[1.0],[2.0]]]" => (1, 2, 1)
# "[[[1.0,2.0]]]" => (2, 1, 1)
# length of innermost array is 1st dim
function discover_dims(style, x)
    @assert arraylike(style, x)
    len = applylength(x)
    ret = (
        applyeach(x) do _, v
            return arraylike(style, v) ? EarlyReturn(discover_dims(style, v)) : EarlyReturn(())
        end
    )::EarlyReturn
    return (ret.value..., len)
end

struct MultiDimClosure{S,A}
    style::S
    arr::A
    dims::Vector{Int}
    cur_dim::Base.RefValue{Int}
end

function (f::MultiDimClosure{S,A})(i::Int, val) where {S,A}
    f.dims[f.cur_dim[]] = i
    if arraylike(f.style, val)
        f.cur_dim[] -= 1
        st = applyeach(f, f.style, val)
        f.cur_dim[] += 1
    else
        val, st = make(f.style, eltype(f.arr), val)
        setindex!(f.arr, val, f.dims...)
    end
    return st
end

struct MultiDimValFunc{S,A}
    style::S
    arr::A
    dims::Vector{Int}
end

(f::MultiDimValFunc{S,A})(x) where {S,A} = setindex!(f.arr, x, f.dims...)

"""
    StructUtils.make(T, source) -> T
    StructUtils.make(T, source, style) -> T
    StructUtils.make(style, T, source) -> Tuple{T, Any}
    StructUtils.make!(style, x::T, source)

Construct a struct of type `T` from `source` using the given `style`. The `source` can be any
type of object, and the `style` can be any `StructStyle` subtype (default `StructUtils.DefaultStyle()`).

`make` will use any knowledge of `noarg`, `arraylike`, or `dictlike` in order to
determine how to construct an instance of `T`. The fallback for structs is to rely on
the automatic "all argument" constructor that structs have defined by default (e.g. `T(fields...)`).

`make` calls `applyeach` on the `source` object, where the key-value pairs
from `source` will be used in constructing `T`.

The 3rd definition takes a `style` argument, allowing for overloads of non-owned types `T`.
The main difference between this and the 2nd definition is that the 3rd definition allows for
the `make` function to return a tuple of the constructed struct and any side-effect state
derived from making the struct.

The 4th definition allows passing in an already-constructed instance of `T` (`x`),
which must be mutable, and source key-value pairs will be applied as
to `x` as source keys are matched to struct field names.

For structs, `fieldtags` will be accounted for and certain tags can be used
to influence the construction of the struct.
"""
function make end

function make(::Type{T}, source, style::StructStyle=DefaultStyle()) where {T}
    x, _ = make(style, T, source)
    return x
end

if isdefined(Base, :delete) && applicable(Base.delete, (a=1,), :a)
    const _delete = Base.delete
else
    Base.@constprop :aggressive function delete(a::NamedTuple{an}, field::Symbol) where {an}
        names = Base.diff_names(an, (field,))
        NamedTuple{names}(a)
    end
    const _delete = delete
end

function make(style::StructStyle, T::Type, source, tags)
    if haskey(tags, :choosetype)
        return make(style, tags.choosetype(source), source, _delete(tags, :choosetype))
    end
    if T !== Any
        if T >: Missing && T !== Missing
            if nulllike(style, source)
                return make(style, Missing, source, tags)
            else
                return make(style, nonmissingtype(T), source, tags)
            end
        elseif T >: Nothing && T !== Nothing
            if nulllike(style, source)
                return make(style, Nothing, source, tags)
            else
                return make(style, Base.nonnothingtype(T), source, tags)
            end
        end
    end
    if T <: Tuple || dictlike(style, T) || arraylike(style, T) || noarg(style, T) || structlike(style, T)
        return make(style, T, source)
    else
        return lift(style, T, source, tags)
    end
end

function make(style::StructStyle, T::Type, source)
    # start with some hard-coded Union cases
    if T !== Any
        if T >: Missing && T !== Missing
            if nulllike(style, source)
                return make(style, Missing, source)
            else
                return make(style, nonmissingtype(T), source)
            end
        elseif T >: Nothing && T !== Nothing
            if nulllike(style, source)
                return make(style, Nothing, source)
            else
                return make(style, Base.nonnothingtype(T), source)
            end
        end
    end
    if T <: Tuple
        return maketuple(style, T, source)
    elseif dictlike(style, T)
        return makedict(style, T, source)
    elseif arraylike(style, T)
        return makearray(style, T, source)
    elseif noarg(style, T)
        return makenoarg(style, T, source)
    elseif structlike(style, T)
        return makestruct(style, T, source)
    else
        return lift(style, T, source)
    end
end

if VERSION < v"1.11"
    mem(n) = Vector{Any}(undef, n)
else
    mem(n) = Memory{Any}(undef, n)
end

macro _t(i)
    esc(:(@inbounds(vals[$i]::fieldtype(T, $i))))
end

@generated function _tuple(::Type{T}, vals) where {T}
    n = fieldcount(T)
    ex = Expr(:block)
    push!(ex.args, :(Base.@_inline_meta))
    push!(ex.args, Expr(:tuple, [:(@_t($i)) for i = 1:n]...))
    return ex
end

struct TupleClosure{T,A,S}
    vals::A
    style::S
    i::Ptr{Int}
end

function (f::TupleClosure{T,A,S})(k, v) where {T,A,S}
    st = _foreach(T) do i
        if typeof(k) == Int
            if k == i
                intval, intst = make(f.style, fieldtype(T, i), v)
                @inbounds f.vals[i] = intval
                return EarlyReturn(intst)
            end
        else
            j = unsafe_load(f.i)
            if j == i
                unsafe_store!(f.i, i + 1)
                elseval, elsest = make(f.style, fieldtype(T, i), v)
                @inbounds f.vals[i] = elseval
                return EarlyReturn(elsest)
            end
        end
    end
    return st === nothing ? defaultstate(f.style) : st
end

function maketuple(style, ::Type{T}, source) where {T}
    vals = mem(fieldcount(T))
    ref = Ref(1)
    GC.@preserve ref begin
        i = Base.unsafe_convert(Ptr{Int}, ref)
        st = applyeach(style, TupleClosure{T,typeof(vals),typeof(style)}(vals, style, i), source)
        return _tuple(T, vals), st
    end
end

struct DictClosure{T,S}
    dict::T
    style::S
end

function (f::DictClosure{T,S})(k, v) where {T,S}
    val, st = make(f.style, _valtype(f.dict), v)
    addkeyval!(f.dict, liftkey(f.style, _keytype(f.dict), k), val)
    return st
end

makedict(style, ::Type{T}, source) where {T} = makedict(style, initialize(style, T, source), source)

function makedict(style, dict::T, source) where {T}
    st = applyeach(style, DictClosure(dict, style), source)
    return dict, st
end

struct ArrayClosure{T,S}
    arr::T
    style::S
end

function (f::ArrayClosure{T,S})(_, v) where {T,S}
    val, st = make(f.style, eltype(f.arr), v)
    push!(f.arr, val)
    return st
end

makearray(style, ::Type{T}, source) where {T} = @inline makearray(style, initialize(style, T, source), source)

function makearray(style, x::T, source) where {T}
    if !(T <: AbstractSet) && ndims(T) > 1
        # multidimensional arrays
        n = ndims(T)
        st = applyeach(style, MultiDimClosure(style, x, ones(Int, n), Ref(n)), source)
        return x, st
    else
        st = applyeach(style, ArrayClosure(x, style), source)
        return x, st
    end
end

@generated function fieldnamestrings(::Type{T}) where {T}
    :($(Tuple(String(fieldname(T, i)) for i in 1:fieldcount(T))))
end

@generated function fieldnamesymbols(::Type{T}) where {T}
    :($(Tuple(fieldname(T, i) for i in 1:fieldcount(T))))
end

struct StructClosure{T,A,S,FS,FSS}
    vals::A # Memory{Any} for structs, T for mutable structs
    style::S
    fsyms::FS
    fstrs::FSS
end

StructClosure{T}(vals::A, style::S, fsyms::FS, fstrs::FSS) where {T,A,S,FS,FSS} = StructClosure{T,A,S,FS,FSS}(vals, style, fsyms, fstrs)

if VERSION < v"1.11"
    setval!(vals::Vector{Any}, x, i) = @inbounds vals[i] = x
else
    setval!(vals::Memory{Any}, x, i) = @inbounds vals[i] = x
end

setval!(vals::T, x, i) where {T} = _setfield!(vals, i, x)

function findfield(::Type{T}, k, v, f) where {T}
    st = _foreach(T) do i
        if typeof(k) == Symbol
            fn = f.fsyms[i]
            ftags = fieldtags(f.style, T, fn)
            field = get(ftags, :name, fn)
            if k == field
                symval, symst = make(f.style, fieldtype(T, i), v, ftags)
                setval!(f.vals, symval, i)
                return EarlyReturn(symst)
            end
        elseif typeof(k) == Int
            if k == i
                ftags = fieldtags(f.style, T, f.fsyms[i])
                intval, intst = make(f.style, fieldtype(T, i), v, ftags)
                setval!(f.vals, intval, i)
                return EarlyReturn(intst)
            end
        else
            fn = f.fsyms[i]
            fstr = f.fstrs[i]
            ftags = fieldtags(f.style, T, fn)
            field = get(ftags, :name, fstr)
            if keyeq(k, field)
                strval, strst = make(f.style, fieldtype(T, i), v, ftags)
                setval!(f.vals, strval, i)
                return EarlyReturn(strst)
            end
        end
    end
    return st === nothing ? defaultstate(f.style) : st
end

(f::StructClosure{T,A,S,FS,FSS})(k, v) where {T,A,S,FS,FSS} = findfield(T, k, v, f)

@inline makenoarg(style, ::Type{T}, source) where {T} = makenoarg(style, initialize(style, T, source), source)

function makenoarg(style, y::T, source) where {T}
    fsyms = fieldnamesymbols(T)
    fstrs = fieldnamestrings(T)
    st = applyeach(style, StructClosure{T}(y, style, fsyms, fstrs), source)
    return y, st
end

macro _v(i)
    esc(:(isassigned(vals, $i) ? @inbounds(vals[$i])::fieldtype(T, $i) : fielddefault(style, T, @inbounds(fsyms[$i]))::fieldtype(T, $i)))
end

@generated function _construct(::Type{T}, vals, style, fsyms) where {T}
    n = fieldcount(T)
    ex = Expr(:block)
    push!(ex.args, :(Base.@_inline_meta))
    push!(ex.args, Expr(:call, Any[:T, [:(@_v($i)) for i = 1:n]...]...))
    return ex
end

function makestruct(style, ::Type{T}, source) where {T}
    vals = mem(fieldcount(T))
    fsyms = fieldnamesymbols(T)
    fstrs = fieldnamestrings(T)
    st = applyeach(style, StructClosure{T}(vals, style, fsyms, fstrs), source)
    if T <: NamedTuple
        return T(_tuple(T, vals)), st
    else
        return _construct(T, vals, style, fsyms), st
    end
end

make!(x::T, source; style::StructStyle=DefaultStyle()) where {T} = make!(style, x, source)
make!(::Type{T}, source; style::StructStyle=DefaultStyle()) where {T} = make!(style, T, source)

function make!(style::StructStyle, x::T, source) where {T}
    if dictlike(style, x)
        _, st = makedict(style, x, source)
    elseif arraylike(style, x)
        _, st = makearray(style, x, source)
    elseif noarg(style, x)
        _, st = makenoarg(style, x, source)
    else
        throw(ArgumentError("Type `$T` does not support in-place `make!`"))
    end
    return st
end

function make!(style::StructStyle, ::Type{T}, source) where {T}
    x = initialize(style, T, source)
    make!(style, x, source)
    return x
end

@doc (@doc make) make!

"""
  StructUtils.reset!(x::T)

If `T` was defined with default values via `@defaults`, `@tags`, `@kwarg`, or `@noarg`,
`reset!` will reset the fields of `x` to their default values.
`T` must be a mutable struct type.
"""
function reset!(x::T; style::StructStyle=DefaultStyle()) where {T}
    if @generated
        N = fieldcount(T)
        ex = quote
            defs = fielddefaults(style, T)
        end
        for i = 1:N
            fname = Meta.quot(fieldname(T, i))
            push!(ex.args, quote
                if haskey(defs, $fname)
                    _setfield!(x, $i, defs[$fname])
                end
            end)
        end
        push!(ex.args, :(return x))
        return ex
    else
        defs = fielddefaults(style, T)
        for i = 1:fieldcount(T)
            fname = fieldname(T, i)
            if haskey(defs, fname)
                _setfield!(x, i, defs[fname])
            end
        end
        return x
    end
end

include("selectors.jl")

"""
    StructUtils.@choosetype T func
    StructUtils.@choosetype style T func

Convenience macro for defining a `StructUtils.make!` overload for an abstract type `T` where
`func` is a function that "chooses" a concrete type `S` at runtime. `func` can be one of two forms:
  * `source -> S`
  * `(source, tags) -> S)`

That is, it either takes just the `source` object that is passed to `make` and must choose a concrete
type `S`, or it can take both the `source` and a set of fieldtags that may be present for the field
of a type being "made".

The 2nd definition also takes a `style` argument, allowing for overloads of non-owned types `T`.

Example:

```julia
abstract type Vehicle end

struct Car <: Vehicle
    make::String
    model::String
    seatingCapacity::Int
    topSpeed::Float64
end

struct Truck <: Vehicle
    make::String
    model::String
    payloadCapacity::Float64
end

StructUtils.@choosetype Vehicle x -> x["type"] == "car" ? Car : x["type"] == "truck" ? Truck : throw(ArgumentError("Unknown vehicle type: \$(x["type"])"))

x = StructUtils.make(Vehicle, Dict("type" => "car", "make" => "Toyota", "model" => "Corolla", "seatingCapacity" => 4, "topSpeed" => 120.5))
@test x == Car("Toyota", "Corolla", 4, 120.5)
```
"""
macro choosetype(T, ex)
    esc(quote
        function StructUtils.make(st::StructUtils.StructStyle, ::Type{$T}, source, tags)
            func = $(ex)
            StructUtils.make(st, applicable(func, source, tags) ? func(source, tags) : func(source), source, tags)
        end
        function StructUtils.make(st::StructUtils.StructStyle, ::Type{$T}, source)
            func = $(ex)
            StructUtils.make(st, func(source), source)
        end
    end)
end

macro choosetype(style, T, ex)
    esc(quote
        function StructUtils.make(st::$(style), ::Type{$T}, source, tags)
            func = $(ex)
            StructUtils.make(st, applicable(func, source, tags) ? func(source, tags) : func(source), source, tags)
        end
        function StructUtils.make(st::$(style), ::Type{$T}, source)
            func = $(ex)
            StructUtils.make(st, func(source), source)
        end
    end)
end

end
