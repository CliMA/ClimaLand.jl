module Extents

struct Fix2kw{F,B,KW}
    f::F
    b::B
    kw::KW
end
Fix2kw(f, b; kw...) = Fix2kw{typeof(b),typeof(kw)}(b, kw)
(f2::Fix2kw)(a; kw...) = f2.f(a, f2.b; f2.kw..., kw...)

export Extent, extent, bounds

## DO NOT export anything else ##

const ORDER_DOC = """
The order of dimensions is ignored in all cases.
"""

"""
    Extent

    Extent(; kw...)
    Extent(bounds::NamedTuple)

A wrapper for a `NamedTuple` of tuples holding
the lower and upper bounds for each dimension of the object.

`keys(extent)` will return the dimension name Symbols,
in the order the dimensions are used in the object.

`values(extent)` will return a tuple of tuples: `(lowerbound, upperbound)` for each
dimension.

# Examples
```julia-repl
julia> ext = Extent(X = (1.0, 2.0), Y = (3.0, 4.0))
Extent(X = (1.0, 2.0), Y = (3.0, 4.0))

julia> keys(ext)
(:X, :Y)

julia> values(ext)
((1.0, 2.0), (3.0, 4.0))
```
"""
struct Extent{K,V}
    bounds::NamedTuple{K,V}
    function Extent{K,V}(bounds::NamedTuple{K,V}) where {K,V}
        bounds = map(b -> promote(b...), bounds)
        new{K,typeof(values(bounds))}(bounds)
    end
end
Extent(; kw...) = Extent(values(kw))
# Allow vectors or tuples
Extent{K}(vals) where {K} = Extent{K}(NamedTuple{K}(vals))
# Subset K2 to K1
Extent{K1}(vals::NamedTuple{K2,V}) where {K1,K2,V} = Extent(NamedTuple{K1}(vals))
Extent(vals::NamedTuple{K,V}) where {K,V} = Extent{K,V}(vals)

bounds(ext::Extent) = getfield(ext, :bounds)

@inline function Base.getproperty(ext::Extent, key::Symbol)
    haskey(bounds(ext), key) || _ext_no_key(key)
    getproperty(bounds(ext), key)
end
Base.propertynames(ext::Extent) = propertynames(getfield(ext, :bounds))

@inline Base.getindex(ext::Extent, keys::NTuple{<:Any,Symbol}) = Extent{keys}(bounds(ext))
@inline Base.getindex(ext::Extent, keys::AbstractVector{Symbol}) = ext[Tuple(keys)]
@inline function Base.getindex(ext::Extent, key::Symbol)
    haskey(bounds(ext), key) || _ext_no_key(key)
    getindex(bounds(ext), key)
end
Base.@propagate_inbounds function Base.getindex(ext::Extent, i::Int)
    haskey(bounds(ext), i) || _ext_no_key(i)
    getindex(bounds(ext), i)
end
Base.haskey(ext::Extent, x) = haskey(bounds(ext), x)
Base.keys(ext::Extent) = keys(bounds(ext))
Base.values(ext::Extent) = values(bounds(ext))
Base.length(ext::Extent) = length(bounds(ext))
Base.iterate(ext::Extent, args...) = iterate(bounds(ext), args...)
Base.map(f, ext::Extent) = Extent(map(f, bounds(ext)))
Base.NamedTuple(ext::Extent) = bounds(ext)

function Base.isapprox(a::Extent{K1}, b::Extent{K2}; kw...) where {K1,K2}
    _check_keys_match(a, b) || return false
    values_match = map(K1) do k
        bounds_a = a[k]
        bounds_b = b[k]
        if isnothing(bounds_a) && isnothing(bounds_b) 
            true
        else
            map(bounds_a, bounds_b) do val_a, val_b
                isapprox(val_a, val_b; kw...)
            end |> all
        end
    end
    return all(values_match)
end

function Base.:(==)(a::Extent{K1}, b::Extent{K2}) where {K1,K2}
    _check_keys_match(a, b) || return false
    values_match = map(K1) do k
        bounds_a = a[k]
        bounds_b = b[k]
        isnothing(bounds_a) && isnothing(bounds_b) || bounds_a == bounds_b
    end
    return all(values_match)
end

function Base.show(io::IO, extent::Extent)
    print(io, "Extent")
    show(io, bounds(extent))
end

"""
    extent(x)

Returns an [`Extent`](@ref), holding the bounds for each dimension of the object.
"""
function extent end

extent(extent) = nothing
extent(extent::Extent) = extent

"""
    union(ext1::Extent, ext2::Extent; strict=false)

Get the union of two extents, e.g. the combined extent of both objects
for all dimensions.

$ORDER_DOC
"""
function union(ext1::Extent, ext2::Extent; strict=false)
    _maybe_check_keys_match(ext1, ext2, strict) || return nothing
    keys = _shared_keys(ext1, ext2)
    if length(keys) == 0
        return nothing
    else
        values = map(keys) do k
            k = _unwrap(k)
            k_exts = (ext1[k], ext2[k])
            a = min(map(first, k_exts)...)
            b = max(map(last, k_exts)...)
            (a, b)
        end
        return Extent{map(_unwrap, keys)}(values)
    end
end
union(a::Extent, ::Nothing; strict=false) = strict ? nothing : a
union(::Nothing, b::Extent; strict=false) = strict ? nothing : b
union(::Nothing, ::Nothing; kw...) = nothing
union(a, b; kw...) = union(extent(a), extent(b))
union(a, b, c, args...; kw...) = union(union(a, b), c, args...)

"""
    intersection(ext1::Extent, ext2::Extent; strict=false)

Get the intersection of two extents as another `Extent`, e.g.
the area covered by the shared dimensions for both extents.

If there is no intersection for any shared dimension, `nothing` will be returned.

$ORDER_DOC
"""
function intersection(a::Extent, b::Extent; strict=false)
    _maybe_check_keys_match(a, b, strict) || return nothing
    intersects(a, b) || return nothing
    keys = _shared_keys(a, b)
    values = map(keys) do k
        # Get a symbol from `Val{:k}`
        k = _unwrap(k)
        # Acces the k symbol of `a` and `b`
        k_exts = (a[k], b[k])
        maxs = max(map(first, k_exts)...)
        mins = min(map(last, k_exts)...)
        (maxs, mins)
    end
    return Extent{map(_unwrap, keys)}(values)
end
intersection(a::Extent, b::Nothing; kw...) = nothing
intersection(a::Nothing, b::Extent; kw...) = nothing
intersection(a::Nothing, b::Nothing; kw...) = nothing
intersection(a, b; kw...) = intersection(extent(a), extent(b); kw...)
intersection(a, b, c, args...; kw...) = intersection(intersection(a, b), c, args...; kw...)

"""
    buffer(ext::Extent, buff::NamedTuple)

buffer `Extent` by corresponding name-pair values supplied in `buff` NamedTuple.

# Examples

```julia-repl
julia> ext = Extent(X = (1.0, 2.0), Y = (3.0, 4.0))
Extent(X = (1.0, 2.0), Y = (3.0, 4.0))

julia> ext_buffered = Extents.buffer(ext, (X=1, Y=3))
Extent(X = (0.0, 3.0), Y = (0.0, 7.0))
```
"""
function buffer(ext::Extent{K}, buff::NamedTuple) where {K}
    bounds = map(map(Val, K)) do k
        if haskey(buff, _unwrap(k))
            map(+, ext[_unwrap(k)], (-buff[_unwrap(k)], +buff[_unwrap(k)]))
        else
            ext[_unwrap(k)]
        end
    end
    Extent{K}(bounds)
end
buffer(ext::Nothing, buff) = nothing

"""
    grow(ext::Extent, x)
    grow(ext::Extent; kw...)

Grow the bounds of the extent by `x`, as a fraction of the current size of the extent.

If `x` is a `Tuple` the lower and upper bounds are grow by those amounts.

Keyword arguments or a `NamedTuple` also be passed, with the same or a subset 
of the keys as `ext`. This can hold `Real` or `Tuple{Real,Real}` values for 
each named dimension.

## Examples

```jldoctest
julia> ext = Extent(X = (1.0, 1.8), Y = (3.0, 5.0))
Extent(X = (1.0, 1.8), Y = (3.0, 5.0))

julia> Extents.grow(ext, 0.5) 
Extent(X = (0.6, 2.2), Y = (2.0, 6.0))

````
"""
function grow(ext::Extent, x::Union{Real,Tuple{<:Real,<:Real}})
    map(ext) do bs
        _grow(bs, x)
    end
end
function grow(ext::Extent{K1}, xs::NamedTuple{K2}) where {K1,K2}
    map(K1) do k
        if k in K2
            _grow(ext[k], xs[k])
        else
            ext[k]
        end
    end |> Extent{K1}
end
grow(ext::Extent; kw...) = grow(ext, NamedTuple(kw))

function _grow(bs, x::Real)
    amount = (bs[2] - bs[1]) * x
    (bs[1] - amount, bs[2] + amount)
end
function _grow(bs, x::Tuple{Real,Real})
    range = (bs[2] - bs[1]) 
    (bs[1] - x * grow[1], bs[2] + range * x[2])
end

# DE_9IM predicates

const STRICT_DOC = """
Dimensions that are not shared are ignored by default with `strict=false`.
When `strict=true`, any unshared dimensions cause the function to return `nothng`.
"""

const DE_9IM_DOC = """
Conforms to the DE-9IM spatial predicates standard
https://en.wikipedia.org/wiki/DE-9IM
"""

"""
    contains(a, b; strict=false)
    contains(b; strict=false)(a)

Extent `a` contains extent `b` if no points of `b` lie in the exterior of `a`, 
and at least one point of the interior of `b` lies in the interior of `a`.
If `b` has no interior points it is not contained in `a`.

Identical to [`within`](@ref) with argument order reversed.

$STRICT_DOC

If there are no common dimensions, `false` is returned.

$ORDER_DOC

$DE_9IM_DOC
"""
contains(a::Extent, b::Extent; strict=false) = _do_bounds(all, _contain, a, b, strict)

# Must contain interior points, not just boundary
_contain(a::Tuple, b::Tuple) = _cover(a, b) && _hasinterior(b)

"""
    within(a, b; strict=false)
    within(b; strict=false)(a)

Extent `a` is within extent `b` if no points of `a` lie in the exterior of `b`, 
and at least one point of the interior of `a` lies in the interior of `b`.
If `a` has no interior points it is not contained in `b`.

Identical to [`contains`](@ref) with argument order reversed.

$STRICT_DOC

If there are no common dimensions, `false` is returned.

$ORDER_DOC

$DE_9IM_DOC 
"""
within(a, b; kw...) = contains(b, a; kw...) # swapped order of `contains`

"""
    intersects(a, b; strict=false)
    intersects(b; strict=false)(a)

`a` intersects `b` if `a` and `b` have at least one point in common
(the inverse of [`disjoint`](@ref)).

Returns `true` if the extents of all common dimensions share some values,
including just the edge values of their range.

$STRICT_DOC

If there are no common dimensions with `strict=false`, `false` is returned.

$ORDER_DOC

$DE_9IM_DOC 
"""
intersects(a::Extent, b::Extent; strict=false) = _do_bounds(all, _intersect, a, b, strict)

_intersect((min_a, max_a)::Tuple, (min_b, max_b)::Tuple) = 
    (min_a <= min_b && max_a >= min_b) || (min_b <= min_a && max_b >= min_a)

"""
    disjoint(a, b; strict=false)
    disjoint(b; strict=false)(a)

Extents `a` and `b` are disjoint if they have no point in common
(the inverse of [`intersects`](@ref)).

Returns `false` if the extents of all common dimensions share some values,
including just the edge values of their range.

$STRICT_DOC

If there are no common dimensions when `strict=false`, `true` is returned.

$ORDER_DOC

$DE_9IM_DOC 
"""
disjoint(a, b; kw...) = !intersects(a, b; kw...)

"""
    touches(a, b; strict=false)
    touches(b; strict=false)(a)

Extents `a` and `b` have at least one point in common, but their interiors do not intersect. 

Returns `true` if the extents of any common dimensions share boundaries.

$STRICT_DOC

If there are no common dimensions with `strict=false`, `false` is returned.

$ORDER_DOC

$DE_9IM_DOC 
"""
function touches(a::Extent, b::Extent; strict=false)
    if intersects(a, b)
        # At least one bound must just touch
        _do_bounds(any, _touch, a, b, strict)
    else
        false
    end
end

_touch((min_a, max_a)::Tuple, (min_b, max_b)::Tuple) = (min_a == max_b || max_a == min_b)


"""
    covers(a, b; strict=false)
    covers(b; strict=false)(a)

At least one point of extent `b` lies in extent `a`, 
and no point of `b` lies in the exterior of `a`.

Every point of `b` is a point in the interior or boundary of `a`. 

Identical to [`coveredby`](@ref) with argument order reversed.

$STRICT_DOC

If there are no common dimensions with `strict=false`, `false` is returned.

$ORDER_DOC

$DE_9IM_DOC 
"""
covers(a::Extent, b::Extent; strict=false) = _do_bounds(all, _cover, a, b, strict)

_cover((min_a, max_a)::Tuple, (min_b, max_b)::Tuple) = (min_a <= min_b && max_a >= max_b)

"""
    coveredby(a, b; strict=false)
    coveredby(b; strict=false)(a)

At least one point of extent `a` lies in extent `b`, 
and no point of `a` lies in the exterior of `b`.

Every point of `a` is a point in the interior or boundary of `b`. 

Identical to [`covers`](@ref) with argument order reversed.

$STRICT_DOC

If there are no common dimensions with `strict=false`, `false` is returned.

$ORDER_DOC

$DE_9IM_DOC 
"""
coveredby(a, b; kw...) = covers(b, a; kw...) # swapped order of `covers`


"""
    overlaps(a, b; strict=false)
    overlaps(b; strict=false)(a)

Extent `a` overlaps extent `b`: they have some but not all points in common, 
they have the same dimension, and the intersection of the interiors
of the two geometries has the same dimension as the geometries themselves.

Returns `true` if the extents of common dimensions overlap.

$STRICT_DOC

If there are no common dimensions with `strict=false`, `false` is returned.

$ORDER_DOC

$DE_9IM_DOC 
"""
function overlaps(a::Extent, b::Extent; strict=false)
    if intersects(a, b; strict)
        # Bounds can't just touch, they must share interior points
        !_do_bounds(any, _touch, a, b, strict)
    else
        false
    end
end

"""
    equals(a, b; strict=false)
    equals(b; strict=false)(a)

Extents `a` and `b` are topologically equal: their interiors intersect and no 
part of the interior or boundary of one intersects the exterior of the other.

$STRICT_DOC

If there are no common dimensions with `strict=false`, `false` is returned.

$ORDER_DOC

$DE_9IM_DOC 
"""
equals(a::Extent, b::Extent; strict=false) = _do_bounds(all, _equal, a, b, strict)

_equal(a::Tuple, b::Tuple) = a == b

# Handle `nothing` bounds for all methods
for f in (:_intersect, :_cover, :_contain, :_touch, :_equal)
    @eval begin
        $f(::Nothing, ::Tuple) = nothing 
        $f(::Tuple, ::Nothing) = nothing 
        $f(::Nothing, ::Nothing) = nothing 
    end
end

# Handle objects with `extent` methods and `nothing` extents
for f in (:intersects, :covers, :contains, :touches, :equals, :overlaps)
    @eval begin
        $f(a, b; kw...) = $f(extent(a), extent(b); kw...)
        $f(a::Extent, b::Nothing; kw...) = false
        $f(a::Nothing, b::Extent; kw...) = false
        $f(a::Nothing, b::Nothing; kw...) = false
    end
end

for f in (:intersects, :disjoint, :covers, :coveredby, :contains, :within, :touches, :equals, :overlaps)
    @eval begin 
        $f(b; kw...) = Fix2kw($f, extent(b), kw)
    end
end


# Internal utils

_maybe_keys_match(ext1, ext2, strict) = !strict || _keys_match(ext1, ext2)

# Keys

_maybe_check_keys_match(a, b, strict) = !strict || _check_keys_match(a, b)

function _check_keys_match(::Extent{K1}, ::Extent{K2}) where {K1,K2}
    length(K1) == length(K2) || return false
    keys_match = map(K2) do k
        k in K1
    end |> all
end

# _shared_keys uses a static `Val{k}` instead of a `Symbol` to
# represent keys, because constant propagation fails through `reduce`
# meaning most of the time of `union` or `intersect` is doing the `Symbol` lookup.
# So we help the compiler out a little by doing manual constant propagation.
# We know K1 and K2 at compile time, and wrapping them in `Val{k}() maintains
# that through reduce. This makes union/intersect 15x faster, at ~10ns.
function _shared_keys(ext1::Extent{K1}, ext2::Extent{K2}) where {K1,K2}
    reduce(K1; init=()) do acc, k
        k in K2 ? (acc..., Val{k}()) : acc
    end
end

@noinline _ext_no_key(key) = throw(ErrorException("Extent has no field $key"))

_unwrap(::Val{X}) where {X} = X


# Bounds comparisons

# compare all bounds and reduce the result to a Bool or `nothing`
# Running `compare` and `boolreduce` should be the only runtime costs
function _do_bounds(boolreduce::Function, compare::Function, a::Extent, b::Extent, strict::Bool)
    _maybe_check_keys_match(a, b, strict) || return false
    keys = _shared_keys(a, b)
    if length(keys) == 0
        # There are no shared dimensions. Maybe this should return `nothing`?
        # But we need to handle it otherwise `all` returns `true` for empty tuples
        return false 
    else
        maybe_comparisons = map(keys) do k
            compare(a[_unwrap(k)], b[_unwrap(k)])
        end
        comparisons = _skipnothing(maybe_comparisons...)
        if length(comparisons) == 0
            return nothing
        else
            return boolreduce(comparisons)
        end
    end
end

_skipnothing(v1, vals...) = (v1, _skipnothing(vals...)...)
_skipnothing(::Nothing, vals...) = _skipnothing(vals...)
_skipnothing() = ()

_hasinterior(ex::Extent) = all(map(_hasinterior, bounds(ex)))
_hasinterior((min, max)::Tuple) = min != max

end