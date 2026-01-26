"""
    Selection syntax

Special "selection syntax" is provided that allows easy querying of objects/arrays that implement `StructUtils.applyeach` using a syntax similar to XPath or CSS selectors,
applied using common Julia syntax.

This syntax mainly uses various forms of `getindex` to select elements of an object or array.
Supported syntax includes:
  * `x["key"]` / `x.key` / `x[:key]` / `x[1]` - select the value associated for a key in object `x` (key can be a String, Symbol, or Integer for an array)
  * `x[:]` - select all values in object or array `x`, returned as a `Selectors.List`, which is a custom array type that supports the selection syntax
  * `x.key` - when `x` is a `List`, select the value for `key` in each element of the `List` (like a broadcasted `getindex`)
  * `x[~, key]` - recursively select all values in object or array `x` that have `key`
  * `x[~, :]` - recursively select all values in object or array `x`, returned as a flattened `List`
  * `x[:, (k, v) -> Bool]` - apply a key-value function `f` to each key-value/index-value in object or array `x`, and return a `List` of all values for which `f` returns `true`
"""
module Selectors

using ..StructUtils

export List, @selectors

"""
    List(...)

A custom array wrapper that supports the Selectors selection syntax.
"""
struct List{T} <: AbstractVector{T}
    items::Vector{T}
end

items(x::List) = getfield(x, :items)
Base.getindex(x::List) = map(getindex, items(x))
List(T=Any) = List(T[])
Base.size(x::List) = size(items(x))
Base.eltype(::List{T}) where {T} = T
Base.isassigned(x::List, args::Integer...) = isassigned(items(x), args...)

Base.push!(x::List, item) = push!(items(x), item)
Base.append!(x::List, items_to_append) = append!(items(x), items_to_append)

StructUtils.arraylike(::StructUtils.StructStyle, ::List) = true

function StructUtils.applyeach(::StructUtils.StructStyle, f, x::List)
    # note that there should *never* be #undef
    # values in a list, since we only ever initialize empty
    # then push!/append! to it
    for (i, v) in enumerate(items(x))
        ret = f(i, v)
        ret isa StructUtils.EarlyReturn && return ret
    end
    return
end

const KeyInd = Union{AbstractString, Symbol}
const Inds = Union{AbstractVector{<:KeyInd}, NTuple{N, <:KeyInd} where {N},
    AbstractVector{<:Integer}, NTuple{N, <:Integer} where {N}}

function _getindex_array(x, key::Union{KeyInd, Integer})
    values = List()
    StructUtils.applyeach(x) do _, item
        if StructUtils.structlike(StructUtils.DefaultStyle(), item)
            push!(values, _getindex(item, key))
        end
    end
    return values
end

function _getindex(x, key::Union{KeyInd, Integer})
    if StructUtils.arraylike(StructUtils.DefaultStyle(), x) && key isa KeyInd
        # indexing an array with a key, so we check
        # each element if it's an object and if the
        # object has the key
        # like a broadcasted getindex over x
        return _getindex_array(x, key)
    elseif StructUtils.structlike(StructUtils.DefaultStyle(), x) || StructUtils.arraylike(StructUtils.DefaultStyle(), x)
        # indexing object w/ key or array w/ index
        # returns a single value
        ret = StructUtils.applyeach(x) do k, v
            StructUtils.keyeq(k, key) && return StructUtils.EarlyReturn(v)
            return
        end
        ret isa StructUtils.EarlyReturn || throw(KeyError(key))
        return ret.value
    else
        noselection(x)
    end
end

# return all values of an object or elements of an array as a List
function _getindex(x, ::Colon)
    selectioncheck(x)
    values = List()
    StructUtils.applyeach(x) do _, v
        push!(values, v)
        return
    end
    return values
end

# a list is already a list of all its elements
_getindex(x::List, ::Colon) = x

# indexing object or array w/ a list of keys/indexes
function _getindex(x, inds::Inds)
    selectioncheck(x)
    values = List()
    StructUtils.applyeach(x) do k, v
        i = findfirst(StructUtils.keyeq(k), inds)
        i !== nothing && push!(values, v)
        return
    end
    return values
end

# return all values of an object or elements of an array as a List
# that satisfy a key-value function
function _getindex(x, S::Union{typeof(~), Colon}, f::Base.Callable)
    selectioncheck(x)
    values = List()
    StructUtils.applyeach(x) do k, v
        f(k, v) && push!(values, v)
        if S == ~
            if StructUtils.structlike(StructUtils.DefaultStyle(), v)
                ret = _getindex(v, ~, f)
                append!(values, ret)
            elseif StructUtils.arraylike(StructUtils.DefaultStyle(), v)
                ret = _getindex(v, ~, f)
                append!(values, ret)
            end
        end
        return
    end
    return values
end

# recursively return all values of an object or elements of an array as a List (:)
# as a single flattened List; or all properties that match key
function _getindex(x, ::typeof(~), key::Union{KeyInd, Colon})
    values = List()
    return _getindex(x, ~, key, values, key === Colon())
end

function _getindex(x, ::typeof(~), key::Union{KeyInd, Colon}, values, keymatched)
    if StructUtils.structlike(StructUtils.DefaultStyle(), x)
        StructUtils.applyeach(x) do k, v
            if key === Colon()
                _getindex(v, ~, key, values, true)
            elseif StructUtils.keyeq(k, key)
                # stop recursion if we find a match
                push!(values, v)
            else
                _getindex(v, ~, key, values, false)
            end
        end
    elseif StructUtils.arraylike(StructUtils.DefaultStyle(), x)
        StructUtils.applyeach(x) do _, item
            _getindex(item, ~, key, values, true)
        end
    elseif keymatched
        push!(values, x)
    end
    return values
end

# _get, like Base.get for objects
function _get(x, key::Union{KeyInd, Integer}, default)
    if StructUtils.arraylike(StructUtils.DefaultStyle(), x) && key isa KeyInd
        # indexing an array with a key, so we check
        # each element if it's an object and if the
        # object has the key
        # like a broadcasted getindex over x
        values = List()
        StructUtils.applyeach(x) do _, item
            if StructUtils.structlike(StructUtils.DefaultStyle(), item)
                # if array elements are objects, we do a broadcasted getproperty with `key`
                # should we try-catch and ignore KeyErrors?
                push!(values, _getindex(item, key))
            else
                # non-objects are just ignored
            end
            return
        end
        return values
    elseif StructUtils.structlike(StructUtils.DefaultStyle(), x) || StructUtils.arraylike(StructUtils.DefaultStyle(), x)
        # indexing object w/ key or array w/ index
        # returns a single value
        ret = StructUtils.applyeach(x) do k, v
            StructUtils.keyeq(k, key) && return StructUtils.EarlyReturn(v)
            return
        end
        ret isa StructUtils.EarlyReturn || return default
        return ret.value
    else
        noselection(x)
    end
end

selectioncheck(x) = StructUtils.structlike(StructUtils.DefaultStyle(), x) || StructUtils.arraylike(StructUtils.DefaultStyle(), x) || noselection(x)
@noinline noselection(x) = throw(ArgumentError("Selection syntax not defined for this object of type: `$(typeof(x))`"))

_symbol(x::AbstractString) = Symbol(x)
_symbol(x) = convert(Symbol, x)

# build up propertynames by iterating over each key-value pair
function _propertynames(x)
    selectioncheck(x)
    nms = Symbol[]
    StructUtils.applyeach(x) do k, _
        push!(nms, _symbol(k))
        return
    end
    return nms
end

# convenience macro for defining high-level getindex/getproperty methods
macro selectors(T)
    esc(quote
        Base.getindex(x::$T, arg) = StructUtils.Selectors._getindex(x, arg)
        Base.getindex(x::$T, ::Colon, arg) = StructUtils.Selectors._getindex(x, :, arg)
        Base.getindex(x::$T, ::typeof(~), arg) = StructUtils.Selectors._getindex(x, ~, arg)
        Base.getindex(x::$T, ::typeof(~), key, val) = StructUtils.Selectors._getindex(x, ~, key, val)
        Base.get(x::$T, key::Symbol, def) = StructUtils.Selectors._get(x, key, def)
        Base.getproperty(x::$T, key::Symbol) = StructUtils.Selectors._getindex(x, key)
        Base.propertynames(x::$T) = StructUtils.Selectors._propertynames(x)
        Base.hasproperty(x::$T, key::Symbol) = key in propertynames(x)
        Base.length(x::$T) = StructUtils.applylength(x)
    end)
end

@selectors List

end # module