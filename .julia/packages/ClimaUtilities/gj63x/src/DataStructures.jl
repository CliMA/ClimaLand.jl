"""
    DataStructures

The `DataStructures` module implements useful data structures.

Currently, we have implemented a Least Recently Used (LRU) cache. The cache is
implemented as a dictionary with a maximum size. When the cache is full, the
least recently used item is removed to make space for the new item.

Within ClimaUtilities, the LRU cache is used by the `FileReaders` module to
store files that are currently open, and by the `DataHandler` module to store
regridded fields.
"""
module DataStructures

import Base: Callable

export LRUCache

struct LRUCache{K, V} <: AbstractDict{K, V}
    """The cache itself, containing key-value pairs of information."""
    cache::Dict{K, V}

    """The maximum number of key-value pairs in the cache."""
    max_size::Int

    """A list of keys, ordered by their last access time, which serves as a
    priority queue for the cache. Note that another data structure could be
    more efficient here, but a built-in vector is sufficient given that
    the cache is not expected to be very large."""
    priority::Vector{K}
end

"""
    LRUCache{K, V}(; max_size::Int = 128) where {K, V}

Construct an empty `LRUCache` with a maximum size of `max_size`.
"""
function LRUCache{K, V}(; max_size::Int = 128) where {K, V}
    return LRUCache{K, V}(Dict{K, V}(), max_size, Vector{K}())
end

"""
    Base.setindex!(cache::LRUCache{K, V}, value::V, key::K)

Store the mapping from `key` to `value` in `cache`. Then, return `cache`.
"""
function Base.setindex!(cache::LRUCache{K, V}, value::V, key::K) where {K, V}
    _update_priority!(cache, key)
    _enforce_size!(cache)
    cache.cache[key] = value
    return cache
end

"""
    Base.getindex(cache::LRUCache{K, V}, key::K)

Return the mapping for `key` in `cache`.
"""
function Base.getindex(cache::LRUCache{K, V}, key::K) where {K, V}
    value = cache.cache[key]
    _update_priority!(cache, key)
    return value
end

"""
    Base.get!(cache::LRUCache{K, V}, key::K, default::V) where {K, V}

Get the value associated with `key` in `cache`. If the key is not in the
cache, add it with the value `default`. In any case, update the key's priority
and make sure the cache doesn't exceed its maximum size.
"""
function Base.get!(cache::LRUCache{K, V}, key::K, default::V) where {K, V}
    _update_priority!(cache, key)
    _enforce_size!(cache)
    return get!(cache.cache, key, default)
end

"""
    Base.get!(default::Callable, cache::LRUCache{K, V}, key::K) where {K, V}

Get the value associated with `key` in `cache`. If the key is not in the
cache, add it with the value `default()`. In any case, update the key's priority
and make sure the cache doesn't exceed its maximum size.

This method is intended to be used with `do` block syntax.
"""
function Base.get!(
    default::Callable,
    cache::LRUCache{K, V},
    key::K,
) where {K, V}
    _update_priority!(cache, key)
    _enforce_size!(cache)
    return get!(default, cache.cache, key)
end

"""
    Base.length(cache::LRUCache{K, V})

Return the number of elements in `cache`.
"""
function Base.length(cache::LRUCache{K, V}) where {K, V}
    return length(cache.cache)
end

"""
    ==(cache1::LRUCache{K1, V2}, cache2::LRUCache{K2, V2})

Return whether the two caches are identical in content and priority.
"""
function Base.:(==)(
    cache1::LRUCache{K1, V1},
    cache2::LRUCache{K2, V2},
) where {K1, V1, K2, V2}
    return cache1.max_size == cache2.max_size &&
           cache1.cache == cache2.cache &&
           cache1.priority == cache2.priority
end

"""
    iterate(cache::LRUCache{K, V} [, state])

Advance the iterator to obtain the next element. If no elements remain, nothing
should be returned. Otherwise, a 2-tuple of the next element and the new
iteration state should be returned.
"""
function Base.iterate(cache::LRUCache{K, V}) where {K, V}
    return iterate(cache.cache)
end

function Base.iterate(cache::LRUCache{K, V}, state) where {K, V}
    return iterate(cache.cache, state)
end

"""
    Base.get(cache::LRUCache{K, V}, key::K, default::V)

Return the value associated with `key` in `cache`, or if `key` is not in `cache`,
return `default`. 
"""
function Base.get(cache::LRUCache{K, V}, key::K, default::V) where {K, V}
    haskey(cache, key) && return cache[key]
    return default
end

"""
    Base.get(default::Callable, cache::LRUCache{K, V}, key::K)

Return the value associated with `key` in `cache`, or if the key is not in the
cache, return `f()`.
"""
function Base.get(f::Callable, cache::LRUCache{K, V}, key::K) where {K, V}
    haskey(cache, key) && return cache[key]
    return f()
end

"""
    Base.haskey(cache::LRUCache{K, V}, key::K)
    
Return whether `key` has a mapping in `cache`. 
"""
function Base.haskey(cache::LRUCache{K, V}, key::K) where {K, V}
    return haskey(cache.cache, key)
end

"""
    Base.copy(cache::LRUCache{K, V})

Return a shallow copy of `cache`.
"""
function Base.copy(cache::LRUCache{K, V}) where {K, V}
    return LRUCache{K, V}(
        copy(cache.cache),
        cache.max_size,
        copy(cache.priority),
    )
end

"""
    Base.deepcopy(cache::LRUCache{K, V})

Return a deep copy of `cache`.
"""
function Base.deepcopy(cache::LRUCache{K, V}) where {K, V}
    return LRUCache{K, V}(
        deepcopy(cache.cache),
        cache.max_size,
        deepcopy(cache.priority),
    )
end

"""
    Base.empty(cache::LRUCache{K, V}, key_type::DataType=K, value_type::DataType=V)

Create an empty `LRUCache` with keys of type `key_type` and values of type `value_type`.

The second and third arguments are optional and default to the input's `keytype` and `valuetype`. 
If only one of the two type is specified, it is assumed to be the `value_type`, and `key_type` is 
assumed to the key type of `cache`.
"""
function Base.empty(
    cache::LRUCache{K, V},
    key_type::DataType,
    value_type::DataType,
) where {K, V}
    return LRUCache{key_type, value_type}(max_size = cache.max_size)
end

function Base.empty(
    cache::LRUCache{K, V},
    value_type::DataType = V,
) where {K, V}
    return LRUCache{K, value_type}(max_size = cache.max_size)
end

"""
    Base.empty!(cache::LRUCache{K, V})

Remove all key/value mappings from the input and return the emptied input.
"""
function Base.empty!(cache::LRUCache{K, V}) where {K, V}
    empty!(cache.cache)
    empty!(cache.priority)
    return cache
end

"""
    Base.pop!(cache::LRUCache{K, V}, key::K, [default::V])

Delete the mapping for `key` in `cache` and return the associated `value`, or 
if the `key` does not exist, return `default` or throw an error if `default` is not specified
"""
function Base.pop!(cache::LRUCache{K, V}, key::K, default::V) where {K, V}
    value = pop!(cache.cache, key, default)
    filter!(k -> k != key, cache.priority)
    return value
end

function Base.pop!(cache::LRUCache{K, V}, key::K) where {K, V}
    value = pop!(cache.cache, key)
    filter!(k -> k != key, cache.priority)
    return value
end

"""
    Base.delete!(cache::LRUCache{K, V}, key::K)

Delete mapping for `key`, if it is in`cache`, and return `cache`
"""
function Base.delete!(cache::LRUCache{K, V}, key::K) where {K, V}
    value = delete!(cache.cache, key)
    filter!(k -> k != key, cache.priority)
    return cache
end

"""
   Base.merge(cache_1::LRUCache{K1, V1}, cache_2::LRUCache{K2, V2})

Merge not implemented for LRUCache
"""
function Base.merge(
    cache_1::LRUCache{K1, V1},
    cache_2::LRUCache{K2, V2},
) where {K1, V1, K2, V2}
    error("Merge not implemented for LRUCache")
    return nothing
end

"""
    _update_priority!(cache, key)

Update the priority of `key` in `cache` to reflect its most recent access.
"""
function _update_priority!(cache, key)
    # Remove the key if it's already in the priority queue
    filter!(k -> k != key, cache.priority)
    # Add the key to the end of the priority queue
    push!(cache.priority, key)
    return
end

"""
    _enforce_size!(cache::LRUCache{K, V})

Remove the least recently used items from `cache` until its size is
less than or equal to `cache.max_size`.
"""
function _enforce_size!(cache::LRUCache{K, V}) where {K, V}
    if length(cache.priority) > cache.max_size
        key = popfirst!(cache.priority)
        delete!(cache.cache, key)
    end
end


end # module
