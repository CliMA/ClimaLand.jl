# HashArrayMappedTries.jl

A [HashArrayMappedTrie](https://en.wikipedia.org/wiki/Hash_array_mapped_trie) or
HAMT for short, is a data-structure that can be used efficient persistent hash tables.

## Usage

```julia
dict = HAMT{Int, Int}()
dict[1] = 1
delete!(dict, 1)
```

### Persitency

```julia
dict = HAMT{Int, Int}()
dict = insert(dict, 1, 1)
dict = delete(dict, 1)
```

## Robustness against hash collisions

The HAMT is robust to hash collision as long as they are not collection.
As an example of a devious hash take.

```
mutable struct CollidingHash
end
Base.hash(::CollidingHash, h::UInt) = hash(UInt(0), h)

ch1 = CollidingHash()
ch2 = CollidingHash()
```

For all `h` `hash(ch1, h) == hash(ch2, h)`. `Base.Dict` is robust under those
hashes as well.

```
dict = Dict{CollidingHash, Nothing}()
dict[CollidingHash()] = nothing
dict[CollidingHash()] = nothing
display(dict)

# Dict{CollidingHash, Nothing} with 2 entries:
#  CollidingHash() => nothing
#  CollidingHash() => nothing
```

Whereas HAMT.

```
dict = HAMT{CollidingHash, Nothing}()
dict[CollidingHash()] = nothing
dict[CollidingHash()] = nothing
ERROR: Perfect hash collision detected
```
