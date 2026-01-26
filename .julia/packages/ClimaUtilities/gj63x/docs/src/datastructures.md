# `DataStructures`

The `DataStructures` module implements helpful data structures to be used by
other ClimaUtilities.jl modules or external packages.

## `LRUCache`

`DataStructures` implements an `LRUCache`, which is used in both `DataHandlingExt`
and `NCFileReaderExt`. `LRUCache` can be used to store pieces of information that
need to be accessed multiple times and may be expensive to compute, such as
regridded fields or loaded files. Instead of recomputing the values each time
they're needed, the previously-computed information is saved in the cache and
retrieved directly from it. To prevent the cache from growing so large that it
takes up significant memory, the least-recently-used (LRU) scheme maintains
a maximum cache size: every time adding an element to the cache would
lead its size to grow larger thant the maximum allowed, the element that was
accessed the least recently is deleted first.

!!! note

    All the methods supported by dictionaries are currently implemented.
    If you need one that is not implemented, please open an issue or a pull request.

### Example
In many ways, `LRUCache`s behave like Julia dictionaries. To use one, we
first need to initialize the empty cache specifying types of keys and values
and optionally the maximum allowed size:

```julia
import ClimaUtilities.DataStructures
cache = DataStructures.LRUCache{Int, Int}(; max_size = 128)
```

Once we have the `cache`, we can access and insert elements with `get!`.
`get!` retrieves the value associated to the key if available, otherwise
it inserts a new key with a provided default. The default can be passed as
third argument or can be the return statement of a function provided as
first (as typically done in `do` blocks). Continuing the example above

```julia
println("The value associated to 1 is: ", get!(cache, 1, 10))
other_value = get!(cache, 2) do
    return 20
end
```

This cache can be used to implement some recursive function more efficiently,
as in the famous case of the factorial:
```julia
function factorial(n)
    n in (0, 1) && return 1
    return n * get!(cache, n - 1, factorial(n - 1))
end
# The first time we call factorial it will compute everything
@time factorial(10)
# The second time it will only compute the last element
@time factorial(11)
```

## API

```@docs
ClimaUtilities.DataStructures.LRUCache
Base.get!
Base.get
Base.haskey
Base.copy
Base.deepcopy
Base.empty
Base.empty!
Base.pop!
Base.delete!
Base.getindex
Base.setindex!
Base.length
```
