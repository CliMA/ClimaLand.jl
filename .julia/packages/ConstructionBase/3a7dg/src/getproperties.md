    getproperties(obj)::NamedTuple
    getproperties(obj::Tuple)::Tuple

Return the properties of `obj` as a `NamedTuple`. Since `Tuple` don't have symbolic properties,
`getproperties` is the identity function on tuples.

# Examples
```jldoctest
julia> using ConstructionBase

julia> struct S
           a
           b
           c
       end

julia> s = S(1, 2, 3)
S(1, 2, 3)

julia> getproperties(s)
(a = 1, b = 2, c = 3)

julia> getproperties((10,20))
(10, 20)
```

## Specification

`getproperties` belongs to [the semantic level](@ref the-semantic-level).
`getproperties` guarantees a couple of invariants. When overloading it, the user is responsible for ensuring them:

1. `getproperties` should be consistent with `Base.propertynames`, `Base.getproperty`, `Base.setproperty!`. 
    Semantically it should be equivalent to:
    ```julia
    function getproperties(obj)
        fnames = propertynames(obj)
        NamedTuple{fnames}(getproperty.(Ref(obj), fnames))
    end
    ```
2. `getproperties` is defined in relation to `setproperties` so that:
   ```julia
   obj == setproperties(obj, getproperties(obj))
   ```
   The only exception from this semantics is that undefined properties may be avoided 
   in the return value of `getproperties`.

# Implementation

`getproperties` is defined by default for all objects. It should be very rare that a custom type `MyType`, has to implement `getproperties(obj::MyType)`. Reasons to do so are undefined fields or performance considerations.
