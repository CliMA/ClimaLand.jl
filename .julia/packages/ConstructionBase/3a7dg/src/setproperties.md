    setproperties(obj, patch::NamedTuple)


Return a copy of `obj` with properties updated according to `patch`.

# Examples
```jldoctest
julia> using ConstructionBase

julia> struct S
           a
           b
           c
       end

julia> s = S(1,2,3)
S(1, 2, 3)

julia> setproperties(s, (a=10,c=4))
S(10, 2, 4)

julia> setproperties((a=1,c=2,b=3), (a=10,c=4))
(a = 10, c = 4, b = 3)
```

There is also a convenience method, which builds the `patch` argument from
keywords:

    setproperties(obj; kw...)

# Examples
```jldoctest
julia> using ConstructionBase

julia> struct S
           a
           b
           c
       end

julia> o = S(10, 2, 4)
S(10, 2, 4)

julia> setproperties(o, a="A", c="cc")
S("A", 2, "cc")
```

## Specification

`setproperties` belongs to [the semantic level](@ref the-semantic-level). If satisfies the following invariants:

1. Purity: `setproperties` is supposed to have no side effects. In particular `setproperties(obj, patch::NamedTuple)` may not mutate `obj`.
2. Relation to `propertynames` and `fieldnames`: `setproperties` relates to `propertynames` and `getproperty`, not to `fieldnames` and `getfield`.
   This means that any subset `p₁, p₂, ..., pₙ` of `propertynames(obj)` is a valid set of properties, with respect to which the lens laws below must hold.
3. `setproperties` is defined in relation to `getproperties` so that:
   ```julia
   obj == setproperties(obj, getproperties(obj))
   ```
4. `setproperties` should satisfy the lens laws:

For any valid set of properties `p₁, p₂, ..., pₙ`, following equalities must hold:

* You get what you set.

```julia
let obj′ = setproperties(obj, ($p₁=v₁, $p₂=v₂, ..., $pₙ=vₙ))
    @assert obj′.$p₁ == v₁
    @assert obj′.$p₂ == v₂
    ...
    @assert obj′.$pₙ == vₙ
end
```

* Setting what was already there changes nothing:

```julia
@assert setproperties(obj, ($p₁=obj.$p₁, $p₂=obj.$p₂, ..., $pₙ=obj.$pₙ)) == obj
```

* The last set wins:
```julia
let obj′ = setproperties(obj, ($p₁=v₁, $p₂=v₂, ..., $pₙ=vₙ)),
    obj′′ = setproperties(obj′, ($p₁=w₁, $p₂=w₂, ..., $pₙ=wₙ))
    @assert obj′′.$p₁ == w₁
    @assert obj′′.$p₂ == w₂
    ...
    @assert obj′′.$pₙ == wₙ
end
```

# Implementation

For a custom type `MyType`, a method `setproperties(obj::MyType, patch::NamedTuple)`
may be defined. When doing so it is important to ensure compliance with the specification.

* Prefer to overload [`constructorof`](@ref) whenever makes sense (e.g., no `getproperty`
  method is defined).  Default `setproperties` is defined in terms of `constructorof` and `getproperties`.

* If `getproperty` is customized, it may be a good idea to define `setproperties`.

!!! warning
    The signature `setproperties(obj::MyType; kw...)` should never be overloaded.
    Instead `setproperties(obj::MyType, patch::NamedTuple)` should be overloaded.

