# CompositionsBase.jl: exports `∘`, `⨟`, `compose`, and `opcompose`

## API

    f ∘ g
    g ⨟ f
    compose(f, g)
    opcompose(g, f)

Composition of morphisms.  `∘` is the operator defined in `Base`.
CompositionsBase.jl defines the opposite composition operator `⨟` as

    ⨟(fs...) = ∘(reverse(fs)...)

and also the ASCII aliases `compose` and `opcompose`.

As `⨟`, `compose`, and `opcompose` are all defined in terms of `∘`,
single-argument call is the identity function.

### Examples
```julia
julia> using CompositionsBase

julia> tuple ∘ inv === compose(tuple, inv) === inv ⨟ tuple === opcompose(inv, tuple)
true

julia> ∘(tuple) === compose(tuple) === ⨟(tuple) === opcompose(tuple) === tuple
true
```
