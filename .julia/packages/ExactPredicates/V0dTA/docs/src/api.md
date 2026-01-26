# API

The pseudo type `:: 2` indicate `Tuple{Float64,Float64}` or any type `T` such
that `ExactPredicates.coord(::T)` or `Base.Tuple(::T)` outputs a
`Tuple{Float64,Float64}`. Similarly, `:: 3` indicates
`Tuple{Float64,Float64,Float64}` or any type convertible to it `coord` or
`Tuple`.

## Planar predicates

```@autodocs
Modules = [ExactPredicates]
Pages = ["plane.jl"]
```

## Spatial predicates


```@autodocs
Modules = [ExactPredicates]
Pages = ["space.jl"]
```

