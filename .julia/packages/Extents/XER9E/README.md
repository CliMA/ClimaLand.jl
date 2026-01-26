# Extents

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://rafaqz.github.io/Extents.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://rafaqz.github.io/Extents.jl/dev)
[![Build Status](https://github.com/rafaqz/Extents.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/rafaqz/Extents.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/rafaqz/Extents.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/rafaqz/Extents.jl)

Extents.jl is a small package that defines an `Extent` object that can be used by the
different Julia spatial data packages. `Extent` is a wrapper for a NamedTuple of tuples
holding the lower and upper bounds for each dimension of a object. It is used in
[GeoInterface.jl](https://github.com/JuliaGeo/GeoInterface.jl/), as the required return type
for the `extent` function, and in [DimensionalData.jl](https://github.com/rafaqz/DimensionalData.jl)
and [Rasters.jl](https://github.com/rafaqz/Rasters.jl) to subset arrays with named dimensions.

# Quick start

Install the package:
```julia
]
add Extents
```

Then to use it:
```julia
using Extents
```

Create Extents:
```julia
extent1 = Extent(X = (1.0, 2.0), Y = (3.0, 4.0));
extent2 = Extent(X = (1.5, 2.5), Y = (3.0, 4.0));
extent3 = Extent(X = (-1.0, 5.0), Y = (2.0, 5.0));
```

Find dimensions:
```julia
keys(extent1)
```
```
(:X, :Y)
```

Get extent for a single dimension by name:
```julia
extent1.X
```
```julia
(1.0, 2.0)
```

Get underlying NamedTuple:
```julia
bounds(extent1)
```
```julia
(X = (1.0, 2.0), Y = (3.0, 4.0))
```

Get intersection of two extents:
```julia
Extents.intersection(extent1, extent2)
```
```julia
Extent(X = (1.5, 2.0), Y = (3.0, 4.0))
```

Get union of two extents:
```julia
Extents.union(extent1, extent2)
```
```julia
Extent(X = (1.0, 2.5), Y = (3.0, 4.0))
```

Use reduce() to operate over collections:
```julia
extents = [extent1, extent2, extent3]
```
```julia
3-element Vector{Extent{(:X, :Y), Tuple{Tuple{Float64, Float64}, Tuple{Float64, Float64}}}}:
 Extent(X = (1.0, 2.0), Y = (3.0, 4.0))
 Extent(X = (1.5, 2.5), Y = (3.0, 4.0))
 Extent(X = (-1.0, 5.0), Y = (2.0, 5.0))
```
```julia
reduce(Extents.intersection, extents)
Extent(X = (1.5, 2.0), Y = (3.0, 4.0))
```
```julia
reduce(Extents.union, extents)
Extent(X = (-1.0, 5.0), Y = (2.0, 5.0))
```

Extents.jl also defines spatial predicates following the [DE-9IM](https://en.wikipedia.org/wiki/DE-9IM) standard.

```julia
Extents.intersects(extent1, extent2)
true
```
```julia
Extents.disjoint(extent1, extent2)
false
```
```julia
Extents.touches(extent1, extent2)
false
```
```julia
Extents.overlaps(extent1, extent2)
true
```

See [the docs](https://rafaqz.github.io/Extents.jl/stable) for all available methods.