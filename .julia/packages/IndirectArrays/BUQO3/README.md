# IndirectArrays

[![Build Status](https://travis-ci.org/JuliaArrays/IndirectArrays.jl.svg?branch=master)](https://travis-ci.org/JuliaArrays/IndirectArrays.jl)

[![codecov.io](http://codecov.io/github/JuliaArrays/IndirectArrays.jl/coverage.svg?branch=master)](http://codecov.io/github/JuliaArrays/IndirectArrays.jl?branch=master)

An `IndirectArray` is one that encodes data using a combination of an
`index` and a `value` table. Each element is assigned its own index, which
is used to retrieve the value from the `value` table.  Concretely, if
`A` is an `IndirectArray`, then `A[i,j...] = value[index[i,j,...]]`.

Among other uses, `IndirectArrays` can represent
[indexed images](https://en.wikipedia.org/wiki/Indexed_color),
sometimes called "colormap images" or "paletted images."

## Installation

```jl
Pkg.add("IndirectArrays")
```

## Usage

For example:

```
using IndirectArrays, Colors

colors = distinguishable_colors(6)
index = rand(1:6, 32, 32)
A = IndirectArray(index, colors)
```

![random image](randimage.png)

which has only 6 colors in it.

The `value` array can be of any type; it does not have to be color information.

## Related packages

- [CategoricalArrays](https://github.com/nalimilan/CategoricalArrays.jl) offers an even more flexible interface for dealing with arrays in which values are looked up in an index.
- [PooledArrays](https://github.com/JuliaComputing/PooledArrays.jl)
