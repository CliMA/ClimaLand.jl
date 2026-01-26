# GilbertCurves.jl

[![Build Status](https://github.com/CliMA/GilbertCurves.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/CliMA/GilbertCurves.jl/actions/workflows/CI.yml?query=branch%3Amain)

This is a Julia implementation of the generalized Hilbert ("gilbert") space-filling curve algorithm, by Jakub Červený (https://github.com/jakubcerveny/gilbert).
It provides space-filling curves for ectangular domains of arbitrary (non-power-of-two) sizes. Currently only 2D domains are supported, but it could be extended to 3D.


# Usage

Currently it exports one function, `gilbertindices` which returns a vector of
`CartesianIndex{2}` objects corresponding to their order on the curve:
```julia
julia> using GilbertCurves

julia> list = gilbertindices((5,5))
25-element Vector{CartesianIndex{2}}:
 CartesianIndex(1, 1)
 CartesianIndex(2, 1)
 CartesianIndex(2, 2)
 CartesianIndex(1, 2)
 CartesianIndex(1, 3)
 CartesianIndex(1, 4)
 CartesianIndex(1, 5)
 ⋮
 CartesianIndex(5, 3)
 CartesianIndex(5, 2)
 CartesianIndex(4, 2)
 CartesianIndex(3, 2)
 CartesianIndex(3, 1)
 CartesianIndex(4, 1)
 CartesianIndex(5, 1)
```

Two non-exported functions are also provided. `GilbertCurves.linearindices` takes the output of
`gilbertindices`, returning an integer-valued matrix of the gilbert indices of each component.

```julia
julia> GilbertCurves.linearindices(list)
5×5 Matrix{Int64}:
  1   4   5   6   7
  2   3  10   9   8
 23  22  11  12  13
 24  21  18  17  14
 25  20  19  16  15
```

`GilbertCurves.gilbertorder` constructs a vector containing the elements of a matrix in the
gilbert curve order.

```julia
julia> GilbertCurves.gilbertorder(reshape(1:9,3,3))
9-element Vector{Int64}:
 1
 4
 7
 8
 9
 6
 5
 2
 3
```

# Example

```julia
julia> using Plots

julia> list = gilbertindices((67,29));

julia> plot([c[1] for c in list], [c[2] for c in list], line_z=1:length(list), legend=false)
```
![Gilbert curve on 67 x 29 elements](https://raw.githubusercontent.com/CliMA/GilbertCurves.jl/master/img/67x29.png)

# Notes

The algorithm is not able to avoid non-diagonal moves in the case when the
larger dimension is odd and the smaller is even.

```julia
julia> list = gilbertindices((15,12));

julia> plot([c[1] for c in list], [c[2] for c in list], line_z=1:length(list), legend=false)
```

![Gilbert curve on 15 x 12 elements](https://raw.githubusercontent.com/CliMA/GilbertCurves.jl/master/img/15x12.png)
