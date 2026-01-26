# PaddedViews

[![][action-img]][action-url]
[![][pkgeval-img]][pkgeval-url]
[![][codecov-img]][codecov-url]

## Summary

PaddedViews provides a simple wrapper type, `PaddedView`, to add
"virtual" padding to any array without copying data. Edge values not
specified by the array are assigned a `fillvalue`.  Multiple arrays
may be "promoted" to have common indices using the `paddedviews`
function.

`PaddedView` arrays are read-only, meaning that you cannot assign
values to them. The original array may be extracted using `A =
parent(P)`, where `P` is a `PaddedView`.

## Examples

For padding a single array:

```julia
julia> a = collect(reshape(1:9, 3, 3))
3×3 Array{Int64,2}:
 1  4  7
 2  5  8
 3  6  9

julia> PaddedView(-1, a, (4, 5))
4×5 PaddedView(-1, ::Array{Int64,2}, (Base.OneTo(4), Base.OneTo(5))) with eltype Int64:
  1   4   7  -1  -1
  2   5   8  -1  -1
  3   6   9  -1  -1
 -1  -1  -1  -1  -1

 julia> PaddedView(-1, a, (1:5,1:5), (2:4,2:4))
 5×5 PaddedView(-1, OffsetArray(::Array{Int64,2}, 2:4, 2:4), (1:5, 1:5)) with eltype Int64 with indices 1:5×1:5:
 -1  -1  -1  -1  -1
 -1   1   4   7  -1
 -1   2   5   8  -1
 -1   3   6   9  -1
 -1  -1  -1  -1  -1

 julia> PaddedView(-1, a, (0:4, 0:4))
 5×5 PaddedView(-1, ::Array{Int64,2}, (0:4, 0:4)) with eltype Int64 with indices 0:4×0:4:
  -1  -1  -1  -1  -1
  -1   1   4   7  -1
  -1   2   5   8  -1
  -1   3   6   9  -1
  -1  -1  -1  -1  -1

julia> PaddedView(-1, a, (5,5), (2,2))
5×5 PaddedView(-1, OffsetArray(::Array{Int64,2}, 2:4, 2:4), (Base.OneTo(5), Base.OneTo(5))) with eltype Int64:
 -1  -1  -1  -1  -1
 -1   1   4   7  -1
 -1   2   5   8  -1
 -1   3   6   9  -1
 -1  -1  -1  -1  -1
```

For padding multiple arrays to have common indices:

```julia
julia> a1 = reshape([1, 2, 3], 3, 1)
3×1 Array{Int64,2}:
 1
 2
 3

julia> a2 = [4 5 6]
1×3 Array{Int64,2}:
 4  5  6

julia> a1p, a2p = paddedviews(-1, a1, a2);

julia> a1p
3×3 PaddedView(-1, ::Array{Int64,2}, (Base.OneTo(3), Base.OneTo(3))) with eltype Int64:
 1  -1  -1
 2  -1  -1
 3  -1  -1

julia> a2p
3×3 PaddedView(-1, ::Array{Int64,2}, (Base.OneTo(3), Base.OneTo(3))) with eltype Int64:
  4   5   6
 -1  -1  -1
 -1  -1  -1
```

If you want original arrays in the center of padded results:

```julia
julia> a1 = reshape([1, 2, 3], 3, 1)
3×1 Array{Int64,2}:
 1
 2
 3

julia> a2 = [4 5 6]
1×3 Array{Int64,2}:
 4  5  6

julia> a1p, a2p = sym_paddedviews(-1, a1, a2);

julia> a1p
3×3 PaddedView(-1, ::Array{Int64,2}, (1:3, 0:2)) with eltype Int64 with indices 1:3×0:2:
 -1  1  -1
 -1  2  -1
 -1  3  -1

julia> a2p
3×3 PaddedView(-1, ::Array{Int64,2}, (0:2, 1:3)) with eltype Int64 with indices 0:2×1:3:
 -1  -1  -1
  4   5   6
 -1  -1  -1
```


We can use `dims` keyword to specify the direction of padding for the output arrays as shown above. They can be integers or tuple of integers.

```julia
julia> a1p, a2p = paddedviews(-1, a1, a2; dims=1);

julia> a1p
3×1 PaddedView(-1, ::Matrix{Int64}, (Base.OneTo(3), Base.OneTo(1))) with eltype Int64:
 1
 2
 3

julia> a2p
3×3 PaddedView(-1, ::Matrix{Int64}, (Base.OneTo(3), Base.OneTo(3))) with eltype Int64:
  4   5   6
 -1  -1  -1
 -1  -1  -1

julia> a1p, a2p = sym_paddedviews(-1, a1, a2; dims=2);

julia> a1p
3×3 PaddedView(-1, ::Matrix{Int64}, (1:3, 0:2)) with eltype Int64 with indices 1:3×0:2:
 -1  1  -1
 -1  2  -1
 -1  3  -1

julia> a2p
1×3 PaddedView(-1, ::Matrix{Int64}, (1:1, 1:3)) with eltype Int64 with indices 1:1×1:3:
 4  5  6
```


<!-- badges and urls -->

[pkgeval-img]: https://juliaci.github.io/NanosoldierReports/pkgeval_badges/P/PaddedViews.svg
[pkgeval-url]: https://juliaci.github.io/NanosoldierReports/pkgeval_badges/report.html
[action-img]: https://github.com/JuliaArrays/PaddedViews.jl/workflows/CI/badge.svg
[action-url]: https://github.com/JuliaArrays/PaddedViews.jl/actions
[codecov-img]: https://codecov.io/github/JuliaArrays/PaddedViews.jl/coverage.svg?branch=master
[codecov-url]: https://codecov.io/github/JuliaArrays/PaddedViews.jl?branch=master
