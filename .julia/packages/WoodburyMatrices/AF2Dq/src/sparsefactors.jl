# Used by Interpolations.jl

"""
Given a type `T`, an integer `n` and `m` tuples of the form (i, j, v), build
sparse matrices `rows`, `vals`, `cols` such that the product
`out = rows * vals * cols` is equivalent to:

```julia
out = zeros(T, n, n)

for (i, j, v) in args
    out[i, j] = v
end
```

The first two components (`i` and `j`) of each tuple should be integers
whereas the third component should be of type `T`

Example:


```
julia> r, v, c = WoodburyMatrices.sparse_factors(Float64, 3,
                                                 (1, 1, 2.0),
                                                 (2, 2, 3.0),
                                                 (3, 3, 4.0));

julia> r*c*v - Diagonal([2.0, 3.0, 4.0])
3x3 sparse matrix with 0 Float64 entries:
```

"""
function sparse_factors(::Type{T}, n::Int, args::Tuple{Int, Int, Any}...) where {T}
    m = length(args)
    rows = spzeros(T, n, m)
    cols = spzeros(T, m, n)
    vals = zeros(T, m, m)

    ix = 1
    for (i, (row, col, val)) in enumerate(args)
        rows[row, ix] = 1
        cols[ix, col] = 1
        vals[ix, ix] = val
        ix += 1
    end

    rows, vals, cols
end
