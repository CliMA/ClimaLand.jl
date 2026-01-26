# TensorCore.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaMath.github.io/TensorCore.jl/stable)
[![CI](https://github.com/JuliaMath/TensorCore.jl/workflows/CI/badge.svg)](https://github.com/JuliaMath/TensorCore.jl/actions?query=workflow%3ACI)
[![Codecov](https://codecov.io/gh/JuliaMath/TensorCore.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaMath/TensorCore.jl)

This package is intended as a lightweight foundation for tensor operations across the Julia ecosystem.
Currently it exports three operations:
* `hadamard` elementwise multiplication, with unicode operator `⊙`,
* `tensor` product preserves all dimensions, operator `⊗`, and
* `boxdot` contracts neighbouring dimensions, named after the unicode `⊡`.

```julia
julia> using TensorCore

julia> A = [1 2 3; 4 5 6]; B = [7 8 9; 0 10 20];

julia> A ⊙ B  # hadamard(A, B)
2×3 Matrix{Int64}:
 7  16   27
 0  50  120

julia> V = [1, 10];

julia> C = A ⊗ V  # tensor(A, V)
2×3×2 Array{Int64, 3}:
[:, :, 1] =
 1  2  3
 4  5  6

[:, :, 2] =
 10  20  30
 40  50  60

julia> summary(A ⊗ B)
"2×3×2×3 Array{Int64, 4}"

julia> C ⊡ V  # boxdot(C, V)
2×3 Matrix{Int64}:
 101  202  303
 404  505  606

julia> summary(C ⊡ rand(2,5,7))
"2×3×5×7 Array{Float64, 4}"
 ```
