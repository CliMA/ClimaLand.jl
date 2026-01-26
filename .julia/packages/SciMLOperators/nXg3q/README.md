# `SciMLOperators.jl`

*Unified operator interface for `SciML.ai` and beyond*

[![Join the chat at https://julialang.zulipchat.com #sciml-bridged](https://img.shields.io/static/v1?label=Zulip&message=chat&color=9558b2&labelColor=389826)](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged)
[![Global Docs](https://img.shields.io/badge/docs-SciML-blue.svg)](https://docs.sciml.ai/SciMLOperators/stable)

[![codecov](https://codecov.io/gh/SciML/SciMLOperators.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/SciML/SciMLOperators.jl)
[![Build Status](https://github.com/SciML/SciMLOperators.jl/workflows/CI/badge.svg)](https://github.com/SciML/SciMLOperators.jl/actions?query=workflow%3ACI)

[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor%27s%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)

`SciMLOperators` is a package for managing linear, nonlinear,
time-dependent, and parameter dependent operators acting on vectors,
(or column-vectors of matrices). We provide wrappers for matrix-free
operators, fast tensor-product evaluations, pre-cached mutating
evaluations, as well as `Zygote`-compatible non-mutating evaluations.

The lazily implemented operator algebra allows the user to update the
operator state by passing in an update function that accepts arbitrary
parameter objects. Further, our operators behave like `AbstractMatrix` types
thanks to  overloads defined for methods in `Base`, and `LinearAlgebra`.

Therefore, an `AbstractSciMLOperator` can be passed to `LinearSolve.jl`,
or `NonlinearSolve.jl` as a linear/nonlinear operator, or to
`OrdinaryDiffEq.jl` as an `ODEFunction`. Examples of usage within the
`SciML` ecosystem are provided in the documentation.

## Installation

`SciMLOperators.jl` is a registered package and can be installed via

```
julia> import Pkg
julia> Pkg.add("SciMLOperators")
```

## Examples

Let `M`, `D`, `F` be matrix-based, diagonal-matrix-based, and function-based
`SciMLOperators` respectively.

Let `M`, `D`, `F` be matrix-based, diagonal-matrix-based, and function-based
`SciMLOperators` respectively.

```@example operator_algebra
using SciMLOperators, LinearAlgebra
N = 4
function f(v, u, p, t)
    u .* v
end
function f(w, v, u, p, t)
    w .= u .* v
end

u = rand(4)
p = nothing # parameter struct
t = 0.0     # time

M = MatrixOperator(rand(N, N))
D = DiagonalOperator(rand(N))
F = FunctionOperator(f, zeros(N), zeros(N); u, p, t)
```

Then, the following codes just work.

```@example operator_algebra
L1 = 2M + 3F + LinearAlgebra.I + rand(N, N)
L2 = D * F * M'
L3 = kron(M, D, F)
L4 = lu(M) \ D
L5 = [M; D]' * [M F; F D] * [F; D]
```

Each `L#` can be applied to `AbstractVector`s of appropriate sizes:

```@example operator_algebra
v = rand(N)
w = L1(v, u, p, t) # == L1 * v

v_kron = rand(N^3)
w_kron = L3(v_kron, u, p, t) # == L3 * v_kron
```

For mutating operator evaluations, call `cache_operator` to generate an
in-place cache, so the operation is nonallocating.

```@example operator_algebra
α, β = rand(2)

# allocate cache
L2 = cache_operator(L2, u)
L4 = cache_operator(L4, u)

# allocation-free evaluation
L2(w, v, u, p, t) # == mul!(w, L2, v)
L4(w, v, u, p, t, α, β) # == mul!(w, L4, v, α, β)
```

## Roadmap

  - [ ] [Complete integration with `SciML` ecosystem](https://github.com/SciML/SciMLOperators.jl/issues/142)
  - [ ] [Block-matrices](https://github.com/SciML/SciMLOperators.jl/issues/161)
  - [x] [Benchmark and speed-up tensorbproduct evaluations](https://github.com/SciML/SciMLOperators.jl/issues/58)
  - [ ] [Fast tensor-sum (`kronsum`) evaluation](https://github.com/SciML/SciMLOperators.jl/issues/53)
  - [ ] [Fully flesh out operator array algebra](https://github.com/SciML/SciMLOperators.jl/issues/62)
  - [ ] [Operator fusion/matrix chain multiplication at constant `(u, p, t)`-slices](https://github.com/SciML/SciMLOperators.jl/issues/51)
