"""
$(README)
"""
module SciMLOperators

using DocStringExtensions

using LinearAlgebra

import ArrayInterface
import MacroTools: @forward
import Accessors: @reset

# overload
import Base: show
import Base: zero, one, oneunit
import Base: +, -, *, /, \, ∘, ==, conj, exp, kron
import Base: iszero, inv, adjoint, transpose, size, convert
import LinearAlgebra: mul!, ldiv!, lmul!, rmul!, factorize
import LinearAlgebra: Matrix, Diagonal

# Used for downstream checking
const isv1 = true

"""
$(TYPEDEF)

Subtypes of `AbstractSciMLOperator` represent linear, nonlinear,
time-dependent operators acting on vectors, or matrix column-vectors.
A lazy operator algebra is also defined for `AbstractSciMLOperator`s.

## Mathematical Notation

An `AbstractSciMLOperator` ``L`` is an operator which is used to represent
the following type of equation:

```math
w = L(u,p,t)[v]
```

where `L[v]` is the operator application of ``L`` on the vector ``v``.

## Interface

An `AbstractSciMLOperator` can be called  like a function in the following ways:

  - `L(v, u, p, t)` - Out-of-place application where `v` is the action vector and `u` is the update vector
  - `L(w, v, u, p, t)` - In-place application where `w` is the destination, `v` is the action vector, and `u` is the update vector
  - `L(w, v, u, p, t, α, β)` - In-place application with scaling: `w = α*(L*v) + β*w`

Operator state can be updated separately from application:

  - `update_coefficients!(L, u, p, t)` for in-place operator update
  - `L = update_coefficients(L, u, p, t)` for out-of-place operator update

SciMLOperators also overloads `Base.*`, `LinearAlgebra.mul!`,
`LinearAlgebra.ldiv!` for operator evaluation without updating operator state.
An `AbstractSciMLOperator` behaves like a matrix in these methods.
Allocation-free methods, suffixed with a `!` often need cache arrays.
To precache an `AbstractSciMLOperator`, call the function
`L = cache_operator(L, input_vector)`.

## Overloaded Actions

The behavior of a `SciMLOperator` is
indistinguishable from an `AbstractMatrix`. These operators can be
passed to linear solver packages, and even to ordinary differential
equation solvers. The list of overloads to the `AbstractMatrix`
interface includes, but is not limited to, the following:

  - `Base: size, zero, one, +, -, *, /, \\, ∘, inv, adjoint, transpose, convert`
  - `LinearAlgebra: mul!, ldiv!, lmul!, rmul!, factorize, issymmetric, ishermitian, isposdef`
  - `SparseArrays: sparse, issparse`

## Multidimensional arrays and batching

SciMLOperator can also be applied to `AbstractMatrix` subtypes where
operator-evaluation is done column-wise.

```julia
K = 10
u_mat = rand(N, K)

v_mat = F(u_mat, p, t) # == mul!(v_mat, F, u_mat)
size(v_mat) == (N, K) # true
```

`L` can also be applied to `AbstractArray`s that are not
`AbstractVecOrMat`s so long as their size in the first dimension is appropriate
for matrix-multiplication. Internally, `SciMLOperator`s reshapes an
`N`-dimensional array to an `AbstractMatrix`, and applies the operator via
matrix-multiplication.

## Operator update

This package can also be used to write state-dependent, time-dependent, and
parameter-dependent operators, whose state can be updated per
a user-defined function.
The updates can be done in-place, i.e. by mutating the object,
or out-of-place, i.e. in a non-mutating, `Zygote`-compatible way.

For example,

```julia
u = rand(N)
p = rand(N)
t = rand()

# out-of-place update
mat_update_func = (A, u, p, t) -> t * (p * u')
sca_update_func = (a, u, p, t) -> t * sum(p)

M = MatrixOperator(zero(N, N); update_func = mat_update_func)
α = ScalarOperator(zero(Float64); update_func = sca_update_func)

L = α * M
L = cache_operator(L, v)

# L is initialized with zero state
L * v == zeros(N) # true

# update operator state with `(u, p, t)`
L = update_coefficients(L, u, p, t)
# and multiply
L * v != zeros(N) # true

# updates state and evaluates L*v at (u, p, t)
L(v, u, p, t) != zeros(N) # true
```

The out-of-place evaluation function `L(v, u, p, t)` calls
`update_coefficients` under the hood, which recursively calls
the `update_func` for each component `SciMLOperator`.
Therefore, the out-of-place evaluation function is equivalent to
calling `update_coefficients` followed by `Base.*`. Notice that
the out-of-place evaluation does not return the updated operator.

On the other hand, the in-place evaluation function, `L(w, v, u, p, t)`,
mutates `L`, and is equivalent to calling `update_coefficients!`
followed by `mul!`. The in-place update behavior works the same way,
with a few `<!>`s appended here and there. For example,

```julia
w = rand(N)
v = rand(N)
u = rand(N)
p = rand(N)
t = rand()

# in-place update
_A = rand(N, N)
_d = rand(N)
mat_update_func! = (A, u, p, t) -> (copy!(A, _A); lmul!(t, A); nothing)
diag_update_func! = (diag, u, p, t) -> copy!(diag, N)

M = MatrixOperator(zero(N, N); update_func! = mat_update_func!)
D = DiagonalOperator(zero(N); update_func! = diag_update_func!)

L = D * M
L = cache_operator(L, v)

# L is initialized with zero state
L * v == zeros(N) # true

# update L in-place
update_coefficients!(L, v, p, t)
# and multiply
mul!(w, v, u, p, t) != zero(N) # true

# updates L in-place, and evaluates w=L*v at (u, p, t)
L(w, v, u, p, t) != zero(N) # true
```

The update behavior makes this package flexible enough to be used
in `OrdinaryDiffEq`. As the parameter object `p` is often reserved
for sensitivity computation via automatic-differentiation, a user may
prefer to pass in state information via other arguments. For that
reason, we allow update functions with arbitrary keyword arguments.

```julia
mat_update_func = (A, u, p, t; scale = 0.0) -> scale * (p * u')

M = MatrixOperator(zero(N, N); update_func = mat_update_func,
    accepted_kwargs = (:state,))

M(v, u, p, t) == zeros(N) # true
M(v, u, p, t; scale = 1.0) != zero(N)
```
"""
abstract type AbstractSciMLOperator{T} end

"""
$(TYPEDEF)

An `AbstractSciMLScalarOperator` represents a linear scaling operation
that may be applied to `Number`, `AbstractVecOrMat` subtypes. Addition,
multiplication, division of `AbstractSciMLScalarOperator`s is defined
lazily so operator state may be updated.
"""
abstract type AbstractSciMLScalarOperator{T} <: AbstractSciMLOperator{T} end

include("utils.jl")
include("interface.jl")
include("left.jl")

include("scalar.jl")
include("matrix.jl")
include("basic.jl")
include("batch.jl")
include("func.jl")
include("tensor.jl")

export
       IdentityOperator,
       NullOperator,
       ScalarOperator,
       MatrixOperator,
       DiagonalOperator,
       InvertibleOperator,
       AffineOperator,
       AddVector,
       FunctionOperator,
       TensorProductOperator

export update_coefficients!,
       update_coefficients, isconstant,
       iscached,
       cache_operator, issquare,
       islinear,
       concretize,
       isconvertible, has_adjoint,
       has_expmv,
       has_expmv!,
       has_exp,
       has_mul,
       has_mul!,
       has_ldiv,
       has_ldiv!

end # module
