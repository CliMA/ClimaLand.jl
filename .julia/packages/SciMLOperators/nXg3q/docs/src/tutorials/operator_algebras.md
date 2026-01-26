# [Demonstration of Operator Algebras and Kron](@id operator_algebras)

Let `M`, `D`, `F` be matrix-based, diagonal-matrix-based, and function-based
`SciMLOperators` respectively. Here are some examples of composing operators
in order to build more complex objects and using their operations.

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
