using LinearAlgebra

const n = 8
a, b = legendre_coefs(Float64, n)

Tn = SymTridiagonal(copy(a), b[2:n])
F = eigen(Tn)

v = zeros(n)
special_eigenproblem!(a, b, v, 30)

idx = sortperm(a)
@test a[idx] ≈ F.values
@test abs.(v[idx]) ≈ abs.(F.vectors[1,:])
