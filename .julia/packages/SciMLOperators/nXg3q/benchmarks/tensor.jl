using SciMLOperators, LinearAlgebra, BenchmarkTools
using SciMLOperators: IdentityOperator, ⊗

N = 12
K = 100
Id = IdentityOperator(N)
A = rand(N, N)
B = rand(N, N)
C = rand(N, N)

println("#===============================#")
println("2D Tensor Products")
println("#===============================#")

println("⊗(A, B)")

u = rand(N^2, K)
v = rand(N^2, K)

T = ⊗(A, B)
T = cache_operator(T, u)

@btime *($T, $u)
@btime mul!($v, $T, $u)

println("⊗(I, B)")

u = rand(N^2, K)
v = rand(N^2, K)

T = ⊗(Id, B)
T = cache_operator(T, u)

@btime *($T, $u)
@btime mul!($v, $T, $u)

println("⊗(A, I)")

u = rand(N^2, K)
v = rand(N^2, K)

T = ⊗(A, Id)
T = cache_operator(T, u)

@btime *($T, $u)
@btime mul!($v, $T, $u)

println("#===============================#")
println("3D Tensor Products")
println("#===============================#")

println("⊗(⊗(A, B), C)")

u = rand(N^3, K)
v = rand(N^3, K)

T = ⊗(⊗(A, B), C)
T = cache_operator(T, u)

@btime *($T, $u)
@btime mul!($v, $T, $u); #

println("⊗(A, ⊗(B, C))")

u = rand(N^3, K)
v = rand(N^3, K)

T = ⊗(A, ⊗(B, C))
T = cache_operator(T, u)

@btime *($T, $u)
@btime mul!($v, $T, $u); #

println("#===============================#")
#
