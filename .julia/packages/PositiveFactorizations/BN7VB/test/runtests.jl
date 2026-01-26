using PositiveFactorizations
using LinearAlgebra, Test

import ForwardDiff, ReverseDiff

@testset "PositiveFactorizations" begin
for pivot in (Val{false}, Val{true})
    A = [1 0; 0 1]
    F = cholesky(Positive, A, pivot)
    @test Matrix(F) ≈ A
    F, d = ldlt(Positive, A, pivot)
    @test Matrix(F) ≈ A
    @test d == [1,1]

    A = [1 0; 0 -1]
    F = cholesky(Positive, A, pivot)
    @test Matrix(F) ≈ Matrix(1.0I,2,2)
    F, d = ldlt(Positive, A, pivot)
    @test Matrix(F) ≈ Matrix(1.0I,2,2)
    @test d == [1,-1]

    A = [-1 0.5; 0.5 4]
    target = pivot == Val{false} ? [1 -0.5; -0.5 4.5] : [1.125 0.5; 0.5 4]
    dtarget = pivot == Val{false} ? [-1,1] : [1,-1]
    F = cholesky(Positive, A, pivot)
    @test Matrix(F) ≈ target
    F = cholesky(Positive, 10*A, pivot)
    @test Matrix(F) ≈ 10*target
    F, d = ldlt(Positive, A, pivot)
    @test Matrix(F) ≈ target
    @test d == dtarget

    A = [0 1; 1 0]
    F = cholesky(Positive, A, pivot)
    @test Matrix(F) ≈ Matrix(1.0I,2,2)
    F, d = ldlt(Positive, A, pivot)
    @test Matrix(F) ≈ Matrix(1.0I,2,2)
    @test d == [0,0]

    A = rand(201,200); A = A'*A
    F = cholesky(Positive, A, pivot)
    @test Matrix(F) ≈ A
    F, d = ldlt(Positive, A, pivot)
    @test Matrix(F) ≈ A
    @test all(d .== 1)

    # factorization of (not too small) BigFloat matrices passes
    a = BigFloat.(1:15); A = a * a'
    F = cholesky(Positive, A, pivot)
    F, d = ldlt(Positive, A, pivot)

    # Differentiability test
    g = ForwardDiff.gradient((x) -> det(cholesky(Positive, Matrix(Hermitian(Diagonal(x))), pivot)), [ 2.0, 3.0 ])
    @test g ≈ [ 3.0, 2.0 ]

    g = ReverseDiff.gradient((x) -> det(cholesky(Positive, Matrix(Hermitian(Diagonal(x))), pivot)), [ 2.0, 3.0 ])
    @test g ≈ [ 3.0, 2.0 ]

    # Extra 'collect' are needed for ReverseDiff
    vec_to_hermitian = (v) -> begin A = I - 2 * v * collect(v'); A = collect(A') * A end;

    v = rand(10)
    g1 = ForwardDiff.gradient((x) -> det(cholesky(Positive, vec_to_hermitian(x), pivot)), v)
    g2 = ForwardDiff.gradient((x) -> det(cholesky(vec_to_hermitian(x))), v)
    @test g1 ≈ g2
    g1 = ReverseDiff.gradient((x) -> det(cholesky(Positive, vec_to_hermitian(x), pivot)), v)
    g2 = ReverseDiff.gradient((x) -> det(cholesky(vec_to_hermitian(x))), v)
    @test g1 ≈ g2
end

A = [1 0; 0 -2]
F = eigen(Positive, A)

# TODO: Use this when we drop v0.4 support
# @test Matrix(F) ≈ abs.(A)
absA = abs.(A)
absA = convert(Array{Int}, absA) # v0.4 fix
@test Matrix(F) ≈ absA

A = [1 0; 0 0]
F = eigen(Positive, A)
@test Matrix(F) ≈ Matrix(1.0I,2,2)

# Test whether necessary matrix operations are supported for SubArrays
n = PositiveFactorizations.default_blocksize(Float64)
B = rand(n+3,n+4); C = rand(size(B)...); A = B'*B - C'*C
ldlt!(Positive, A)
end   # @testset
