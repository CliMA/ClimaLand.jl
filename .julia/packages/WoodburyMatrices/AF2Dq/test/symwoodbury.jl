using WoodburyMatrices
using Test

@testset "SymWoodbury" begin
seed!(123)
n = 5

for elty in (Float32, Float64, ComplexF32, ComplexF64, Int)

    elty = Float64

    a = rand(n); B = rand(n,2); D = Symmetric(rand(2,2)); v = rand(n)

    if elty == Int
        v = rand(1:100, n)
        a = rand(1:100, n)
        B = rand(1:100, 2, n)
        D = Symmetric(rand(1:100, 2, 2))
    else
        v = convert(Vector{elty}, v)
        a = convert(Vector{elty}, a)
        B = convert(Matrix{elty}, B)
        D = convert(Matrix{elty}, D)
    end

    ε = eps(abs2(float(one(elty))))
    A = Diagonal(a)

    for W in (SymWoodbury(A, B, D), SymWoodbury(A, B, D; allocatetmp=true), SymWoodbury(A, B[:,1][:], 2.))

        @test issymmetric(W)
        F = Matrix(W)
        @test (2*W)*v ≈ 2*(W*v)
        @test W'*v ≈ W*v
        @test (W'W)*v ≈ Matrix(W)*(Matrix(W)*v)
        @test (W*W)*v ≈ Matrix(W)*(Matrix(W)*v)
        @test (W*W')*v ≈ Matrix(W)*(Matrix(W)*v)
        @test W[1:3,1:3]*v[1:3] ≈ Matrix(W)[1:3,1:3]*v[1:3]
        @test sparse(W) ≈ Matrix(W)
        @test W === W'
        @test W*Matrix(1.0I, n, n) ≈ Matrix(W)
        @test W'*Matrix(1.0I, n, n) ≈ Matrix(W)

        Z = randn(n,n)
        @test Matrix(W*Z) ≈ Matrix(W)*Z

        R = rand(n,n)

        for v = (rand(n, 1), view(rand(n,1), 1:n), view(rand(n,2),1:n,1:2))
            @test (2*W)*v ≈ 2*(W*v)
            @test (W*2)*v ≈ 2*(W*v)
            @test (W/2)*v ≈ (W*v)/2
            @test (W'W)*v ≈ Matrix(W)*(Matrix(W)*v)
            @test (W*W)*v ≈ Matrix(W)*(Matrix(W)*v)
            @test (W*W')*v ≈ Matrix(W)*(Matrix(W)*v)
            @test W[1:3,1:3]*v[1:3] ≈ Matrix(W)[1:3,1:3]*v[1:3]
            @test Matrix(WoodburyMatrices.conjm(W, R)) ≈ R*Matrix(W)*R'
            @test Matrix(copy(W)'W)*v ≈ Matrix(W)*(Matrix(W)*v)
            @test Matrix(W + A) ≈ Matrix(W)+Matrix(A)
            @test Matrix(A + W) ≈ Matrix(W)+Matrix(A)
        end

        v = rand(n,1)
        W2 = convert(Woodbury, W)
        @test Matrix(W2) ≈ Matrix(W)

        if elty != Int
            @test inv(W)*v ≈ inv(Matrix(W))*v
            @test W\v ≈ inv(Matrix(W))*v
            @test WoodburyMatrices.partialInv(W)[1] ≈ inv(W).B
            @test WoodburyMatrices.partialInv(W)[2] ≈ inv(W).D
            @test det(W) ≈ det(Matrix(W))
        end

    end

end


for elty in (Float32, Float64, ComplexF32, ComplexF64, Int)

    elty = Float64

    a1 = rand(n); B1 = rand(n,2); D1 = Symmetric(rand(2,2)); v = rand(n)
    a2 = rand(n); B2 = rand(n,2); D2 = Symmetric(rand(2,2));

    if elty == Int
        v = rand(1:100, n)

        a1 = rand(1:100, n)
        B1 = rand(1:100, 2, n)
        D1 = Symmetric(rand(1:100, 2, 2))

        a2 = rand(1:100, n)
        B2 = rand(1:100, 2, n)
        D2 = Symmetric(rand(1:100, 2, 2))
    else
        v = convert(Vector{elty}, v)

        a1 = convert(Vector{elty}, a1)
        B1 = convert(Matrix{elty}, B1)
        D1 = convert(Matrix{elty}, D1)

        a2 = convert(Vector{elty}, a2)
        B2 = convert(Matrix{elty}, B2)
        D2 = convert(Matrix{elty}, D2)
    end

    ε = eps(abs2(float(one(elty))))

    A1 = Diagonal(a1)
    A2 = Matrix(Diagonal(a2))

    W1 = SymWoodbury(A1, B1, D1)
    W2 = SymWoodbury(A2, B2, D2)

    W1r = SymWoodbury(A1, B1[:,1][:], 2.)
    W2r = SymWoodbury(A2, B2[:,1][:], 3.)

    for (W1, W2) = ((W1,W2), (W1r, W2), (W1, W2r), (W1r,W2r))
        @test (W1 + W2)*v ≈ (Matrix(W1) + Matrix(W2))*v
        @test (Matrix(W1) + W2)*v ≈ (Matrix(W1) + Matrix(W2))*v
        @test (W1 + 2*Matrix(Diagonal((a1))))*v ≈ (Matrix(W1) + Matrix(2*Matrix(Diagonal((a1)))))*v
        @test_throws MethodError W1*W2
    end

end

# Sparse U and D

A = Diagonal((rand(n)))
B = sprandn(n,2,1.)
D = sprandn(2,2,1.); D = (D + D')/2
W = SymWoodbury(A, B, D)
v = randn(n)
vdiag = Diagonal(v)
V = randn(n,1)

@test size(W) == (n,n)
@test size(W,1) == n
@test size(W,2) == n

@test inv(W)*v ≈ inv(Matrix(W))*v
@test (2*W)*v ≈ 2*(W*v)
@test (W'W)*v ≈ Matrix(W)*(Matrix(W)*v)
@test (W*W)*v ≈ Matrix(W)*(Matrix(W)*v)
@test (W*W')*v ≈ Matrix(W)*(Matrix(W)*v)

@test inv(W)*vdiag ≈ W\vdiag
@test W\vdiag ≈ W\Matrix(vdiag)
@test inv(W)*vdiag ≈ inv(Matrix(W))*vdiag

@test inv(W)*V ≈ inv(Matrix(W))*V
@test (2*W)*V ≈ 2*(W*V)
@test (W'W)*V ≈ Matrix(W)*(Matrix(W)*V)
@test (W*W)*V ≈ Matrix(W)*(Matrix(W)*V)
@test (W*W')*V ≈ Matrix(W)*(Matrix(W)*V)

@test Matrix(2*W) ≈ 2*Matrix(W)
@test Matrix(W*2) ≈ 2*Matrix(W)
R = Symmetric(rand(size(W)...))
# Not sure when this got fixed, but check that failures are for a "known" reason
canadd = try (W + R; true;) catch err; (@test err isa MethodError && err.f == ldiv!; false;) end
if canadd
    @test Matrix(W + R) ≈ Matrix(W) + R
    @test Matrix(R + W) ≈ Matrix(W) + R
else
    @test_broken Matrix(W + R) ≈ Matrix(W) + R
    @test_broken Matrix(R + W) ≈ Matrix(W) + R
end
Wm = SymWoodbury(A, Matrix(B), D)
@test Matrix(Wm + R) ≈ Matrix(Wm) + R
@test Matrix(R + Wm) ≈ Matrix(Wm) + R

@test transpose(W)*v ≈ transpose(Matrix(W))*v

# Factorization for A
A = SymTridiagonal(rand(5).+2, rand(4))
B = rand(5)
D = 2
x = rand(5)
W1 = SymWoodbury(A, B, D)
@test W1 \ x ≈ Matrix(W1) \ x
for AA in (ldlt(A), cholesky(Matrix(A)), bunchkaufman(Matrix(A), true))
    W2 = SymWoodbury(AA, B, D)
    @test W2 \ x ≈ W1 \ x
end
A = Symmetric(rand(5, 5))
W1 = SymWoodbury(A, B, D)
@test W1 \ x ≈ Matrix(W1) \ x
for AA in (bunchkaufman(Matrix(A), true),)
    W2 = SymWoodbury(AA, B, D)
    @test W2 \ x ≈ W1 \ x
end
A = sparse(A)
W1 = SymWoodbury(A, B, D)
@test W1 \ x ≈ Matrix(W1) \ x

# Mismatched sizes
@test_throws DimensionMismatch SymWoodbury(rand(5,5),rand(5,2),rand(2,3))
@test_throws DimensionMismatch SymWoodbury(rand(5,5),rand(5,2),rand(3,3))
@test_throws DimensionMismatch SymWoodbury(rand(5,5),rand(3),1.)

# Asymmetric
@test_throws ArgumentError SymWoodbury(rand(5,5),rand(5,2),rand(2,2))

# Display
iob = IOBuffer()
show(iob, MIME("text/plain"), W)
str = String(take!(iob))
@test occursin("D:", str)
show(iob, W)
str = String(take!(iob))
@test occursin("SymWoodbury{Float64}", str)
@test occursin("D=", str)

# logdet
# make sure all matrices are PSD because I don't want complex numbers
W = SymWoodbury([randpsd(50) for _ in 1:3]...)
@test logdet(W) ≈ log(det(W)) ≈ logdet(Array(W))
@test all(logabsdet(W) .≈ logabsdet(Array(W)))

@testset "Diagonal (issue #35)" begin
    A = Diagonal( Float64[1,2,3,4] )
    V = [5  6  7  8;
         9 10 11 12]
    U = V'
    D = Matrix( I, 2, 2 )

    W = SymWoodbury( A, U, D )
    @test W \ [13, 14, 15, 16] ≈ Matrix(W) \ [13, 14, 15, 16]
end

@testset "Empty B" begin
    A = Diagonal( Float64[1,2,3,4] )
    B = Matrix{Float64}(undef, 4, 0)
    D = Diagonal(Float64[])

    W = SymWoodbury( A, B, D)
    @test W \ [13, 14, 15, 16] ≈ Matrix(W) \ [13, 14, 15, 16]
    @test det(W) ≈ det(Matrix(W))
    @test logdet(W) ≈ logdet(Matrix(W))
    @test all(logabsdet(W) .≈ logabsdet(Matrix(W)))
end

end # @testset "SymWoodbury"
