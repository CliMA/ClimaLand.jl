@testset "complex" begin
  @testset "factorizable" begin
    # this matrix possesses an LDLᵀ factorization without pivoting
    #! format: off
    A = [
      1.7     0       0   0   0       0   0       0       0.13+im 0
      0       1.0     0   0   0.02+im 0   0       0       0       0.01+im
      0       0       1.5 0   0       0   0       0       0       0
      0       0       0   1.1 0       0   0       0       0       0
      0       0.02-im 0   0   2.6     0   0.16+im 0.09+im 0.52+im 0.53+im
      0       0       0   0   0       1.2 0       0       0       0
      0       0       0   0   0.16-im 0   1.3     0       0       0.56+im
      0       0       0   0   0.09-im 0   0       1.6     0.11+im 0
      0.13-im 0       0   0   0.52-im 0   0       0.11-im 1.4     0
      0       0.01-im 0   0   0.53-im 0   0.56-im 0       0       3.1
    ]
    #! format: on
    b = [
      0.287 + 0.9im,
      0.22 + 1.5im,
      0.45,
      0.44,
      2.486 + 3.2im,
      0.72,
      1.55 + 0.5im,
      1.424 + 0.4im,
      1.621 - 1.4im,
      3.759 - 1.4im,
    ]
    ϵ = sqrt(eps(real(eltype(A))))

    LDLT = ldl(A)
    x = LDLT \ b

    r = A * x - b
    @test norm(r) ≤ ϵ * norm(b)

    y = collect(0.1:0.1:1)
    @test norm(x - y) ≤ ϵ * norm(y)

    x2 = copy(b)
    ldiv!(LDLT, x2)

    r2 = A * x2 - b
    @test norm(r2) ≤ ϵ * norm(b)

    @test norm(x2 - y) ≤ ϵ * norm(y)

    @test norm(A[LDLT.P, LDLT.P] - (LDLT.L + I) * LDLT.D * (LDLT.L + I)') ≤ ϵ

    # test properties
    @test eltype(LDLT) == eltype(A)
    @test nnz(LDLT) == nnz(LDLT.L) + length(LDLT.d)
    @test size(LDLT) == size(A)
    @test typeof(LDLT.D) <: Diagonal
    @test propertynames(LDLT) == (:L, :D, :P)
  end

  @testset "factorizable_upper" begin
    # Using only the upper triangle tests
    #! format: off
    A = [
      1.7     0       0   0   0       0   0       0       0.13+im 0
      0       1.0     0   0   0.02+im 0   0       0       0       0.01+im
      0       0       1.5 0   0       0   0       0       0       0
      0       0       0   1.1 0       0   0       0       0       0
      0       0.02-im 0   0   2.6     0   0.16+im 0.09+im 0.52+im 0.53+im
      0       0       0   0   0       1.2 0       0       0       0
      0       0       0   0   0.16-im 0   1.3     0       0       0.56+im
      0       0       0   0   0.09-im 0   0       1.6     0.11+im 0
      0.13-im 0       0   0   0.52-im 0   0       0.11-im 1.4     0
      0       0.01-im 0   0   0.53-im 0   0.56-im 0       0       3.1
    ]
    #! format: on
    b = [
      0.287 + 0.9im,
      0.22 + 1.5im,
      0.45,
      0.44,
      2.486 + 3.2im,
      0.72,
      1.55 + 0.5im,
      1.424 + 0.4im,
      1.621 - 1.4im,
      3.759 - 1.4im,
    ]
    ϵ = sqrt(eps(real(eltype(A))))
    LDLT = ldl(A)

    A_upper = Hermitian(triu(A), :U)
    LDLT_upper = ldl(A_upper)
    x = LDLT_upper \ b

    y = collect(0.1:0.1:1)
    @test norm(x - y) ≤ ϵ * norm(y)

    r = A * x - b
    @test norm(r) ≤ ϵ * norm(b)

    x2 = copy(b)
    ldiv!(LDLT, x2)
    @test norm(x2 - y) ≤ ϵ * norm(y)

    r2 = A * x2 - b
    @test norm(r2) ≤ ϵ * norm(b)

    @test nnz(LDLT_upper) == nnz(LDLT_upper.L) + length(LDLT_upper.d)

    # test with separate analyze/factorize phases
    S = ldl_analyze(A)
    ldl_factorize!(A, S)
    x2 = copy(b)
    ldiv!(LDLT, x2)
    @test norm(x2 - y) ≤ ϵ * norm(y)
    r2 = A * x2 - b
    @test norm(r2) ≤ ϵ * norm(b)

    # test with a permutation of a different int type
    p = Int32.(collect(size(A, 1):-1:1))
    LDLT_upper = ldl(A_upper, P = p)
    x = LDLT_upper \ b

    y = collect(0.1:0.1:1)
    @test norm(x - y) ≤ ϵ * norm(y)

    r = A * x - b
    @test norm(r) ≤ ϵ * norm(b)

    x2 = copy(b)
    ldiv!(LDLT, x2)
    @test norm(x2 - y) ≤ ϵ * norm(y)

    r2 = A * x2 - b
    @test norm(r2) ≤ ϵ * norm(b)

    # Tests with multiple right-hand sides
    B = one(ComplexF64) * [i + j for j = 1:10, i = 0:3]
    X = A \ B
    Y = similar(B)
    ldiv!(Y, LDLT, B)
    @test norm(Y - X) ≤ ϵ * norm(X)
  end
end
