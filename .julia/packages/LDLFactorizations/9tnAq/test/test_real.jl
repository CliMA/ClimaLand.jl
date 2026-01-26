@testset "real" begin
  @testset "factorizable" begin
    # this matrix possesses an LDLᵀ factorization without pivoting
    A = [
      1.7 0 0 0 0 0 0 0 0.13 0
      0 1.0 0 0 0.02 0 0 0 0 0.01
      0 0 1.5 0 0 0 0 0 0 0
      0 0 0 1.1 0 0 0 0 0 0
      0 0.02 0 0 2.6 0 0.16 0.09 0.52 0.53
      0 0 0 0 0 1.2 0 0 0 0
      0 0 0 0 0.16 0 1.3 0 0 0.56
      0 0 0 0 0.09 0 0 1.6 0.11 0
      0.13 0 0 0 0.52 0 0 0.11 1.4 0
      0 0.01 0 0 0.53 0 0.56 0 0 3.1
    ]
    b = [0.287, 0.22, 0.45, 0.44, 2.486, 0.72, 1.55, 1.424, 1.621, 3.759]
    ϵ = sqrt(eps(eltype(A)))

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

  @testset "not_factorizable" begin
    # this matrix does not possess an LDLᵀ factorization without pivoting
    A = [
      0 1
      1 1
    ]
    S = ldl(A, P = [1, 2])
    @test !factorized(S)

    A = Symmetric(sparse(triu(A)))
    S = ldl_analyze(A)
    ldl_factorize!(A, S)
    @test !factorized(S)
  end

  @testset "sparse" begin
    for Ti in (Int32, Int), Tf in (Float32, Float64, BigFloat)
      A = sparse(Ti[1, 2, 1, 2], Ti[1, 1, 2, 2], Tf[10, 2, 2, 5])
      z = ones(Tf, 2)
      b = A * z
      LDLT = ldl(A)
      x = LDLT \ b
      @test norm(x - z) ≤ sqrt(eps(Tf)) * norm(z)
      r = A * x - b
      @test norm(r) ≤ sqrt(eps(Tf)) * norm(b)

      x2 = copy(b)
      ldiv!(LDLT, x2)
      @test norm(x2 - z) ≤ sqrt(eps(Tf)) * norm(z)
      r2 = A * x2 - b
      @test norm(r2) ≤ sqrt(eps(Tf)) * norm(b)

      y = similar(b)
      ldiv!(y, LDLT, b)
      @test norm(y - z) ≤ sqrt(eps(Tf)) * norm(z)
      r2 = A * y - b
      @test norm(r2) ≤ sqrt(eps(Tf)) * norm(b)

      @test nnz(LDLT) == nnz(LDLT.L) + length(LDLT.d)

      # test with separate analyze/factorize phases
      # A = Symmetric(sparse(triu(A)))
      S = ldl_analyze(A)
      ldl_factorize!(A, S)
      x2 = copy(b)
      ldiv!(LDLT, x2)
      @test norm(x2 - z) ≤ sqrt(eps(Tf)) * norm(z)
      r2 = A * x2 - b
      @test norm(r2) ≤ sqrt(eps(Tf)) * norm(b)
    end
  end

  @testset "factorizable_upper" begin
    # Using only the upper triangle tests
    A = [
      1.7 0 0 0 0 0 0 0 0.13 0
      0 1.0 0 0 0.02 0 0 0 0 0.01
      0 0 1.5 0 0 0 0 0 0 0
      0 0 0 1.1 0 0 0 0 0 0
      0 0.02 0 0 2.6 0 0.16 0.09 0.52 0.53
      0 0 0 0 0 1.2 0 0 0 0
      0 0 0 0 0.16 0 1.3 0 0 0.56
      0 0 0 0 0.09 0 0 1.6 0.11 0
      0.13 0 0 0 0.52 0 0 0.11 1.4 0
      0 0.01 0 0 0.53 0 0.56 0 0 3.1
    ]
    b = [0.287, 0.22, 0.45, 0.44, 2.486, 0.72, 1.55, 1.424, 1.621, 3.759]
    ϵ = sqrt(eps(eltype(A)))
    LDLT = ldl(A)

    A_upper = Symmetric(triu(A), :U)
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
    # A = Symmetric(sparse(triu(A)))
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
    B = 1.0 * [i + j for j = 1:10, i = 0:3]
    X = A \ B
    Y = similar(B)
    ldiv!(Y, LDLT, B)
    @test norm(Y - X) ≤ ϵ * norm(X)
  end

  @testset "not_factorizable_upper" begin
    # this matrix does not possess an LDLᵀ factorization without pivoting
    A = triu([
      0 1
      1 1
    ])
    S = ldl(A, P = [1, 2])
    @test !factorized(S)

    S = ldl_analyze(Symmetric(A, :U))
    ldl_factorize!(Symmetric(A, :U), S)
    @test !S.__factorized
  end

  @testset "sparse_upper" begin
    for Ti in (Int32, Int), Tf in (Float32, Float64, BigFloat)
      A = sparse(Ti[1, 2, 1, 2], Ti[1, 1, 2, 2], Tf[10, 2, 2, 5])
      A_upper = Symmetric(triu(A), :U)
      b = A * ones(Tf, 2)
      LDLT = ldl(A_upper)
      x = LDLT \ b
      r = A * x - b
      @test norm(r) ≤ sqrt(eps(Tf)) * norm(b)

      x2 = copy(b)
      ldiv!(LDLT, x2)
      r2 = A * x2 - b
      @test norm(r2) ≤ sqrt(eps(Tf)) * norm(b)

      @test nnz(LDLT) == nnz(LDLT.L) + length(LDLT.d)
    end
  end

  @testset "positive_semidefinite" begin
    A = [
      0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 4.0 0.0
      0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 5.0 0.0
      2.0 4.0 5.0 -2 4.0 1.0 2.0 2.0 2.0 0.0
      0.0 0.0 0.0 0.0 1.0 9.0 9.0 1.0 7.0 1.0
      0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0
      1.0 3.0 2.0 1.0 4.0 3.0 1.0 0.0 0.0 7.0
      -3.0 8.0 0.0 0.0 0.0 0.0 -2.0 0.0 0.0 1.0
      0.0 0.0 0.0 5.0 7.0 9.0 0.0 2.0 7.0 1.0
      3.0 2.0 0.0 0.0 0.0 0.0 1.0 3.0 3.0 2.0
      0.0 0.0 0.0 0.0 -3 -4 0.0 0.0 0.0 0.0
    ]
    ϵ = sqrt(eps(eltype(A)))
    M = A * A'  # det(A) = 0 => M positive semidefinite
    b = M * ones(10)
    x = copy(b)
    S = ldl_analyze(Symmetric(triu(M), :U))
    S.r1 = -ϵ
    S.r2 = ϵ
    S.tol = ϵ
    S.n_d = 0
    S = ldl_factorize!(Symmetric(triu(M), :U), S)
    x = ldiv!(S, x)
    r = M * x - b
    @test norm(r) ≤ sqrt(eps()) * norm(b)
  end

  @testset "SQD" begin
    A = [
      0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 4.0 0.0
      0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 5.0 0.0
      2.0 4.0 5.0 -2 4.0 1.0 2.0 2.0 2.0 0.0
      0.0 0.0 0.0 0.0 1.0 9.0 9.0 1.0 7.0 1.0
      0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0
      1.0 3.0 2.0 1.0 4.0 3.0 1.0 0.0 0.0 7.0
      -3.0 8.0 0.0 0.0 0.0 0.0 -2.0 0.0 0.0 1.0
      0.0 0.0 0.0 5.0 7.0 9.0 0.0 2.0 7.0 1.0
      3.0 2.0 0.0 0.0 0.0 0.0 1.0 3.0 3.0 2.0
      0.0 0.0 0.0 0.0 -3 -4 0.0 0.0 0.0 0.0
    ]
    ϵ = sqrt(eps(eltype(A)))
    M = spzeros(20, 20)
    M[1:10, 1:10] = -A * A'
    M[11:20, 11:20] = A * A'
    # M = [-A*A'    0
    #        0     A*A'] where A*A' is symmetric positive semidefinite
    b = M * ones(20)
    x = copy(b)
    S = ldl_analyze(Symmetric(triu(M), :U))
    S.r1 = -ϵ
    S.r2 = ϵ
    S.tol = ϵ
    S.n_d = 10
    S = ldl_factorize!(Symmetric(triu(M), :U), S)
    x = ldiv!(S, x)
    r = M * x - b
    @test norm(r) ≤ sqrt(eps()) * norm(b)
  end

  @testset "SQD_semi_dynamic" begin
    A = [
      0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 4.0 0.0
      0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 5.0 0.0
      2.0 4.0 5.0 -2 4.0 1.0 2.0 2.0 2.0 0.0
      0.0 0.0 0.0 0.0 1.0 9.0 9.0 1.0 7.0 1.0
      0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0
      1.0 3.0 2.0 1.0 4.0 3.0 1.0 0.0 0.0 7.0
      -3.0 8.0 0.0 0.0 0.0 0.0 -2.0 0.0 0.0 1.0
      0.0 0.0 0.0 5.0 7.0 9.0 0.0 2.0 7.0 1.0
      3.0 2.0 0.0 0.0 0.0 0.0 1.0 3.0 3.0 2.0
      0.0 0.0 0.0 0.0 -3 -4 0.0 0.0 0.0 0.0
    ]
    ϵ = sqrt(eps(eltype(A)))
    M = spzeros(20, 20)
    M[1:10, 1:10] = A * A' + ϵ * I
    M[11:20, 11:20] = -A * A'
    # M = [-A*A'    0
    #        0     A*A'] where A*A' is symmetric positive semidefinite
    b = M * ones(20)
    x = copy(b)
    S = ldl_analyze(Symmetric(triu(M), :U))
    S.r1 = zero(eltype(A))
    S.r2 = -ϵ
    S.tol = ϵ
    S.n_d = 10
    S = ldl_factorize!(Symmetric(triu(M), :U), S)
    x = ldiv!(S, x)
    r = M * x - b
    @test norm(r) ≤ sqrt(eps()) * norm(b)
  end

  @testset "Test booleans and allocations" begin
    A = spdiagm(0 => 2 * ones(10), 1 => -ones(9))
    ϵ = √eps()
    M1 = spzeros(20, 20)
    M1[1:10, 1:10] = A
    M1[11:20, 11:20] = spdiagm(0 => -ϵ * ones(10))
    M2 = spzeros(20, 20)
    M2[1:10, 1:10] = A
    M2[11:20, 11:20] = spdiagm(0 => -100ϵ * ones(10))
    M3 = spzeros(20, 20)
    M3[1:10, 1:10] = A
    M3[11:20, 11:20] = spdiagm(0 => zeros(10))

    S = ldl_analyze(M1)
    @test !factorized(S)
    _allocs1 = @allocated ldl_factorize!(M1, S)
    @test S.d[11:20] ≈ -ϵ * ones(10)
    @test factorized(S)
    _allocs2 = @allocated ldl_factorize!(M2, S)
    @test S.d[11:20] ≈ -100ϵ * ones(10)
    @test _allocs1 == _allocs2
    ldl_factorize!(M3, S)
    @test !factorized(S)
    b = ones(20)
    @test_throws LDLFactorizations.SQDException ldiv!(S, b)
    @test_throws LDLFactorizations.SQDException lmul!(S, b)
    B = ones(20, 2)
    @test_throws LDLFactorizations.SQDException ldiv!(S, B)
  end

  @testset "ldl_mul!" begin
    A0 = [
      0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 4.0 0.0
      0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 5.0 0.0
      2.0 4.0 5.0 -2 4.0 1.0 2.0 2.0 2.0 0.0
      0.0 0.0 0.0 0.0 1.0 9.0 9.0 1.0 7.0 1.0
      0.0 -7.0 0.0 4.0 0.0 0.0 0.0 0.0 1.0 0.0
      1.0 3.0 2.0 1.0 4.0 3.0 1.0 0.0 0.0 7.0
      -3.0 8.0 0.0 0.0 0.0 0.0 -2.0 0.0 0.0 1.0
      0.0 0.0 0.0 5.0 7.0 9.0 0.0 2.0 7.0 1.0
      3.0 2.0 0.0 0.0 0.0 0.0 1.0 3.0 3.0 2.0
      0.0 0.0 0.0 0.0 -3 -4 0.0 0.0 0.0 0.0
    ]
    A = A0 * A0'
    b = A * ones(10)
    T = eltype(A)
    ϵ = sqrt(eps(T))
    S = ldl_analyze(Symmetric(triu(A), :U))
    ldl_factorize!(Symmetric(triu(A), :U), S)
    x = copy(b)
    ldiv!(S, x)
    r1 = copy(x)
    lmul!(S, r1)
    @test norm(r1 - A * x) ≤ sqrt(ϵ) * norm(r1)
    r1 .-= b  # r1 = LDL*x-b
    @test norm(r1) ≤ sqrt(ϵ) * norm(b)
    mul!(r1, S, x) # same test using mul!
    @test norm(r1 - A * x) ≤ sqrt(ϵ) * norm(r1)
    r1 .-= b
    @test norm(r1) ≤ sqrt(ϵ) * norm(b)

    LDL1 = spzeros(T, size(A)...)
    LDL1[diagind(LDL1)] .= one(T)
    LDL1 = lmul!(S, LDL1)
    r2 = similar(b)
    mul!(r2, LDL1, x)
    @test norm(r2 - A * x) ≤ sqrt(ϵ) * norm(r2)
    r2 .-= b # r2 = LDL*I*x-b
    @test norm(r2) ≤ sqrt(ϵ) * norm(b)
    LDL2 = spzeros(T, size(A)...)
    LDL2[diagind(LDL2)] .= one(T)
    mul!(LDL2, S, copy(LDL2))
    mul!(r2, LDL2, x) # same test using mul!
    @test norm(r2 - A * x) ≤ sqrt(ϵ) * norm(r2)
    r2 .-= b
    @test norm(r2) ≤ sqrt(ϵ) * norm(b)
  end

  @testset "factorization precision" begin
    A = [
      1.7 0 0 0 0 0 0 0 0.13 0
      0 1.0 0 0 0.02 0 0 0 0 0.01
      0 0 1.5 0 0 0 0 0 0 0
      0 0 0 1.1 0 0 0 0 0 0
      0 0.02 0 0 2.6 0 0.16 0.09 0.52 0.53
      0 0 0 0 0 1.2 0 0 0 0
      0 0 0 0 0.16 0 1.3 0 0 0.56
      0 0 0 0 0.09 0 0 1.6 0.11 0
      0.13 0 0 0 0.52 0 0 0.11 1.4 0
      0 0.01 0 0 0.53 0 0.56 0 0 3.1
    ]
    b = [0.287, 0.22, 0.45, 0.44, 2.486, 0.72, 1.55, 1.424, 1.621, 3.759]
    ϵ = sqrt(eps(eltype(A)))
    m, n = size(A)

    LDL32 = ldl_analyze(A, Float32)
    @test typeof(LDL32) == LDLFactorizations.LDLFactorization{Float32, Int, Int, Int}
    @test eltype(LDL32) == Float32
    ldl_factorize!(A, LDL32)
    @test LDL32.__factorized

    LDL32 = ldl_analyze(Symmetric(triu(A), :U), Float32)
    @test typeof(LDL32) == LDLFactorizations.LDLFactorization{Float32, Int, Int, Int}
    @test eltype(LDL32) == Float32
    ldl_factorize!(Symmetric(triu(A), :U), LDL32)
    @test LDL32.__factorized

    LDL32 = ldl(A, Float32)
    @test typeof(LDL32) == LDLFactorizations.LDLFactorization{Float32, Int, Int, Int}
    @test eltype(LDL32) == Float32
    @test LDL32.__factorized

    LDL32 = ldl(Symmetric(triu(A), :U), Float32)
    @test typeof(LDL32) == LDLFactorizations.LDLFactorization{Float32, Int, Int, Int}
    @test eltype(LDL32) == Float32
    @test LDL32.__factorized

    b, c = rand(n), rand(n)
    d = LDL32 \ b
    @test eltype(d) == Float64
    ldiv!(c, LDL32, b)
    ldiv!(LDL32, b)
    @test eltype(b) == Float64
    @test all(b .== d .== c)

    mul!(c, LDL32, b)
    lmul!(LDL32, b)
    @test eltype(c) == eltype(b) == Float64
    @test all(b .== c)
  end

  @testset "Operations on LDL struct" begin
    A = rand(2, 2)
    b = A * ones(2)
    A = Symmetric(triu(A), :U)
    S = ldl_analyze(A)
    ldl_factorize!(A, S)
    x1 = S \ (-b)
    x2 = -S \ b
    @test x1 == x2
  end
end
