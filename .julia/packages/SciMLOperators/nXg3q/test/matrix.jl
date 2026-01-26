using SciMLOperators, LinearAlgebra
using Random
using Test

using SciMLOperators: InvertibleOperator, InvertedOperator, ⊗, AbstractSciMLOperator
using FFTW

Random.seed!(0)
N = 8
K = 19

@testset "MatrixOperator, InvertibleOperator" begin
    # Vectors for testing
    u = rand(N, K)  # Both update and action vector
    v = rand(N, K)  # Output/action vector
    w = zeros(N, K)  # Output vector

    p = nothing
    t = 0
    α = rand()
    β = rand()

    A = rand(N, N)
    At = A'

    AA = MatrixOperator(A)
    AAt = AA'

    @test AA isa MatrixOperator
    @test AAt isa MatrixOperator

    @test isconstant(AA)
    @test isconstant(AAt)

    @test issquare(AA)
    @test islinear(AA)

    FF = factorize(AA)
    FFt = FF'

    @test FF isa InvertibleOperator
    @test FFt isa InvertibleOperator

    @test isconstant(FF)
    @test isconstant(FFt)

    @test_throws MethodError resize!(AA, N)

    @test eachindex(A) === eachindex(AA)
    @test eachindex(A') === eachindex(AAt) === eachindex(MatrixOperator(At))

    @test A ≈ convert(AbstractMatrix, AA) ≈ convert(AbstractMatrix, FF)
    @test At ≈ convert(AbstractMatrix, AAt) ≈ convert(AbstractMatrix, FFt)

    @test A ≈ Matrix(AA) ≈ Matrix(FF)
    @test At ≈ Matrix(AAt) ≈ Matrix(FFt)

    # Test with new interface - same vector for update and action
    @test A * u ≈ AA(u, u, p, t)
    @test At * u ≈ AAt(u, u, p, t)

    # Test with different vectors for update and action
    @test A * v ≈ AA(v, u, p, t)
    @test At * v ≈ AAt(v, u, p, t)

    @test A \ u ≈ AA \ u ≈ FF \ u
    @test At \ u ≈ AAt \ u ≈ FFt \ u

    # Test in-place operations
    copy!(w, zeros(N, K))
    AA(w, v, u, p, t)
    @test w ≈ A * v

    # Test in-place with scaling
    copy!(w, rand(N, K))
    orig_w = copy(w)
    AA(w, v, u, p, t, α, β)
    @test w ≈ α * (A * v) + β * orig_w
end

@testset "InvertibleOperator test" begin
    # Vectors for testing
    u = rand(N, K)  # Update vector
    v = rand(N, K)  # Action vector
    w = zeros(N, K)  # Output vector

    p = nothing
    t = 0
    α = rand()
    β = rand()

    d = rand(N, K)
    D = DiagonalOperator(d)
    Di = DiagonalOperator(inv.(d)) |> InvertedOperator

    L = InvertibleOperator(D, Di)
    L = cache_operator(L, u)

    @test iscached(L)

    # Test with new interface
    @test L(v, u, p, t) ≈ d .* v
    @test L \ u ≈ d .\ u

    # Test in-place operations
    copy!(w, zeros(N, K))
    L(w, v, u, p, t)
    @test w ≈ d .* v

    # Test in-place with scaling
    copy!(w, rand(N, K))
    orig_w = copy(w)
    L(w, v, u, p, t, α, β)
    @test w ≈ α * (d .* v) + β * orig_w

    # Test division operations
    copy!(w, zeros(N, K))
    ldiv!(w, L, u)
    @test w ≈ d .\ u

    # Existing test for in-place ldiv!
    v = copy(u)
    ldiv!(L, v)
    @test v ≈ d .\ u
end

@testset "MatrixOperator update test" begin
    # Vectors for testing
    u = rand(N, K)  # Update vector
    v = rand(N, K)  # Action vector
    w = zeros(N, K)  # Output vector

    p = rand(N)
    t = rand()
    α = rand()
    β = rand()

    L = MatrixOperator(zeros(N, N);
        update_func = (A, u, p, t) -> p * p',
        update_func! = (A, u, p, t) -> A .= p * p')

    @test !isconstant(L)

    # Expected matrix after update
    A = p * p'

    # Test with new interface - same vector for update and action
    @test L(u, u, p, t) ≈ A * u

    # Test with different vectors for update and action
    @test L(v, u, p, t) ≈ A * v

    # Test in-place operation
    copy!(w, zeros(N, K))
    L(w, v, u, p, t)
    @test w ≈ A * v

    # Test in-place with scaling
    copy!(w, rand(N, K))
    orig_w = copy(w)
    L(w, v, u, p, t, α, β)
    @test w ≈ α * (A * v) + β * orig_w

    A = [-2.0 1 0 0 0
         1 -2 1 0 0
         0 1 -2 1 0
         0 0 1 -2 1
         0 0 0 1 -2]
    v = [3.0, 2.0, 1.0, 2.0, 3.0]
    opA = MatrixOperator(A)

    function update_function!(B, u, p, t)
        dt = p
        B .= A .* u + dt*I
    end

    u = Array(1:1.0:5);
    p = 0.1;
    t = 0.0
    opB = MatrixOperator(copy(A); update_func! = update_function!)

    function Bfunc!(w, v, u, p, t)
        dt = p
        w[1] = -(2*u[1]-dt)*v[1] + v[2]*u[1]
        for i in 2:4
            w[i] = v[i - 1]*u[i] - (2*u[i]-dt)*v[i] + v[i + 1]*u[i]
        end
        w[5] = v[4]*u[5] - (2*u[5]-dt)*v[5]
        nothing
    end

    function Bfunc!(v, u, p, t)
        w = zeros(5)
        Bfunc!(w, v, u, p, t)
        w
    end

    mfopB = FunctionOperator(Bfunc!, zeros(5), zeros(5); u, p, t, isconstant = false)

    @test iszero(opB(v, Array(2:1.0:6), 0.5, nothing) -
                 mfopB(v, Array(2:1.0:6), 0.5, nothing))
end

@testset "DiagonalOperator update test" begin
    # Vectors for testing
    u = rand(N, K)  # Update vector
    v = rand(N, K)  # Action vector
    w = zeros(N, K)  # Output vector

    p = rand(N)
    t = rand()
    α = rand()
    β = rand()

    D = DiagonalOperator(zeros(N);
        update_func = (diag, u, p, t) -> p * t,
        update_func! = (diag, u, p, t) -> diag .= p * t)

    @test !isconstant(D)
    @test issquare(D)
    @test islinear(D)

    # Expected result after update
    expected = (p * t) .* v

    # Test with new interface - different vectors for update and action
    @test D(v, u, p, t) ≈ expected

    # Test in-place operation
    copy!(w, zeros(N, K))
    D(w, v, u, p, t)
    @test w ≈ expected

    # Test in-place with scaling
    copy!(w, rand(N, K))
    orig_w = copy(w)
    D(w, v, u, p, t, α, β)
    @test w ≈ α * expected + β * orig_w
end

@testset "Batched Diagonal Operator" begin
    # Vectors for testing
    u = rand(N, K)  # Update vector
    v = rand(N, K)  # Action vector
    w = zeros(N, K)  # Output vector

    d = rand(N, K)
    α = rand()
    β = rand()
    p = nothing
    t = 0.0

    L = DiagonalOperator(d)
    @test isconstant(L)
    @test_throws MethodError resize!(L, N)

    @test issquare(L)
    @test islinear(L)

    # Test with new interface
    @test L(v, u, p, t) ≈ d .* v

    # Test in-place operation
    copy!(w, zeros(N, K))
    L(w, v, u, p, t)
    @test w ≈ d .* v

    # Test in-place with scaling
    copy!(w, rand(N, K))
    orig_w = copy(w)
    L(w, v, u, p, t, α, β)
    @test w ≈ α * (d .* v) + β * orig_w

    # Test division operations
    @test L \ u ≈ d .\ u

    copy!(w, zeros(N, K))
    ldiv!(w, L, u)
    @test w ≈ d .\ u

    # Existing test for in-place ldiv!
    v = copy(u)
    ldiv!(L, u)
    @test u ≈ d .\ v
end

@testset "Batched DiagonalOperator update test" begin
    # Vectors for testing
    u = rand(N, K)  # Update vector
    v = rand(N, K)  # Action vector
    w = zeros(N, K)  # Output vector

    d = zeros(N, K)
    p = rand(N, K)
    t = rand()

    D = DiagonalOperator(d;
        update_func = (diag, u, p, t) -> p * t,
        update_func! = (diag, u, p, t) -> diag .= p * t)

    @test !isconstant(D)
    @test issquare(D)
    @test islinear(D)

    # Expected result after update
    expected = (p * t) .* v

    # Test with new interface - different vectors for update and action
    @test D(v, u, p, t) ≈ expected

    # Test in-place operation
    copy!(w, zeros(N, K))
    D(w, v, u, p, t)
    @test w ≈ expected
end

@testset "AffineOperator" begin
    # Vectors for testing
    u = rand(N, K)  # Update vector
    v = rand(N, K)  # Action vector
    w = zeros(N, K)  # Output vector

    A = rand(N, N)
    B = rand(N, N)
    D = Diagonal(A)
    b = rand(N, K)
    α = rand()
    β = rand()
    p = nothing
    t = 0.0

    L = AffineOperator(MatrixOperator(A), MatrixOperator(B), b)
    @test isconstant(L)
    @test issquare(L)
    @test !islinear(L)

    # Test with new interface
    expected = A * v + B * b
    @test L(v, u, p, t) ≈ expected

    # Test in-place operation
    copy!(w, zeros(N, K))
    L(w, v, u, p, t)
    @test w ≈ expected

    # Test in-place with scaling
    copy!(w, rand(N, K))
    orig_w = copy(w)
    L(w, v, u, p, t, α, β)
    @test w ≈ α * expected + β * orig_w

    L = AffineOperator(MatrixOperator(D), MatrixOperator(B), b)
    @test issquare(L)
    @test !islinear(L)

    # Test division operations
    @test L \ u ≈ D \ (u - B * b)

    copy!(w, zeros(N, K))
    ldiv!(w, L, u)
    @test w ≈ D \ (u - B * b)

    # Existing test for in-place ldiv!
    v = copy(u)
    ldiv!(L, u)
    @test u ≈ D \ (v - B * b)

    L = AddVector(b)
    @test issquare(L)
    @test !islinear(L)

    # Test with new interface
    expected = v + b
    @test L(v, u, p, t) ≈ expected
    @test L \ u ≈ u - b

    # Test in-place operation
    copy!(w, zeros(N, K))
    L(w, v, u, p, t)
    @test w ≈ expected

    # Test in-place with scaling
    copy!(w, rand(N, K))
    orig_w = copy(w)
    L(w, v, u, p, t, α, β)
    @test w ≈ α * expected + β * orig_w

    # Test division operations
    copy!(w, zeros(N, K))
    ldiv!(w, L, u)
    @test w ≈ u - b

    # Existing test for in-place ldiv!
    v = copy(u)
    ldiv!(L, u)
    @test u ≈ v - b

    L = AddVector(MatrixOperator(B), b)
    @test issquare(L)
    @test !islinear(L)

    # Test with new interface
    expected = v + B * b
    @test L(v, u, p, t) ≈ expected
    @test L \ u ≈ u - B * b

    # Test in-place operation
    copy!(w, zeros(N, K))
    L(w, v, u, p, t)
    @test w ≈ expected

    # Test in-place with scaling
    copy!(w, rand(N, K))
    orig_w = copy(w)
    L(w, v, u, p, t, α, β)
    @test w ≈ α * expected + β * orig_w

    # Test division operations
    copy!(w, zeros(N, K))
    ldiv!(w, L, u)
    @test w ≈ u - B * b

    # Existing test for in-place ldiv!
    v = copy(u)
    ldiv!(L, u)
    @test u ≈ v - B * b
end

@testset "AffineOperator update test" begin
    # Vectors for testing
    u = rand(N, K)  # Update vector
    v = rand(N, K)  # Action vector
    w = zeros(N, K)  # Output vector

    A = rand(N, N)
    B = rand(N, N)
    p = rand(N, K)
    t = rand()
    α = rand()
    β = rand()

    L = AffineOperator(A, B, zeros(N, K);
        update_func = (b, u, p, t) -> p * t,
        update_func! = (b, u, p, t) -> b .= p * t)

    @test !isconstant(L)

    # Expected updated bias and result
    b = p * t
    expected = A * v + B * b

    # Test with new interface - different vectors for update and action
    @test L(v, u, p, t) ≈ expected

    # Test in-place operation
    copy!(w, zeros(N, K))
    L(w, v, u, p, t)
    @test w ≈ expected

    # Test in-place with scaling
    copy!(w, rand(N, K))
    orig_w = copy(w)
    L(w, v, u, p, t, α, β)
    @test w ≈ α * expected + β * orig_w
end

@testset "TensorProductOperator, square = $square" for square in [false, true]
    m1, n1 = 3, 5
    m2, n2 = 7, 11
    m3, n3 = 13, 17

    if square
        n1, n2, n3 = m1, m2, m3
    end

    A = rand(m1, n1)
    B = rand(m2, n2)
    C = rand(m3, n3)
    α = rand()
    β = rand()
    p = nothing
    t = 0.0

    AB = kron(A, B)
    ABC = kron(A, B, C)

    # test Base.kron overload
    # ensure kron(mat, mat) is not a TensorProductOperator
    @test !isa(AB, AbstractSciMLOperator)
    @test !isa(ABC, AbstractSciMLOperator)

    # test Base.kron overload
    _A = rand(N, N)
    @test kron(_A, MatrixOperator(_A)) isa TensorProductOperator
    @test kron(MatrixOperator(_A), _A) isa TensorProductOperator

    @test kron(MatrixOperator(_A), MatrixOperator(_A)) isa TensorProductOperator

    # Inputs/Update vectors
    u2 = rand(n1 * n2, K)
    u3 = rand(n1 * n2 * n3, K)

    # Action vectors (same as update vectors initially)
    v2 = copy(u2)
    v3 = copy(u3)

    # Output vectors
    w2 = zeros(m1 * m2, K)
    w3 = zeros(m1 * m2 * m3, K)

    opAB = TensorProductOperator(A, B)
    opABC = TensorProductOperator(A, B, C)

    @test opAB isa TensorProductOperator
    @test opABC isa TensorProductOperator

    @test isconstant(opAB)
    @test isconstant(opABC)

    @test islinear(opAB)
    @test islinear(opABC)

    if square
        @test issquare(opAB)
        @test issquare(opABC)
    else
        @test !issquare(opAB)
        @test !issquare(opABC)
    end

    @test AB ≈ convert(AbstractMatrix, opAB)
    @test ABC ≈ convert(AbstractMatrix, opABC)

    # factorization tests
    opAB_F = factorize(opAB)
    opABC_F = factorize(opABC)

    @test isconstant(opAB_F)
    @test isconstant(opABC_F)

    @test opAB_F isa TensorProductOperator
    @test opABC_F isa TensorProductOperator

    @test AB ≈ convert(AbstractMatrix, opAB_F)
    @test ABC ≈ convert(AbstractMatrix, opABC_F)

    # Test with new interface
    @test AB * v2 ≈ opAB(v2, u2, p, t)
    @test ABC * v3 ≈ opABC(v3, u3, p, t)

    @test AB \ w2 ≈ opAB \ w2
    @test AB \ w2 ≈ opAB_F \ w2
    @test ABC \ w3 ≈ opABC \ w3
    @test ABC \ w3 ≈ opABC_F \ w3

    @test !iscached(opAB)
    @test !iscached(opABC)

    @test !iscached(opAB_F)
    @test !iscached(opABC_F)

    opAB = cache_operator(opAB, u2)
    opABC = cache_operator(opABC, u3)

    opAB_F = cache_operator(opAB_F, u2)
    opABC_F = cache_operator(opABC_F, u3)

    @test iscached(opAB)
    @test iscached(opABC)

    @test iscached(opAB_F)
    @test iscached(opABC_F)

    N2 = n1 * n2
    N3 = n1 * n2 * n3
    M2 = m1 * m2
    M3 = m1 * m2 * m3

    # Test in-place operations with new interface
    v2 = rand(n1 * n2, K)  # Action vector
    w2 = zeros(M2, K)      # Output vector
    opAB(w2, v2, u2, p, t)
    @test w2 ≈ AB * v2

    v3 = rand(n1 * n2 * n3, K)  # Action vector
    w3 = zeros(M3, K)           # Output vector
    opABC(w3, v3, u3, p, t)
    @test w3 ≈ ABC * v3

    # Test in-place with scaling
    v2 = rand(n1 * n2, K)    # Action vector
    w2 = rand(M2, K)         # Output vector
    orig_w2 = copy(w2)
    opAB(w2, v2, u2, p, t, α, β)
    @test w2 ≈ α * AB * v2 + β * orig_w2

    v3 = rand(n1 * n2 * n3, K)  # Action vector
    w3 = rand(M3, K)            # Output vector
    orig_w3 = copy(w3)
    opABC(w3, v3, u3, p, t, α, β)
    @test w3 ≈ α * ABC * v3 + β * orig_w3

    if square
        # Test division operations with new interface
        v2 = rand(M2, K)     # Action vector (size of output space)
        u2 = rand(N2, K)     # Update vector (size of input space)
        w2 = zeros(N2, K)    # Output vector (size of input space)

        # ldiv! with new interface
        ldiv!(w2, opAB_F, v2)
        @test w2 ≈ AB \ v2

        v3 = rand(M3, K)     # Action vector
        u3 = rand(N3, K)     # Update vector
        w3 = zeros(N3, K)    # Output vector
        ldiv!(w3, opABC_F, v3)
        @test w3 ≈ ABC \ v3

        # In-place ldiv! (original style)
        v2 = rand(M2, K)
        u2 = copy(v2)
        ldiv!(opAB_F, v2)
        @test v2 ≈ AB \ u2

        v3 = rand(M3, K)
        u3 = copy(v3)
        ldiv!(opABC_F, v3)
        @test v3 ≈ ABC \ u3
    else # TODO
        v2 = rand(M2, K)     # Action vector
        u2 = rand(N2, K)     # Update vector
        w2 = zeros(N2, K)    # Output vector

        if VERSION < v"1.9-"
            @test_broken ldiv!(w2, opAB_F, v2) ≈ AB \ v2
        else
            @test ldiv!(w2, opAB_F, v2) ≈ AB \ v2
        end

        v3 = rand(M3, K)     # Action vector
        u3 = rand(N3, K)     # Update vector
        w3 = zeros(N3, K)    # Output vector
        @test_broken ldiv!(w3, opABC_F, v3) ≈ ABC \ v3 # errors
    end
end
