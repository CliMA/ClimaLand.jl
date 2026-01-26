#
using SciMLOperators, LinearAlgebra
using Random

using SciMLOperators: InvertibleOperator, ⊗
using FFTW

Random.seed!(0)
N = 8
K = 12

@testset "FunctionOperator FFTW test" begin
    n = 256
    L = 2π

    dx = L / n
    x = (-L / 2):dx:(L / 2 - dx) |> Array

    k = rfftfreq(n, 2π * n / L) |> Array
    m = length(k)
    P = plan_rfft(x)

    fwd(v, u, p, t) = P * v
    bwd(v, u, p, t) = P \ v

    fwd(w, v, u, p, t) = mul!(w, P, v)
    bwd(w, v, u, p, t) = ldiv!(w, P, v)

    ftr = FunctionOperator(fwd, x, im * k;
        T = ComplexF64, op_adjoint = bwd,
        op_inverse = bwd,
        op_adjoint_inverse = fwd, islinear = true)

    @test size(ftr) == (length(k), length(x))

    # derivative test
    ik = im * DiagonalOperator(k)
    Dx = ftr \ ik * ftr
    Dx = cache_operator(Dx, x)
    D2x = cache_operator(Dx * Dx, x)

    v = @. sin(5x)cos(7x)
    w = @. 5cos(5x)cos(7x) - 7sin(5x)sin(7x)
    w2x = @. 5(-5sin(5x)cos(7x) - 7cos(5x)sin(7x)) +
             -7(5cos(5x)sin(7x) + 7sin(5x)cos(7x))

    @test ≈(Dx * v, w; atol = 1e-8)
    @test ≈(D2x * v, w2x; atol = 1e-8)

    w2 = zero(w)
    @test ≈(mul!(w2, D2x, v), w2x; atol = 1e-8)
    w2 = zero(w)
    @test ≈(mul!(w2, Dx, v), w; atol = 1e-8)

    itr = inv(ftr)
    ftt = ftr'
    itt = itr'

    @test itr isa FunctionOperator
    @test ftt isa FunctionOperator
    @test itt isa FunctionOperator

    @test size(ftr) == (m, n)
    @test size(itr) == (n, m)
    @test size(ftt) == (n, m)
    @test size(itt) == (m, n)

    @test ftt.op == ftr.op_adjoint
    @test ftt.op_adjoint == ftr.op
    @test ftt.op_inverse == ftr.op_adjoint_inverse
    @test ftt.op_adjoint_inverse == ftr.op_inverse

    @test itr.op == ftr.op_inverse
    @test itr.op_adjoint == ftr.op_adjoint_inverse
    @test itr.op_inverse == ftr.op
    @test itr.op_adjoint_inverse == ftr.op_adjoint

    @test itt.op == ftr.op_adjoint_inverse
    @test itt.op_adjoint == ftr.op_inverse
    @test itt.op_inverse == ftr.op_adjoint
    @test itt.op_adjoint_inverse == ftr.op
end

@testset "Operator Algebra" begin
    N2 = N * N

    v = rand(N2, K)
    u = rand()
    p = rand()
    t = rand()

    A = rand(N, N)

    # Introduce update function for B
    B = MatrixOperator(zeros(N, N); update_func! = (A, u, p, t) -> (A .= p))

    # FunctionOp
    _C = rand(N, N) |> Symmetric
    f(v, u, p, t) = _C * v
    f(w, v, u, p, t) = mul!(w, _C, v)
    C = FunctionOperator(f, zeros(N); batch = true, issymmetric = true, p = p, u = u)

    # Introduce update function for D dependent on kwarg "matrix"
    D = MatrixOperator(zeros(N, N);
        update_func! = (A, u, p, t; matrix) -> (A .= p * t * matrix),
        accepted_kwargs = Val((:matrix,)))

    matrix = rand(N, N)
    diag = rand(N2)
    α = rand()
    β = rand()

    T1 = ⊗(A, B)
    T2 = ⊗(C, D)

    D1 = DiagonalOperator(zeros(N2); update_func! = (d, u, p, t) -> d .= p)
    D2 = DiagonalOperator(
        zeros(N2); update_func! = (d, u, p, t; diag) -> d .= p * t * diag,
        accepted_kwargs = Val((:diag,)))

    TT = [T1, T2]
    DD = Diagonal([D1, D2])

    op = TT' * DD * TT
    op = cache_operator(op, v)

    # Update operator
    @test_nowarn update_coefficients!(op, u, p, t; diag, matrix)

    # Form dense operator manually 
    dense_T1 = kron(A, p * ones(N, N))
    dense_T2 = kron(_C, (p * t) .* matrix)
    dense_DD = Diagonal(vcat(p * ones(N2), p * t * diag))
    dense_op = hcat(dense_T1', dense_T2') * dense_DD * vcat(dense_T1, dense_T2)

    # Test correctness of op
    @test op * v ≈ dense_op * v

    # Test consistency with three-arg mul!
    w = rand(N2, K)
    @test mul!(w, op, v) ≈ op * v

    # Test consistency with in-place five-arg mul!
    w = rand(N2, K)
    w2 = copy(w)
    @test mul!(w, op, v, α, β) ≈ α * (op * v) + β * w2

    # Create a fresh operator for each test
    op_fresh = TT' * DD * TT
    op_fresh = cache_operator(op_fresh, v)
    # Use in-place update directly in test
    result1 = similar(v)
    mul!(result1, op_fresh, v)
    update_coefficients!(op_fresh, u, p, t; diag, matrix)
    @test result1 ≈ dense_op * v
end

@testset "Resize! test" begin
    M1 = 4
    M2 = 12

    v = rand(N)
    v1 = rand(M1)
    v2 = rand(M2)

    f(v, u, p, t) = 2 * v
    f(w, v, u, p, t) = (copy!(w, v); lmul!(2, w))

    fi(v, u, p, t) = 0.5 * v
    fi(w, v, u, p, t) = (copy!(w, v); lmul!(0.5, w))

    F = FunctionOperator(f, v, v; islinear = true, op_inverse = fi, issymmetric = true)

    multest(L, v) = @test mul!(zero(v), L, v) ≈ L * v

    function multest(L::SciMLOperators.AdjointOperator, v)
        @test mul!(adjoint(zero(v)), adjoint(v), L) ≈ adjoint(v) * L
    end

    function multest(L::SciMLOperators.TransposedOperator, v)
        @test mul!(transpose(zero(v)), transpose(v), L) ≈ transpose(v) * L
    end

    function multest(L::SciMLOperators.InvertedOperator, v)
        @test ldiv!(zero(v), L, v) ≈ L \ v
    end

    for (L, LT) in ((F, FunctionOperator),
        (F + F, SciMLOperators.AddedOperator),
        (F * 2, SciMLOperators.ScaledOperator),
        (F ∘ F, SciMLOperators.ComposedOperator),
        (AffineOperator(F, F, v), AffineOperator),
        (SciMLOperators.AdjointOperator(F), SciMLOperators.AdjointOperator),
        (SciMLOperators.TransposedOperator(F), SciMLOperators.TransposedOperator),
        (SciMLOperators.InvertedOperator(F), SciMLOperators.InvertedOperator),
        (SciMLOperators.InvertibleOperator(F, F), SciMLOperators.InvertibleOperator))
        L = deepcopy(L)
        L = cache_operator(L, v)

        @test L isa LT
        @test size(L) == (N, N)
        multest(L, v)

        resize!(L, M1)
        @test size(L) == (M1, M1)
        multest(L, v1)

        resize!(L, M2)
        @test size(L) == (M2, M2)
        multest(L, v2)
    end

    # InvertedOperator
    # AffineOperator
    # FunctionOperator
end
#
