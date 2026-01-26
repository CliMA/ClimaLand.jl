using SciMLOperators
using SciMLOperators: AbstractSciMLScalarOperator,
                      ComposedScalarOperator,
                      AddedScalarOperator,
                      InvertedScalarOperator,
                      IdentityOperator,
                      AddedOperator,
                      ScaledOperator

using LinearAlgebra, Random, Test

Random.seed!(0)
N = 8
K = 12

@testset "ScalarOperator Basic Operations" begin
    a = rand()
    b = rand()
    x = rand()
    α = ScalarOperator(x)
    u = rand(N, K)

    @test α isa ScalarOperator
    @test iscached(α)
    @test issquare(α)
    @test islinear(α)

    @test convert(Float32, α) isa Float32
    @test convert(ScalarOperator, a) isa ScalarOperator

    @test size(α) == ()
    @test isconstant(α)

    # Original lmul!/rmul! tests
    v = copy(u)
    @test lmul!(α, u) ≈ v * x
    v = copy(u)
    @test rmul!(u, α) ≈ x * v

    # Original mul!/ldiv! tests
    v = rand(N, K)
    @test mul!(v, α, u) ≈ u * x
    v = rand(N, K)
    w = copy(v)
    @test mul!(v, α, u, a, b) ≈ a * (x * u) + b * w

    v = rand(N, K)
    @test ldiv!(v, α, u) ≈ u / x
    w = copy(u)
    @test ldiv!(α, u) ≈ w / x

    # Original axpy! test
    X = rand(N, K)
    Y = rand(N, K)
    Z = copy(Y)
    a_scalar = rand()
    aa = ScalarOperator(a_scalar)
    @test axpy!(aa, X, Y) ≈ a_scalar * X + Z

    # Tests with the new interface
    v = copy(u)  # Action vector
    w = zeros(N, K)  # Output vector

    # Test with new interface
    result = α(v, u, nothing, 0.0)
    @test result ≈ v * x

    # Test in-place operations
    α(w, v, u, nothing, 0.0)
    @test w ≈ v * x

    # Test in-place operations with scaling
    orig_w = rand(N, K)
    copy!(w, orig_w)
    α(w, v, u, nothing, 0.0, a, b)
    @test w ≈ a * (x * v) + b * orig_w

    # Operations with UniformScalar
    α = ScalarOperator(1)
    @test α * I == α
    @test α + I isa AddedScalarOperator
end

@testset "ScalarOperator Combinations" begin
    x = rand()
    α = ScalarOperator(x)
    u = rand(N, K)  # Update vector
    v = rand(N, K)  # Action vector

    # Test scalar operator combinations
    β = α + α
    @test β isa AddedScalarOperator
    @test β(v, u, nothing, 0.0) ≈ x * v + x * v
    @test β * u ≈ x * u + x * u  # Original style test
    @inferred convert(Float32, β)
    @test convert(Number, β) ≈ x + x

    β = α * α
    @test β isa ComposedScalarOperator
    @test β(v, u, nothing, 0.0) ≈ x * x * v
    @test β * u ≈ x * x * u  # Original style test
    @inferred convert(Float32, β)
    @test convert(Number, β) ≈ x * x

    β = inv(α)
    @test β isa InvertedScalarOperator
    @test β(v, u, nothing, 0.0) ≈ (1 / x) * v
    @test β * u ≈ (1 / x) * u  # Original style test
    @inferred convert(Float32, β)
    @test convert(Number, β) ≈ 1 / x

    β = α * inv(α)
    @test β isa ComposedScalarOperator
    @test β * u ≈ u
    @inferred convert(Float32, β)
    @test convert(Number, β) ≈ true

    β = α / α
    @test β isa ComposedScalarOperator
    @test β * u ≈ u
    @inferred convert(Float32, β)
    @test convert(Number, β) ≈ true

    # Test combination with other operators
    for op in (MatrixOperator(rand(N, N)), SciMLOperators.IdentityOperator(N))
        @test α + op isa SciMLOperators.AddedOperator
        @test (α + op) * u ≈ x * u + op * u

        L = α + op
        @test L isa SciMLOperators.AddedOperator
        @test L(v, u, nothing, 0.0) ≈ x * v + op * v

        @test α * op isa SciMLOperators.ScaledOperator
        @test (α * op) * u ≈ x * (op * u)

        L = α * op
        @test L isa SciMLOperators.ScaledOperator
        @test L(v, u, nothing, 0.0) ≈ x * (op * v)

        # Division tests from original
        @test all(map(T -> (T isa SciMLOperators.ScaledOperator),
            (α / op, op / α, op \ α, α \ op)))
        @test (α / op) * u ≈ (op \ α) * u ≈ α * (op \ u)
        @test (op / α) * u ≈ (α \ op) * u ≈ 1 / α * op * u
    end

    # Test for ComposedScalarOperator nesting (from original)
    α_new = ScalarOperator(rand())
    L = α_new * (α_new * α_new) * α_new
    @test L isa ComposedScalarOperator
    for op in L.ops
        @test !isa(op, ComposedScalarOperator)
    end
end

@testset "ScalarOperator scalar argument test" begin
    a = rand()
    u = rand()  # Update scalar
    v = rand()  # Action scalar
    p = nothing
    t = 0.0

    α = ScalarOperator(a)
    @test_throws MethodError α(u, p, t) ≈ u * a  # Original style
    @test α(v, u, p, t) ≈ v * a  # New interface
    @test_throws ArgumentError α(v, u, p, t, 1, 2)  # Keep error test
end

@testset "ScalarOperator update test" begin
    u = ones(N, K)  # Update vector
    v = ones(N, K)  # Action vector
    w = zeros(N, K)  # Output vector
    p = 2.0
    t = 4.0

    c = rand()
    d = rand()

    α = ScalarOperator(0.0; update_func = (a, u, p, t) -> p)
    β = ScalarOperator(0.0; update_func = (a, u, p, t) -> t)

    @test !isconstant(α)
    @test !isconstant(β)

    @test convert(Float32, α) isa Float32
    @test convert(Float32, β) isa Float32

    @test convert(Number, α) ≈ 0.0
    @test convert(Number, β) ≈ 0.0

    # Test update_coefficients
    update_coefficients!(α, u, p, t)
    update_coefficients!(β, u, p, t)

    @test convert(Number, α) ≈ p
    @test convert(Number, β) ≈ t

    # Original style tests
    @test_throws MethodError α(u, p, t) ≈ p * u
    @test_throws MethodError β(u, p, t) ≈ t * u

    # Tests with new interface
    @test α(v, u, p, t) ≈ p * v
    @test β(v, u, p, t) ≈ t * v

    # Test in-place with scaling
    orig_w = rand(N, K)
    copy!(w, orig_w)
    α(w, v, u, p, t, c, d)
    @test w ≈ c * p * v + d * orig_w

    # Retain original test with random vectors
    v_rand = rand(N, K)
    @test α(v_rand, u, p, t) ≈ p * v_rand
    v_rand = rand(N, K)
    w_rand = copy(v_rand)
    @test_broken α(v_rand, u, p, t, c, d) ≈ c * p * u + d * w_rand

    # Test operator combinations
    num = α + 2 / β * 3 - 4
    val = p + 2 / t * 3 - 4

    @test convert(Number, num) ≈ val

    # Test with keyword arguments
    γ = ScalarOperator(0.0; update_func = (args...; dtgamma) -> dtgamma,
        accepted_kwargs = Val((:dtgamma,)))

    dtgamma = rand()
    # Original tests
    @test_throws MethodError γ(u, p, t; dtgamma) ≈ dtgamma * u

    # New interface tests
    @test γ(v, u, p, t; dtgamma) ≈ dtgamma * v

    # In-place test with keywords
    w_test = zeros(N, K)
    γ(w_test, v, u, p, t; dtgamma)
    @test w_test ≈ dtgamma * v

    γ_added = γ + α
    # Original tests
    @test_throws MethodError γ_added(u, p, t; dtgamma) ≈ (dtgamma + p) * u

    # New interface tests
    @test γ_added(v, u, p, t; dtgamma) ≈ (dtgamma + p) * v

    # In-place test with keywords for combined operator
    w_test = zeros(N, K)
    γ_added(w_test, v, u, p, t; dtgamma)
    @test w_test ≈ (dtgamma + p) * v

    # In-place test with scaling and keywords
    w_test = rand(N, K)
    w_orig = copy(w_test)
    γ_added(w_test, v, u, p, t, c, d; dtgamma)
    @test w_test ≈ c * (dtgamma + p) * v + d * w_orig
end
