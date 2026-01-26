using SciMLOperators, LinearAlgebra
using Random
using Test

using SciMLOperators: ⊗

Random.seed!(0)
N = 8
K = 12
NK = N * K

@testset "(Unbatched) FunctionOperator ND array" begin
    N1, N2, N3 = 3, 4, 5
    M1, M2, M3 = 4, 5, 6

    u = nothing
    p = nothing
    t = 0.0
    α = rand()
    β = rand()

    for (sz_in, sz_out) in (((N1, N2, N3), (N1, N2, N3)), # equal size
        ((N1, N2, N3), (M1, M2, M3)))
        N = prod(sz_in)
        M = prod(sz_out)

        A = rand(M, N)
        u = nothing
        v = rand(sz_in...) # action vector 
        w = rand(sz_out...) # output vector for in-place tests

        _mul(A, v) = reshape(A * vec(v), sz_out)
        f(v, u, p, t) = _mul(A, v)
        f(w, v, u, p, t) = (mul!(vec(w), A, vec(v)); w)

        kw = (;) # FunctionOp kwargs

        if sz_in == sz_out
            F = lu(A)
            _div(A, v) = reshape(A \ vec(v), sz_in)
            fi(v, u, p, t) = _div(A, v)
            fi(w, v, u, p, t) = (ldiv!(vec(w), F, vec(v)); w)

            kw = (; op_inverse = fi)
        end

        L = FunctionOperator(f, v, w; kw...)
        L = cache_operator(L, v)

        # test with ND-arrays and new interface
        @test _mul(A, v) ≈ L(v, u, p, t) ≈ L * v ≈ mul!(zero(w), L, v)
        @test α * _mul(A, v) + β * w ≈ mul!(copy(w), L, v, α, β)

        # Test with different update and action vectors
        action_vec = rand(sz_in...)
        @test _mul(A, action_vec) ≈ L(action_vec, u, p, t)

        if sz_in == sz_out
            @test _div(A, w) ≈ L \ w ≈ ldiv!(zero(v), L, w) ≈ ldiv!(L, copy(w))
        end

        # test with vec(Array)
        @test vec(_mul(A, v)) ≈ L(vec(v), u, p, t) ≈ L * vec(v) ≈
              mul!(vec(zero(w)), L, vec(v))
        @test vec(α * _mul(A, v) + β * w) ≈ mul!(vec(copy(w)), L, vec(v), α, β)

        if sz_in == sz_out
            @test vec(_div(A, w)) ≈ L \ vec(w) ≈ ldiv!(vec(zero(v)), L, vec(w)) ≈
                  ldiv!(L, vec(copy(w)))
        end

        # Test in-place with different update and action vectors
        output_vec = zeros(sz_out...)
        L(output_vec, action_vec, u, p, t)
        @test output_vec ≈ _mul(A, action_vec)

        # Test in-place with scaling
        output_vec = rand(sz_out...)
        orig_output = copy(output_vec)
        L(output_vec, action_vec, u, p, t, α, β)
        @test output_vec ≈ α * _mul(A, action_vec) + β * orig_output

        @test_throws DimensionMismatch mul!(vec(w), L, v)
        @test_throws DimensionMismatch mul!(w, L, vec(v))
    end
end

@testset "(Unbatched) FunctionOperator" begin
    v = rand(N, K) # action vector
    w = zeros(N, K) # Output vector
    u = nothing
    p = nothing
    t = 0.0
    α = rand()
    β = rand()

    _mul(A, v) = reshape(A * vec(v), N, K)
    _div(A, v) = reshape(A \ vec(v), N, K)

    A = rand(NK, NK) |> Symmetric
    F = lu(A)
    Ai = inv(A)

    f1(v, u, p, t) = _mul(A, v)
    f1i(v, u, p, t) = _div(A, v)

    f2(w, v, u, p, t) = (mul!(vec(w), A, vec(v)); w)
    f2(w, v, u, p, t, α, β) = (mul!(vec(w), A, vec(v), α, β); w)
    f2i(w, v, u, p, t) = (ldiv!(vec(w), F, vec(v)); w)
    f2i(w, v, u, p, t, α, β) = (mul!(vec(w), Ai, vec(v), α, β); w)

    # out of place
    op1 = FunctionOperator(f1, v; op_inverse = f1i, ifcache = false, islinear = true,
        opnorm = true,
        issymmetric = true,
        ishermitian = true,
        isposdef = true)

    # in place
    op2 = FunctionOperator(f2, v, w; op_inverse = f2i, ifcache = false, islinear = true,
        opnorm = true,
        issymmetric = true,
        ishermitian = true,
        isposdef = true)

    # Test traits
    @test issquare(op1)
    @test issquare(op2)
    @test islinear(op1)
    @test islinear(op2)
    @test op1' === op1

    # Test operator properties
    @test size(op1) == (NK, NK)
    @test has_adjoint(op1)
    @test has_mul(op1)
    @test has_mul!(op1)
    @test has_ldiv(op1)
    @test !has_ldiv!(op1)

    @test size(op2) == (NK, NK)
    @test has_adjoint(op2)
    @test has_mul(op2)
    @test has_mul!(op2)
    @test has_ldiv(op2)
    @test has_ldiv!(op2)

    @test !iscached(op1)
    @test !iscached(op2)
    @test !op1.traits.has_mul5
    @test op2.traits.has_mul5

    # 5-arg mul! (w/o cache)
    v = rand(N, K)
    w = copy(v)
    @test α * _mul(A, v) + β * w ≈ mul!(w, op2, v, α, β)

    # Create test vectors for new interface
    action_vec = rand(N, K)  # Action vector
    result_vec = zeros(N, K) # Result vector

    # Cache operators
    op1 = cache_operator(op1, v)
    op2 = cache_operator(op2, v)

    @test iscached(op1)
    @test iscached(op2)

    # Test standard operator operations (from original test)
    w = rand(N, K)
    @test _mul(A, v) ≈ op1 * v ≈ mul!(w, op2, v) ≈ mul!(w, op1, v)
    w = rand(N, K)
    @test _mul(A, v) ≈ op1(v, u, p, t) ≈ op2(v, u, p, t)
    v = rand(N, K)
    w = copy(v)
    @test α * _mul(A, v) + β * w ≈ mul!(w, op2, v, α, β)

    w = rand(N, K)
    @test _div(A, w) ≈ op1 \ w ≈ ldiv!(v, op2, w)
    w = copy(v)
    @test _div(A, w) ≈ ldiv!(op2, w)

    # Test with new interface - out of place
    @test _mul(A, action_vec) ≈ op1(action_vec, u, p, t)

    # Test with new interface - in place
    op2(result_vec, action_vec, u, p, t)
    @test result_vec ≈ _mul(A, action_vec)

    # Test in-place with scaling
    result_vec = rand(N, K)
    orig_result = copy(result_vec)
    op2(result_vec, action_vec, u, p, t, α, β)
    @test result_vec ≈ α * _mul(A, action_vec) + β * orig_result

    # Test inverse operations with new interface
    inv_result = zeros(N, K)
    @test _div(A, action_vec) ≈ op1 \ action_vec
    ldiv!(inv_result, op2, action_vec)
    @test inv_result ≈ _div(A, action_vec)
end

@testset "Batched FunctionOperator" begin
    v = rand(N, K)
    u = nothing
    p = nothing
    t = 0.0
    α = rand()
    β = rand()

    A = rand(N, N) |> Symmetric
    F = lu(A)
    Ai = inv(A)

    f1(v, u, p, t) = A * v
    f1i(v, u, p, t) = A \ v

    f2(w, v, u, p, t) = mul!(w, A, v)
    f2(w, v, u, p, t, α, β) = mul!(w, A, v, α, β)
    f2i(w, v, u, p, t) = ldiv!(w, F, v)
    f2i(w, v, u, p, t, α, β) = mul!(w, Ai, v, α, β)

    # out of place
    op1 = FunctionOperator(f1, v, A * v; op_inverse = f1i, ifcache = false,
        batch = true,
        islinear = true,
        opnorm = true,
        issymmetric = true,
        ishermitian = true,
        isposdef = true)

    # in place
    op2 = FunctionOperator(f2, v, A * v; op_inverse = f2i, ifcache = false,
        batch = true,
        islinear = true,
        opnorm = true,
        issymmetric = true,
        ishermitian = true,
        isposdef = true)

    @test issquare(op1)
    @test issquare(op2)

    @test islinear(op1)
    @test islinear(op2)

    @test op1' === op1

    @test size(op1) == (N, N)
    @test has_adjoint(op1)
    @test has_mul(op1)
    @test has_mul!(op1)
    @test has_ldiv(op1)
    @test !has_ldiv!(op1)

    @test size(op2) == (N, N)
    @test has_adjoint(op2)
    @test has_mul(op2)
    @test has_mul!(op2)
    @test has_ldiv(op2)
    @test has_ldiv!(op2)

    @test !iscached(op1)
    @test !iscached(op2)

    @test !op1.traits.has_mul5
    @test op2.traits.has_mul5

    # 5-arg mul! (w/o cache)
    v = rand(N, K)
    w = copy(v)
    @test α * *(A, v) + β * w ≈ mul!(w, op2, v, α, β)

    op1 = cache_operator(op1, v)
    op2 = cache_operator(op2, v)

    @test iscached(op1)
    @test iscached(op2)

    v = rand(N, K)
    @test *(A, v) ≈ op1 * v ≈ mul!(w, op2, v)

    # Test with new interface
    v = rand(N, K)
    @test *(A, v) ≈ op1(w, v, u, p, t) ≈ op2(w, v, u, p, t)

    v = rand(N, K)
    w = copy(v)
    @test α * *(A, v) + β * w ≈ mul!(w, op2, v, α, β)

    # Test old style calls
    w = rand(N, K)
    @test \(A, w) ≈ op1 \ w ≈ ldiv!(v, op2, w)
    w = copy(v)
    @test \(A, w) ≈ ldiv!(op2, w)

    # Test new interface ldiv
    w = rand(N, K)
    ldiv_result = zeros(N, K)
    ldiv!(ldiv_result, op2, w)
    @test ldiv_result ≈ A \ w
end

@testset "FunctionOperator update test" begin
    u = rand(N, N)  # Update vector
    v = rand(N, K)  # Action vector
    w = zeros(N, K) # Result vector
    p = rand(N)
    t = rand()
    scale = rand()

    # Accept a kwarg "scale" in operator action
    f(w, v, u, p, t; scale = 1.0) = mul!(w, Diagonal(u * p * t * scale), v)
    f(v, u, p, t; scale = 1.0) = Diagonal(u * p * t * scale) * v

    # Test with both tuple and Val forms of accepted_kwargs
    for acc_kw in ((:scale,), Val((:scale,)))
        # Function operator with keyword arguments
        L = FunctionOperator(f, v, w;
            u = u,
            p = zero(p),
            t = zero(t),
            batch = true,
            accepted_kwargs = acc_kw,
            scale = 1.0)

        @test_throws ArgumentError FunctionOperator(
            f, v, w; u = u, p = zero(p), t = zero(t), batch = true,
            accepted_kwargs = acc_kw)

        @test size(L) == (N, N)

        # Expected result with scaling
        A = Diagonal(u * p * t * scale)
        expected = A * v
        ans = u * p .* t .* scale

        # Test with new interface
        @test L(v, u, p, t; scale) ≈ expected

        # Test in-place with new interface
        copy!(w, zeros(N, K))
        L(w, v, u, p, t; scale)
        @test w ≈ expected

        # Test in-place with scaling
        copy!(w, rand(N, K))
        orig_w = copy(w)
        α_val = rand()
        β_val = rand()
        L(w, v, u, p, t, α_val, β_val; scale)
        @test w ≈ α_val * expected + β_val * orig_w

        # Test that outputs aren't accidentally mutated
        v1 = rand(N, K)
        v2 = rand(N, K)
        w1 = rand(N, K)
        w2 = rand(N, K)

        # Expected results with different vectors
        result1 = A * v1
        result2 = A * v2

        # Test output consistency
        w1 = zeros(N, K)
        w2 = zeros(N, K)

        L(w1, v1, u, p, t; scale)
        L(w2, v2, u, p, t; scale)

        @test w1 ≈ result1
        @test w2 ≈ result2

        # Test matrix-vector multiplication
        w1 = L * v1
        @test w1 ≈ A * v1
        w2 = L * v2
        @test w2 ≈ A * v2
        @test w1 ≈ A * v1  # Check v1 hasn't changed
        @test w1 + w2 ≈ A * (v1 + v2)

        # Test in-place matrix-vector multiplication
        v1 .= 0.0
        v2 .= 0.0

        mul!(w1, L, v1)
        @test w1 ≈ A * v1
        mul!(w2, L, v2)
        @test w2 ≈ A * v2
        @test w1 ≈ A * v1
        @test w1 + w2 ≈ A * (v1 + v2)

        # Test scaling
        v1 = rand(N, K)
        w1 = copy(v1)
        v2 = rand(N, K)
        w2 = copy(v2)
        a1, a2, b1, b2 = rand(4)

        res = copy(w1)
        mul!(res, L, v1, a1, b1)
        @test res ≈ a1 * A * v1 + b1 * w1
        res2 = copy(w2)
        mul!(res2, L, v2, a2, b2)
        @test res2 ≈ a2 * A * v2 + b2 * w2
        @test res ≈ a1 * A * v1 + b1 * w1
        @test res + res2 ≈ (a1 * A * v1 + b1 * w1) + (a2 * A * v2 + b2 * w2)
    end
end
