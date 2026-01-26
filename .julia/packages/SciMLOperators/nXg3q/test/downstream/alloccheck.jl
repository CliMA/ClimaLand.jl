using SciMLOperators, Random, SparseArrays, Test, LinearAlgebra
using SciMLOperators: IdentityOperator,
                      NullOperator,
                      ScaledOperator,
                      AddedOperator,
                      ComposedOperator,
                      cache_operator

function apply_op!(H, w, v, u, p, t)
    H(w, v, u, p, t)
    return nothing
end

test_apply_noalloc(H, w, v, u, p, t) = @test (@allocations apply_op!(H, w, v, u, p, t)) == 0

@testset "Allocations Check" begin
    Random.seed!(0)
    N = 8
    K = 12
    A = rand(N, N) |> MatrixOperator
    B = rand(N, N) |> MatrixOperator
    u = rand(N, K)       # Update vector
    v = rand(N, K)       # Action vector
    w = zeros(N, K)      # Output vector
    p = ()
    t = 0
    op = AddedOperator(A, B)

    apply_op!(op, w, v, u, p, t) # Warm up
    if VERSION < v"1.10" || VERSION >= v"1.11"
        test_apply_noalloc(op, w, v, u, p, t)
    else
        # Julia 1.10 has a known allocation issue with AddedOperator
        @test (@allocations apply_op!(op, w, v, u, p, t)) == 2
    end

    for T in (Float32, Float64, ComplexF32, ComplexF64)
        N = 100
        A1_sparse = MatrixOperator(sprand(T, N, N, 5 / N))
        A2_sparse = MatrixOperator(sprand(T, N, N, 5 / N))
        A3_sparse = MatrixOperator(sprand(T, N, N, 5 / N))

        A1_dense = MatrixOperator(rand(T, N, N))
        A2_dense = MatrixOperator(rand(T, N, N))
        A3_dense = MatrixOperator(rand(T, N, N))

        coeff1(a, u, p, t) = sin(p.ω * t)
        coeff2(a, u, p, t) = cos(p.ω * t)
        coeff3(a, u, p, t) = sin(p.ω * t) * cos(p.ω * t)

        c1 = ScalarOperator(rand(T), coeff1)
        c2 = ScalarOperator(rand(T), coeff2)
        c3 = ScalarOperator(rand(T), coeff3)

        H_sparse = c1 * A1_sparse + c2 * A2_sparse + c3 * A3_sparse
        H_dense = c1 * A1_dense + c2 * A2_dense + c3 * A3_dense

        u = rand(T, N)
        v = rand(T, N)
        w = similar(u)
        p = (ω = 0.1,)
        t = 0.1

        apply_op!(H_sparse, w, v, u, p, t) # Warm up
        apply_op!(H_dense, w, v, u, p, t) # Warm up
        test_apply_noalloc(H_sparse, w, v, u, p, t)
        test_apply_noalloc(H_dense, w, v, u, p, t)
    end

    # Test ComposedOperator allocations (PR #316)
    # Before the fix, tuple splatting caused many allocations.
    # After the fix, we should have minimal allocations (Julia 1.11 has 1, earlier versions have 0).
    @testset "ComposedOperator minimal allocations" begin
        N = 100

        # Create operators for composition
        A1 = MatrixOperator(rand(N, N))
        A2 = MatrixOperator(rand(N, N))
        A3 = MatrixOperator(rand(N, N))

        # Create ComposedOperator
        L = A1 * A2 * A3

        # Set up cache
        v = rand(N)
        w = similar(v)
        L = cache_operator(L, v)

        u = rand(N)
        p = nothing
        t = 0.0

        # Warm up
        mul!(w, L, v)
        L(w, v, u, p, t)

        # Test mul! - should have minimal allocations
        # Julia 1.11 has a known minor allocation issue (1 allocation)
        # Earlier versions should have 0 allocations
        allocs_mul = @allocations mul!(w, L, v)
        @test allocs_mul <= 1

        # Test operator call - should have minimal allocations
        allocs_call = @allocations L(w, v, u, p, t)
        @test allocs_call <= 1

        # Test with matrices
        K = 5
        V = rand(N, K)
        W = similar(V)
        L_mat = cache_operator(A1 * A2 * A3, V)

        # Warm up
        mul!(W, L_mat, V)
        L_mat(W, V, u, p, t)

        # Test with matrices - should have minimal allocations
        allocs_mul_mat = @allocations mul!(W, L_mat, V)
        @test allocs_mul_mat <= 1

        allocs_call_mat = @allocations L_mat(W, V, u, p, t)
        @test allocs_call_mat <= 1
    end

    # Test accepted_kwargs allocations (PR #313)
    # With Val(tuple), kwarg filtering should be compile-time with minimal allocations
    @testset "accepted_kwargs with Val" begin
        N = 50

        # Create a MatrixOperator with accepted_kwargs using Val for compile-time filtering
        J = rand(N, N)

        update_func! = (M, u, p, t; dtgamma = 1.0) -> begin
            M .= dtgamma .* J
            nothing
        end

        op = MatrixOperator(
            copy(J);
            update_func! = update_func!,
            accepted_kwargs = Val((:dtgamma,))  # Use Val for compile-time filtering
        )

        u = rand(N)
        p = nothing
        t = 0.0

        # Warm up
        update_coefficients!(op, u, p, t; dtgamma = 0.5)

        # Test that update_coefficients! with accepted_kwargs has minimal allocations
        # The Val approach significantly reduces allocations compared to plain tuples
        allocs_update = @allocations update_coefficients!(op, u, p, t; dtgamma = 0.5)
        @test allocs_update <= 6  # Some allocations may occur due to Julia version/kwarg handling

        # Test with different dtgamma values - should have similar behavior
        allocs_update2 = @allocations update_coefficients!(op, u, p, t; dtgamma = 1.0)
        @test allocs_update2 <= 6

        allocs_update3 = @allocations update_coefficients!(op, u, p, t; dtgamma = 2.0)
        @test allocs_update3 <= 6

        # Test operator application after update
        v = rand(N)
        w = similar(v)
        op(w, v, u, p, t; dtgamma = 0.5)  # Warm up
        allocs_call = @allocations op(w, v, u, p, t; dtgamma = 0.5)
        @test allocs_call <= 6
    end
end
