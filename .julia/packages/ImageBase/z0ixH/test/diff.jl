using ImageBase.FiniteDiff
# TODO(johnnychen94): remove this after we delete `Imagebase.fdiff` and `ImageBase.fdiff!` entrypoints
using ImageBase.FiniteDiff: fdiff, fdiff!

@testset "fdiff" begin
    # Base.diff doesn't promote integer to float
    @test ImageBase.FiniteDiff.maybe_floattype(Int) == Int
    @test ImageBase.FiniteDiff.maybe_floattype(N0f8) == Float32
    @test ImageBase.FiniteDiff.maybe_floattype(RGB{N0f8}) == RGB{Float32}

    @testset "API" begin
        # fdiff! works the same as fdiff
        mat_in = rand(3, 3, 3)
        mat_out = similar(mat_in)
        fdiff!(mat_out, mat_in, dims = 2)
        @test mat_out == fdiff(mat_in, dims = 2)

        mat_in = rand(3, 3, 3)
        mat_out = similar(mat_in)
        fdiff!(mat_out, mat_in, dims = 3, rev=true)
        @test mat_out == fdiff(mat_in, dims = 3, rev=true)

        # default `dims` are only available when input is a vector
        v = rand(9)
        @test fdiff(v) == fdiff(v; dims=1)
        m = rand(3, 3)
        @test_throws UndefKeywordError fdiff(m)
    end

    @testset "NumericalTests" begin
        a = reshape(collect(1:9), 3, 3)
        b_fd_1 = [1 1 1; 1 1 1; -2 -2 -2]
        b_fd_2 = [3 3 -6; 3 3 -6; 3 3 -6]
        b_bd_1 = [-2 -2 -2; 1 1 1; 1 1 1]
        b_bd_2 = [-6 3 3; -6 3 3; -6 3 3]
        out = similar(a)

        @test fdiff(a, dims = 1) == b_fd_1
        @test fdiff(a, dims = 2) == b_fd_2
        @test fdiff(a, dims = 1, rev=true) == b_bd_1
        @test fdiff(a, dims = 2, rev=true) == b_bd_2
        fdiff!(out, a, dims = 1)
        @test out == b_fd_1
        fdiff!(out, a, dims = 2)
        @test out == b_fd_2
        fdiff!(out, a, dims = 1, rev=true)
        @test out == b_bd_1
        fdiff!(out, a, dims = 2, rev=true)
        @test out == b_bd_2

        @testset "Boundary" begin
            # These originally lives in Images and we use it to check the compatibility
            forwarddiffy(u::AbstractMatrix) = [u[2:end,:]; u[end:end,:]] - u
            forwarddiffx(u::AbstractMatrix) = [u[:,2:end] u[:,end:end]] - u
            backdiffy(u::AbstractMatrix) = u - [u[1:1,:]; u[1:end-1,:]]
            backdiffx(u::AbstractMatrix) = u - [u[:,1:1] u[:,1:end-1]]

            X = rand(3, 3)
            @test forwarddiffx(X) == fdiff(X, dims=2, boundary=:zero)
            @test forwarddiffy(X) == fdiff(X, dims=1, boundary=:zero)
            @test backdiffx(X) == fdiff(X, dims=2, rev=true, boundary=:zero)
            @test backdiffy(X) == fdiff(X, dims=1, rev=true, boundary=:zero)
        end

        if VERSION >= v"1.1"
            # dims keyword for diff only exists after Julia 1.1
            # check numerical results with base implementation
            drop_last_slice(X, dims) = collect(StackView(collect(eachslice(X, dims=dims))[1:end-1]..., dims=dims))
            function drop_last_slice(X::AbstractVector, dims)
                @assert dims==1
                X[1:end-1]
            end
            drop_first_slice(X, dims) = collect(StackView(collect(eachslice(X, dims=dims))[2:end]..., dims=dims))
            function drop_first_slice(X::AbstractVector, dims)
                @assert dims==1
                X[2:end]
            end

            for N in 1:3
                sz = ntuple(_->5, N)
                A = rand(sz...)
                A_out = similar(A)

                for dims = 1:N
                    out_base = diff(A; dims=dims)
                    out = fdiff(A; dims=dims)
                    @test out_base == drop_last_slice(out, dims)

                    out_base = .-reverse(diff(reverse(A; dims=dims); dims=dims); dims=dims)
                    out = fdiff(A; dims=dims, rev=true)
                    @test out_base == drop_first_slice(out, dims)
                end
            end
        end
    end

    @testset "OffsetArrays" begin
        A = OffsetArray(rand(3, 3), -1, -1)
        A_out = fdiff(A, dims=1)
        @test axes(A_out) == (0:2, 0:2)
        @test A_out.parent == fdiff(parent(A), dims=1)

        A = OffsetArray(rand(3, 3), -1, -1)
        A_out = fdiff(A, dims=1, rev=true)
        @test axes(A_out) == (0:2, 0:2)
        @test A_out.parent == fdiff(parent(A), dims=1, rev=true)
    end

    @testset "FixedPoint" begin
        A = rand(N0f8, 6, 6)
        @test fdiff(A, dims=1) == fdiff(float.(A), dims=1)
    end
end

@testset "fgradient" begin
    for T in generate_test_types([N0f8, Float32], [Gray, RGB])
        for sz in [(7, ), (7, 5), (7, 5, 3)]
            A = rand(T, sz...)

            ∇A = fgradient(A)
            for i in 1:length(sz)
                @test ∇A[i] == fdiff(A, dims=i)
            end
            ∇A = fgradient(A; boundary=:zero)
            for i in 1:length(sz)
                @test ∇A[i] == fdiff(A, dims=i, boundary=:zero)
            end
            ∇A = map(similar, ∇A)
            @test fgradient!(∇A, A) == fgradient(A)
            @test fgradient!(∇A, A; boundary=:zero) == fgradient(A; boundary=:zero)

            ∇tA = fgradient(A, adjoint=true)
            for i in 1:length(sz)
                @test ∇tA[i] == .-fdiff(A, dims=i, rev=true)
            end
            ∇tA = fgradient(A, adjoint=true, boundary=:zero)
            for i in 1:length(sz)
                @test ∇tA[i] == .-fdiff(A, dims=i, rev=true, boundary=:zero)
            end
            ∇tA = map(similar, ∇tA)
            @test fgradient!(∇tA, A, adjoint=true) == fgradient(A, adjoint=true)
            @test fgradient!(∇tA, A, adjoint=true, boundary=:zero) == fgradient(A, adjoint=true, boundary=:zero)
        end
    end
end

@testset "fdiv/flaplacian" begin
    ref_laplacian(X) = imfilter(X, Kernel.Laplacian(ntuple(x->true, ndims(X))), "circular")

    X = [
        5  3  8  1  2  2  3
        5  5  1  3  3  1  1
        1  8  2  9  1  2  7
        7  3  4  5  8  1  5
        1  4  1  8  8  9  7
    ]
    ΔX_ref = [
        -8   10  -26   17    6    7    3
        -8   -3   14    2   -5    4   12
        23  -21   14  -25   18    2  -19
       -18   11   -5    9  -17   20    2
        19   -8   20  -17   -5  -18  -10
    ]
    ΔX = ref_laplacian(X)
    # Base.diff doesn't promote Int to floats so we should probably do the same for laplacian
    @test eltype(ΔX) == Int
    @test ΔX_ref == ΔX

    for T in generate_test_types([N0f8, Float32], [Gray, RGB])
        for sz in [(7,), (7, 7), (7, 7, 7)]
            A = rand(T, sz...)
            ∇A = fgradient(A)
            out = fdiv(∇A)
            @test out ≈ ref_laplacian(A)

            fill!(out, zero(T))
            fdiv!(out, ∇A...)
            @test out == fdiv(∇A)
        end
    end

    for T in generate_test_types([N0f8, Float32], [Gray, RGB])
        for sz in [(7,), (7, 7), (7, 7, 7)]
            A = rand(T, sz...)
            out = flaplacian(A)
            @test out ≈ ref_laplacian(A)

            ∇A = fgradient(A)
            foreach((out, ∇A...)) do x
                fill!(x, zero(T))
            end
            flaplacian!(out, ∇A, A)
            @test out == flaplacian(A)
            @test ∇A == fgradient(A)
        end
    end
end
