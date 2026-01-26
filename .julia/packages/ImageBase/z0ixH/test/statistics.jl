using ImageBase
using ImageBase: varmult
using Statistics
using Test

@testset "Reductions" begin
    _abs(x::Colorant) = mapreducec(abs, +, 0, x)

    @testset "sumfinite, meanfinite, varfinite" begin
        for T in generate_test_types([N0f8, Float32], [Gray, RGB])
            A = rand(T, 5, 5)
            s12 = sum(A, dims=(1,2))
            @test eltype(s12) <: Union{T, float(T), float64(T)}

            @test sumfinite(A) ≈ sum(A)
            @test sumfinite(A, dims=1) ≈ sum(A, dims=1)
            @test sumfinite(A, dims=(1, 2)) ≈ sum(A, dims=(1, 2))

            @test meanfinite(A) ≈ mean(A)
            @test meanfinite(A, dims=1) ≈ mean(A, dims=1)
            @test meanfinite(A, dims=(1, 2)) ≈ mean(A, dims=(1, 2))

            @test varfinite(A) ≈ varmult(⋅, A)
            @test varfinite(A, dims=1) ≈ varmult(⋅, A, dims=1)
            @test varfinite(A, dims=(1, 2)) ≈ varmult(⋅, A, dims=(1, 2))

            # test NaN/Inf
            if eltype(T) != N0f8
                A = rand(T, 5, 5) .- 0.5 .* oneunit(T)
                A[1] = Inf
                @test sum(A) ≈ A[1]
                @test sum(abs, A) ≈ A[1]
                @test sumfinite(A) ≈ sum(A[2:end])
                @test sumfinite(abs, A) ≈ sum(abs, A[2:end])
                A[1] = NaN
                @test isnan(sum(A))
                @test isnan(sum(abs, A))
                @test sumfinite(A) ≈ sum(A[2:end])
                @test sumfinite(abs, A) ≈ sum(abs, A[2:end])

                A = rand(T, 5, 5) .- 0.5 .* oneunit(T)
                A[1] = Inf
                @test mean(A) ≈ A[1]
                @test mean(abs, A) ≈ A[1]
                @test meanfinite(A) ≈ mean(A[2:end])
                @test meanfinite(abs, A) ≈ mean(abs, A[2:end])
                A[1] = NaN
                @test isnan(mean(A))
                @test isnan(mean(abs, A))
                @test meanfinite(A) ≈ mean(A[2:end])
                @test meanfinite(abs, A) ≈ mean(abs, A[2:end])

                A = rand(T, 5, 5)
                A[1] = Inf
                @test isnan(varmult(⋅, A))
                @test varfinite(A) ≈ varmult(⋅, A[2:end])
                A[1] = NaN
                @test isnan(varmult(⋅, A))
                @test varfinite(A) ≈ varmult(⋅, A[2:end])
            end
        end

        A = [NaN, 1, 2, 3]
        @test meanfinite(A, dims=1) ≈ [2]
        @test varfinite(A, dims=1) ≈ [1]

        A = [NaN NaN 1;
            1 2 3]
        vf = varfinite(A, dims=2)
        @test isnan(vf[1])

        A = [NaN 1 2 3;
            NaN 6 5 4]
        mf = meanfinite(A, dims=1)
        vf = varfinite(A, dims=1)
        @test isnan(mf[1])
        @test mf[2:end] ≈ [3.5,3.5,3.5]
        @test isnan(vf[1])
        @test vf[2:end] ≈ [12.5,4.5,0.5]

        @test meanfinite(A, dims=2) ≈ reshape([2, 5], 2, 1)
        @test varfinite(A, dims=2) ≈ reshape([1, 1], 2, 1)

        @test meanfinite(A, dims=(1,2)) ≈ [3.5]
        @test varfinite(A, dims=(1,2)) ≈ [3.5]

        # Ensure we're consistant with our decision to `abs2` in ColorVectorSpace
        # See also: https://github.com/JuliaGraphics/ColorVectorSpace.jl/blob/master/README.md#abs-and-abs2
        A = rand(Gray, 4, 4)
        @test varfinite(A) ≈ varfinite(RGB.(A))
        A[1] = Inf
        @test varfinite(A) ≈ varfinite(RGB.(A))
        A[1] = NaN
        @test varfinite(A) ≈ varfinite(RGB.(A))
    end

    @testset "minfinite, maxfinite, maxabsfinite" begin
        for T in generate_test_types([N0f8, Float32], [Gray, ])
            A = rand(T, 5, 5)
            @test @inferred(minimum_finite(A)) == minimum(A)
            @test @inferred(maximum_finite(A)) == maximum(A)
            @test minimum_finite(A; dims=1) == minimum(A; dims=1)
            @test maximum_finite(A; dims=1) == maximum(A; dims=1)
            @test_broken @inferred maximum_finite(A; dims=1)
            @test_broken @inferred minimum_finite(A; dims=1)

            @test maximum_finite(abs2, A) == maximum(abs2, A)
            @test_broken @inferred maximum(abs2, A)
            @test_broken @inferred minimum(abs2, A)

            if eltype(T) != N0f8
                A = rand(T, 5, 5) .- 0.5 * rand(T, 5, 5)
                A[1] = Inf

                @test @inferred(minimum_finite(A)) == minimum(A[2:end])
                @test @inferred(maximum_finite(A)) == maximum(A[2:end])
                @test minimum_finite(abs, A) == minimum(abs, A[2:end])
                @test maximum_finite(abs2, A) == maximum(abs2, A[2:end])
            end
        end

        # minimum_finite and maximum_finite for RGB are processed per channel
        A = rand(RGB{Float32}, 5, 5)
        @test minimum_finite(A) == RGB(minimum(channelview(A), dims=(2, 3))...)
        @test maximum_finite(A) == RGB(maximum(channelview(A), dims=(2, 3))...)

        # Container of abstract type
        A = Any[1, 2, Gray(0.3)]
        @test minimum_finite(A) == 0.3
        @test maximum_finite(A) == 2.0
        @test Base.return_types(minimum_finite, (typeof(A),)) == Base.return_types(minimum, (typeof(A), ))
    end

end
