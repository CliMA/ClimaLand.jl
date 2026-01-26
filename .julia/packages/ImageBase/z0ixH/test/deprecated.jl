@testset "deprecation" begin
    @testset "restrict" begin
        A = rand(N0f8, 4, 5, 3)
        @test restrict(A, [1, 2]) == restrict(A, (1, 2))
    end

    @testset "statistics" begin
        A = rand(Float32, 4, 4) .- 0.5
        @test minfinite(A, dims=1) == minimum_finite(A, dims=1)
        @test minfinite(A) == minimum_finite(A)
        @test maxfinite(A, dims=1) == maximum_finite(A, dims=1)
        @test maxfinite(A) == maximum_finite(A)

        @test maxabsfinite(A) == maximum_finite(abs, A)
        @test maxabsfinite(A, dims=1) == maximum_finite(abs, A, dims=1)
    end

    @testset "fdiff entrypoints" begin
        A = rand(Float32, 5)
        @test ImageBase.fdiff(A, rev=true) == ImageBase.FiniteDiff.fdiff(A, rev=true)
        out = similar(A)
        @test ImageBase.fdiff!(out, A, rev=false) == ImageBase.FiniteDiff.fdiff!(out, A, rev=false)
    end
end
