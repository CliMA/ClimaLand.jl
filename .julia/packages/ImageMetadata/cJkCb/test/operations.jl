using ImageMetadata, ImageCore, Test

@testset "operations" begin
    function checkmeta(A, B)
        @test isa(A, ImageMeta)
        @test A == B
        if isa(B, ImageMeta)
            @test properties(A) == properties(B)
        end
        nothing
    end
    for A in (rand(Bool, 3, 5), rand(3, 5),
              rand(Gray{N0f8}, 3, 5), rand(RGB{N0f8}, 3, 5))
        M = ImageMeta(A)
        M2 = similar(M)
        checkmeta(-M, -A)
        checkmeta(M .+ zero(eltype(M)), M)
        checkmeta(zero(eltype(M)) .+ M, M)
        checkmeta(M .- zero(eltype(M)), M)
        checkmeta(zero(eltype(M)) .- M, -M)
        B = falses(size(M))
        if !(eltype(A) <: RGB)
            checkmeta(M + B, M)
            checkmeta(M .+ B, M)
            checkmeta((B .+ M .+ B), M)
            checkmeta(B + M, M)
            checkmeta(B .+ M, M)
            @test_throws ErrorException M + M2
            @test_throws ErrorException M .+ M2
            checkmeta(A + M, A+A)
            checkmeta(A .+ M, A+A)
            checkmeta(M + A, A+A)
            checkmeta(M .+ A, A+A)
            checkmeta(M - B, M)
            checkmeta(M .- B, M)
            checkmeta(B - M, -M)
            checkmeta(B .- M, -M)
            @test_throws ErrorException M - M2
            @test_throws ErrorException M .- M2
            checkmeta(A - M, 0*M)
            checkmeta(A .- M, 0*M)
            checkmeta(M - A, 0*M)
            checkmeta(M .- A, 0*M)
            checkmeta(A .- M .+ A, M)
            checkmeta(M .- A .+ M, M)
        end
        checkmeta(M*2, 2*M)
        checkmeta(2 .* M, 2*M)
        checkmeta(M .* 2, 2*M)
        checkmeta(M/2, M/2)
        checkmeta(M./2, M/2)
        checkmeta(M.*B, 0*M)
        checkmeta(B.*M, 0*M)
        if !(eltype(A) <: RGB)
            @test_throws ErrorException M.*M2
            checkmeta(A.*M, A.*A)
            checkmeta(M.*A, A.*A)
            checkmeta(M.^1, M)
        end
        B1 = trues(size(M))
        checkmeta(M./B1, M)
        @test_throws Union{MethodError,ErrorException} M./M2
        if !(eltype(A) <: RGB)
            checkmeta(M .+ false, M)
            checkmeta(false .+ M, M)
            checkmeta(M .+ 0.0, M)
            checkmeta(0.0 .+ M, M)
            if eltype(A) == Float64
                checkmeta(M./A, ones(size(M)))
            end
            @test (M .< 0.5) == (A .< 0.5)
            @test (M .> 0.5) == (A .> 0.5)
            @test (M .< A) == B
            @test (M .> A) == B
            @test (M .== 0.5) == (A .== 0.5)
            @test (M .== A) == B1
        end
    end
end

M = ImageMeta([1,2,3,4])
@test minimum(M) == 1
@test maximum(M) == 4
Mp = M'
@test ndims(Mp) == 2
Ms = dropdims(Mp, dims=1)
@test Ms == M

img = convert(ImageMeta{Gray{N0f16}}, [0.01164 0.01118; 0.01036 0.01187])
@test all((1 .- img) .== (1 .- img))

nothing
