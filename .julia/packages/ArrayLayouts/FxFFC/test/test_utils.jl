module TestUtils

using ArrayLayouts, LinearAlgebra, FillArrays, Test

@testset "_copy_oftype" begin
    # Copied from original implementation in FillArrays
    @testset "FillArrays" begin
        m = Eye(10)
        D = Diagonal(Fill(2,10))

        @test ArrayLayouts._copy_oftype(m, eltype(m)) ≡ m
        @test ArrayLayouts._copy_oftype(m, Int) ≡ Eye{Int}(10)
        @test ArrayLayouts._copy_oftype(D, eltype(D)) ≡ D
        @test ArrayLayouts._copy_oftype(D, Float64) ≡ Diagonal(Fill(2.0,10))

        # test that _copy_oftype does, in fact, copy the array
        D2 = Diagonal([1,1])
        @test ArrayLayouts._copy_oftype(D2, Float64) isa Diagonal{Float64}
        @test ArrayLayouts._copy_oftype(D2, eltype(D2)) == D2
        @test ArrayLayouts._copy_oftype(D2, eltype(D2)) !== D2
    end

    @testset "general" begin
        for T in (Float32, Float64)
            u = T(1):T(10)
            v = collect(u)
            for S in (Float32, Float64)
                # Usually an actual copy
                @test ArrayLayouts._copy_oftype(u, S) isa AbstractRange{S}
                @test ArrayLayouts._copy_oftype(u, S) == u
                @test (ArrayLayouts._copy_oftype(u, S) === u) === (T === S && copy(u) === u)

                # Always an actual copy
                @test ArrayLayouts._copy_oftype(v, S) isa Array{S}
                @test ArrayLayouts._copy_oftype(v, S) == v
                @test ArrayLayouts._copy_oftype(v, S) !== v
            end
        end

        A = [1 3; 2 4]
        ArrayLayouts._fill_lmul!(2.0, A)
        @test A == 2 * [1 3; 2 4]
        ArrayLayouts._fill_lmul!(0, A)
        @test all(==(0), A)
    end
end

end
