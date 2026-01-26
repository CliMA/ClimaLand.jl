using WoodburyMatrices
using Test

@testset "Woodbury" begin
    d  = 2 .+ rand(5)
    dl = -rand(4)
    du = -rand(4)
    M = Tridiagonal(dl, d, du)
    F = lu(M)
    U = sprand(5,2,0.2)
    V = sprand(2,5,0.2)
    C = rand(2,2)
    W = Woodbury(F, U, C, V)

    src = rand(5, 8)
    @test W\src ≈ AxisAlgorithms.A_ldiv_B_md(W, src, 1)
    src = rand(5, 5, 5, 5)
    for dim = 1:4
        dest1 = mapslices(x->W\x, copy(src), dims=dim)
        dest2 = similar(src)
        AxisAlgorithms.A_ldiv_B_md!(dest2, W, src, dim)
        @test dest1 ≈ dest2
    end
end
