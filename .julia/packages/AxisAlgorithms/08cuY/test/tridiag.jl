using Test, OffsetArrays

@testset "Tridiag" begin
    d  = 2 .+ rand(5)
    dl = rand(4)
    du = rand(4)
    M = Tridiagonal(dl, d, du)
    F = lu(M)
    src = rand(5,5,5)
    for dim = 1:3
        dest1 = mapslices(x->ldiv!(F, x), copy(src), dims=dim)
        dest2 = similar(src)
        AxisAlgorithms.A_ldiv_B_md!(dest2, F, src, dim)
        @test dest1 ≈ dest2
    end

    src = OffsetArray(src, -2:2, 1:5, 0:4)
    dest1 = mapslices(x->ldiv!(F, x), copy(src), dims=2)
    dest2 = similar(src)
    AxisAlgorithms.A_ldiv_B_md!(dest2, F, src, 2)
    @test dest1 ≈ dest2
end
