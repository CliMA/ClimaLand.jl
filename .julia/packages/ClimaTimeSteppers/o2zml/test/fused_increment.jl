using Test
import Base.Broadcast: broadcasted, materialize
using StaticArrays
using ClimaTimeSteppers: SparseCoeffs, fused_increment, fused_increment!, assign_fused_increment!, zero_coeff
using Random

mat(args...) = materialize(args...)
function dummy_coeffs(S)
    Random.seed!(1234)
    coeffs = rand(S...)
    for I in eachindex(coeffs)
        rand() < 0.5 && continue
        coeffs[I] = 0
    end
    return coeffs
end

function dummy_coeffs_example(S)
    Random.seed!(1234)
    coeffs = rand(S...)
    coeffs[1] = 0
    coeffs[5] = 0
    coeffs[6] = 0
    coeffs[9] = 0
    coeffs[10] = 0
    coeffs[11] = 0
    coeffs[13] = 0
    coeffs[14] = 0
    coeffs[15] = 0
    coeffs[16] = 0
    return SArray{Tuple{S...}}(coeffs)
end

@testset "Test indices" begin
    S = (3, 3)
    coeffs = dummy_coeffs(S)
    mask = BitArray(iszero.(coeffs))
    TMC = typeof(SparseCoeffs(coeffs))
    for i in 1:S[1], j in 1:S[2]
        @test zero_coeff(TMC, i, j) == mask[i, j]
    end
    S = (3,)
    coeffs = dummy_coeffs(S)
    mask = BitArray(iszero.(coeffs))
    TMC = typeof(SparseCoeffs(coeffs))
    for i in 1:S[1]
        @test zero_coeff(TMC, i) == mask[i]
    end
end

import Random
@testset "increment 2D" begin
    FT = Float64
    U = FT[1, 2, 3]
    u = FT[1, 2, 3]
    tend = ntuple(i -> u .* i, 3)
    coeffs = dummy_coeffs((3, 3))
    coeffs .= 0
    sc = SparseCoeffs(coeffs)
    dt = 0.5
    # edge case: zero coeffs
    @test fused_increment(u, dt, sc, tend, Val(3)) == u
    @test fused_increment!(u, dt, sc, tend, Val(3)) == nothing

    FT = Float64
    u = FT[1, 2, 3]
    tend = ntuple(i -> u .* (i + 3), 3)
    coeffs = dummy_coeffs((3, 3))
    coeffs .= 1
    sc = SparseCoeffs(coeffs)
    dt = 0.5

    @test fused_increment(u, dt, sc, tend, Val(1)) == u

    bc2 = broadcasted(+, u, broadcasted(*, dt * coeffs[2, 1], tend[1]))
    @test fused_increment(u, dt, sc, tend, Val(2)) == bc2
    @test mat(fused_increment(u, dt, sc, tend, Val(2))) == @. u + dt * coeffs[2, 1] * tend[1]

    bc3 = broadcasted(+, u, broadcasted(*, dt * coeffs[3, 1], tend[1]), broadcasted(*, dt * coeffs[3, 2], tend[2]))
    @test mat(fused_increment(u, dt, sc, tend, Val(3))) == mat(bc3)

    @test materialize(bc2) == @. u + dt * coeffs[2, 1] * tend[1]

    assign_fused_increment!(U, u, dt, sc, tend, Val(2))
    @test U == @. u + dt * coeffs[2, 1] * tend[1]

    FT = Float64
    u = FT[1, 2, 3]
    tend = ntuple(i -> u .* (i + 1), 3)
    coeffs = dummy_coeffs_example((4, 4))
    sc = SparseCoeffs(coeffs)
    dt = 0.5

    bcb = fused_increment(u, dt, sc, tend, Val(4))
    fused_increment!(u, dt, sc, tend, Val(4))
end

@testset "increment 1D" begin
    FT = Float64
    U = FT[1, 2, 3]
    u = FT[1, 2, 3]
    tend = ntuple(i -> u .* i, 3)
    coeffs = dummy_coeffs((3,))
    coeffs .= 0
    sc = SparseCoeffs(coeffs)
    dt = 0.5
    # Edge case (zero coeffs)
    @test fused_increment(u, dt, sc, tend, Val(1)) == u
    @test fused_increment!(u, dt, sc, tend, Val(1)) == nothing

    FT = Float64
    u = FT[1, 2, 3]
    tend = ntuple(i -> u .* i, 3)
    coeffs = dummy_coeffs((3,))
    coeffs .= 1
    sc = SparseCoeffs(coeffs)
    dt = 0.5

    bc2 = broadcasted(+, u, broadcasted(*, dt * coeffs[1], tend[1]))
    @test fused_increment(u, dt, sc, tend, Val(1)) == bc2

    bc3 = broadcasted(+, u, broadcasted(*, dt * coeffs[1], tend[1]), broadcasted(*, dt * coeffs[2], tend[2]))
    @test fused_increment(u, dt, sc, tend, Val(2)) == bc3

    @test Base.Broadcast.materialize(bc2) == @. u + dt * coeffs[1] * tend[1]

    assign_fused_increment!(U, u, dt, sc, tend, Val(1))
    @test U == @. u + dt * coeffs[1] * tend[1]
end

@testset "increment 1D mask" begin
    FT = Float64
    u = FT[1, 2, 3]
    tend = ntuple(i -> u .* i, 3)
    coeffs = dummy_coeffs((3,))
    coeffs .= 1
    coeffs[2] = 0
    sc = SparseCoeffs(coeffs)
    dt = 0.5

    bc3 = broadcasted(+, u, broadcasted(*, dt * coeffs[1], tend[1]))
    @test fused_increment(u, dt, sc, tend, Val(2)) == bc3
end

@testset "increment 2D mask" begin
    FT = Float64
    u = FT[1, 2, 3]
    tend = ntuple(i -> u .* i, 3)
    coeffs = dummy_coeffs((3, 3))
    coeffs .= 1
    coeffs[3, 2] = 0
    sc = SparseCoeffs(coeffs)
    dt = 0.5

    bc3 = broadcasted(+, u, broadcasted(*, dt * coeffs[3, 1], tend[1]))
    @test fused_increment(u, dt, sc, tend, Val(3)) == bc3
end
