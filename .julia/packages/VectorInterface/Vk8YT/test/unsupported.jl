module Unsupported
using VectorInterface
using Test
using TestExtras

include("simplevec.jl")

deepcollect(x) = vcat(map(deepcollect, x)...)
deepcollect(x::Number) = x

x = SimpleVec(randn(4, 4))
y = SimpleVec(randn(4, 4))

@testset "scalartype" begin
    s = @constinferred scalartype(x)
    @test s == Float64
    @test_throws ArgumentError scalartype(())
    @test_throws ArgumentError scalartype((randn(Float64, 2), randn(ComplexF64, 3)))
end

@testset "zerovector" begin
    z = @constinferred zerovector(x)
    @test all(iszero, deepcollect(z))
    @test all(deepcollect(z) .== zero(scalartype(x)))
    z1 = @constinferred zerovector!!(deepcopy(x))
    @test all(deepcollect(z1) .== zero(scalartype(x)))

    z2 = @constinferred zerovector(x, ComplexF64)
    @test all(deepcollect(z2) .== zero(ComplexF64))
    z3 = @constinferred zerovector!!(deepcopy(x), ComplexF64)
    @test all(deepcollect(z3) .== zero(ComplexF64))

    xy = [deepcopy(x), deepcopy(y)]
    z4 = @constinferred zerovector!(xy)
    @test all(deepcollect(z4) .== zero(scalartype(xy)))
    @test_throws MethodError zerovector!(xy, ComplexF64)
end

@testset "scale" begin
    α = randn()
    z = @constinferred scale(x, α)
    @test all(deepcollect(z) .== α .* deepcollect(x))
    z = @constinferred scale(x, One())
    @test all(deepcollect(z) .== deepcollect(x))
    z2 = @constinferred scale!!(deepcopy(x), α)
    @test deepcollect(z2) ≈ (α .* deepcollect(x))
    z2 = @constinferred scale!!(deepcopy(x), One())
    @test all(deepcollect(z2) .== deepcollect(x))
    xcopy = deepcopy(x)
    z2 = @constinferred scale!!(deepcopy(y), xcopy, α)
    @test deepcollect(z2) ≈ (α .* deepcollect(x))
    @test all(deepcollect(xcopy) .== deepcollect(x))
    z3 = @constinferred scale!(deepcopy(x), α)
    @test deepcollect(z3) ≈ (α .* deepcollect(x))
    z3 = @constinferred scale!(deepcopy(x), One())
    @test all(deepcollect(z3) .== deepcollect(x))

    α = randn(ComplexF64)
    z4 = @constinferred scale(x, α)
    @test deepcollect(z4) ≈ (α .* deepcollect(x))
    xcopy = deepcopy(x)
    z5 = @constinferred scale!!(xcopy, α)
    @test deepcollect(z5) ≈ (α .* deepcollect(x))
    @test all(deepcollect(xcopy) .== deepcollect(x))

    α = randn(ComplexF64)
    xcopy = deepcopy(x)
    z6 = @constinferred scale!!(zerovector(x), xcopy, α)
    @test deepcollect(z6) ≈ (α .* deepcollect(x))
    @test all(deepcollect(xcopy) .== deepcollect(x))
    xz = @constinferred zerovector(x, ComplexF64)
    z6 = @constinferred scale!!(xz, xcopy, α)
    @test deepcollect(z6) ≈ (α .* deepcollect(x))
    @test all(deepcollect(xcopy) .== deepcollect(x))
end

@testset "add" begin
    α, β = randn(2)
    z = @constinferred add(y, x)
    @test all(deepcollect(z) ≈ deepcollect(x) .+ deepcollect(y))
    z = @constinferred add(y, x, α)
    @test all(deepcollect(z) ≈ muladd.(deepcollect(x), α, deepcollect(y)))
    z = @constinferred add(y, x, α, β)
    # for some reason, on Julia 1.6, but only in test mode
    # there is a small floating point discrepancy, which makes the following test fail:
    # @test all(deepcollect(z) .== muladd.(deepcollect(x), α, deepcollect(y) .* β))
    @test deepcollect(z) ≈ muladd.(deepcollect(x), α, deepcollect(y) .* β)

    z = @constinferred add(y, x, One(), β)
    @test deepcollect(z) ≈ muladd.(deepcollect(y), β, deepcollect(x))
    z = @constinferred add(y, x, α, One())
    @test deepcollect(z) ≈ muladd.(deepcollect(x), α, deepcollect(y))
    z = @constinferred add(y, x, One(), One())
    @test all(deepcollect(z) .== deepcollect(x) .+ deepcollect(y))
    z = @constinferred add(y, x, One())
    @test all(deepcollect(z) .== deepcollect(x) .+ deepcollect(y))

    α, β = randn(2)
    z2 = @constinferred add!!(deepcopy(y), deepcopy(x))
    @test deepcollect(z2) ≈ (deepcollect(x) .+ deepcollect(y))
    z2 = @constinferred add!!(deepcopy(y), deepcopy(x), α)
    @test deepcollect(z2) ≈ (muladd.(deepcollect(x), α, deepcollect(y)))
    z2 = @constinferred add!!(deepcopy(y), deepcopy(x), α, β)
    @test deepcollect(z2) ≈ (muladd.(deepcollect(x), α, deepcollect(y) .* β))

    z2 = @constinferred add!!(deepcopy(y), deepcopy(x), One(), β)
    @test deepcollect(z2) ≈ muladd.(deepcollect(y), β, deepcollect(x))
    z2 = @constinferred add!!(deepcopy(y), deepcopy(x), α, One())
    @test deepcollect(z2) ≈ muladd.(deepcollect(x), α, deepcollect(y))
    z2 = @constinferred add!!(deepcopy(y), deepcopy(x), One(), One())
    @test all(deepcollect(z2) .== deepcollect(x) .+ deepcollect(y))
    z2 = @constinferred add!!(deepcopy(y), deepcopy(x), One())
    @test all(deepcollect(z2) .== deepcollect(x) .+ deepcollect(y))

    α, β = randn(2)
    z3 = @constinferred add!(deepcopy(y), deepcopy(x))
    @test deepcollect(z3) ≈ (deepcollect(x) .+ deepcollect(y))
    z3 = @constinferred add!(deepcopy(y), deepcopy(x), α)
    @test deepcollect(z3) ≈ (muladd.(deepcollect(x), α, deepcollect(y)))
    z3 = @constinferred add!(deepcopy(y), deepcopy(x), α, β)
    @test deepcollect(z3) ≈ (muladd.(deepcollect(x), α, deepcollect(y) .* β))

    z3 = @constinferred add!!(deepcopy(y), deepcopy(x), One(), β)
    @test deepcollect(z3) ≈ muladd.(deepcollect(y), β, deepcollect(x))
    z3 = @constinferred add!!(deepcopy(y), deepcopy(x), α, One())
    @test deepcollect(z3) ≈ muladd.(deepcollect(x), α, deepcollect(y))
    z3 = @constinferred add!!(deepcopy(y), deepcopy(x), One(), One())
    @test all(deepcollect(z3) .== deepcollect(x) .+ deepcollect(y))
    z3 = @constinferred add!!(deepcopy(y), deepcopy(x), One())
    @test all(deepcollect(z3) .== deepcollect(x) .+ deepcollect(y))

    α, β = randn(ComplexF64, 2)
    z4 = @constinferred add(y, x, α)
    @test deepcollect(z4) ≈ (muladd.(deepcollect(x), α, deepcollect(y)))
    z4 = @constinferred add(y, x, α, β)
    @test deepcollect(z4) ≈ (muladd.(deepcollect(x), α, deepcollect(y) .* β))

    α, β = randn(ComplexF64, 2)
    z5 = @constinferred add!!(deepcopy(y), deepcopy(x), α)
    @test deepcollect(z5) ≈ (muladd.(deepcollect(x), α, deepcollect(y)))
    z5 = @constinferred add!!(deepcopy(y), deepcopy(x), α, β)
    @test deepcollect(z5) ≈ (muladd.(deepcollect(x), α, deepcollect(y) .* β))
end

@testset "inner" begin
    s = @constinferred inner(x, y)
    @test s ≈ inner(deepcollect(x), deepcollect(y))

    α, β = randn(ComplexF64, 2)
    s2 = @constinferred inner(scale(x, α), scale(y, β))
    @test s2 ≈ inner(α * deepcollect(x), β * deepcollect(y))
end

end
