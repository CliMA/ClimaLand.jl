module Simple
using VectorInterface
using Test
using TestExtras

deepcollect(x) = vcat(map(deepcollect, x)...)
deepcollect(x::Number) = x

x = randn(3, 3, 3)
y = randn(3, 3, 3)

@testset "scalartype" begin
    s = @constinferred scalartype(x)
    @test s == Float64
    @test_throws ArgumentError scalartype(())
    @test_throws ArgumentError scalartype((randn(Float64, 2), randn(ComplexF64, 3)))
end

@testset "zerovector" begin
    z = @constinferred zerovector(x)
    @test all(iszero, deepcollect(z))
    @test all(deepcollect(z) .=== zero(scalartype(x)))
    z1 = @constinferred zerovector!!(deepcopy(x))
    @test all(deepcollect(z1) .=== zero(scalartype(x)))
    z2 = @constinferred zerovector!(deepcopy(x))
    @test all(deepcollect(z2) .=== zero(scalartype(x)))

    z3 = @constinferred zerovector(x, ComplexF64)
    @test all(deepcollect(z3) .=== zero(ComplexF64))
    z4 = @constinferred zerovector!!(deepcopy(x), ComplexF64)
    @test all(deepcollect(z4) .=== zero(ComplexF64))
    @test_throws MethodError zerovector!(deepcopy(x), ComplexF64)
end

@testset "scale" begin
    α = randn()
    z = @constinferred scale(x, α)
    @test all(deepcollect(z) .== α .* deepcollect(x))

    z2 = @constinferred scale!!(deepcopy(x), α)
    @test deepcollect(z2) ≈ (α .* deepcollect(x))
    xcopy = deepcopy(x)
    z2 = @constinferred scale!!(deepcopy(y), xcopy, α)
    @test deepcollect(z2) ≈ (α .* deepcollect(x))
    @test all(deepcollect(xcopy) .== deepcollect(x))

    z3 = @constinferred scale!(deepcopy(x), α)
    @test deepcollect(z3) ≈ (α .* deepcollect(x))
    xcopy = deepcopy(x)
    z3 = @constinferred scale!(zerovector(x), xcopy, α)
    @test deepcollect(z3) ≈ (α .* deepcollect(x))
    @test all(deepcollect(xcopy) .== deepcollect(x))

    α = randn(ComplexF64)
    z4 = @constinferred scale(x, α)
    @test deepcollect(z4) ≈ (α .* deepcollect(x))
    xcopy = deepcopy(x)
    z5 = @constinferred scale!!(xcopy, α)
    @test deepcollect(z5) ≈ (α .* deepcollect(x))
    @test all(deepcollect(xcopy) .== deepcollect(x))
    @test_throws InexactError scale!(xcopy, α)

    α = randn(ComplexF64)
    xcopy = deepcopy(x)
    z6 = @constinferred scale!!(zerovector(x), xcopy, α)
    @test deepcollect(z6) ≈ (α .* deepcollect(x))
    @test all(deepcollect(xcopy) .== deepcollect(x))
    @test_throws InexactError scale!(zerovector(x), xcopy, α)

    xz = @constinferred zerovector(x, ComplexF64)
    z6 = @constinferred scale!!(xz, xcopy, α)
    @test deepcollect(z6) ≈ (α .* deepcollect(x))
    @test all(deepcollect(xcopy) .== deepcollect(x))

    z7 = @constinferred scale!(xz, xcopy, α)
    @test deepcollect(z7) ≈ (α .* deepcollect(x))
    @test all(deepcollect(xcopy) .== deepcollect(x))

    ycomplex = zerovector(y, ComplexF64)
    α = randn(Float64)
    xcopy = deepcopy(x)
    z8 = @constinferred scale!!(ycomplex, xcopy, α)
    @test z8 === ycomplex
    @test all(deepcollect(z8) .== α .* deepcollect(xcopy))
end

@testset "add" begin
    α, β = randn(2)
    z = @constinferred add(y, x)
    @test all(deepcollect(z) .== deepcollect(x) .+ deepcollect(y))
    z = @constinferred add(y, x, α)
    # for some reason, on some Julia versions on some platforms, but only in test mode
    # there is a small floating point discrepancy, which makes the following test fail:
    # @test all(deepcollect(z) .== muladd.(deepcollect(x), α, deepcollect(y)))
    @test deepcollect(z) ≈ muladd.(deepcollect(x), α, deepcollect(y))
    z = @constinferred add(y, x, α, β)
    # for some reason, on some Julia versions on some platforms, but only in test mode
    # there is a small floating point discrepancy, which makes the following test fail:
    # @test all(deepcollect(z) .== muladd.(deepcollect(x), α, deepcollect(y) .* β))
    @test deepcollect(z) ≈ muladd.(deepcollect(x), α, deepcollect(y) .* β)

    α, β = randn(2)
    xcopy = deepcopy(x)
    z2 = @constinferred add!!(deepcopy(y), xcopy)
    @test deepcollect(z2) ≈ (deepcollect(x) .+ deepcollect(y))
    @test all(deepcollect(xcopy) .== deepcollect(x))
    z2 = @constinferred add!!(deepcopy(y), xcopy, α)
    @test deepcollect(z2) ≈ (muladd.(deepcollect(x), α, deepcollect(y)))
    @test all(deepcollect(xcopy) .== deepcollect(x))
    z2 = @constinferred add!!(deepcopy(y), xcopy, α, β)
    @test deepcollect(z2) ≈ (muladd.(deepcollect(x), α, deepcollect(y) .* β))
    @test all(deepcollect(xcopy) .== deepcollect(x))

    α, β = randn(2)
    z3 = @constinferred add!(deepcopy(y), xcopy)
    @test deepcollect(z3) ≈ (deepcollect(y) .+ deepcollect(x))
    @test all(deepcollect(xcopy) .== deepcollect(x))
    z3 = @constinferred add!(deepcopy(y), xcopy, α)
    @test all(deepcollect(xcopy) .== deepcollect(x))
    @test deepcollect(z3) ≈ (muladd.(deepcollect(x), α, deepcollect(y)))
    z3 = @constinferred add!(deepcopy(y), xcopy, α, β)
    @test deepcollect(z3) ≈ (muladd.(deepcollect(x), α, deepcollect(y) .* β))
    @test all(deepcollect(xcopy) .== deepcollect(x))

    α, β = randn(ComplexF64, 2)
    z4 = @constinferred add(y, x, α)
    @test deepcollect(z4) ≈ (muladd.(deepcollect(x), α, deepcollect(y)))
    z4 = @constinferred add(y, x, α, β)
    @test deepcollect(z4) ≈ (muladd.(deepcollect(x), α, deepcollect(y) .* β))

    α, β = randn(ComplexF64, 2)
    z5 = @constinferred add!!(deepcopy(y), xcopy, α)
    @test deepcollect(z5) ≈ (muladd.(deepcollect(x), α, deepcollect(y)))
    @test all(deepcollect(xcopy) .== deepcollect(x))
    z5 = @constinferred add!!(deepcopy(y), xcopy, α, β)
    @test deepcollect(z5) ≈ (muladd.(deepcollect(x), α, deepcollect(y) .* β))
    @test all(deepcollect(xcopy) .== deepcollect(x))

    α, β = randn(ComplexF64, 2)
    z5 = @constinferred add!!(deepcopy(y), xcopy, α)
    @test deepcollect(z5) ≈ (muladd.(deepcollect(x), α, deepcollect(y)))
    @test all(deepcollect(xcopy) .== deepcollect(x))
    z5 = @constinferred add!!(deepcopy(y), xcopy, α, β)
    @test deepcollect(z5) ≈ (muladd.(deepcollect(x), α, deepcollect(y) .* β))
    @test all(deepcollect(xcopy) .== deepcollect(x))

    α, β = randn(ComplexF64, 2)
    @test_throws InexactError add!(deepcopy(y), xcopy, α)
    @test_throws InexactError add!(deepcopy(y), xcopy, α, β)
end

@testset "inner" begin
    s = @constinferred inner(x, y)
    @test s ≈ inner(deepcollect(x), deepcollect(y))

    α, β = randn(ComplexF64, 2)
    s2 = @constinferred inner(scale(x, α), scale(y, β))
    @test s2 ≈ inner(α * deepcollect(x), β * deepcollect(y))
end

end
