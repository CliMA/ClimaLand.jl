using NaNMath
using Test

@test isnan(NaNMath.log(-10))
@test isnan(NaNMath.log1p(-100))
@test isnan(NaNMath.pow(-1.5,2.3))
@test isnan(NaNMath.pow(-1.5f0,2.3f0))
@test isnan(NaNMath.pow(-1.5,2.3f0))
@test isnan(NaNMath.pow(-1.5f0,2.3))
@test NaNMath.pow(-1.5f0,2) isa Float32
@test NaNMath.pow(-1.5f0,2//1) isa Float32
@test NaNMath.pow(-1.5f0,2.3f0) isa Float32
@test NaNMath.pow(-1.5f0,2.3) isa Float64
@test NaNMath.pow(-1.5,2) isa Float64
@test NaNMath.pow(-1.5,2//1) isa Float64
@test NaNMath.pow(-1.5,2.3f0) isa Float64
@test NaNMath.pow(-1.5,2.3) isa Float64
@test NaNMath.pow(-1,2) === 1
@test NaNMath.pow(2,2) === 4
@test NaNMath.pow(1.0, 1.0+im) === 1.0 + 0.0im
@test NaNMath.pow(1.0+im, 1) === 1.0 + 1.0im
@test NaNMath.pow(1.0+im, 1.0) === 1.0 + 1.0im
@test isnan(NaNMath.sqrt(-5))
@test NaNMath.sqrt(5) == Base.sqrt(5)
@test isnan(NaNMath.sqrt(-3.2f0)) && NaNMath.sqrt(-3.2f0) isa Float32
@test isnan(NaNMath.sqrt(-BigFloat(7.0))) && NaNMath.sqrt(-BigFloat(7.0)) isa BigFloat 
@test isnan(NaNMath.sqrt(-7)) && NaNMath.sqrt(-7) isa Float64 
@inferred NaNMath.sqrt(5)
@inferred NaNMath.sqrt(5.0)
@inferred NaNMath.sqrt(5.0f0)
@inferred NaNMath.sqrt(-5)
@inferred NaNMath.sqrt(-5.0)
@inferred NaNMath.sqrt(-5.0f0)
@test NaNMath.sum([1., 2., NaN]) == 3.0
@test NaNMath.sum([1. 2.; NaN 1.]) == 4.0
@test isnan(NaNMath.sum([NaN, NaN]))
@test NaNMath.sum(Float64[]) == 0.0
@test NaNMath.sum([1f0, 2f0, NaN32]) === 3.0f0
@test NaNMath.maximum([1., 2., NaN]) == 2.0
@test NaNMath.maximum([1. 2.; NaN 1.]) == 2.0
@test NaNMath.minimum([1., 2., NaN]) == 1.0
@test NaNMath.minimum([1. 2.; NaN 1.]) == 1.0
@test NaNMath.extrema([1., 2., NaN]) == (1.0, 2.0)
@test NaNMath.extrema([2., 1., NaN]) == (1.0, 2.0)
@test NaNMath.extrema([1. 2.; NaN 1.]) == (1.0, 2.0)
@test NaNMath.extrema([2. 1.; 1. NaN]) == (1.0, 2.0)
@test NaNMath.extrema([NaN, -1., NaN]) == (-1.0, -1.0)
@test NaNMath.mean([1., 2., NaN]) == 1.5
@test NaNMath.mean([1. 2.; NaN 3.]) == 2.0
@test NaNMath.var([1., 2., NaN]) == 0.5
@test NaNMath.std([1., 2., NaN]) == 0.7071067811865476

@test NaNMath.median([1.]) == 1.
@test NaNMath.median([1., NaN]) == 1.
@test NaNMath.median([NaN, 1., 3.]) == 2.
@test NaNMath.median([1., 3., 2., NaN]) == 2.
@test NaNMath.median([NaN, 1, 3]) == 2.
@test NaNMath.median([1, 2, NaN]) == 1.5
@test NaNMath.median([1 2; NaN NaN]) == 1.5
@test NaNMath.median([NaN 2; 1 NaN]) == 1.5
@test isnan(NaNMath.median(Float64[]))
@test isnan(NaNMath.median(Float32[]))
@test isnan(NaNMath.median([NaN]))

@test NaNMath.min(1, 2) == 1
@test NaNMath.min(1.0, 2.0) == 1.0
@test NaNMath.min(1, 2.0) == 1.0
@test NaNMath.min(BigFloat(1.0), 2.0) == BigFloat(1.0)
@test NaNMath.min(BigFloat(1.0), BigFloat(2.0)) == BigFloat(1.0)
@test NaNMath.min(NaN, 1) == 1.0
@test NaNMath.min(NaN32, 1) == 1.0f0
@test isnan(NaNMath.min(NaN, NaN))
@test isnan(NaNMath.min(NaN))
@test NaNMath.min(NaN, NaN, 0.0, 1.0) == 0.0

@test NaNMath.max(1, 2) == 2
@test NaNMath.max(1.0, 2.0) == 2.0
@test NaNMath.max(1, 2.0) == 2.0
@test NaNMath.max(BigFloat(1.0), 2.0) == BigFloat(2.0)
@test NaNMath.max(BigFloat(1.0), BigFloat(2.0)) == BigFloat(2.0)
@test NaNMath.max(NaN, 1) == 1.0
@test NaNMath.max(NaN32, 1) == 1.0f0
@test isnan(NaNMath.max(NaN, NaN))
@test isnan(NaNMath.max(NaN))
@test NaNMath.max(NaN, NaN, 0.0, 1.0) == 1.0

## Based on https://github.com/sethaxen/NaNMath.jl/blob/41b3e7edd9dd4cb6c2873abf6e0d90acf43138ec/test/runtests.jl
@testset "findmin/findmax" begin
    if VERSION ≥ v"1.7"
        xvals = [
            [1., 2., 3., 3., 1.],
            (1., 2., 3., 3., .1),
            (1f0, 2f0, 3f0, -1f0),
            (x=1.0, y=3f0, z=-4.0, w=-2f0),
            Dict(:a => 1.0, :b => 1.0, :d => 3.0, :c => 2.0),
        ]
        @testset for x in xvals
            @test NaNMath.findmin(x) === findmin(x)
            @test NaNMath.findmax(x) === findmax(x)
            @test NaNMath.findmin(identity, x) === findmin(identity, x)
            @test NaNMath.findmax(identity, x) === findmax(identity, x)
            @test NaNMath.findmin(sin, x) === findmin(sin, x)
            @test NaNMath.findmax(sin, x) === findmax(sin, x)
        end
    end

    x = [7, 7, NaN, 1, 1, NaN]
    @test NaNMath.findmin(x) === (1.0, 4)
    @test NaNMath.findmax(x) === (7.0, 1)
    @test NaNMath.findmin(identity, x) === (1.0, 4)
    @test NaNMath.findmax(identity, x) === (7.0, 1)
    @test NaNMath.findmin(-, x) === (-7.0, 1)
    @test NaNMath.findmax(-, x) === (-1.0, 4)

    x = [NaN, NaN]
    @test NaNMath.findmin(x) === (NaN, 1)
    @test NaNMath.findmax(x) === (NaN, 1)
    @test NaNMath.findmin(identity, x) === (NaN, 1)
    @test NaNMath.findmax(identity, x) === (NaN, 1)
    @test NaNMath.findmin(sin, x) === (NaN, 1)
    @test NaNMath.findmax(sin, x) === (NaN, 1)

    x = Dict(:a => 1.0, :b => 1 + 2im, :d => 3.0, :c => 2.0)
    @test_throws ErrorException NaNMath.findmin(x)
    @test_throws ErrorException NaNMath.findmax(x)

    x = [3, missing, NaN, -1]
    @test_throws ErrorException NaNMath.findmin(x)

    x = Dict('a' => 1.0, missing => NaN, 'c' => 2.0)
    @test NaNMath.findmin(x) === (1.0, 'a')
    @test NaNMath.findmax(x) === (2.0, 'c')

    x = Dict(:x => 3.0, :w => 2f0, :y => -1.0, :z => NaN)
    @test NaNMath.findmin(x) === (-1.0, :y)
    @test NaNMath.findmax(x) === (3.0, :x)
    @test NaNMath.findmin(identity, x) === (-1.0, :y)
    @test NaNMath.findmax(identity, x) === (3.0, :x)
    @test NaNMath.findmin(-, x) === (-3.0, :x)
    @test NaNMath.findmax(-, x) === (1.0, :y)
    @test NaNMath.findmin(exp, x) === (exp(-1.0), :y)
    @test NaNMath.findmax(exp, x) === (exp(3.0), :x)

    x = (x=1.0, y=NaN, z=NaN, w=-2.0)
    @test NaNMath.findmin(x) === (-2.0, :w)
    @test NaNMath.findmax(x) === (1.0, :x)
    @test NaNMath.findmin(-,x) === (-1.0, :x)
    @test NaNMath.findmax(-,x) === (2.0, :w)

    x = [2.0 3.0; 2.0 -1.0]
    @test NaNMath.findmin(x) === (-1.0, CartesianIndex(2, 2))
    @test NaNMath.findmax(x) === (3.0, CartesianIndex(1, 2))
    @test NaNMath.findmin(exp,x) === (exp(-1), CartesianIndex(2, 2))
    @test NaNMath.findmax(exp,x) === (exp(3.0), CartesianIndex(1, 2))
end

@testset "argmin/argmax" begin
    if VERSION ≥ v"1.7"
        xvals = [
            [1., 2., 4., 3., 1.],
            (1., 2., 4., 3., .1),
            (1f0, 2f0, 3f0, -1f0),
            (x=1.0, y=3f0, z=-4.0, w=-2f0),
            Dict(:a => 1.0, :b => 1.0, :d => 3.0, :c => 2.0),
        ]    
        @testset for x in xvals
            @test NaNMath.argmin(x) === argmin(x)
            @test NaNMath.argmax(x) === argmax(x)
            x isa Dict || @test NaNMath.argmin(identity, x) === argmin(identity, x)
            x isa Dict || @test NaNMath.argmax(identity, x) === argmax(identity, x)
            x isa Dict || @test NaNMath.argmin(sin, x) === argmin(sin, x)
            x isa Dict || @test NaNMath.argmax(sin, x) === argmax(sin, x)
        end
    end
    x = [7, 7, NaN, 1, 1, NaN]
    @test NaNMath.argmin(x) === 4
    @test NaNMath.argmax(x) === 1
    @test NaNMath.argmin(identity, x) === 1.0
    @test NaNMath.argmax(identity, x) === 7.0
    @test NaNMath.argmin(-, x) === 7.0
    @test NaNMath.argmax(-, x) === 1.0

    x = [NaN, NaN]
    @test NaNMath.argmin(x) === 1
    @test NaNMath.argmax(x) === 1
    @test NaNMath.argmin(identity, x) === NaN
    @test NaNMath.argmax(identity, x) === NaN
    @test NaNMath.argmin(-, x) === NaN
    @test NaNMath.argmax(-, x) === NaN

    x = [3, missing, NaN, -1]
    @test_throws ErrorException NaNMath.argmin(x)
    @test_throws ErrorException NaNMath.argmax(x)

    x = Dict('a' => 1.0, missing => NaN, 'c' => 2.0)
    @test NaNMath.argmin(x) === 'a'
    @test NaNMath.argmax(x) === 'c'

    x = Dict(:v => NaN, :w => 2.1f0, :x => 3.1, :z => -1.0, :y => NaN)
    @test NaNMath.argmin(x) === :z
    @test NaNMath.argmax(x) === :x
    @test NaNMath.argmin(-, x) === 3.1
    @test NaNMath.argmax(-, x) === -1.0
    @test NaNMath.argmin(exp, x) === -1.0
    @test NaNMath.argmax(exp, x) === 3.1

    x = (x=1.1, y=NaN, z=NaN, w=-2.3)
    @test NaNMath.argmin(x) === :w
    @test NaNMath.argmax(x) === :x
    @test NaNMath.argmin(exp, x) === -2.3
    @test NaNMath.argmax(exp, x) === 1.1

    x = [2.0 3.0; 2.0 -1.0]
    @test NaNMath.argmin(x) === CartesianIndex(2, 2)
    @test NaNMath.argmax(x) === CartesianIndex(1, 2)
    @test NaNMath.argmin(exp,x) === -1.0
    @test NaNMath.argmax(exp,x) === 3.0
end

# Test forwarding
x = 1 + 2im
for f in (:sin, :cos, :tan, :asin, :acos, :acosh, :atanh, :log, :log2, :log10,
          :log1p, :sqrt)
    @test @eval (NaNMath.$f)(x) == $f(x)
end
@test_throws MethodError NaNMath.lgamma(x)

struct A end
Base.isless(::A, ::A) = false
y = A()
for f in (:max, :min)
    @test @eval (NaNMath.$f)(y, y) == $f(y, y)
end
@test NaNMath.pow(x, x) == ^(x, x)

@testset "acosh" begin
    for T in (Float16, Float32, Float64)
        y = @inferred(NaNMath.acosh(T(0.5)))
        @test y isa T
        @test isnan(y)
        y = NaNMath.acosh(T(2.1))
        @test y isa T
        @test !isnan(y)
        @test y === acosh(T(2.1))
    end
end
