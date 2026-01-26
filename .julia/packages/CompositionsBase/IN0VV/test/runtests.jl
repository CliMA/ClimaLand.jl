module TestCompositionsBase
import Compat
using CompositionsBase
using CompositionsBase: @var_str
using Test

@testset "binary" begin
    @test tuple ∘ inv ===
          compose(tuple, inv) ===
          var"⨟"(inv, tuple) ===
          opcompose(inv, tuple)
end

@testset "unary" begin
    @test ∘(tuple) === compose(tuple) === var"⨟"(tuple) === opcompose(tuple) === tuple
end

if VERSION >= v"1.5.0-DEV.302"  # for ⨟
    include("test_julia15.jl")
end
if VERSION >= v"1.6"
    @testset "decompose, deopcompose" begin
        @test decompose(sin) === (sin,)
        @test decompose(sin∘cos) === (sin,cos)
        @test decompose(sin∘cos∘ tan) === (sin,cos,tan)
        @test decompose((sin∘cos)∘ tan) === (sin,cos,tan)
        @test decompose(sin∘(cos∘ tan)) === (sin,cos,tan)
        @test decompose((sqrt∘sin)∘(cos∘ tan)) === (sqrt,sin,cos,tan)

        @test deopcompose((sqrt∘sin)∘(cos∘ tan)) === (tan,cos,sin,sqrt)
        @inferred deopcompose((sqrt∘sin)∘(cos∘ tan))
        @inferred decompose((sqrt∘sin)∘(cos∘ tan))
    end
end
if VERSION >= v"1.9-"
    using InverseFunctions

    @testset "inverses" begin
        InverseFunctions.test_inverse(decompose, sin ∘ tan ∘ cos; compare= ==)
        InverseFunctions.test_inverse(deopcompose, sin ∘ tan ∘ cos; compare= ==)
        InverseFunctions.test_inverse(splat(compose), (sin, tan, cos); compare= ==)
        InverseFunctions.test_inverse(splat(opcompose), (sin, tan, cos); compare= ==)
    end
end

end  # module
