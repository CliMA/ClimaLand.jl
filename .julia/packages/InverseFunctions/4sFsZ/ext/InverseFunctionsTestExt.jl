module InverseFunctionsTestExt

using Test: @test, @testset
using InverseFunctions: InverseFunctions, inverse

function InverseFunctions.test_inverse(f, x; compare=isapprox, kwargs...)
    @testset "test_inverse: $f with input $x" begin
        y = f(x)
        inverse_f = inverse(f)
        @test compare(inverse_f(y), x; kwargs...)
        inverse_inverse_f = inverse(inverse_f)
        @test compare(inverse_inverse_f(x), y; kwargs...)
    end
    return nothing
end

end
