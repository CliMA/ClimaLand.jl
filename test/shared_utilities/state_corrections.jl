using Test
import ClimaComms
ClimaComms.@import_required_backends
using ClimaLand
using ClimaLand: LogLinearFactor
using ClimaLand.Domains: Point
import ClimaCore.Fields

for FT in (Float32, Float64)
    @testset "LogLinearFactor, FT = $FT" begin
        # Zero weights and bias: identity
        f0 = LogLinearFactor{FT}(; weights = (0, 0))
        @test f0(FT(1.3), FT(-2)) == FT(1)
        @test eltype(f0) == FT

        # Standardization and bias
        f = LogLinearFactor{FT}(;
            weights = [0.5, -1.0],
            center = [2.0, 1.0],
            scale = [2.0, 0.5],
            bias = 0.1,
        )
        x1, x2 = FT(4), FT(1.5)
        z = 0.1 + 0.5 * (4 - 2) / 2 - 1.0 * (1.5 - 1) / 0.5
        @test abs(z) < log(3)
        @test f(x1, x2) ≈ FT(exp(z))
        @test f(FT(2), FT(1)) ≈ FT(exp(0.1))

        # Bounds: |log f| ≤ log_fmax
        fb = LogLinearFactor{FT}(; weights = (10,), log_fmax = log(3))
        @test fb(FT(100)) ≈ FT(3)
        @test fb(FT(-100)) ≈ FT(1 / 3)
        @test fb(FT(0)) == FT(1)

        # Bad scales are rejected
        @test_throws AssertionError LogLinearFactor{FT}(;
            weights = (1,),
            scale = (0,),
        )
        @test_throws AssertionError LogLinearFactor{FT}(;
            weights = (1, 1),
            center = (0,),
        )

        # Broadcasting over fields
        domain = Point(; z_sfc = FT(0))
        a = Fields.ones(domain.space.surface) .* FT(4)
        b = Fields.ones(domain.space.surface) .* FT(1.5)
        g = @. f(a, b)
        @test all(parent(g) .≈ FT(exp(z)))
        g2 = @. min(1, a / 10 * f0(a, b))
        @test all(parent(g2) .≈ FT(0.4))
    end
end
