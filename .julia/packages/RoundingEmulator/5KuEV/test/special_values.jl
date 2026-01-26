for T in (Float32, Float64)
    @testset "$(T): Signed zero" begin
        @test isequal(add_up(zero(T), zero(T)), zero(T))
        @test isequal(add_up(zero(T), -zero(T)), zero(T))
        @test isequal(add_up(-zero(T), zero(T)), zero(T))
        @test isequal(add_up(-zero(T), -zero(T)), -zero(T))

        @test isequal(add_down(zero(T), zero(T)), zero(T))
        @test isequal(add_down(zero(T), -zero(T)), -zero(T))
        @test isequal(add_down(-zero(T), zero(T)), -zero(T))
        @test isequal(add_down(-zero(T), -zero(T)), -zero(T))

        @test isequal(sub_up(zero(T), zero(T)), zero(T))
        @test isequal(sub_up(zero(T), -zero(T)), zero(T))
        @test isequal(sub_up(-zero(T), zero(T)), -zero(T))
        @test isequal(sub_up(-zero(T), -zero(T)), zero(T))

        @test isequal(sub_down(zero(T), zero(T)), -zero(T))
        @test isequal(sub_down(zero(T), -zero(T)), zero(T))
        @test isequal(sub_down(-zero(T), zero(T)), -zero(T))
        @test isequal(sub_down(-zero(T), -zero(T)), -zero(T))

        @test isequal(mul_up(zero(T), zero(T)), zero(T))
        @test isequal(mul_up(zero(T), -zero(T)), -zero(T))
        @test isequal(mul_up(-zero(T), zero(T)), -zero(T))
        @test isequal(mul_up(-zero(T), -zero(T)), zero(T))

        @test isequal(mul_down(zero(T), zero(T)), zero(T))
        @test isequal(mul_down(zero(T), -zero(T)), -zero(T))
        @test isequal(mul_down(-zero(T), zero(T)), -zero(T))
        @test isequal(mul_down(-zero(T), -zero(T)), zero(T))

        @test isequal(div_up(zero(T), zero(T)), T(NaN))
        @test isequal(div_up(zero(T), -zero(T)), T(NaN))
        @test isequal(div_up(-zero(T), zero(T)), T(NaN))
        @test isequal(div_up(-zero(T), -zero(T)), T(NaN))

        @test isequal(div_down(zero(T), zero(T)), T(NaN))
        @test isequal(div_down(zero(T), -zero(T)), T(NaN))
        @test isequal(div_down(-zero(T), zero(T)), T(NaN))
        @test isequal(div_down(-zero(T), -zero(T)), T(NaN))

        @test isequal(sqrt_up(zero(T)), zero(T))
        @test isequal(sqrt_down(-zero(T)), -zero(T))
    end
end

@testset "Corner cases" begin
    # TODO
    # Add tests for Float32

    @testset "twosum intermediate overflow" begin
        # http://verifiedby.me/adiary/09
        a = 3.5630624444874539e+307
        b = -floatmax(Float64)
        x = a + b
        @test isfinite(x)
        tmp = x - a
        @test isinf(tmp)
        
        @test isequal(add_up(a, b), -1.4413868904135702e308)
        @test isequal(add_down(a, b), -1.4413868904135704e308) 
    end

    @testset "twoprod intermediate overflow" begin
        # http://verifiedby.me/adiary/09
        function split(a)
            tmp = a * (2.0^27 + 1.0)
            x = tmp - (tmp - a)
            y = a - x
            x, y
        end
        a = 6.929001713869936e+236
        b = 2.5944475251952003e+71
        x = a * b
        @test isfinite(x)
        a1, _ = split(a)
        b1, _ = split(a)
        tmp = a1 * b1
        @test isinf(tmp)

        @test isequal(mul_up(a, b), floatmax(Float64))
        @test isequal(mul_down(a, b), prevfloat(floatmax(Float64)))
    end
end
