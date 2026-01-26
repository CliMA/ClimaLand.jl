using IrrationalConstants
using Documenter
using Test

const ALLCONSTANTS = filter!(
    x -> x isa IrrationalConstants.IrrationalConstant,
    map(Base.Fix1(getproperty, IrrationalConstants), names(IrrationalConstants)),
)

@testset "k*pi" begin
    @test isapprox(2*pi, twoπ)
    @test isapprox(4*pi, fourπ)
    @test isapprox(pi/2, halfπ)
    @test isapprox(pi/4, quartπ)
end

@testset "k/pi" begin
    @test isapprox(1/pi, invπ)
    @test isapprox(2/pi, twoinvπ)
    @test isapprox(4/pi, fourinvπ)
end

@testset "1/(k*pi)" begin  
    @test isapprox(1/(2pi), inv2π)
    @test isapprox(1/(4pi), inv4π)
end

@testset "sqrt" begin
    @test isapprox(sqrt(2), sqrt2)
    @test isapprox(sqrt(3), sqrt3)
    @test isapprox(sqrt(pi), sqrtπ)
    @test isapprox(sqrt(2pi), sqrt2π)
    @test isapprox(sqrt(4pi), sqrt4π)
    @test isapprox(sqrt(pi/2), sqrthalfπ)
    @test isapprox(sqrt(1/2), invsqrt2)
    @test isapprox(sqrt(1/(pi)), invsqrtπ)
    @test isapprox(sqrt(1/(2pi)), invsqrt2π)
end

@testset "log" begin
    @test isapprox(log(1/2), loghalf)
    @test isapprox(log(2), logtwo)
    @test isapprox(log(10), logten)
    @test isapprox(log(pi), logπ)
    @test isapprox(log(2pi), log2π)
    @test isapprox(log(4pi), log4π)
end

@testset "type system" begin
    @test twoπ === IrrationalConstants.Twoπ()
    @test twoπ isa IrrationalConstants.IrrationalConstant
    @test twoπ isa AbstractIrrational
end

@testset "hash" begin
    for i in ALLCONSTANTS, j in ALLCONSTANTS
        @test isequal(i==j, hash(i)==hash(j))
    end
end

@testset "doctests" begin
    doctest(IrrationalConstants; manual=false)
end

# copied from https://github.com/JuliaLang/julia/blob/cf5ae0369ceae078cf6a29d7aa34f48a5a53531e/test/numbers.jl
# and adapted to irrationals in this package

@testset "IrrationalConstants zero and one" begin
    for i in ALLCONSTANTS
        @test one(i) === true
        @test zero(i) === false
        @test one(typeof(i)) === true
        @test zero(typeof(i)) === false
    end
end

@testset "IrrationalConstants iszero, isfinite, isinteger, and isone" begin
    for i in ALLCONSTANTS
        @test !iszero(i)
        @test !isone(i)
        @test !isinteger(i)
        @test isfinite(i)
    end
end

@testset "IrrationalConstants promote_type" begin
    for T in (Float16, Float32, Float64)
        for i in ALLCONSTANTS
            @test T(2.0) * i ≈ T(2.0) * T(i)
            @test T(2.0) * i isa T
        end
    end
end

@testset "IrrationalConstants compared with IrrationalConstants" begin
    for i in ALLCONSTANTS, j in ALLCONSTANTS
        @test isequal(i==j, Float64(i)==Float64(j))
        @test isequal(i!=j, Float64(i)!=Float64(j))
        @test isequal(i<=j, Float64(i)<=Float64(j))
        @test isequal(i>=j, Float64(i)>=Float64(j))
        @test isequal(i<j, Float64(i)<Float64(j))
        @test isequal(i>j, Float64(i)>Float64(j))
    end
end

@testset "IrrationalConstants Inverses, JuliaLang/Julia Issue #30882" begin
    @test @inferred(inv(twoπ)) ≈ 0.15915494309189535
end

@testset "IrrationalConstants compared with Rationals and Floats" begin
    for i in ALLCONSTANTS
        @test Float64(i, RoundDown) < i
        @test Float64(i, RoundUp) > i
        @test !(Float64(i, RoundDown) > i)
        @test !(Float64(i, RoundUp) < i)
        @test Float64(i, RoundDown) <= i
        @test Float64(i, RoundUp) >= i
        @test Float64(i, RoundDown) != i
        @test Float64(i, RoundUp) != i

        @test Float32(i, RoundDown) < i
        @test Float32(i, RoundUp) > i
        @test !(Float32(i, RoundDown) > i)
        @test !(Float32(i, RoundUp) < i)

        @test prevfloat(big(i)) < i
        @test nextfloat(big(i)) > i
        @test !(prevfloat(big(i)) > i)
        @test !(nextfloat(big(i)) < i)
    end

    @test 5293386250278608690//842468587426513207 < twoπ
    @test !(5293386250278608690//842468587426513207 > twoπ)
    @test 5293386250278608690//842468587426513207 != twoπ
end
IrrationalConstants.@irrational i46051 4863.185427757 1548big(pi)
@testset "IrrationalConstant printing" begin
    @test sprint(show, "text/plain", twoπ) == "twoπ = 6.2831853071795..."
    @test sprint(show, "text/plain", twoπ, context=:compact => true) == "twoπ"
    @test sprint(show, twoπ) == "twoπ"
    # JuliaLang/Julia issue #46051
    @test sprint(show, "text/plain", i46051) == "i46051 = 4863.185427757..."
end

@testset "IrrationalConstant/Bool multiplication" begin
    @test false*twoπ === 0.0
    @test twoπ*false === 0.0
    @test true*twoπ === Float64(twoπ)
    @test twoπ*true === Float64(twoπ)
end

# JuliaLang/Julia issue #26324
@testset "irrational promotion" begin
    @test twoπ*ComplexF32(2) isa ComplexF32
    @test twoπ/ComplexF32(2) isa ComplexF32
    @test log(twoπ, ComplexF32(2)) isa ComplexF32
end

# issue #23
@testset "rounding irrationals" begin
    # without rounding mode
    @test @inferred(round(twoπ)) == 6.0
    @test @inferred(round(sqrt2)) == 1.0
    @test @inferred(round(sqrt3)) == 2.0
    @test @inferred(round(loghalf)) == -1.0

    # with rounding modes
    for mode in (RoundDown, RoundToZero, RoundNearest, RoundNearestTiesAway, RoundNearestTiesUp)
        @test @inferred(round(twoπ, mode)) == 6.0
        @test @inferred(round(sqrt2, mode)) == 1.0
    end
    @test @inferred(round(sqrt3, RoundUp)) == 2.0
    for mode in (RoundUp, RoundToZero)
        @test @inferred(round(loghalf, mode)) == 0.0
    end
end

@testset "trigonometric functions" begin
    # 2π, 4π
    for (n, x) in ((2, twoπ), (4, fourπ))
        @test sin(x) === sinpi(n) === sin(0.0)
        @test cos(x) === cospi(n) === cos(0.0)
    end

    # halfπ, quartπ
    for (r, x) in ((big"0.5", halfπ), (big"0.25", quartπ))
        @test sin(x) === Float64(sinpi(r))
        @test cos(x) === Float64(cospi(r))
    end

    # Check consistency of definitions
    for x in (twoπ, fourπ, halfπ, quartπ)
        @test sincos(x) === (sin(x), cos(x))
        @test tan(x) === sin(x) / cos(x)
    end

    # Check `csc`, `sec`, and `cot`
    for x in (twoπ, fourπ, halfπ)
        # These are defined automatically via sin, cos, and tan
        @test csc(x) === 1 / sin(x)
        @test sec(x) === 1 / cos(x)
        @test cot(x) === csc(x) / sec(x)
    end
    @test csc(quartπ) === Float64(csc(big(quartπ)))
    @test sec(quartπ) === Float64(sec(big(quartπ)))
    @test cot(quartπ) === Float64(cot(big(quartπ)))
end

# Ref https://github.com/JuliaLang/julia/pull/46054
IrrationalConstants.@irrational irrational_1548_pi 4863.185427757 1548big(pi)
IrrationalConstants.@irrational irrational_inv_1548_pi 1/big(irrational_1548_pi)
@testset "IrrationalConstants.@irrational" begin
    @test irrational_1548_pi ≈ 1548big(pi)
    @test Float64(irrational_1548_pi) == 1548π
    @test irrational_1548_pi ≈ 1548pi
    @test irrational_1548_pi != 1548pi
    @test irrational_inv_1548_pi ≈ inv(1548big(pi))
    @test Float64(irrational_inv_1548_pi) == 1/(1548π)
    @test irrational_inv_1548_pi ≈ inv(1548pi)
    @test irrational_inv_1548_pi != inv(1548pi)
end

# Ref https://github.com/JuliaLang/julia/pull/50894
@testset "irrational special values" begin
    for v ∈ ALLCONSTANTS
        @test v === typemin(v) === typemax(v)
    end
end

# Ref https://github.com/JuliaLang/julia/pull/55911
@testset "logtwo to `BigFloat` with `setrounding`" begin
    function irrational_to_big_float(c::AbstractIrrational)
        BigFloat(c)
    end

    function irrational_to_big_float_with_rounding_mode(c::AbstractIrrational, rm::RoundingMode)
        f = () -> irrational_to_big_float(c)
        setrounding(f, BigFloat, rm)
    end

    function irrational_to_big_float_with_rounding_mode_and_precision(c::AbstractIrrational, rm::RoundingMode, prec::Int)
        f = () -> irrational_to_big_float_with_rounding_mode(c, rm)
        setprecision(f, BigFloat, prec)
    end

    # logtwo is the only constant defined based on an MPFR constant (similar to π, γ, catalan)
    c = logtwo
    for p ∈ 1:40
        @test (
            irrational_to_big_float_with_rounding_mode_and_precision(c, RoundDown, p) < c <
            irrational_to_big_float_with_rounding_mode_and_precision(c, RoundUp, p)
        )
    end
end

# issues #37 and #40
@testset "slow comparisons" begin
    @test iszero(@allocated(3.0 <= invsqrt2))
end

# issues #43
@testset "macro error" begin
    msg = "The name of the irrational constant (Myπ) and its type (Myπ) cannot be the same. Please choose a different name for the constant or specify a different type name as the last argument to the macro."
    @test_throws ArgumentError(msg) @macroexpand(IrrationalConstants.@irrational Myπ big(π))
    @test_throws ArgumentError(msg) @macroexpand(IrrationalConstants.@irrational Myπ 1.0 big(π))
    @test_throws ArgumentError(msg) @macroexpand(IrrationalConstants.@irrational Myπ big(π) Myπ)
    @test_throws ArgumentError(msg) @macroexpand(IrrationalConstants.@irrational Myπ 1.0 big(π) Myπ)
end

# test that defining a type that already exists throws an error
module TestTypeCollision
    using IrrationalConstants
    using Test
    struct MyExistingType end

    @testset "type collision" begin
        msg1 = "Type `MyExistingType` of irrational constant `myExistingType` is already defined in module `Main.TestTypeCollision`."
        @test_throws ArgumentError(msg1) @macroexpand(IrrationalConstants.@irrational myExistingType big(π))
        msg2 = "Type `MyExistingType` of irrational constant `myconst` is already defined in module `Main.TestTypeCollision`."
        @test_throws ArgumentError(msg2) @macroexpand(IrrationalConstants.@irrational myconst big(π) MyExistingType)
        @test_throws ArgumentError(msg2) @macroexpand(IrrationalConstants.@irrational myconst 1.0 big(π) MyExistingType)
    end
end
