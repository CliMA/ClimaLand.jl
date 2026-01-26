using LinearAlgebra, Statistics, SpecialFunctions
using ColorVectorSpace, Colors, FixedPointNumbers

using Test

using ColorVectorSpace: _mul

n8sum(x,y) = Float64(N0f8(x)) + Float64(N0f8(y))

macro test_colortype_approx_eq(a, b)
    :(test_colortype_approx_eq($(esc(a)), $(esc(b)), $(string(a)), $(string(b))))
end

function test_colortype_approx_eq(a::Colorant, b::Colorant, astr, bstr)
    @test typeof(a) == typeof(b)
    n = length(fieldnames(typeof(a)))
    for i = 1:n
        @test getfield(a, i) ≈ getfield(b,i)
    end
end

struct RatRGB <: AbstractRGB{Rational{Int}}
    r::Rational{Int}
    g::Rational{Int}
    b::Rational{Int}
end
ColorTypes.red(c::RatRGB)   = c.r
ColorTypes.green(c::RatRGB) = c.g
ColorTypes.blue(c::RatRGB)  = c.b

struct RGBA32 <: AbstractRGBA{RGB24, N0f8}
    color::UInt32
    RGBA32(c::UInt32, ::Type{Val{true}}) = new(c)
end
function RGBA32(r, g, b, alpha=1N0f8)
    u32 = reinterpret(UInt32, ARGB32(r, g, b, alpha))
    RGBA32((u32 << 0x8) | (u32 >> 0x18), Val{true})
end
ColorTypes.red(  c::RGBA32) = reinterpret(N0f8, (c.color >> 0x18) % UInt8)
ColorTypes.green(c::RGBA32) = reinterpret(N0f8, (c.color >> 0x10) % UInt8)
ColorTypes.blue( c::RGBA32) = reinterpret(N0f8, (c.color >> 0x08) % UInt8)
ColorTypes.alpha(c::RGBA32) = reinterpret(N0f8, c.color % UInt8)
ColorTypes.comp4(c::RGBA32) = alpha(c)

struct GrayA32 <: AbstractGrayA{Gray24, N0f8}
    color::UInt32
    GrayA32(c::UInt32, ::Type{Val{true}}) = new(c)
end
GrayA32(g, alpha=1N0f8) = GrayA32(bswap(reinterpret(UInt32, AGray32(g, alpha))), Val{true})
ColorTypes.gray( c::GrayA32) = reinterpret(N0f8, (c.color >> 0x18) % UInt8)
ColorTypes.alpha(c::GrayA32) = reinterpret(N0f8, c.color % UInt8)
ColorTypes.comp2(c::RGBA32) = alpha(c)

@testset "Colortypes" begin
    @testset "ambiguities" begin
        @test isempty(detect_ambiguities(ColorVectorSpace))
    end

    @testset "MathTypes" begin
        @test Gray{Float32} <: ColorVectorSpace.MathTypes{Float32} <: ColorVectorSpace.MathTypes
        @test AGray{Float32} <: ColorVectorSpace.MathTypes{Float32}
        @test GrayA{Float32} <: ColorVectorSpace.MathTypes{Float32}
        @test RGB{Float32} <: ColorVectorSpace.MathTypes{Float32}
        @test RGBA{Float32} <: ColorVectorSpace.MathTypes{Float32}
        @test ARGB{Float32} <: ColorVectorSpace.MathTypes{Float32}
        @test !(HSV{Float32} <: ColorVectorSpace.MathTypes)
        @test !(HSVA{Float32} <: ColorVectorSpace.MathTypes)
        @test !(AHSV{Float32} <: ColorVectorSpace.MathTypes)
        @test !(HSV <: ColorVectorSpace.MathTypes)
        @test !(HSVA <: ColorVectorSpace.MathTypes)
        @test !(AHSV <: ColorVectorSpace.MathTypes)
        @test AbstractGray <: ColorVectorSpace.MathTypes
        @test AbstractRGB <: ColorVectorSpace.MathTypes
        @test AbstractGray{Float32} <: ColorVectorSpace.MathTypes{Float32}
        @test AbstractRGB{Float32} <: ColorVectorSpace.MathTypes{Float32}
        @test_broken TransparentGray <: ColorVectorSpace.MathTypes
        @test_broken TransparentRGB <: ColorVectorSpace.MathTypes
    end

    @testset "convert" begin
        for x in (0.5, 0.5f0, NaN, NaN32, N0f8(0.5))
            @test @inferred(convert(Gray{typeof(x)}, x))  === @inferred(convert(Gray, x))  === Gray(x)
            @test @inferred(convert(RGB{typeof(x)}, x))   === @inferred(convert(RGB, x))   === RGB(x, x, x)
            # These should be fixed by a future release of ColorTypes
            @test @inferred(convert(AGray{typeof(x)}, x)) === @inferred(convert(AGray, x)) === AGray(x, 1)
            @test @inferred(convert(ARGB{typeof(x)}, x))  === @inferred(convert(ARGB, x))  === ARGB(x, x, x, 1)
            @test @inferred(convert(GrayA{typeof(x)}, x)) === @inferred(convert(GrayA, x)) === GrayA(x, 1)
            @test @inferred(convert(RGBA{typeof(x)}, x))  === @inferred(convert(RGBA, x))  === RGBA(x, x, x, 1)
        end
    end

    @testset "nan" begin
        function make_checked_nan(::Type{T}) where T
            x = nan(T)
            isa(x, T) && mapreducec(isnan, &, true, x)
        end
        for S in (Float32, Float64)
            @test make_checked_nan(S)
            @test make_checked_nan(Gray{S})
            @test make_checked_nan(AGray{S})
            @test make_checked_nan(GrayA{S})
            @test make_checked_nan(RGB{S})
            @test make_checked_nan(ARGB{S})
            @test make_checked_nan(RGBA{S})
        end
    end

    @testset "traits" begin
        @test floattype(Gray{N0f8}) === Gray{float(N0f8)}
    end

    @testset "_mul" begin
        n0f8p = ((x, y) for x = N0f8(0):eps(N0f8):N0f8(1), y = N0f8(0):eps(N0f8):N0f8(1))
        @test all(((x, y),) -> _mul(x, y) === N0f8(Float32(x) * Float32(y)), n0f8p)
        n0f16a = N0f16(0):eps(N0f16):N0f16(1)
        n0f16b = reinterpret(N0f16, 0xBADb)
        @test all(a -> _mul(a, n0f16b) === N0f16(Float64(a) * 0xBADb / 0xFFFF), n0f16a)
        @test _mul(0.6N0f8, 0.2N0f16) === _mul(0.2N0f16, 0.6N0f8) === 0.12N0f16
    end

    @testset "Arithmetic with Gray" begin
        cf = Gray{Float32}(0.1)
        @test @inferred(+cf) === cf
        @test @inferred(-cf) === Gray(-0.1f0)
        @test @inferred(one(cf)*cf) === cf
        @test oneunit(cf) === Gray(1.0f0)
        ccmp = Gray{Float32}(0.2)
        @test @inferred(2*cf) === cf*2 === 2.0f0*cf === cf*2.0f0 === ccmp
        @test @inferred(ccmp/2) === cf
        @test @inferred(cf*cf) === Gray{Float32}(0.1f0*0.1f0)
        @test @inferred(Gray{N0f32}(0.5)*Gray(0.5f0)) === Gray(Float64(N0f32(0.5)) * 0.5)
        @test @inferred(cf^2 ) === Gray{Float32}(0.1f0*0.1f0)
        @test @inferred(cf^3.0f0) === Gray{Float32}(0.1f0^3.0f0)
        @test @inferred(2.0*cf) === cf*2.0 === Gray(2.0*0.1f0)
        cf64 = Gray(0.2)
        @test cf / cf64 === Gray(0.1f0/0.2)
        @test_throws MethodError cf ÷ cf
        @test cf + 0.1     === 0.1 + cf        === Gray(Float64(0.1f0) + 0.1)
        @test cf64 - 0.1f0 === -(0.1f0 - cf64) === Gray( 0.2 - Float64(0.1f0))
        @test abs2(ccmp) === abs2(gray(ccmp))
        @test norm(cf) == norm(cf, 2) == norm(gray(cf))
        @test norm(cf, 1)   == norm(gray(cf), 1)
        @test norm(cf, Inf) == norm(gray(cf), Inf)
        @test @inferred(abs(cf)) === Gray(0.1f0)
        cu = Gray{N0f8}(0.1)
        @test @inferred(2*cu) === cu*2 === Gray(2*gray(cu))
        @test @inferred(2.0f0*cu) === cu*2.0f0 === Gray(2.0f0*gray(cu))
        f = N0f8(0.5)
        @test @inferred(gray(f*cu)) === gray(cu*f) ===f*gray(cu)
        @test @inferred(cf/2.0f0) === Gray{Float32}(0.05)
        @test @inferred(cu/2) === Gray(cu.val/2)
        @test @inferred(cu/0.5f0) === Gray(cu.val/0.5f0)
        @test @inferred(cu/0.1N0f8) === Gray{N0f8}(1) # issue #154
        @test_colortype_approx_eq @inferred(cf/0.6N0f8) Gray{Float32}(1/6) # issue #154
        @test @inferred(Gray{Q0f7}(0.25) / 0.5Q0f7) === Gray{Q0f7}(0.5) # issue #154
        @test @inferred(1+0im / Gray(1)) === @inferred(1 / Gray(1)) === Gray{Float32}(1)
        @test @inferred(cf+cf) === ccmp
        @test isfinite(cf)
        @test isfinite(Gray(true))
        @test !isinf(cf)
        @test !isinf(Gray(f))
        @test !isnan(cf)
        @test !isfinite(Gray(NaN))
        @test !isinf(Gray(NaN))
        @test isnan(Gray(NaN))
        @test !isfinite(Gray(Inf))
        @test Gray(Inf) |> isinf
        @test !isnan(Gray(Inf))
        @test abs(Gray(0.1)) === Gray(0.1)
        @test eps(Gray{N0f8}) === Gray(eps(N0f8))  # #282
        @test atan(Gray(0.1), Gray(0.2)) == atan(0.1, 0.2)
        @test hypot(Gray(0.2), Gray(0.3)) === hypot(0.2, 0.3)
        # Multiplication
        @test cf ⋅ cf   === gray(cf)^2
        @test cf ⋅ cf64 === gray(cf) * gray(cf64)
        @test cf ⊙ cf   === Gray(gray(cf)^2)
        @test cf ⊙ cf64 === Gray(gray(cf) * gray(cf64))
        @test cf ⊗ cf   === Gray(gray(cf)^2)
        @test cf ⊗ cf64 === Gray(gray(cf) * gray(cf64))

        acu = Gray{N0f8}[cu]
        acf = Gray{Float32}[cf]
        @test @inferred(acu./trues(1)) == acu
        @test typeof(acu./trues(1)) == Vector{typeof(cu/true)}
        @test @inferred(ones(Int, 1)./acu) == [1/cu]
        @test typeof(ones(Int, 1)./acu) == Vector{typeof(1/cu)}
        @test @inferred(acu./acu) == [1]
        @test typeof(acu./acu) == Vector{typeof(cu/cu)}
        @test typeof(acu+acf) == Vector{Gray{Float32}}
        @test typeof(acu-acf) == Vector{Gray{Float32}}
        @test typeof(acu.+acf) == Vector{Gray{Float32}}
        @test typeof(acu.-acf) == Vector{Gray{Float32}}
        @test typeof(acu.+cf) == Vector{Gray{Float32}}
        @test typeof(acu.-cf) == Vector{Gray{Float32}}
        @test typeof(2*acf) == Vector{Gray{Float32}}
        @test typeof(2 .* acf) == Vector{Gray{Float32}}
        @test typeof(0x02*acu) == Vector{Gray{Float32}}
        @test typeof(acu/2) == Vector{Gray{typeof(N0f8(0.5)/2)}}
        @test typeof(acf.^2) == Vector{Gray{Float32}}
        @test (acu./Gray{N0f8}(0.5))[1] == gray(acu[1])/N0f8(0.5)
        @test (acf./Gray{Float32}(2))[1] ≈ 0.05f0
        @test (acu/2)[1] == Gray(gray(acu[1])/2)
        @test (acf/2)[1] ≈ Gray{Float32}(0.05f0)

        @test gray(0.8) === 0.8

        a = Gray{N0f8}[0.8,0.7]
        @test a == a
        @test a === a
        @test isapprox(a, a)
        @test sum(a) == Gray(n8sum(0.8,0.7))
        @test sum(a[1:1]) == a[1]
        @test abs( varmult(*, a) - (a[1]-a[2])^2 / 2 ) <= 0.001
        ab = Gray[true, true]
        @test sum(ab) === Gray(n8sum(true, true))
        @test prod(ab) === Gray(true)
        ab = Gray[true]
        @test sum(ab) === Gray(n8sum(true, false))
        @test prod(ab) === Gray(true)
        @test sum(Gray{Bool}[]) === Gray(n8sum(false, false))
        @test prod(Gray{Bool}[]) === one(Gray{Bool})

        @test real(Gray{Float32}) <: Real
        @test zero(ColorTypes.Gray) == 0
        @test oneunit(ColorTypes.Gray) == 1

        @test typeof(float(Gray{N0f16}(0.5))) <: AbstractFloat
        @test quantile( Gray{N0f16}[0.0,0.5,1.0], 0.1) ≈ 0.1 atol=eps(N0f16)
        @test middle(Gray(0.2)) === Gray(0.2)
        @test middle(Gray(0.2), Gray(0.4)) === Gray((0.2+0.4)/2)

        # issue #56
        @test Gray24(0.8)*N0f8(0.5) === Gray{N0f8}(0.4)
        @test Gray24(0.8)*0.5 === Gray(0.4)
        @test Gray24(0.8)/2   === Gray(0.5f0*N0f8(0.8))
        @test Gray24(0.8)/2.0 === Gray(0.4)

        # issue #133
        @test Gray24(0.2) + Gray24(0.4) === Gray24(0.6)
        @test Gray24(1)   - Gray24(0.2) === Gray24(0.8)
        @test Gray24(1)   * Gray24(0.2) === Gray24(0.2)
    end

    @testset "Comparisons with Gray" begin
        g1 = Gray{N0f8}(0.2)
        g2 = Gray{N0f8}(0.3)
        @test isless(g1, g2)
        @test !(isless(g2, g1))
        @test g1 < g2
        @test !(g2 < g1)
        @test isless(g1, 0.5)
        @test !(isless(0.5, g1))
        @test g1 < 0.5
        @test !(0.5 < g1)
        @test @inferred(max(g1, g2)) === g2
        @test @inferred(max(g1, Gray(0.3))) === Gray(0.3)
        @test max(g1, 0.1) === max(0.1, g1) === Float64(gray(g1))
        @test (@inferred(min(g1, g2)) ) == g1
        @test min(g1, 0.1) === min(0.1, g1) === 0.1
        a = Gray{Float64}(0.9999999999999999)
        b = Gray{Float64}(1.0)

        @test (Gray(0.3) < Gray(NaN)) == (0.3 < NaN)
        @test (Gray(NaN) < Gray(0.3)) == (NaN < 0.3)
        @test isless(Gray(0.3), Gray(NaN)) == isless(0.3, NaN)
        @test isless(Gray(NaN), Gray(0.3)) == isless(NaN, 0.3)
        @test isless(Gray(0.3), NaN) == isless(0.3, NaN)
        @test isless(Gray(NaN), 0.3) == isless(NaN, 0.3)
        @test isless(0.3, Gray(NaN)) == isless(0.3, NaN)
        @test isless(NaN, Gray(0.3)) == isless(NaN, 0.3)

        @test isapprox(a, b)
        a = Gray{Float64}(0.99)
        @test !(isapprox(a, b, rtol = 0.01))
        @test isapprox(a, b, rtol = 0.1)
    end

    @testset "Unary operations with Gray" begin
        ntested = 0
        for g in (Gray(0.4), Gray{N0f8}(0.4))
            @test @inferred(zero(g)) === typeof(g)(0)
            @test @inferred(oneunit(g)) === typeof(g)(1)
            SFE = if isdefined(Base, :get_extension)
                Base.get_extension(ColorVectorSpace, :SpecialFunctionsExt)
            else
                ColorVectorSpace.SpecialFunctionsExt
            end
            for mod in (
                ColorVectorSpace,
                SFE,
                (; unaryops = (:trunc, :floor, :round, :ceil, :eps, :bswap)),
            )
                for op in mod.unaryops
                    op ∈ (:frexp, :exponent, :modf, :logfactorial) && continue
                    op === :~ && eltype(g) === Float64 && continue
                    op === :significand && eltype(g) === N0f8 && continue
                    try
                        v = @eval $op(gray($g))  # if this fails, don't bother with the next test
                        @test @eval($op($g)) === Gray(v)
                        ntested += 1
                    catch ex
                        @test ex isa Union{DomainError,ArgumentError}
                    end
                end
            end
        end
        @test ntested > 130
        @test logabsgamma(Gray(0.2)) == (Gray(logabsgamma(0.2)[1]), 1)
        for g in (Gray{N0f8}(0.4), Gray{N0f8}(0.6))
            for op in (:trunc, :floor, :round, :ceil)
                v = @eval $op(Bool, gray($g))
                @test @eval($op(Bool, $g)) === v
            end
        end
        for g in (Gray(1.4), Gray(1.6))
            for op in (:trunc, :floor, :round, :ceil)
                v = @eval $op(Int, gray($g))
                @test @eval($op(Int, $g)) === v
            end
        end
        for (g1, g2) in ((Gray(0.4), Gray(0.3)), (Gray(N0f8(0.4)), Gray(N0f8(0.3))))
            for op in (:mod, :rem, :mod1)
                v = @eval $op(gray($g1), gray($g2))
                @test @eval($op($g1, $g2)) === Gray(v)
            end
        end
        u = N0f8(0.4)
        @test ~Gray(u) == Gray(~u)
        @test -Gray(u) == Gray(-u)
    end

    @testset "Arithmetic with TransparentGray" begin
        p1 = GrayA{Float32}(Gray(0.8), 0.2)
        @test @inferred(zero(p1)) === GrayA{Float32}(0,0)
        @test @inferred(oneunit(p1)) === GrayA{Float32}(1,1)
        @test +p1 == p1
        @test -p1 == GrayA(-0.8f0, -0.2f0)
        p2 = GrayA{Float32}(Gray(0.6), 0.3)
        @test_colortype_approx_eq p1+p2 GrayA{Float32}(Gray(1.4),0.5)
        @test_colortype_approx_eq (p1+p2)/2 GrayA{Float32}(Gray(0.7),0.25)
        @test_colortype_approx_eq 0.4f0*p1+0.6f0*p2 GrayA{Float32}(Gray(0.68),0.26)
        @test_colortype_approx_eq ([p1]+[p2])[1] GrayA{Float32}(Gray(1.4),0.5)
        @test_colortype_approx_eq ([p1].+[p2])[1] GrayA{Float32}(Gray(1.4),0.5)
        @test_colortype_approx_eq ([p1].+p2)[1] GrayA{Float32}(Gray(1.4),0.5)
        @test_colortype_approx_eq ([p1]-[p2])[1] GrayA{Float32}(Gray(0.2),-0.1)
        @test_colortype_approx_eq ([p1].-[p2])[1] GrayA{Float32}(Gray(0.2),-0.1)
        @test_colortype_approx_eq ([p1].-p2)[1] GrayA{Float32}(Gray(0.2),-0.1)
        @test_colortype_approx_eq ([p1]/2)[1] GrayA{Float32}(Gray(0.4),0.1)
        @test_colortype_approx_eq (0.4f0*[p1]+0.6f0*[p2])[1] GrayA{Float32}(Gray(0.68),0.26)

        cf = AGray{Float32}(0.8, 0.2)
        cu = AGray{N0f8}(0.8, 0.2)
        @test @inferred(cu / 0.8N0f8) === AGray{N0f8}(1, 0.25) # issue #154
        @test_colortype_approx_eq @inferred(cf / 0.8N0f8) AGray{Float32}(1, 0.25) # issue #154
        @test @inferred(AGray{Q0f7}(0.25, 0.125) / 0.5Q0f7) === AGray{Q0f7}(0.5, 0.25) # issue #154

        a = GrayA{N0f8}[GrayA(0.8,0.7), GrayA(0.5,0.2)]
        @test sum(a) == GrayA(n8sum(0.8,0.5), n8sum(0.7,0.2))
        @test isapprox(a, a)
        a = AGray{Float64}(1.0, 0.9999999999999999)
        b = AGray{Float64}(1.0, 1.0)

        @test a ≈ b
        a = AGray{Float64}(1.0, 0.99)
        @test !isapprox(a, b, rtol = 0.01)
        @test isapprox(a, b, rtol = 0.1)

        # issue #56
        @test AGray32(0.8,0.2)*N0f8(0.5) === AGray{N0f8}(0.4,0.1)
        @test AGray32(0.8,0.2)*0.5 === AGray(0.4,0.1)
        @test AGray32(0.8,0.2)/2   === AGray(0.5f0*N0f8(0.8),0.5f0*N0f8(0.2))
        @test AGray32(0.8,0.2)/2.0 === AGray(0.4,0.1)

        # issue #133
        @test AGray32(1, 0.4) - AGray32(0.2, 0.2) === AGray32(0.8, 0.2)

        # issue #146
        @test @inferred(GrayA32(0.8,0.2)*N0f8(0.5)) === GrayA{N0f8}(0.4,0.1)
        @test @inferred(GrayA32(0.8,0.2)*0.5) === GrayA(0.4,0.1)
        @test @inferred(GrayA32(0.8,0.2)/2)   === GrayA(0.5f0*N0f8(0.8),0.5f0*N0f8(0.2))
        @test @inferred(GrayA32(0.8,0.2)/2.0) === GrayA(0.4,0.1)
        @test @inferred(GrayA32(1, 0.4) - GrayA32(0.2, 0.2)) === GrayA32(0.8, 0.2)

        # Multiplication
        @test_throws MethodError cf * cf
        @test_throws MethodError cf ⋅ cf
        @test_throws MethodError cf ⊗ cf
        cf64 = mapc(Float64, cf)
        @test @inferred(cf ⊙ cf)   === AGray{Float32}(0.8f0^2, 0.2f0^2)
        @test @inferred(cf ⊙ cf64) === AGray{Float64}(0.8f0*(0.8f0*1.0), 0.2f0*(0.2f0*1.0))
    end

    @testset "Arithemtic with RGB" begin
        cf = RGB{Float32}(0.1,0.2,0.3)
        @test @inferred(zero(cf)) === RGB{Float32}(0,0,0)
        @test @inferred(oneunit(cf)) === RGB{Float32}(1,1,1)
        @test +cf == cf
        @test -cf == RGB(-0.1f0, -0.2f0, -0.3f0)
        ccmp = RGB{Float32}(0.2,0.4,0.6)
        @test 2*cf == ccmp
        @test cf*2 == ccmp
        @test ccmp/2 == cf
        @test 2.0f0*cf == ccmp
        @test eltype(2.0*cf) == Float64
        cu = RGB{N0f8}(0.1,0.2,0.3)
        @test_colortype_approx_eq 2*cu RGB(2*cu.r, 2*cu.g, 2*cu.b)
        @test_colortype_approx_eq 2.0f0*cu RGB(2.0f0*cu.r, 2.0f0*cu.g, 2.0f0*cu.b)
        f = N0f8(0.5)
        @test (f*cu).r ≈ f*cu.r
        @test cf/2.0f0 == RGB{Float32}(0.05,0.1,0.15)
        @test cu/2 ≈ RGB(cu.r/2,cu.g/2,cu.b/2)
        @test cu/0.5f0 ≈ RGB(cu.r/0.5f0, cu.g/0.5f0, cu.b/0.5f0)
        @test @inferred(cu/0.4N0f8) === RGB{N0f8}(26/102, 51/102, 76/102)
        @test_colortype_approx_eq @inferred(cf / 0.4N0f8) RGB{Float32}(0.25, 0.5, 0.75)
        cq0f7 = RGB{Q0f7}(0.125, 0.25, 0.375)
        @test @inferred(cq0f7 / 0.5Q0f7) === RGB{Q0f7}(0.25, 0.5, 0.75) # issue #154
        @test cf+cf == ccmp
        @test cu * 1//2 == mapc(x->Float64(Rational(x)/2), cu)
        @test_colortype_approx_eq (cf.*[0.8f0])[1] RGB{Float32}(0.8*0.1,0.8*0.2,0.8*0.3)
        @test_colortype_approx_eq ([0.8f0].*cf)[1] RGB{Float32}(0.8*0.1,0.8*0.2,0.8*0.3)
        @test isfinite(cf)
        @test !isinf(cf)
        @test !isnan(cf)
        @test !isfinite(RGB(NaN, 1, 0.5))
        @test !isinf(RGB(NaN, 1, 0.5))
        @test isnan(RGB(NaN, 1, 0.5))
        @test !isfinite(RGB(1, Inf, 0.5))
        @test isinf(RGB(1, Inf, 0.5))
        @test !isnan(RGB(1, Inf, 0.5))
        @test abs(RGB(0.1,0.2,0.3)) == RGB(0.1,0.2,0.3)
        cv = RGB(0.1,0.2,0.3)
        @test ColorVectorSpace.Future.abs2(cv) == cv ⋅ cv
        @test_logs (:warn, r"is now consistent") ColorVectorSpace.Future.abs2(cv)
        @test abs2(cv) == cv ⋅ cv
        @test abs2(cv) ≈ norm(cv)^2
        @test_throws MethodError sum(abs2, RGB(0.1,0.2,0.3))
        @test norm(RGB(0.1,0.2,0.3)) ≈ sqrt(0.14)/sqrt(3)

        acu = RGB{N0f8}[cu]
        acf = RGB{Float32}[cf]
        @test typeof(acu+acf) == Vector{RGB{Float32}}
        @test typeof(acu-acf) == Vector{RGB{Float32}}
        @test typeof(acu.+acf) == Vector{RGB{Float32}}
        @test typeof(acu.-acf) == Vector{RGB{Float32}}
        @test typeof(acu.+cf) == Vector{RGB{Float32}}
        @test typeof(acu.-cf) == Vector{RGB{Float32}}
        @test typeof(2*acf) == Vector{RGB{Float32}}
        @test typeof(convert(UInt8, 2)*acu) == Vector{RGB{Float32}}
        @test typeof(acu/2) == Vector{RGB{typeof(N0f8(0.5)/2)}}
        rcu = rand(RGB{N0f8}, 3, 5)
        @test @inferred(rcu./trues(3, 5)) == rcu
        @test typeof(rcu./trues(3, 5)) == Matrix{typeof(cu/true)}

        a = RGB{N0f8}[RGB(1,0,0), RGB(1,0.8,0)]
        @test sum(a) == RGB(2.0,0.8,0)
        @test sum(typeof(a)()) == RGB(0.0,0.0,0)
        ab = [RGB(true, false, true), RGB(true, false, false)]
        @test sum(ab) === RGB(n8sum(true, true), n8sum(false, false), n8sum(true, false))
        ab = [RGB(true, false, true)]
        @test sum(ab) === RGB(n8sum(true, false), n8sum(false, false), n8sum(true, false))

        @test isapprox(a, a)
        a = RGB{Float64}(1.0, 1.0, 0.9999999999999999)
        b = RGB{Float64}(1.0, 1.0, 1.0)
        @test isapprox(a, b)
        a = RGB{Float64}(1.0, 1.0, 0.99)
        @test !(isapprox(a, b, rtol = 0.01))
        @test isapprox(a, b, rtol = 0.1)
        # issue #56
        @test RGB24(1,0,0)*N0f8(0.5) === RGB{N0f8}(0.5,0,0)
        @test RGB24(1,0,0)*0.5 === RGB(0.5,0,0)
        @test RGB24(1,0,0)/2   === RGB(0.5f0,0,0)
        @test RGB24(1,0,0)/2.0 === RGB(0.5,0,0)
        # issue #133
        @test RGB24(1, 0, 0) + RGB24(0, 0, 1) === RGB24(1, 0, 1)
        # issue #166
        @test 0.6N0f8 * RGB{N0f16}(0.2) === RGB{N0f16}(0.12, 0.12, 0.12)

        # Multiplication
        @test_throws MethodError cf*cf
        cf64 = mapc(Float64, cf)
        @test cf ⋅ cf   === (red(cf)^2 + green(cf)^2 + blue(cf)^2)/3
        @test cf ⋅ cf64 === (red(cf)*red(cf64) + green(cf)*green(cf64) + blue(cf)*blue(cf64))/3
        @test cf ⊙ cf   === RGB(red(cf)^2, green(cf)^2, blue(cf)^2)
        @test cf ⊙ cf64 === RGB(red(cf)*red(cf64), green(cf)*green(cf64), blue(cf)*blue(cf64))
        c2 = rand(RGB{Float64})
        rr = cf ⊗ c2
        matrix_rr = Matrix(rr)
        @test matrix_rr  == [red(cf)*red(c2)   red(cf)*green(c2)   red(cf)*blue(c2);
                             green(cf)*red(c2) green(cf)*green(c2) green(cf)*blue(c2);
                             blue(cf)*red(c2)  blue(cf)*green(c2)  blue(cf)*blue(c2)]
        @test RGBRGB(matrix_rr) == rr
        @test RGBRGB{N0f8}(zeros(3, 3)) == zero(RGBRGB{N0f8})
        @test_throws MethodError RGBRGB(rand(9))
        @test_throws DimensionMismatch RGBRGB(rand(9, 1))
        @test +rr === rr
        @test -rr === RGBRGB(-rr.rr, -rr.gr, -rr.br, -rr.rg, -rr.gg, -rr.bg, -rr.rb, -rr.gb, -rr.bb)
        @test rr + rr == 2*rr == rr*2
        @test rr - rr == zero(rr)
        io = IOBuffer()
        print(io, N0f8)
        Tstr = String(take!(io))
        cfn = RGB{N0f8}(0.1, 0.2, 0.3)
        show(io, cfn ⊗ cfn)
        nsuf = string(0.0N0f8)[4:end] # FixedPointNumbers issue #241
        @test String(take!(io)) == "RGBRGB{$Tstr}([0.012$nsuf 0.02$nsuf 0.031$nsuf; " *
                                                  "0.02$nsuf 0.039$nsuf 0.059$nsuf; " *
                                                  "0.031$nsuf 0.059$nsuf 0.09$nsuf])"
        show(io, "text/plain", cfn ⊗ cfn)
        spstr = Base.VERSION >= v"1.5" ? "" : " "
        @test String(take!(io)) == "RGBRGB{$Tstr}:\n 0.012  0.02   0.031\n 0.02   0.039  0.059\n 0.031  0.059  0.09$spstr"
        show(IOContext(io, :compact => false), "text/plain", cf64 ⊗ cf64)
        regex = Regex(raw"^RGBRGB{Float64}:\n" *
                      raw" 0.010\d+\s+0.020\d+\s+0.030\d+\s*\n" *
                      raw" 0.020\d+\s+0.040\d+\s+0.060\d+\s*\n" *
                      raw" 0.030\d+\s+0.060\d+\s+0.090\d+\s*$")
        @test occursin(regex, String(take!(io)))
    end

    @testset "Arithemtic with TransparentRGB" begin
        cf = RGBA{Float32}(0.1,0.2,0.3,0.4)
        @test @inferred(zero(cf)) === RGBA{Float32}(0,0,0,0)
        @test @inferred(oneunit(cf)) === RGBA{Float32}(1,1,1,1)
        @test +cf == cf
        @test -cf == RGBA(-0.1f0, -0.2f0, -0.3f0, -0.4f0)
        ccmp = RGBA{Float32}(0.2,0.4,0.6,0.8)
        @test 2*cf == ccmp
        @test cf*2 == ccmp
        @test ccmp/2 == cf
        @test 2.0f0*cf == ccmp
        @test eltype(2.0*cf) == Float64
        cu = RGBA{N0f8}(0.1,0.2,0.3,0.4)
        @test_colortype_approx_eq 2*cu RGBA(2*cu.r, 2*cu.g, 2*cu.b, 2*cu.alpha)
        @test_colortype_approx_eq 2.0f0*cu RGBA(2.0f0*cu.r, 2.0f0*cu.g, 2.0f0*cu.b, 2.0f0*cu.alpha)
        f = N0f8(0.5)
        @test (f*cu).r ≈ f*cu.r
        @test cf/2.0f0 == RGBA{Float32}(0.05,0.1,0.15,0.2)
        @test cu/2 == RGBA(cu.r/2,cu.g/2,cu.b/2,cu.alpha/2)
        @test cu/0.5f0 == RGBA(cu.r/0.5f0, cu.g/0.5f0, cu.b/0.5f0, cu.alpha/0.5f0)
        @test @inferred(cu/0.4N0f8) === RGBA{N0f8}(26/102, 51/102, 76/102, 102/102)
        @test_colortype_approx_eq @inferred(cf / 0.4N0f8) RGBA{Float32}(0.25, 0.5, 0.75, 1)
        cq0f7 = RGBA{Q0f7}(0.125, 0.25, 0.375, 0.4375)
        @test @inferred(cq0f7 / 0.5Q0f7) === RGBA{Q0f7}(0.25, 0.5, 0.75, 0.875) # issue #154
        @test cf+cf == ccmp
        @test_colortype_approx_eq (cf.*[0.8f0])[1] RGBA{Float32}(0.8*0.1,0.8*0.2,0.8*0.3,0.8*0.4)
        @test_colortype_approx_eq ([0.8f0].*cf)[1] RGBA{Float32}(0.8*0.1,0.8*0.2,0.8*0.3,0.8*0.4)
        @test isfinite(cf)
        @test !isinf(cf)
        @test !isnan(cf)
        @test isnan(RGBA(NaN, 1, 0.5, 0.8))
        @test !isinf(RGBA(NaN, 1, 0.5))
        @test isnan(RGBA(NaN, 1, 0.5))
        @test !isfinite(RGBA(1, Inf, 0.5))
        @test RGBA(1, Inf, 0.5) |> isinf
        @test !isnan(RGBA(1, Inf, 0.5))
        @test !isfinite(RGBA(0.2, 1, 0.5, NaN))
        @test !isinf(RGBA(0.2, 1, 0.5, NaN))
        @test isnan(RGBA(0.2, 1, 0.5, NaN))
        @test !isfinite(RGBA(0.2, 1, 0.5, Inf))
        @test RGBA(0.2, 1, 0.5, Inf) |> isinf
        @test !isnan(RGBA(0.2, 1, 0.5, Inf))
        @test abs(RGBA(0.1,0.2,0.3,0.2)) === RGBA(0.1,0.2,0.3,0.2)

        acu = RGBA{N0f8}[cu]
        acf = RGBA{Float32}[cf]
        @test typeof(acu+acf) == Vector{RGBA{Float32}}
        @test typeof(acu-acf) == Vector{RGBA{Float32}}
        @test typeof(acu.+acf) == Vector{RGBA{Float32}}
        @test typeof(acu.-acf) == Vector{RGBA{Float32}}
        @test typeof(acu.+cf) == Vector{RGBA{Float32}}
        @test typeof(acu.-cf) == Vector{RGBA{Float32}}
        @test typeof(2*acf) == Vector{RGBA{Float32}}
        @test typeof(convert(UInt8, 2)*acu) == Vector{RGBA{Float32}}
        @test typeof(acu/2) == Vector{RGBA{typeof(N0f8(0.5)/2)}}

        a = RGBA{N0f8}[RGBA(1,0,0,0.8), RGBA(0.7,0.8,0,0.9)]
        @test sum(a) == RGBA(n8sum(1,0.7),0.8,0,n8sum(0.8,0.9))
        @test isapprox(a, a)
        a = ARGB{Float64}(1.0, 1.0, 1.0, 0.9999999999999999)
        b = ARGB{Float64}(1.0, 1.0, 1.0, 1.0)

        @test isapprox(a, b)
        a = ARGB{Float64}(1.0, 1.0, 1.0, 0.99)
        @test !(isapprox(a, b, rtol = 0.01))
        @test isapprox(a, b, rtol = 0.1)
        # issue #56
        @test ARGB32(1,0,0,0.8)*N0f8(0.5) === ARGB{N0f8}(0.5,0,0,0.4)
        @test ARGB32(1,0,0,0.8)*0.5 === ARGB(0.5,0,0,0.4)
        @test ARGB32(1,0,0,0.8)/2   === ARGB(0.5f0,0,0,0.5f0*N0f8(0.8))
        @test ARGB32(1,0,0,0.8)/2.0 === ARGB(0.5,0,0,0.4)
        # issue #133
        @test ARGB32(1, 0, 0, 0.2) + ARGB32(0, 0, 1, 0.2) === ARGB32(1, 0, 1, 0.4)

        # issue #146
        @test @inferred(RGBA32(1,0,0,0.8)*N0f8(0.5)) === RGBA{N0f8}(0.5,0,0,0.4)
        @test @inferred(RGBA32(1,0,0,0.8)*0.5) === RGBA(0.5,0,0,0.4)
        @test @inferred(RGBA32(1,0,0,0.8)/2)   === RGBA(0.5f0,0,0,0.5f0*N0f8(0.8))
        @test @inferred(RGBA32(1,0,0,0.8)/2.0) === RGBA(0.5,0,0,0.4)
        @test @inferred(RGBA32(1, 0, 0, 0.2) + RGBA32(0, 0, 1, 0.2)) === RGBA32(1, 0, 1, 0.4)

        # Multiplication
        @test_throws MethodError cf * cf
        @test_throws MethodError cf ⋅ cf
        @test_throws MethodError cf ⊗ cf
        cf64 = mapc(Float64, cf)
        @test @inferred(cf ⊙ cf)   === RGBA{Float32}(0.1f0^2, 0.2f0^2, 0.3f0^2, 0.4f0^2)
        @test @inferred(cf ⊙ cf64) === RGBA{Float64}(0.1f0*(0.1f0*1.0), 0.2f0*(0.2f0*1.0),
                                                     0.3f0*(0.3f0*1.0), 0.4f0*(0.4f0*1.0))
    end

    @testset "Mixed-type arithmetic" begin
        # issue 155
        @test @inferred(Gray(0.2f0) + Gray24(0.2)) === Gray{Float32}(0.2 + 0.2N0f8)
        @test @inferred(RGBX(0, 0, 1) + XRGB(1, 0, 0)) === XRGB{N0f8}(1, 0, 1)
        @test @inferred(BGR(0, 0, 1) + RGB24(1, 0, 0)) === RGB{N0f8}(1, 0, 1)
        @test_throws Exception HSV(100, 0.2, 0.4) + Gray(0.2)

        @test AGray32(0.2, 0.4) + Gray24(0.2) === AGray32(0.4, 0.4N0f8+1N0f8)
        @test AGray32(0.2, 0.4) + Gray(0.2f0) === AGray{Float32}(0.2+0.2N0f8, 0.4N0f8+1)
        @test RGB(1, 0, 0)      + Gray(0.2f0) === RGB{Float32}(1.2, 0.2, 0.2)
        @test RGB(1, 0, 0)      - Gray(0.2f0) === RGB{Float32}(0.8, -0.2, -0.2)
        @test RGB24(1, 0, 0)    + Gray(0.2f0) === RGB{Float32}(1.2, 0.2, 0.2)
        @test RGB24(1, 0, 0)    - Gray(0.2f0) === RGB{Float32}(0.8, -0.2, -0.2)
        @test RGB(1.0f0, 0, 0)  + Gray24(0.2) === RGB{Float32}(1.2, 0.2, 0.2)
        @test RGB(1.0f0, 0, 0)  - Gray24(0.2) === RGB{Float32}(0.8, -0.2, -0.2)
        @test RGB24(1, 0, 0)    + Gray24(0.2) === RGB24(1N0f8+0.2N0f8, 0.2, 0.2)
        @test RGB24(0.4, 0, 0.2)   + AGray32(0.4, 1)   === ARGB32(0.8, 0.4, 0.6, 1N0f8+1N0f8)
        @test RGB24(0.4, 0.6, 0.5) - AGray32(0.4, 0.2) === ARGB32(0, 0.2, 0.1, 0.8)
        @test ARGB32(0.4, 0, 0.2, 0.5) + Gray24(0.4)   === ARGB32(0.8, 0.4, 0.6, 0.5N0f8+1N0f8)
        @test ARGB32(0.4, 0, 0.2, 0.5) + AGray32(0.4, 0.2) === ARGB32(0.8, 0.4, 0.6, 0.5N0f8+0.2N0f8)
        @test ARGB32(0.4, 0, 0.2, 0.5) + RGB(0.4f0, 0, 0) === ARGB{Float32}(0.4N0f8+0.4, 0, 0.2N0f8, 0.5N0f8+1)

        g, rgb = Gray{Float32}(0.2), RGB{Float64}(0.1, 0.2, 0.3)
        ag, argb = AGray{Float64}(0.2, 0.8), ARGB{Float32}(0.1, 0.2, 0.3, 0.4)
        @test g ⋅ rgb === rgb ⋅ g === 0.2f0*(0.1 + 0.2 + 0.3)/3
        @test_throws MethodError g ⋅ ag
        @test_throws MethodError g ⋅ argb
        @test_throws MethodError ag ⋅ rgb
        @test_throws MethodError ag ⋅ argb
        @test_throws MethodError rgb ⋅ argb
        @test g ⊙ rgb === rgb ⊙ g === RGB{Float64}(0.2f0*0.1, 0.2f0*0.2, 0.2f0*0.3)
        @test g ⊙ ag === ag ⊙ g === AGray{Float64}(0.2f0*0.2, 1.0f0*0.8)
        @test g ⊙ argb === argb ⊙ g === ARGB{Float32}(0.2f0*0.1f0, 0.2f0*0.2f0, 0.2f0*0.3f0, 1.0f0*0.4f0)
        @test ag ⊙ rgb === rgb ⊙ ag === ARGB{Float64}(0.2*0.1, 0.2*0.2, 0.2*0.3, 0.8*1.0)
        @test ag ⊙ argb === argb ⊙ ag === ARGB{Float64}(0.2*0.1f0, 0.2*0.2f0, 0.2*0.3f0, 0.8*0.4f0)
        @test rgb ⊙ argb === argb ⊙ rgb === ARGB{Float64}(0.1*0.1f0, 0.2*0.2f0, 0.3*0.3f0, 1.0*0.4f0)
        @test g ⊗ rgb === RGB(g) ⊗ rgb
        @test rgb ⊗ g === rgb ⊗ RGB(g)
        @test Matrix(g ⊗ rgb) == Matrix(rgb ⊗ g)'
        @test_throws MethodError g ⊗ ag
        @test_throws MethodError g ⊗ argb
        @test_throws MethodError ag ⊗ rgb
        @test_throws MethodError ag ⊗ argb
        @test_throws MethodError rgb ⊗ argb
    end

    @testset "Custom RGB arithmetic" begin # see also the `RGBA32` cases above
        cf = RatRGB(1//10, 2//10, 3//10)
        @test cf ⋅ cf   === (Float64(red(cf))^2 + Float64(green(cf))^2 + Float64(blue(cf))^2)/3
    end

    @testset "arithmetic with Bool" begin # issue 148
        cb = Gray{Bool}(1)
        @test @inferred(+cb) === cb
        @test @inferred(-cb) === cb # wrapped around
        @test @inferred(one(cb) * cb) === cb
        @test oneunit(cb) === Gray(true)

        @testset "vs. Bool" begin
            @test_broken @inferred(cb + true) === @inferred(true + cb) === Gray{Float32}(2)
            @test_broken @inferred(cb - true) === Gray{Float32}(0)
            @test_broken @inferred(true - cb) === Gray{Float32}(0)
            @test @inferred(cb + false) === @inferred(false + cb) === Gray{N0f8}(1) # v0.9 behavior
            @test @inferred(cb - true) === Gray{N0f8}(0) # v0.9 behavior
            @test @inferred(true - cb) === Gray{N0f8}(0) # v0.9 behavior
            @test @inferred(cb * true) === @inferred(true * cb) === Gray{Bool}(1)
            @test @inferred(cb / true) === Gray{Float32}(1)
            @test @inferred(cb / false) === Gray{Float32}(Inf32)
            @test @inferred(true / cb) === Gray{Float32}(1)
            @test @inferred(cb^true) === cb
        end
        @testset "vs. Gray{Bool}" begin
            @test @inferred(cb + Gray(true)) === @inferred(Gray(true) + cb) === Gray{Bool}(0) # wrapped around
            @test @inferred(cb - Gray(true)) === Gray{Bool}(0)
            @test @inferred(Gray(false) - cb) === Gray{Bool}(1) # wrapped around
            @test @inferred(cb * Gray(true)) === @inferred(Gray(true) * cb) === Gray{Bool}(1)
            @test @inferred(cb / Gray(true)) === Gray{Float32}(1)
            @test @inferred(cb / Gray(false)) === Gray{Float32}(Inf32)
        end
        @testset "vs. Int" begin
            @test_broken @inferred(cb + 2) === @inferred(2 + cb) === Gray{Float32}(3)
            @test_broken @inferred(cb - 2) === Gray{Float32}(-1)
            @test_broken @inferred(2 - cb) === Gray{Float32}(1)
            @test @inferred(cb + 0) === @inferred(0 + cb) === Gray{N0f8}(1) # v0.9 behavior
            @test @inferred(cb - 1) === Gray{N0f8}(0) # v0.9 behavior
            @test @inferred(2 - cb) === Gray{N0f8}(1) # v0.9 behavior
            @test @inferred(cb * 2) === @inferred(2 * cb) === Gray{Float32}(2)
            @test @inferred(cb / 2) === Gray{Float32}(0.5)
            @test @inferred(2 / cb) === Gray{Float32}(2)
            @test @inferred(cb^1) === cb
        end
        # vs. Float32 and Gray{Float32}
        @testset "vs. $(typeof(x))" for x in (0.5f0, Gray(0.5f0))
            @test @inferred(cb + x) === @inferred(x + cb) === Gray{Float32}(1.5)
            @test @inferred(cb - x) === Gray{Float32}(0.5)
            @test @inferred(x - cb) === Gray{Float32}(-0.5)
            @test @inferred(cb * x) === @inferred(x * cb) === Gray{Float32}(0.5)
            @test @inferred(cb / x) === Gray{Float32}(2)
            @test @inferred(x / cb) === Gray{Float32}(0.5)
            if x isa Real
                @test @inferred(cb^x) === Gray{Float32}(1)
            else
                @test @inferred(x^true) === Gray{Float32}(0.5)
            end
        end
        # vs. N0f8 and Gray{N0f8}
        @testset "vs. $(typeof(x))" for x in (0.6N0f8, Gray(0.6N0f8))
            @test @inferred(cb + x) === @inferred(x + cb) === Gray{N0f8}(1N0f8 + 0.6N0f8)
            @test @inferred(cb - x) === Gray{N0f8}(1N0f8 - 0.6N0f8)
            @test @inferred(x - cb) === Gray{N0f8}(0.6N0f8 - 1N0f8)
            @test @inferred(cb * x) === @inferred(x * cb) === Gray{N0f8}(0.6)
            @test_broken @inferred(cb / x) === Gray{Float32}(1 / 0.6)
            @test_broken @inferred(x / cb) === Gray{Float32}(0.6)
            @test @inferred(cb / oneunit(x)) === Gray{N0f8}(1) # v0.9 behavior
            @test @inferred(x / cb) === Gray{N0f8}(0.6) # v0.9 behavior
            if x isa Gray
                @test_broken @inferred(true / x) === Gray{Float32}(1 / 0.6)
                @test @inferred(true / Gray(1)) === Gray{N0f8}(1.0) # v0.9 behavior
                @test @inferred(x^true) === Gray{N0f8}(0.6)
            end
        end

        @testset "vs. $(typeof(c)) multiplications" for c in (Gray(true), Gray(0.5f0), Gray(0.6N0f8))
            @test @inferred(cb ⋅ c) === @inferred(c ⋅ cb) === gray(c)
            @test @inferred(cb ⊙ c) === @inferred(c ⊙ cb) === c
            @test @inferred(cb ⊗ c) === @inferred(c ⊗ cb) === c
        end
        cf = RGB{Float32}(0.1, 0.2, 0.3)
        @test @inferred(cf + cb) === @inferred(cb + cf) === RGB{Float32}(1.1, 1.2, 1.3)
        @test @inferred(cf - cb) === RGB{Float32}(-0.9, -0.8, -0.7)
        @test @inferred(cb - cf) === RGB{Float32}(0.9, 0.8, 0.7)
        cu = RGB{N0f8}(0.1, 0.2, 0.3)
        @test @inferred(cu + cb) === @inferred(cb + cu) === mapc(v -> v + 1N0f8, cu) # wrapped around
        @test @inferred(cu - cb) === mapc(v -> v - 1N0f8, cu) # wrapped around
        @test @inferred(cb - cu) === mapc(v -> 1N0f8 - v, cu)
        @testset "vs. $(typeof(c))" for c in (cf, cu)
            @test @inferred(c * true) === @inferred(true * c) === c
            if c === cu
                @test_broken @inferred(c / true) === c / 1
                @test @inferred(c / true) == c # v0.9 behavior
            else
                @test @inferred(c / true) === c / 1
            end
            @test @inferred(cb ⋅ c) === @inferred(c ⋅ cb) === Gray(1) ⋅ c
            @test @inferred(cb ⊙ c) === @inferred(c ⊙ cb) === c
            @test @inferred(cb ⊗ c) === Gray(1) ⊗ c
            @test @inferred(c ⊗ cb) === c ⊗ Gray(1)
        end
    end

    @testset "Complement" begin
        @test complement(Gray(0.2)) === Gray(0.8)
        @test complement(AGray(0.2f0, 0.7f0)) === AGray(0.8f0, 0.7f0)
        @test complement(GrayA{N0f8}(0.2, 0.7)) === GrayA{N0f8}(0.8, 0.7)
        @test complement(Gray24(0.2)) === Gray24(0.8)
        @test complement(AGray32(0.2, 0.7)) === AGray32(0.8, 0.7)

        @test complement(RGB(0, 0.3, 1)) === RGB(1, 0.7, 0)
        @test complement(ARGB(0, 0.3f0, 1, 0.7f0)) === ARGB(1, 0.7f0, 0, 0.7f0)
        @test complement(RGBA{N0f8}(0, 0.6, 1, 0.7)) === RGBA{N0f8}(1, 0.4, 0.0, 0.7)
        @test complement(RGB24(0, 0.6, 1)) === RGB24(1, 0.4, 0.0)
        @test complement(ARGB32(0, 0.6, 1, 0.7)) === ARGB32(1, 0.4, 0.0, 0.7)
    end

    @testset "dotc" begin
        @test dotc(0.2, 0.2) == 0.2^2
        @test dotc(Int8(3), Int16(6)) === 18
        @test dotc(0.2, 0.3f0) == 0.2*0.3f0
        @test dotc(N0f8(0.2), N0f8(0.3)) == Float32(N0f8(0.2))*Float32(N0f8(0.3))
        @test dotc(Gray{N0f8}(0.2), Gray24(0.3)) == Float32(N0f8(0.2))*Float32(N0f8(0.3))
        xc, yc = RGB(0.2,0.2,0.2), RGB{N0f8}(0.3,0.3,0.3)
        @test isapprox(dotc(xc, yc) , dotc(convert(Gray, xc), convert(Gray, yc)), atol=1e-6)
        @test dotc(RGB(1,0,0), RGB(0,1,1)) == 0
    end

    @testset "typemin/max" begin
        for T in (Normed{UInt8,8}, Normed{UInt8,6}, Normed{UInt16,16}, Normed{UInt16,14}, Float32, Float64)
            @test typemin(Gray{T}) === Gray{T}(typemin(T))
            @test typemax(Gray{T}) === Gray{T}(typemax(T))
            @test typemin(Gray{T}(0.5)) === Gray{T}(typemin(T))
            @test typemax(Gray{T}(0.5)) === Gray{T}(typemax(T))
            A = maximum(Gray{T}.([1 0 0; 0 1 0]); dims=1)  # see PR#44 discussion
            @test isa(A, Matrix{Gray{T}})
            @test size(A) == (1,3)
        end
    end

    @testset "Colors issue #326" begin
        A = rand(RGB{N0f8}, 2, 2)
        @test @inferred(mean(A)) == mean(map(c->mapc(FixedPointNumbers.Treduce, c), A))
    end

    @testset "Equivalence" begin
        x = 0.4
        g = Gray(x)
        c = RGB(g)
        for p in (0, 1, 2, Inf)
            @test norm(x, p) == norm(g, p) ≈ norm(c, p)
        end
        @test dot(x, x) == dot(g, g) ≈ dot(c, c)
        @test abs2(x) ≈ abs2(g) ≈ abs2(c)
        @test_throws MethodError mapreduce(x->x^2, +, c)   # this risks breaking equivalence & noniterability
    end

    @testset "varmult" begin
        cs = [RGB(0.2, 0.3, 0.4), RGB(0.5, 0.3, 0.2)]
        @test varmult(⋅, cs) ≈ 2*(0.15^2 + 0.1^2)/3    # the /3 is for the 3 color channels, i.e., equivalence
        @test varmult(⋅, cs; corrected=false) ≈ (0.15^2 + 0.1^2)/3
        @test varmult(⋅, cs; mean=RGB(0, 0, 0)) ≈ (0.2^2+0.3^2+0.4^2 + 0.5^2+0.3^2+0.2^2)/3
        @test varmult(⊙, cs) ≈ 2*RGB(0.15^2, 0, 0.1^2)
        @test Matrix(varmult(⊗, cs)) ≈ 2*[0.15^2 0 -0.1*0.15; 0 0 0; -0.1*0.15 0 0.1^2]
        @test stdmult(⋅, cs) ≈ sqrt(2*(0.15^2 + 0.1^2)/3)    # the /3 is for the 3 color channels, i.e., equivalence
        @test stdmult(⋅, cs; corrected=false) ≈ sqrt((0.15^2 + 0.1^2)/3)
        @test stdmult(⋅, cs; mean=RGB(0, 0, 0)) ≈ sqrt((0.2^2+0.3^2+0.4^2 + 0.5^2+0.3^2+0.2^2)/3)
        @test stdmult(⊙, cs) ≈ RGB(sqrt(2*0.15^2), 0, sqrt(2*0.1^2))
        @test_throws DomainError stdmult(⊗, cs)

        cs = [RGB(0.1, 0.2,  0.3)  RGB(0.3, 0.5, 0.3);
              RGB(0.2, 0.21, 0.33) RGB(0.4, 0.51, 0.33);
              RGB(0.3, 0.22, 0.36) RGB(0.5, 0.52, 0.36)]
        v1 = RGB(0.1^2, 0.15^2, 0)
        s2v1 = mapc(sqrt, 2*v1)
        @test varmult(⊙, cs, dims=2) ≈ 2*[v1, v1, v1]
        @test stdmult(⊙, cs, dims=2) ≈ [s2v1, s2v1, s2v1]
        v2 = RGB(0.1^2, 0.01^2, 0.03^2)
        sv2 = mapc(sqrt, v2)
        @test varmult(⊙, cs, dims=1) ≈ [v2 v2]
        @test stdmult(⊙, cs, dims=1) ≈ [sv2 sv2]
        @test var(cs) == varmult(⋅, cs)
    end

    @testset "copy" begin
        g = Gray{N0f8}(0.2)
        @test copy(g) === g
        c = RGB(0.1, 0.2, 0.3)
        @test copy(c) === c
    end

    @testset "ranges" begin
        T = Gray{N0f8}
        for r in (zero(T):eps(T):oneunit(T),
                  zero(T):eps(N0f8):oneunit(T),
                  zero(N0f8):eps(T):oneunit(N0f8),
                  StepRangeLen(zero(T), eps(T), 256),
                  StepRangeLen(zero(T), eps(N0f8), 256),
                  LinRange(zero(T), oneunit(T), 256))
            @test length(r) == 256
            @test step(r) == eps(T)
        end

        r = LinRange(RGB(1, 0, 0), RGB(0, 0, 1), 256)
        @test length(r) == 256
        @test step(r) == RGB(-1.0f0/255, 0, eps(N0f8))
        @test r[2] == RGB(254.0f0/255, 0, 1.0f0/255)
    end

end
