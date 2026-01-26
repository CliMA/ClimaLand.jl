using AdaptivePredicates
using Test
using Supposition
import ExactPredicates: ExactPredicates
using Aqua

@testset "Aqua" begin
    Aqua.test_all(AdaptivePredicates)
end

const AP = AdaptivePredicates
cd("original") do
    include("compile.jl")
end
include("CPredicates.jl")
import .CPredicates as C

include("utils.jl")

@testset "Check exactinit()" begin
    @test AP.InitConstants{Float64}(C.C64_consts...) == AP.IC64
    @test AP.InitConstants{Float32}(C.C32_consts...) == AP.IC32
end

const MACROS = (
    MacroMethod(:Absolute, 1),
    MacroMethod(:Fast_Two_Sum_Tail, 3),
    MacroMethod(:Fast_Two_Sum, 2),
    MacroMethod(:Fast_Two_Diff_Tail, 3),
    MacroMethod(:Fast_Two_Diff, 2),
    MacroMethod(:Two_Sum_Tail, 3),
    MacroMethod(:Two_Sum, 2),
    MacroMethod(:Two_Diff_Tail, 3),
    MacroMethod(:Two_Diff, 2),
    MacroMethod(:Split, 1),
    MacroMethod(:Two_Product_Tail, 3),
    MacroMethod(:Two_Product, 2),
    MacroMethod(:Two_Product_Presplit, 4),
    MacroMethod(:Two_Product_2Presplit, 6),
    MacroMethod(:Square_Tail, 2),
    MacroMethod(:Square, 1),
    MacroMethod(:Two_One_Sum, 3),
    MacroMethod(:Two_One_Diff, 3),
    MacroMethod(:Two_Two_Sum, 4),
    MacroMethod(:Two_Two_Diff, 4),
    MacroMethod(:Four_One_Sum, 5),
    MacroMethod(:Four_Two_Sum, 6),
    MacroMethod(:Four_Four_Sum, 8),
    MacroMethod(:Eight_One_Sum, 9),
    MacroMethod(:Eight_Two_Sum, 10),
    MacroMethod(:Eight_Four_Sum, 12),
    MacroMethod(:Two_One_Product, 3),
    MacroMethod(:Four_One_Product, 5),
    MacroMethod(:Two_Two_Product, 4),
    MacroMethod(:Two_Square, 2)
)
_sum_size(args) = length(args[end])
_grow_size(args) = length(args[end]) + length(args[end-1])
_scale_size(args) = 2length(args[end-1])
const ARITHMETIC = (
    ArithmeticMethod(:grow_expansion, 3, (-1, 1, _grow_size)),
    ArithmeticMethod(:grow_expansion_zeroelim, 3, (-1, 1, _grow_size)),
    ArithmeticMethod(:expansion_sum, 3, (-1, -1, _grow_size)),
    ArithmeticMethod(:expansion_sum_zeroelim1, 3, (-1, -1, _grow_size)),
    ArithmeticMethod(:expansion_sum_zeroelim2, 3, (-1, -1, _grow_size)),
    ArithmeticMethod(:fast_expansion_sum, 3, (-1, -1, _grow_size)),
    ArithmeticMethod(:fast_expansion_sum_zeroelim, 3, (-1, -1, _grow_size)),
    ArithmeticMethod(:linear_expansion_sum, 3, (-1, -1, _grow_size)),
    ArithmeticMethod(:linear_expansion_sum_zeroelim, 3, (-1, -1, _grow_size)),
    ArithmeticMethod(:scale_expansion, 3, (-1, 1, _scale_size)),
    ArithmeticMethod(:scale_expansion_zeroelim, 3, (-1, 1, _scale_size)),
    ArithmeticMethod(:compress, 2, (-1, _sum_size)),
    ArithmeticMethod(:estimate, 1, (-1,))
)
const PREDICATES = (
    PredicateMethod(:orient2, 3, 2),
    PredicateMethod(:orient3, 4, 3),
    PredicateMethod(:incircle, 4, 2),
    PredicateMethod(:insphere, 5, 3)
)

const MACRO_FAILURE = FailedTest[]
const ARITHMETIC_FAILURE = FailedTest[]
const PREDICATE_FAILURE = FailedTest[]

@testset "Exact Equality Tests" begin
    @testset "Macros" begin
        foreach(MACROS) do f
            @repeat test_f(f; failures=MACRO_FAILURE)
        end
    end

    @testset "Arithmetic" begin
        foreach(ARITHMETIC) do f
            @repeat test_f(f; failures=ARITHMETIC_FAILURE)
        end
    end

    @testset "Predicates" begin
        foreach(PREDICATES) do f
            @repeat test_f(f; failures=PREDICATE_FAILURE)
        end
    end
end

check_range(x::Float64) = iszero(x) || exponent(abs(x)) ∈ -142:201 # Ranges are explained in the README
check_range(x::Float32) = iszero(x) || exponent(abs(x)) ∈ -24:24 # Don't know the exact range for Float32, it might be [-3, 24] but ExactPredicates doesn't use Float32 anyway for me to confirm.
@testset "Supposition tests" begin
    for (T, Tn) in ((Float64, 64), (Float32, 32))
        fgen = Data.Floats{T}(infs=false, nans=false)
        cgen = @composed _complex(a=fgen, b=fgen) = a + b * im
        r2gen = @composed _tuple(a=fgen, b=fgen) = (a, b)
        r3gen = @composed _tuple(a=fgen, b=fgen, c=fgen) = (a, b, c)
        igen = Data.Integers(-100, 100)
        i2gen = @composed _ituple(a=igen, b=igen) = (a, b)
        i3gen = @composed _ituple(a=igen, b=igen, c=igen) = (a, b, c)
        # Against C
        ## standard
        @check function _orient2(p=cgen, q=cgen, r=cgen)
            assume!(all(check_range, (p.re, p.im, q.re, q.im, r.re, r.im)))
            ap = orient2(p, q, r)
            c = C.orient2d((p.re, p.im), (q.re, q.im), (r.re, r.im))
            event!("AdaptiveOrient2$(Tn)", ap)
            event!("COrient2$(Tn)", c)
            ap == c
        end

        @check function _orient3(p=r3gen, q=r3gen, r=r3gen, s=r3gen)
            assume!(all(check_range, (p..., q..., r..., s...)))
            ap = orient3(p, q, r, s)
            c = C.orient3d(p, q, r, s)
            event!("AdaptiveOrient3$(Tn)", ap)
            event!("COrient3$(Tn)", c)
            ap == c
        end

        @check function _incircle(p=r2gen, q=r2gen, r=r2gen, s=r2gen)
            assume!(all(check_range, (p..., q..., r..., s...)))
            ap = incircle(p, q, r, s)
            c = C.incircle(p, q, r, s)
            event!("AdaptiveIncircle$(Tn)", ap)
            event!("CIncircle$(Tn)", c)
            ap == c
        end

        @check function _insphere(p=r3gen, q=r3gen, r=r3gen, s=r3gen, t=r3gen)
            assume!(all(check_range, (p..., q..., r..., s..., t...)))
            ap = insphere(p, q, r, s, t)
            c = C.insphere(p, q, r, s, t)
            event!("AdaptiveInsphere$(Tn)", ap)
            event!("CInsphere$(Tn)", c)
            ap == c
        end

        # Against ExactPredicates
        if T ≠ Float32
            @check function _orient2p(p=cgen, q=cgen, r=cgen)
                assume!(all(check_range, (p.re, p.im, q.re, q.im, r.re, r.im)))
                ap = orient2p(p, q, r)
                c = ExactPredicates.orient((p.re, p.im), (q.re, q.im), (r.re, r.im))
                event!("AdaptiveOrient2$(Tn)", ap)
                event!("ExactPredicatesOrient2$(Tn)", c)
                ap == c
            end

            @check function _orient3p(p=r3gen, q=r3gen, r=r3gen, s=r3gen)
                assume!(all(check_range, (p..., q..., r..., s...)))
                ap = orient3p(p, q, r, s)
                c = ExactPredicates.orient(p, q, r, s)
                event!("AdaptiveOrient3$(Tn)", ap)
                event!("ExactPredicatesOrient3$(Tn)", c)
                ap == c
            end

            @check function _incirclep(p=r2gen, q=r2gen, r=r2gen, s=r2gen)
                assume!(all(check_range, (p..., q..., r..., s...)))
                ap = incirclep(p, q, r, s)
                c = ExactPredicates.incircle(p, q, r, s)
                event!("AdaptiveIncircle$(Tn)", ap)
                event!("ExactPredicatesIncircle$(Tn)", c)
                ap == c
            end

            @check function _inspherep(p=r3gen, q=r3gen, r=r3gen, s=r3gen, t=r3gen)
                assume!(all(check_range, (p..., q..., r..., s..., t...)))
                ap = inspherep(p, q, r, s, t)
                c = ExactPredicates.insphere(p, q, r, s, t)
                event!("AdaptiveInsphere$(Tn)", ap)
                event!("ExactPredicatesInsphere$(Tn)", c)
                ap == c
            end

            # integers 

            @check function _intorient2p(p=i2gen, q=i2gen, r=i2gen)
                assume!(all(check_range, (T.(p)..., T.(q)..., T.(r)...)))
                ap = orient2p(T.(p), T.(q), T.(r))
                c = ExactPredicates.orient(T.(p), T.(q), T.(r))
                event!("AdaptiveOrient2Int$(Tn)", ap)
                event!("ExactPredicatesOrient2Int$(Tn)", c)
                ap == c
            end

            @check function _intorient3p(p=i3gen, q=i3gen, r=i3gen, s=i3gen)
                assume!(all(check_range, (T.(p)..., T.(q)..., T.(r)..., T.(s)...)))
                ap = orient3p(T.(p), T.(q), T.(r), T.(s))
                c = ExactPredicates.orient(T.(p), T.(q), T.(r), T.(s))
                event!("AdaptiveOrient3Int$(Tn)", ap)
                event!("ExactPredicatesOrient3Int$(Tn)", c)
                ap == c
            end

            @check function _intincirclep(p=i2gen, q=i2gen, r=i2gen, s=i2gen)
                assume!(all(check_range, (T.(p)..., T.(q)..., T.(r)..., T.(s)...)))
                ap = incirclep(T.(p), T.(q), T.(r), T.(s))
                c = ExactPredicates.incircle(T.(p), T.(q), T.(r), T.(s))
                event!("AdaptiveIncircleInt$(Tn)", ap)
                event!("ExactPredicatesIncircleInt$(Tn)", c)
                ap == c
            end

            @check function _intinspherep(p=i3gen, q=i3gen, r=i3gen, s=i3gen, t=i3gen)
                assume!(all(check_range, (T.(p)..., T.(q)..., T.(r)..., T.(s)..., T.(t)...)))
                ap = inspherep(T.(p), T.(q), T.(r), T.(s), T.(t))
                c = ExactPredicates.insphere(T.(p), T.(q), T.(r), T.(s), T.(t))
                event!("AdaptiveInsphereInt$(Tn)", ap)
                event!("ExactPredicatesInsphereInt$(Tn)", c)
                ap == c
            end
        end
    end
end

@testset "Caches" begin
    function check_non_overlapping(args...)
        ranges = first.(parentindices.(args))
        flag = !issorted(ranges)
        flag && return false
        for i in 1:(length(ranges)-1)
            rᵢ, rᵢ₊₁ = ranges[i], ranges[i+1]
            !isdisjoint(rᵢ, rᵢ₊₁) && return false
        end
        return true
    end
    function check_gaps(args...)
        ranges = first.(parentindices.(args))
        firstdiffs = diff(collect(first.(ranges)))
        boundaries = firstdiffs .% 64
        return all(iszero, boundaries)
    end
    v = zeros(100)
    v1, v2, v3 = view(v, 1:3), view(v, 2:4), view(v, 7:100)
    @test !check_non_overlapping(v1, v2, v3)
    @test !check_gaps(v1, v2, v3)

    function test_cache(f, lengths)
        for T in (Float64, Float32)
            args1 = f(T)
            args2 = f(T, nothing)
            for args in (args1, args2)
                @test check_non_overlapping(args...)
                @test all(lengths .== length.(args))
                @test all(Base.Fix1(===, parent(args[1])), parent.(args))
                @test all(x -> eltype(x) === T, args)
                @test check_gaps(args...)
                @test f(T, args) === args
            end
        end
    end
    test_cache(AP.incircleexact_cache, (48, 48, 96, 96, 96, 96, 192, 192, 384))
    test_cache(AP.incircleslow_cache, (64, 64, 64, 64, 64, 64, 128, 128, 192, 192, 384, 384, 384, 768, 1152))
    test_cache(AP.incircleadapt_cache, (48, 64, 1152, 1152))
    test_cache(AP.insphereslow_cache, (64, 64, 64, 128, 192, 384, 384, 384, 384, 384, 384, 768, 768, 768, 768, 768, 768, 768, 768, 768, 1536, 1536, 1536, 2304, 2304, 2304, 4608, 6912, 6912, 6912, 6912, 13824, 13824, 27648))
    test_cache(AP.insphereexact_cache, (48, 48, 96, 96, 96, 96, 96, 192, 288, 288, 288, 288, 384, 384, 384, 576, 576, 768, 1152, 1152, 1152, 1152, 1152, 2304, 2304, 3456, 5760))
    test_cache(AP.orient3exact_cache, (48, 48, 96))
    test_cache(AP.orient3slow_cache, (64, 64, 64, 128, 192))
    test_cache(AP.orient3adapt_cache, (192, 192))
end


@testset "@check_length" begin
    function test_macro(_val)
        args = (2.0, -2.0, 3.0)
        AP.@check_length test_return args val, len = test_call(_val)
        return val, _val, len
    end
    function test_call(x)
        return 15.0, x > 2 ? 17 : 50
    end
    function test_return(x, y, z)
        return x + y + z
    end
    @test test_macro(5.0) == (15.0, 5.0, 17)
    @test test_macro(1.0) == 3.0
    expr = @macroexpand function test_macro(_val)
        args = (2.0, -2.0, 3.0)
        AP.@check_length test_return args val, len = test_call(_val)
        return val, _val, len
    end
    Base.remove_linenums!(expr)
    @test sprint(show, expr) == ":(function test_macro(_val)\n      args = (2.0, -2.0, 3.0)\n      begin\n          (val, len) = test_call(_val)\n          len < 32 || return test_return(args...)\n      end\n      return (val, _val, len)\n  end)"
end
