using Statistics: mean
using StableRNGs
using PlotUtils
using Random
using Dates
using Test

rng = StableRNG(42)
warn_ticks = :warn, "No strict ticks found"

# ----------------------
# colors

const C = RGBA{Float64}
const C0 = RGBA{PlotUtils.Colors.N0f8}

@testset "colors" begin
    @test plot_color(nothing) == C(0, 0, 0, 0)
    @test plot_color(false) == C(0, 0, 0, 0)
    @test_throws ErrorException plot_color(true)

    @test plot_color(:red) == parse(C, :red)
    @test plot_color("red") == parse(C, "red")
    @test_throws ArgumentError plot_color("notacolor")

    @test plot_color(colorant"red") == C(1, 0, 0, 1)

    grad = cgrad()
    @test typeof(grad) == PlotUtils.ContinuousColorGradient
    @test plot_color(grad) ≡ grad

    # JuliaPlots/Plots.jl/issues/4270
    @test get(grad, BigFloat(1), (-0.0003354626279, 1.0)) isa RGBA{BigFloat}

    grad = cgrad([:red, "blue"])
    @test color_list(grad) == C[colorant"red", colorant"blue"]
    @test grad.values == collect(range(0, stop = 1, length = 2))

    grad = cgrad([:red, "blue"], alpha = 0.5)
    @test C0.(color_list(grad)) == C0[C(1, 0, 0, 0.5), C(0, 0, 1, 0.5)]
    @test grad.values == collect(range(0, stop = 1, length = 2))

    grad = cgrad([:red, :blue], [0, 0.1, 1])
    @test length(color_list(grad)) == 3
    @test grad.values == [0, 0.1, 1]

    Random.seed!(rng, 42)
    cs = plot_color(rand(rng, 10))
    @test typeof(cs) == Vector{C}
    @test length(cs) == 10

    cs = plot_color(rand(rng, 4, 4))
    @test typeof(cs) == Matrix{C}
    @test length(cs) == 16
    @test size(cs) == (4, 4)

    cs = plot_color(rand(rng, 10), 0.5)
    @test typeof(cs) == Vector{C}
    @test length(cs) == 10
    for c ∈ cs
        @test alpha(c) == 0.5
    end

    cs = plot_color(rand(rng, 4, 4), 0.5)
    @test typeof(cs) == Matrix{C}
    @test length(cs) == 16
    @test size(cs) == (4, 4)
    for c ∈ cs
        @test alpha(c) == 0.5
    end
end

# ----------------------
# gradients

@testset "gradients" begin
    grad = cgrad(:inferno)
    @test length(grad) == 256
    @test RGB(grad.colors[1]) ≈ RGB(0.001462, 0.000466, 0.013866)
    @test RGB(grad.colors[end]) ≈ RGB(0.988362, 0.998364, 0.644924)
end

@testset "sampling" begin
    # github.com/MakieOrg/Makie.jl/issues/2635
    cmap = cgrad([:black, :white, :orange], [0, 0.2, 1])
    # sample outside the given values
    @test RGB(get(cmap, 0.15)) ≈ RGB(0.75, 0.75, 0.75)
    @test RGB(get(cmap, 0.5)) ≈ RGB(1.0, 0.86764705, 0.625)
    @test RGB(get(cmap, 0.8)) ≈ RGB(1.0, 0.73529411, 0.25)
end

@testset "reverse" begin
    cmap = reverse(cgrad([:black, :white, :orange], [0, 0.2, 1]))
    # sample outside the given values
    @test RGB(get(cmap, 0.15)) ≈ RGB(1.0, 0.71323529, 0.1875)
    @test RGB(get(cmap, 0.5)) ≈ RGB(1.0, 0.86764705, 0.625)
    @test RGB(get(cmap, 0.75)) ≈ RGB(1.0, 0.97794117, 0.9375)
end

# ----------------------
# ticks

# Copied from Plots.is_uniformly_spaced to avoid dependency on recent version
# on Plots which is not used on Travis.
function is_uniformly_spaced(v; tol = 1e-6)
    dv = diff(v)
    maximum(dv) - minimum(dv) < tol * mean(abs.(dv))
end

function test_ticks(x, y, ticks)
    @test issorted(ticks)
    @test all(x .≤ ticks .≤ y)
    if x < y
        @test length(ticks) ≥ 2
        @test is_uniformly_spaced(ticks)
    end
end

@testset "ticks" begin
    @test optimize_ticks(-1, 2) == ([-1.0, 0.0, 1.0, 2.0], -1.0, 2.0)

    # check if ticks still generate if max - min << abs(min) (i.e. for Float64 ranges)
    @test optimize_ticks(1e11 - 1, 1e11 + 2) == (1e11 .+ (-1:2), 1e11 - 1.0, 1e11 + 2.0)

    @testset "dates" begin
        dt1, dt2 = Dates.value(DateTime(2000)), Dates.value(DateTime(2100))
        @test optimize_datetime_ticks(dt1, dt2) == (
            [63113990400000, 63902908800000, 64691827200000, 65480745600000],
            ["2001-01-01", "2026-01-01", "2051-01-01", "2076-01-01"],
        )
    end

    @testset "small range" begin
        @testset "small range $x, $(i)ϵ" for x ∈ exp10.(-12:12), i ∈ -5:5
            y = x + i * eps(x)
            x, y = minmax(x, y)
            ticks, = optimize_ticks(x, y)
            @test issorted(ticks)
            @test all(x .≤ ticks .≤ y)
            # Fails:
            # @test allunique(ticks)
        end
    end

    @testset "fixed ranges" begin
        @testset "fixed range $x..$y" for (x, y) ∈ [(2, 14), (14, 25), (16, 36), (57, 69)]
            test_ticks(+x, +y, optimize_ticks(+x, +y)[1])
            test_ticks(-y, -x, optimize_ticks(-y, -x)[1])
        end
    end

    @testset "random ranges" begin
        r = [minmax(rand(rng, -100:100, 2)...) .* 10.0^i for _ ∈ 1:10, i ∈ -5:5]
        @testset "random range $x..$y" for (x, y) ∈ r
            test_ticks(x, y, optimize_ticks(x, y)[1])
        end
    end

    @testset "digits" begin
        @testset "digits $((10^n) - 1)*10^$i" for n ∈ 1:9, i ∈ -9:9
            y0 = 10^n
            x0 = y0 - 1
            x, y = (x0, y0) .* 10.0^i
            ticks, = optimize_ticks(x, y)
            test_ticks(x, y, ticks)
        end
    end

    @testset "types" begin
        for T ∈ (Int32, Int64, Float16, Float32, Float64)
            x, y = T(1), T(10)
            ticks, = optimize_ticks(x, y)
            @test eltype(ticks) <: AbstractFloat
            @test eltype(ticks) == (T <: AbstractFloat ? T : float(T))
            test_ticks(x, y, ticks)
        end
    end

    @testset "issues" begin
        @testset "PlotUtils.jl/issues/86" begin
            let x = -1.0, y = 13.0
                test_ticks(x, y, optimize_ticks(x, y, k_min = 4, k_max = 8)[1])
            end
        end

        @testset "Plots.jl/issues/3859" begin
            x, y = extrema([-1.7055509600077687e307, -1.3055509600077687e307, -1.e300])
            test_ticks(x, y, optimize_ticks(x, y, k_min = 4, k_max = 8)[1])
        end

        @testset "PlotUtils.jl/issues/114" begin
            let x = -1.2eps(), y = 1.2eps()
                test_ticks(x, y, optimize_ticks(x, y)[1])
            end
        end

        @testset "PlotUtils.jl/issues/116" begin
            let x = 4.5, y = 5.5
                test_ticks(x, y, optimize_ticks(x, y, scale = :log10)[1])
            end
            let x = 2.5, y = 3.5
                test_ticks(x, y, optimize_ticks(x, y, scale = :log2)[1])
            end
            let x = 0.5, y = 1.5
                test_ticks(x, y, optimize_ticks(x, y, scale = :ln)[1])
            end
        end

        @testset "PlotUtils.jl/issues/129" begin
            # invalid float input
            let x = NaN, y = 1.0
                ticks, = @test_logs warn_ticks optimize_ticks(x, y)
                @test isnan(ticks[1])
                @test ticks[2] ≡ y
                ticks, = @test_logs warn_ticks optimize_ticks(x, y, k_min = 5)
                @test isnan(ticks[1])
                @test ticks[2] ≡ y
            end
            let x = 0.0f0, y = Inf32
                ticks, = @test_logs warn_ticks optimize_ticks(x, y)
                @test ticks[1] ≡ x
                @test isinf(ticks[2])
                ticks, = @test_logs warn_ticks optimize_ticks(x, y, k_min = 5)
                @test ticks[1] ≡ x
                @test isinf(ticks[2])
            end
        end
    end

    @testset "PlotUtils.jl/issues/155" begin
        for n ∈ 1:10
            @test length(palette(:tab10, n)) == n
        end
    end

    @testset "PlotUtils.jl/issues/156" begin
        for c ∈ cgrad([:black, RGBA{Float64}(1, 1, 1, 1)], 5; categorical = true)
            @test 0 ≤ alpha(c) ≤ 1
        end
    end
end

@testset "adapted grid" begin
    let f = sin, int = (0, π)
        xs, fs = adapted_grid(f, int)
        for i ∈ 1:(length(xs) - 1)
            for λ ∈ 0:0.1:1
                # test that `f` is well approximated by a line
                # in the interval `(xs[i], xs[i+1])`
                x = λ * xs[i] + (1 - λ) * xs[i + 1]
                y = λ * fs[i] + (1 - λ) * fs[i + 1]
                @test y ≈ f(x) atol = 1e-2
            end
        end
    end

    let f = sin, int = (2, 2)
        xs, fs = adapted_grid(f, int)
        @test xs == [2]
        @test fs == [f(2)]
    end

    let f = sin, int = (2, 1)
        @test_throws ArgumentError adapted_grid(f, int)
    end

    p(x) = (x < 0.0 || x > 1.0 ? 0.0 : 1.0 + x)
    let f = x -> p(2(x - 3)), int = (-5, 5)  # JuliaPlots/Plots.jl/issues/4106
        xs, fs = adapted_grid(f, int)
        @test count(fs .> 1.5) > 10
        @test fs |> extrema |> collect |> diff |> first > 1.9
    end

    let f = sinc, int = (-40, 40)  # JuliaPlots/Plots.jl/issues/3894
        xs, fs = adapted_grid(sinc, int)
        roots = vcat(int[1]:-1, 1:int[2])
        count_per_extrema = map(1:(length(roots) - 1)) do idx
            left, right = roots[idx:(idx + 1)]
            return count(x -> left < x < right, xs)
        end
        # check that we have at least 5 points for each extrema
        @test all(count_per_extrema .>= 5)
    end
end

@testset "zscale" begin
    Random.seed!(rng, 42)
    bkg = 30 .* randn(rng, 4_096) .+ 1_000
    data = bkg .+ 100 .* randn(rng, 4_096) .+ 2500
    defects = rand(rng, CartesianIndices(bkg), 500)
    data[defects] .= rand(rng, [0, 1e7], 500)
    cmin, cmax = zscale(data)
    # values calculated using IRAF
    @test cmin ≈ 2_784.824 atol = 1e-3
    @test cmax ≈ 4_211.375 atol = 1e-3
    @test cmin > minimum(data)
    @test cmax < maximum(data)

    data = vcat(1:100)
    cmin, cmax = zscale(data)
    @test cmin == 1
    @test cmax == 100

    # make sure output is finite
    data = vcat(0:999, NaN)
    cmin, cmax = zscale(data)
    @test cmin == 0
    @test cmax == 999
end

@testset "allocations" begin  # see PlotUtils.jl/pull/136
    stats = @timed optimize_ticks(0.1123, 100.132)
    @test stats.bytes < 1_000  # ~ 736 (on 1.9)
    @test stats.time < 1e-3  # ~ 0.22ms (on 1.9)
end

if Sys.islinux() && VERSION ≥ v"1.9.0" && isempty(VERSION.prerelease)  # avoid running on `nightly`
    @testset "downstream" begin
        include("downstream.jl")
    end
    @testset "adaptive" begin  # NOTE: must be ran after downstream test (for Plots)
        withenv("GKSwstype" => "nul") do
            include("adaptive_test_functions.jl")
        end
        @test true
    end
end
