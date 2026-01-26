module PlotUtils

using PrecompileTools
using ColorSchemes
using Reexport
using Printf
using Dates

@reexport using Colors
import Base: getindex
import StableRNGs: StableRNG

export ColorGradient,
    ColorPalette,
    cgrad,
    palette,
    color_list,
    cvec,
    rgb_string,
    rgba_string,
    invisible,
    get_color_palette,
    isdark,
    plot_color,
    adapted_grid,
    default_cgrad,
    zscale

include("color_utils.jl")
include("colors.jl")
include("colorschemes.jl")
include("adapted_grid.jl")
include("intervals.jl")

export optimize_ticks, optimize_datetime_ticks

include("ticks.jl")

const _default_colorscheme = generate_colorscheme()

if VERSION ≥ v"1.8.0"
    @compile_workload begin
        for T ∈ (Int, Float64)
            optimize_ticks(-one(T), one(T))
            optimize_ticks(-one(T), one(T); k_min = 2, k_max = 10)
            adapted_grid(sin, (-one(T), one(T)))
            zscale(one(T):10)
            cgrad([:red, :blue], T[0, 1])
        end
        cgrad()
        cgrad([:red, :blue])
        palette(:viridis)
        optimize_datetime_ticks(Dates.value(DateTime(2_000)), Dates.value(DateTime(2_100)))
    end
end

end
