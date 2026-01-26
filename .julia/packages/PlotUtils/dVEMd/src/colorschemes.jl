# An abstract color list is either a Palette or a gradient
abstract type AbstractColorList end

get_colorscheme(acl::AbstractColorList) = acl.colors
color_list(acl::AbstractColorList) = color_list(acl.colors)
plot_color(cl::AbstractColorList) = cl

# Interface
Base.show(io::IO, m::MIME"image/svg+xml", acl::AbstractColorList) =
    show(io, m, color_list(acl))
Base.length(acl::AbstractColorList) = length(color_list(acl))
Base.getindex(acl::AbstractColorList, x) = getindex(color_list(acl), x)
Base.size(acl::AbstractColorList) = size(color_list(acl))
Base.IndexStyle(::Type{<:AbstractColorList}) = IndexLinear()
Base.iterate(acl::AbstractColorList) = iterate(color_list(acl))
Base.iterate(acl::AbstractColorList, s) = iterate(color_list(acl), s)
Base.lastindex(acl::AbstractColorList) = lastindex(color_list(acl))
Base.firstindex(acl::AbstractColorList) = firstindex(color_list(acl))
Base.get(acl::AbstractColorList, args...) = get(get_colorscheme(acl), args...)

## ColorGradient

abstract type ColorGradient <: AbstractColorList end

Base.getindex(cg::ColorGradient, x::Union{AbstractFloat,AbstractVector{<:AbstractFloat}}) =
    get(cg, x)
function Base.get(cg::ColorGradient, v::AbstractArray, rangescale = (0.0, 1.0))
    rangescale ≡ :extrema && (rangescale = extrema(v))
    map(x -> get(cg, x, rangescale), v)
end

Base.:(==)(cg1::ColorGradient, cg2::ColorGradient) =
    color_list(cg1) == color_list(cg2) && cg1.values == cg2.values
Base.hash(cg::ColorGradient) = hash(color_list(cg)) ⊻ hash(cg.values)

## Continuous Color Gradient

struct ContinuousColorGradient <: ColorGradient
    colors::ColorScheme
    values::Vector{Float64}

    function ContinuousColorGradient(colors, values = get_range(colors))
        c, v = prepare_continuous_cgrad_colors(colors, values)
        new(c, v)
    end
end

plot_color(cg::ContinuousColorGradient, α::Number) =
    ContinuousColorGradient(plot_color(cg.colors, α), cg.values)

Base.show(io::IO, m::MIME"image/svg+xml", cg::ContinuousColorGradient) =
    show(io, m, cg[get_range(100)])

function sample_color(cg::ContinuousColorGradient, x::AbstractFloat)
    c, v = cg.colors, cg.values
    if (index = findfirst(==(x), v)) ≡ nothing
        nm1 = length(v) - 1
        i = min(nm1, findlast(<(x), v))
        r = (x - v[i]) / (v[i + 1] - v[i])
        index = (i + r - 1) / nm1
    end
    c[index]
end

function Base.get(cg::ContinuousColorGradient, x::AbstractFloat, rangescale = (0.0, 1.0))
    isfinite(x) || return invisible()
    rangescale = get_rangescale(rangescale)
    allunique(rangescale) || return first(cg.colors)
    x = clamp(x, rangescale...)
    if rangescale != (0.0, 1.0)
        x = ColorSchemes.remap(x, rangescale..., 0, 1)
    end
    sample_color(cg, x)  # specialize for x (boxing issues ?)
end

Base.reverse(cg::ContinuousColorGradient) =
    ContinuousColorGradient(reverse(cg.colors), reverse(1 .- cg.values))

# required for GR's nonuniformcellarray
function ColorSchemes.getinverse(cg::ContinuousColorGradient, c)
    alpha(c) == 0 && return NaN
    z = getinverse(to_rgb(get_colorscheme(cg)), to_rgb(c))
    cr = get_range(cg.colors)
    if (index = findfirst(==(z), cr)) ≢ nothing
        cg.values[index]
    else
        i = min(length(cr) - 1, findlast(<(z), cr))
        ColorSchemes.remap(z, cr[i], cr[i + 1], cg.values[i], cg.values[i + 1])
    end
end

function prepare_continuous_cgrad_colors(c, v)
    v = sort(unique(clamp.(v, 0, 1)))
    0 ∈ v || pushfirst!(v, 0)
    1 ∈ v || push!(v, 1)
    nv = length(v)
    nc = length(c)
    c, v = if nc != nv
        value_range = get_range(nv)
        color_range = get_range(nc)
        values = [0.0]
        for i ∈ 2:nv
            inds = findall(x -> value_range[i - 1] < x < value_range[i], color_range)
            isempty(inds) || append!(
                values,
                ColorSchemes.remap(
                    color_range[inds],
                    value_range[i - 1],
                    value_range[i],
                    v[i - 1],
                    v[i],
                ),
            )
            push!(values, v[i])
        end
        colors = c[sort(unique([value_range; color_range]))]
        colors, values
    else
        color_list(c), v
    end
    ColorScheme(plot_color(c)), v
end

## Categorical Color Gradient

struct CategoricalColorGradient <: ColorGradient
    colors::ColorScheme
    values::Vector{Float64}

    function CategoricalColorGradient(colors, values = get_range(colors))
        c, v = prepare_categorical_cgrad_colors(colors, values)
        new(c, v)
    end
end

plot_color(cg::CategoricalColorGradient, α::Number) =
    CategoricalColorGradient(plot_color(cg.colors, α), cg.values)

Base.show(io::IO, m::MIME"image/svg+xml", cg::CategoricalColorGradient) =
    show(io, m, cg[get_range(100)])

function Base.get(cg::CategoricalColorGradient, x::AbstractFloat, rangescale = (0.0, 1.0))
    isfinite(x) || return invisible()
    rangescale = get_rangescale(rangescale)
    x = clamp.(x, rangescale...)
    if rangescale != (0.0, 1.0)
        x = ColorSchemes.remap(x, rangescale..., 0, 1)
    end
    cg.colors[x == 0 ? 1 : findlast(<(x), cg.values)]
end

Base.reverse(cg::CategoricalColorGradient) =
    CategoricalColorGradient(reverse(cg.colors), reverse(1 .- cg.values))

# required for GR's nonuniformcellarray
function ColorSchemes.getinverse(cg::CategoricalColorGradient, c)
    alpha(c) == 0 && return NaN
    i = findfirst(==(c), color_list(cg))
    (cg.values[i] + cg.values[i + 1]) / 2
end

function prepare_categorical_cgrad_colors(c, v)
    v = sort(unique(clamp.(v, 0, 1)))
    0 ∈ v || pushfirst!(v, 0)
    1 ∈ v || push!(v, 1)
    colors = map(c[get_range(length(v) - 1)]) do col
        RGBA(RGB(col), clamp(alpha(col), 0, 1))
    end
    ColorScheme(plot_color(colors)), v
end

## cgrad

"""
    cgrad(colors, [values]; categorical = nothing, scale = nothing, rev = false, alpha = nothing)

Construct a Colorgradient from `colors` and `values`.

`colors` can be a symbol for ColorSchemes.jl `ColorScheme`s, a `ColorScheme`, a vector of colors, a `ColorGradient` or a `ColorPalette`.
If `values` is an integer, it specifies the numbers of colors chosen equidistantly from the colorscheme specified by colors.
Otherwise vectors are accepted.
For continuous color gradients `values` indicate where between 0 and 1 the colors are positioned.
For categorical color gradients `values` indicate where a color ends and where a new one begins between 0 and 1.
0 and 1 are added to `values` if not already present.

If `rev` is `true` colors are reversed.
`scale` accepts the symbols `:log`, `:log10`, `:log2`, `:ln`, `:exp`, `:exp10` or functions.
If `alpha` is set, it is applied to all colors.
"""
function cgrad(
    colors::ColorScheme,
    values;
    categorical::Union{Nothing,Bool} = nothing,
    scale = nothing,
    rev = false,
    alpha = nothing,
)
    if categorical ≢ nothing && categorical
        colors, values = prepare_categorical_cgrad_colors(colors, values)
    end

    if alpha ≢ nothing
        rgbs = convert.(RGB, colors.colors)
        colors = ColorScheme(RGBA.(rgbs, alpha))
    end
    rev && (colors = reverse(colors))
    values = if scale ∈ (:log, :log10) || scale ≡ log10
        log10.(ColorSchemes.remap(values, 0, 1, 1, 10))
    elseif scale ≡ :log2 || scale ≡ log2
        log2.(ColorSchemes.remap(values, 0, 1, 1, 2))
    elseif scale ≡ :ln || scale ≡ log
        log.(ColorSchemes.remap(values, 0, 1, 1, ℯ))
    elseif scale ∈ (:exp, :exp10) || scale ≡ exp10 || scale ≡ exp
        ColorSchemes.remap(exp10.(values), 1, 10, 0, 1)
    elseif scale isa Function
        v = scale.(values)
        ColorSchemes.remap(v, extrema(v)..., 0, 1)
    else
        values
    end

    if categorical ≢ nothing && categorical
        CategoricalColorGradient(colors, values)
    else
        ContinuousColorGradient(colors, values)
    end
end

function cgrad(
    colors::ColorScheme,
    n::Int = length(colors);
    categorical = nothing,
    kwargs...,
)
    values = get_range(n + (categorical ≢ nothing))
    cgrad(colors, values; categorical = categorical, kwargs...)
end

function cgrad(colors, args...; kwargs...)
    colors ≡ :default && (colors = :inferno)
    cgrad(get_colorscheme(colors), args...; kwargs...)
end

cgrad(; kw...) = cgrad(DEFAULT_COLOR_GRADIENT[]; kw...)

default_cgrad(cg; kw...) = DEFAULT_COLOR_GRADIENT[] = cgrad(cg; kw...)

function get_rangescale(rangescale)
    rangescale ≡ :clamp && return (0.0, 1.0)
    rangescale ≡ :extrema && return extrema(x)
    (rangescale isa NTuple{2,Number}) || error(
        "rangescale ($rangescale) not supported, should be :clamp, :extrema or tuple (minVal, maxVal).  Got $(rangescale).",
    )
    rangescale
end

## ColorPalette

struct ColorPalette <: AbstractColorList
    colors::ColorScheme
end

plot_color(cp::ColorPalette, α::Number) = palette(cp, alpha = α)

Base.reverse(cp::ColorPalette) = ColorPalette(reverse(cp.colors))

"""
    palette(colors, [n]; rev = false, alpha = nothing)

Construct a `ColorPalette`.
Accepts symbols for Colorschemes.jl `ColorScheme`s, `ColorScheme`s, vectors of colors, `ColorGradient`s and `ColorPalettes`.
If 'n' is provided, `n` colors are chosen equidistantly from the colorscheme.
If `rev` is `true` colors are reversed.
"""
function palette(cs; rev = false, alpha = nothing)
    cs = get_colorscheme(cs)
    if alpha ≢ nothing
        rgbs = convert.(RGB, cs.colors)
        cs = ColorScheme(RGBA.(rgbs, alpha))
    end
    rev && (cs = reverse(cs))
    ColorPalette(cs)
end

function palette(cs, n; kwargs...)
    cs = get_colorscheme(cs)
    z = get_range(n)
    palette(cs[z]; kwargs...)
end

## Utils

get_range(n::Int) = range(0, stop = n == 1 ? 0 : 1, length = n)
get_range(cs) = get_range(length(cs))

get_colorscheme(v::AbstractVector{<:Colorant}) = ColorScheme(v)
get_colorscheme(v::AbstractVector) = ColorScheme(parse.(Colorant, v))
function get_colorscheme(sym::Symbol)
    haskey(MISC_COLORSCHEMES, sym) && return MISC_COLORSCHEMES[sym]
    sym = get(COLORSCHEME_ALIASES, sym, sym)
    if sym ≡ :default || sym ≡ :auto
        _default_colorscheme
    elseif haskey(ColorSchemes.colorschemes, sym)
        ColorSchemes.colorschemes[sym]
    else
        error(
            "Unknown ColorScheme `:$sym`. Check https://juliagraphics.github.io/ColorSchemes.jl/stable/ for available ColorSchemes.",
        )
    end
end
get_colorscheme(cs::ColorScheme) = cs

function cvec(cs, n = 10; kw...)
    cg = cgrad(cs; kw...)
    RGBA{Float64}[cg[z] for z ∈ get_range(n)]
end

color_list(c) = get_colorscheme(c).colors
color_list(v::AbstractVector) = plot_color.(v)

get_color_palette(v, n) = palette(v)
get_color_palette(cg::ColorGradient, n) = palette(cg[get_zvalues(n)])

to_rgb(c::Colorant) = RGB(c)
to_rgb(cs::ColorScheme) = ColorScheme(to_rgb.(cs.colors))
to_rgb(cg::ColorGradient) = ColorGradient(to_rgb(cg.colors))
to_rgb(cp::ColorPalette) = ColorPalette(to_rgb(cp.colors))

# allows passing a ColorGradient to rgba_string and get a useful response by picking the first color - introduced because the plotly backend to Plots uses this functionality
rgba_string(cg::T) where {T<:Union{ColorScheme,ColorGradient,ColorPalette}} =
    rgba_string(cg[1])

is_colorscheme(sym) =
    sym ∈ keys(ColorSchemes.colorschemes) ||
    sym ∈ keys(COLORSCHEME_ALIASES) ||
    sym ∈ keys(MISC_COLORSCHEMES)

const DEFAULT_COLOR_GRADIENT = Ref(cgrad(ColorSchemes.colorschemes[:inferno]))

## Compat

const COLORSCHEME_ALIASES = Dict{Symbol,Symbol}(
    # colorbrewer
    :Spectral => :Spectral_11,
    :RdYlGn => :RdYlGn_11,
    :RdBu => :RdBu_11,
    :PiYG => :PiYG_11,
    :PRGn => :PRGn_11,
    :RdYlBu => :RdYlBu_11,
    :BrBG => :BrBG_11,
    :RdGy => :RdGy_11,
    :Set2 => :Set2_8,
    :Accent => :Accent_8,
    :Set1 => :Set1_9,
    :Set3 => :Set3_12,
    :Dark2 => :Dark2_8,
    :Paired => :Paired_12,
    :Pastel2 => :Pastel2_8,
    :Pastel1 => :Pastel1_9,
    :PuOr => :PuOr_11,
    :OrRd => :OrRd_9,
    :PuBu => :PuBu_9,
    :BuPu => :BuPu_9,
    :Oranges => :Oranges_9,
    :BuGn => :BuGn_9,
    :YlOrBr => :YlOrBr_9,
    :YlGn => :YlGn_9,
    :Reds => :Reds_9,
    :RdPu => :RdPu_9,
    :Greens => :Greens_9,
    :YlGnBu => :YlGnBu_9,
    :Purples => :Purples_9,
    :GnBu => :GnBu_9,
    :Greys => :Greys_9,
    :YlOrRd => :YlOrRd_8,
    :PuRd => :PuRd_9,
    :Blues => :Blues_9,
    :PuBuGn => :PuBuGn_9,
    # colorcet
    :bgy => :linear_bgy_10_95_c74_n256,
    :bgyw => :linear_bgyw_15_100_c68_n256,
    :bjy => :diverging_linear_bjy_30_90_c45_n256,
    :bkr => :diverging_bkr_55_10_c35_n256,
    :bky => :diverging_bky_60_10_c30_n256,
    :blues => :linear_blue_5_95_c73_n256,
    :bmw => :linear_bmw_5_95_c86_n256,
    :colorwheel => :cyclic_mygbm_30_95_c78_n256_s25,
    :dimgreys => :linear_grey_10_95_c0_n256,
    :fire => :linear_kryw_0_100_c71_n256,
    :greys => :linear_grey_0_100_c0_n256,
    :gwv => :diverging_gwv_55_95_c39_n256,
    :cinferno => :linear_bmy_10_95_c78_n256,
    :isolum => :isoluminant_cgo_80_c38_n256,
    :kb => :linear_ternary_blue_0_44_c57_n256,
    :kdc => :linear_blue_5_95_c73_n256,
    :kg => :linear_ternary_green_0_46_c42_n256,
    :kgy => :linear_green_5_95_c69_n256,
    :kr => :linear_ternary_red_0_50_c52_n256,
    :diverging => :diverging_bwr_40_95_c42_n256,
    :cyclic1 => :cyclic_mrybm_35_75_c68_n256,
    :cyclic2 => :cyclic_mygbm_30_95_c78_n256,
    :cyclic3 => :cyclic_wrwbw_40_90_c42_n256,
)

## Plots misc colorschemes

const RAINBOW_COLORS = RGBA{Float64}[
    colorant"purple",
    colorant"blue",
    colorant"green",
    colorant"orange",
    colorant"red",
]
const TEST_COLORS = RGBA{Float64}[
    colorant"darkblue",
    colorant"blueviolet",
    colorant"darkcyan",
    colorant"green",
    darken(colorant"yellow", 0.3),
    colorant"orange",
    darken(colorant"red", 0.2),
]

const MISC_COLORSCHEMES = Dict{Symbol,ColorScheme}(
    :blues => ColorScheme(RGBA{Float64}[colorant"lightblue", colorant"darkblue"]),
    :reds => ColorScheme(RGBA{Float64}[colorant"lightpink", colorant"darkred"]),
    :greens => ColorScheme(RGBA{Float64}[colorant"lightgreen", colorant"darkgreen"]),
    :redsblues => ColorScheme(
        RGBA{Float64}[colorant"darkred", RGB(0.8, 0.85, 0.8), colorant"darkblue"],
    ),
    :bluesreds => ColorScheme(
        RGBA{Float64}[colorant"darkblue", RGB(0.8, 0.85, 0.8), colorant"darkred"],
    ),
    :heat => ColorScheme(
        RGBA{Float64}[colorant"lightyellow", colorant"orange", colorant"darkred"],
    ),
    :grays => ColorScheme(RGBA{Float64}[RGB(0.05, 0.05, 0.05), RGB(0.95, 0.95, 0.95)]),
    :rainbow => ColorScheme(RAINBOW_COLORS),
    :lightrainbow => ColorScheme(map(lighten, RAINBOW_COLORS)),
    :darkrainbow => ColorScheme(map(darken, RAINBOW_COLORS)),
    :darktest => ColorScheme(TEST_COLORS),
    :lighttest => ColorScheme(map(c -> lighten(c, 0.3), TEST_COLORS)),
)
