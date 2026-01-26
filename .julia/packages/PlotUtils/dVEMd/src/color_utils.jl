
# --------------------------------------------------------------

# Methods to automatically generate colorschemes for color selection based on
# background color and a short list of seed colors

# here are some magic constants that could be changed if you really want
const _lightness_darkbg = Ref(80.0)
const _lightness_lightbg = Ref(60.0)
const _lch_c_const = Ref(60)

adjust_lch(color, l, c) = convert(RGBA{Float64}, LCHab(l, c, convert(LCHab, color).h))

function lightness_from_background(bgcolor)
    bglight = convert(LCHab, bgcolor).l
    bglight < 50.0 ? _lightness_darkbg[] : _lightness_lightbg[]
end

function generate_colorscheme(
    bgcolor = plot_color(:white);
    color_bases = plot_color([colorant"steelblue", colorant"orangered"]),
    lightness = lightness_from_background(bgcolor),
    chroma = _lch_c_const[],
    n = 17,
)
    seed_colors = vcat(bgcolor, map(c -> adjust_lch(c, lightness, chroma), color_bases))
    seed_colors = convert(Vector{RGB{Float64}}, seed_colors)
    colors = distinguishable_colors(
        n,
        seed_colors,
        lchoices = Float64[lightness],
        cchoices = Float64[chroma],
        hchoices = range(0; stop = 340, length = 20),
    )[2:end]
    ColorScheme(colors)
end

# ----------------------------------------------------------------------------------

function getpctrange(n::Integer)
    n > 0 || error()
    n == 1 && return zeros(1)
    zs = [0.0, 1.0]
    for i ∈ 3:n
        sorted = sort(zs)
        diffs = diff(sorted)
        widestj = 0
        widest = 0.0
        for (j, d) ∈ enumerate(diffs)
            if d > widest
                widest = d
                widestj = j
            end
        end
        push!(zs, sorted[widestj] + 0.5diffs[widestj])
    end
    zs
end

function get_zvalues(n::Integer)
    offsets = getpctrange(ceil(Int, n / 4) + 1) / 4
    offsets = vcat(offsets[1], offsets[3:end])
    zvalues = Float64[]
    for offset ∈ offsets
        append!(zvalues, offset .+ [0.0, 0.5, 0.25, 0.75])
    end
    vcat(zvalues[1], 1.0, zvalues[2:(n - 1)])
end

# ----------------------------------------------------------------------------------

function darken(c, v = 0.1)
    rgba = convert(RGBA, c)
    r = max(0, min(rgba.r - v, 1))
    g = max(0, min(rgba.g - v, 1))
    b = max(0, min(rgba.b - v, 1))
    RGBA(r, g, b, rgba.alpha)
end
lighten(c, v = 0.3) = darken(c, -v)

# borrowed from http://stackoverflow.com/a/1855903:
lightness_level(c::Colorant) = 0.299red(c) + 0.587green(c) + 0.114blue(c)

isdark(c::Colorant)::Bool = lightness_level(c) < 0.5
islight(c::Colorant)::Bool = !isdark(c)

Base.convert(::Type{RGB}, h::Unsigned) =
    RGB([(x & 0x0000FF) / 0xFF for x ∈ (h >> 16, h >> 8, h)]...)

make255(x)::Int = round(Int, 255 * x)

rgb_string(c::Colorant) =
    @sprintf("rgb(%d, %d, %d)", make255(red(c)), make255(green(c)), make255(blue(c)))

rgba_string(c::Colorant) = @sprintf(
    "rgba(%d, %d, %d, %1.3f)",
    make255(red(c)),
    make255(green(c)),
    make255(blue(c)),
    alpha(c)
)
