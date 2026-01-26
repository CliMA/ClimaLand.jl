
invisible() = RGBA{Float64}(0.0, 0.0, 0.0, 0.0)

# the one-arg cases, meant for single colors
plot_color(s::AbstractString) = parse(RGBA{Float64}, s)
plot_color(s::Symbol) =
    try
        parse(RGBA{Float64}, s)
    catch e
        is_colorscheme(s) ? cgrad(s) : rethrow(e)
    end
plot_color(b::Val{true}) = error("plot_color(true) not allowed.")
plot_color(b::Val{false}) = invisible()
plot_color(b::Bool) = plot_color(Val(b))
plot_color(::Nothing) = invisible()
plot_color(c::Colorant) = convert(RGBA{Float64}, c)
# plot_color(cs::AbstractVector) = RGBA{Float64}[plot_color(c) for c ∈ cs]
# plot_color(cs::AbstractArray) = map(plot_color, cs)

# no alpha override
plot_color(x, ::Nothing) = plot_color(x)

# alpha override
plot_color(x, α::Number) = RGBA{Float64}(convert(RGB, plot_color(x)), α)
plot_color(c::Colorant, α::Number) = RGBA{Float64}(red(c), green(c), blue(c), α)
plot_color(s::Symbol, α::Number) = (
    is_colorscheme(s) ? cgrad(s, alpha = α) : RGBA{Float64}(convert(RGB, plot_color(s)), α)
)

function plot_color(cs::AbstractArray)
    a = Array{RGBA{Float64}}(undef, size(cs)...)
    for i ∈ eachindex(cs)
        a[i] = plot_color(cs[i])
    end
    a
end

# plot_color(cs::AbstractVector, α::Number) = RGBA{Float64}[plot_color(c,α) for c ∈ cs]
function plot_color(cs::AbstractArray, α::Number)
    a = Array{RGBA{Float64}}(undef, size(cs)...)
    for i ∈ eachindex(cs)
        a[i] = plot_color(cs[i], α)
    end
    a
end
# map!(c -> plot_color(c,α), a, cs))

# convenience conversions from numeric arrays to gradient values
# note: we need the first version because of dispatch
# function plot_color(zs::AbstractVector{T}) where T<:Number
#     grad = cgrad()
#     zmin, zmax = extrema(zs)
#     RGBA{Float64}[grad[(z-zmin)/(zmax-zmin)] for z ∈ zs]
# end
function plot_color(zs::AbstractArray{T}) where {T<:Number}
    grad = cgrad()
    zmin, zmax = extrema(zs[isfinite.(zs)])
    a = Array{RGBA{Float64}}(undef, size(zs)...)
    for i ∈ eachindex(zs)
        a[i] = grad[(zs[i] - zmin) / (zmax - zmin)]
    end
    a
end

# function plot_color(zs::AbstractVector{T}, α::Number) where T<:Number
#     cs = plot_color(zs)
#     RGBA{Float64}[RGBA{Float64}(convert(RGB, c), α) for c ∈ cs]
# end
function plot_color(zs::AbstractArray{T}, α::Number) where {T<:Number}
    cs = plot_color(zs)
    a = Array{RGBA{Float64}}(undef, size(zs)...)
    for i ∈ eachindex(zs)
        a[i] = RGBA{Float64}(convert(RGB, cs[i]), α)
    end
    a
end
