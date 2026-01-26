# Infrastructure for value transformation

"""
    clamp01(x) -> y

Produce a value `y` that lies between 0 and 1, and equal to `x` when
`x` is already in this range. Equivalent to `clamp(x, 0, 1)` for
numeric values. For colors, this function is applied to each color
channel separately.

See also: [`clamp01!`](@ref), [`clamp01nan`](@ref).
"""
clamp01(x::Union{N0f8,N0f16}) = x
clamp01(x::NumberLike) = clamp(x, zero(x), oneunit(x))
clamp01(c::Union{TransparentGray,AbstractRGB,TransparentRGB}) = mapc(clamp01, c)

"""
    clamp01!(array::AbstractArray)

Restrict values in array to [0, 1], in-place. See also [`clamp01`](@ref).
"""
function clamp01!(img::AbstractArray)
    # Slightly faster than map!(clamp01, img, img)
    @inbounds for i in eachindex(img)
        img[i] = clamp01(img[i])
    end
    img
end

"""
    clamp01nan(x) -> y

Similar to `clamp01`, except that any `NaN` values are changed to 0.

See also: [`clamp01nan!`](@ref), [`clamp01`](@ref).
"""
clamp01nan(x) = clamp01(x)
clamp01nan(x::Union{AbstractFloat,AbstractGray{<:AbstractFloat}}) = ifelse(isnan(x), zero(x), clamp01(x))
clamp01nan(c::Union{TransparentGray,AbstractRGB,TransparentRGB}) = mapc(clamp01nan, c)

"""
    clamp01nan!(array::AbstractArray)

Similar to `clamp01!`, except that any `NaN` values are changed to 0.

See also: [`clamp01!`](@ref), [`clamp01nan`](@ref)
"""
function clamp01nan!(img::GenericImage)
    # slgihtly faster than map!(clamp01nan, img, img)
    @inbounds for i in eachindex(img)
        img[i] = clamp01nan(img[i])
    end
    img
end

"""
    scaleminmax(min, max) -> f
    scaleminmax(T, min, max) -> f

Return a function `f` which maps values less than or equal to `min` to
0, values greater than or equal to `max` to 1, and uses a linear scale
in between. `min` and `max` should be real values.

Optionally specify the return type `T`. If `T` is a colorant (e.g.,
RGB), then scaling is applied to each color channel.

# Examples

## Example 1

```julia
julia> f = scaleminmax(-10, 10)
(::#9) (generic function with 1 method)

julia> f(10)
1.0

julia> f(-10)
0.0

julia> f(5)
0.75
```

## Example 2

```julia
julia> c = RGB(255.0,128.0,0.0)
RGB{Float64}(255.0,128.0,0.0)

julia> f = scaleminmax(RGB, 0, 255)
(::#13) (generic function with 1 method)

julia> f(c)
RGB{Float64}(1.0,0.5019607843137255,0.0)
```

See also: [`takemap`](@ref).
"""
function scaleminmax(min::T, max::T) where T
    @inline function(x)
        xp, minp, maxp = promote(x, min, max)  # improves performance to promote explicitly
        y = clamp(xp, minp, maxp)
        (y-minp)/(maxp-minp)
    end
end
function scaleminmax(::Type{Tout}, min::T, max::T) where {Tout,T}
    @inline function(x)
        xp, minp, maxp = promote(x, min, max)
        y = clamp(xp, minp, maxp)
        smmconvert(Tout, (y-minp)/(maxp-minp))
    end
end
@inline smmconvert(::Type{Tout}, x) where {Tout} = convert(Tout, x)
# since we know the result will be between 0 and 1, we can use rem to save a check
@inline smmconvert(::Type{Tout}, x) where {Tout<:Normed} = rem(x, Tout)

scaleminmax(min, max) = scaleminmax(promote(min, max)...)
scaleminmax(::Type{T}, min, max) where {T} = scaleminmax(T, promote(min, max)...)

# TODO: use triangular dispatch when we can count on Julia 0.6+
function scaleminmax(::Type{C}, min::T, max::T) where {C<:Colorant, T<:Real}
    _scaleminmax(C, eltype(C), min, max)
end
function _scaleminmax(::Type{C}, ::Type{TC}, min::T, max::T) where {C<:Colorant, TC<:Real, T<:Real}
    freal = scaleminmax(TC, min, max)
    @inline function(c)
        C(mapc(freal, c))
    end
end
function _scaleminmax(::Type{C}, ::Type{Any}, min::T, max::T) where {C<:Colorant, T<:Real}
    freal = scaleminmax(min, max)
    @inline function(c)
        C(mapc(freal, c))
    end
end

function takemap(::typeof(scaleminmax), A::AbstractArray{T}) where T<:Real
    min, max = extrema(A)
    scaleminmax(min, max)
end
function takemap(::typeof(scaleminmax), A::AbstractArray{C}) where C<:Colorant
    min, max = extrema(channelview(A))
    scaleminmax(C, min, max)
end
function takemap(::typeof(scaleminmax), ::Type{Tout}, A::AbstractArray{T}) where {Tout,T<:Real}
    min, max = extrema(A)
    scaleminmax(Tout, min, max)
end
function takemap(::typeof(scaleminmax), ::Type{Cout}, A::AbstractArray{C}) where {Cout<:Colorant,C<:Colorant}
    min, max = extrema(channelview(A))
    scaleminmax(Cout, min, max)
end

"""
    scalesigned(maxabs) -> f

Return a function `f` which scales values in the range `[-maxabs,
maxabs]` (clamping values that lie outside this range) to the range
`[-1, 1]`.

See also: [`colorsigned`](@ref).
"""
function scalesigned(maxabs::Real)
    maxabs > 0 || throw(ArgumentError("maxabs must be positive, got $maxabs"))
    x -> clamp(x, -maxabs, maxabs)/maxabs
end

"""
    scalesigned(min, center, max) -> f

Return a function `f` which scales values in the range `[min, center]`
to `[-1,0]` and `[center,max]` to `[0,1]`. Values smaller than
`min`/`max` get clamped to `min`/`max`, respectively.

See also: [`colorsigned`](@ref).
"""
function scalesigned(min::T, center::T, max::T) where T<:Real
    min <= center <= max || throw(ArgumentError("values must be ordered, got $min, $center, $max"))
    sneg, spos = 1/(center-min), 1/(max-center)
    function(x)
        Δy = clamp(x, min, max) - center
        ifelse(Δy < 0, sneg*Δy, spos*Δy)
    end
end
scalesigned(min::Real, center::Real, max::Real) = scalesigned(promote(min, center, max)...)

function takemap(::typeof(scalesigned), A::AbstractArray{T}) where T<:Real
    mn, mx = extrema(A)
    scalesigned(max(abs(mn), abs(mx)))
end

"""
    colorsigned()
    colorsigned(colorneg, colorpos) -> f
    colorsigned(colorneg, colorcenter, colorpos) -> f

Define a function that maps negative values (in the range [-1,0]) to
the linear colormap between `colorneg` and `colorcenter`, and positive
values (in the range [0,1]) to the linear colormap between
`colorcenter` and `colorpos`.

The default colors are:

- `colorcenter`: white
- `colorneg`: green1
- `colorpos`: magenta

See also: [`scalesigned`](@ref).
"""
colorsigned(neg::C, center::C, pos::C) where {C<:Colorant} = function(x)
    y = clamp(x, -oneunit(x), oneunit(x))
    yabs = norm(y)
    C(ifelse(y>zero(y), weighted_color_mean(yabs, pos, center),
                        weighted_color_mean(yabs, neg, center)))
end

function colorsigned(colorneg::C, colorpos::C) where C<:Colorant
    colorsigned(colorneg, C(colorant"white"), colorpos)
end

colorsigned() = colorsigned(colorant"green1", colorant"magenta")


"""
    takemap(f, A) -> fnew
    takemap(f, T, A) -> fnew

Given a value-mapping function `f` and an array `A`, return a
"concrete" mapping function `fnew`. When applied to elements of `A`,
`fnew` should return valid values for storage or display, for example
in the range from 0 to 1 (for grayscale) or valid colorants. `fnew`
may be adapted to the actual values present in `A`, and may not
produce valid values for any inputs not in `A`.

Optionally one can specify the output type `T` that `fnew` should produce.

# Example:
```julia
julia> A = [0, 1, 1000];

julia> f = takemap(scaleminmax, A)
(::#7) (generic function with 1 method)

julia> f.(A)
3-element Array{Float64,1}:
 0.0
 0.001
 1.0
```
"""
takemap

takemap(f, A) = f
takemap(f, ::Type{T}, A) where {T} = x->T(f(x))
