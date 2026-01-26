onfail(body, x) = error("I might have overlooked something: $x")
onfail(body, _::Test.Pass) = nothing
onfail(body, _::Tuple{Test.Fail,T}) where {T} = body()

function is_ci()
    get(ENV, "TRAVIS", "") == "true" ||
    get(ENV, "APPVEYOR", "") in ("true", "True") ||
    get(ENV, "CI", "") in ("true", "True")
end

absdiff(x, y) = x > y ? x - y : y - x

function imdiff(a, b)
    a_view = channelview(a)
    b_view = channelview(b)
    diffscale = max(maxabsfinite(a_view), maxabsfinite(b_view))
    d = sum(absdiff.(a_view, b_view))
    return d / (length(a) * diffscale)
end

_convert(C, T, xs::AbstractArray) =
    collect(colorview(C, map(i -> collect(reinterpret(T, collect(xs)[:, :, i])), 1:size(xs, 3))...))
_convert(C, T, xs::AbstractMatrix) = collect(colorview(C, collect(reinterpret(T, collect(xs)))))

_standardize_grayness(x) = x
_standardize_grayness(x::AbstractArray{<:Gray{Bool}}) = Gray{N0f8}.(x)
_standardize_grayness(x::AbstractArray{<:RGB}) = all(red.(x) .≈ green.(x) .≈ blue.(x)) ? Gray.(red.(x)) : x
_standardize_grayness(x::AbstractArray{<:RGBA}) = all(red.(x) .≈ green.(x) .≈ blue.(x)) ?
    GrayA.(colorview(GrayA, red.(x), alpha.(x))) :
    x


function plotdiffs(p, i)
    """
    Returns a Matrix of Colorants with following structure:
    --------------|---------------|----------------------|---------------------
    Orig Image p  | Orig Image i  | AbsDiff all channels | RelDiff all channels
    --------------|---------------|----------------------|---------------------
    1st channel p | 1st channel i | AbsDiff 1st channel  | RelDiff 1st channel
    --------------|---------------|----------------------|---------------------
    ...
    """
    _p = RGBA{Float64}.(p)
    __p = collect(channelview(_p))
    __p[4, :, :] .= 0.0
    _i = RGBA{Float64}.(i)
    __i = collect(channelview(_i))
    __i[4, :, :] .= 1.0 # make the difference in alpha be 1 so it is visible
    d = collect(colorview(RGBA{Float64}, absdiff.(__p, __i)))
    o = _plotdiffs(p, i)
    hborder = hcat(
        fill(RGBA{Float64}(0, 0, 0, 1), size(o, 1) + size(p, 1), 1),
        fill(RGBA{Float64}(1, 1, 1, 1), size(o, 1) + size(p, 1), 1),
        fill(RGBA{Float64}(0, 0, 0, 1), size(o, 1) + size(p, 1), 1),
    )

    collect(hcat(
        vcat(_p, _plotchannels(p)),
        hborder,
        vcat(_i, _plotchannels(i)),
        hborder,
        vcat(d, o),
        hborder,
        vcat(_rescale(d), _rescale(o)),
    ))
end

vborder(p) = vcat(
    fill(RGBA{Float64}(0, 0, 0, 1), 1, size(p, 2)),
    fill(RGBA{Float64}(1, 1, 1, 1), 1, size(p, 2)),
    fill(RGBA{Float64}(0, 0, 0, 1), 1, size(p, 2)),
)

function _plotdiffs(p, i)
    vcat(
        vborder(p),
        RGBA{Float64}.(absdiff.(red.(i), red.(p)), 0, 0, 1),
        vborder(p),
        RGBA{Float64}.(0, absdiff.(green.(i), green.(p)), 0, 1),
        vborder(p),
        RGBA{Float64}.(0, 0, absdiff.(blue.(i), blue.(p)), 1),
    )
end

function _plotdiffs(p::AbstractArray{<:RGBA}, i::AbstractArray{<:RGBA})
    a = absdiff.(alpha.(i), alpha.(p))
    vcat(
        vborder(p),
        RGBA{Float64}.(absdiff.(red.(i), red.(p)), 0, 0, 1),
        vborder(p),
        RGBA{Float64}.(0, absdiff.(green.(i), green.(p)), 0, 1),
        vborder(p),
        RGBA{Float64}.(0, 0, absdiff.(blue.(i), blue.(p)), 1),
        vborder(p),
        RGBA{Float64}.(a, a, a, 1),
    )
end

function _plotchannels(p)
    vcat(
        vborder(p),
        RGBA{Float64}.(red.(p), 0, 0, 1),
        vborder(p),
        RGBA{Float64}.(0, green.(p), 0, 1),
        vborder(p),
        RGBA{Float64}.(0, 0, blue.(p), 1),
    )
end

function _plotchannels(p::AbstractArray{<:RGBA})
    vcat(
        vborder(p),
        RGBA{Float64}.(red.(p), 0, 0, 1),
        vborder(p),
        RGBA{Float64}.(0, green.(p), 0, 1),
        vborder(p),
        RGBA{Float64}.(0, 0, blue.(p), 1),
        vborder(p),
        RGBA{Float64}.(alpha.(p), alpha.(p), alpha.(p), 1),
    )
end

function _rescale(x::AbstractArray)
    c = channelview(x)
    m = maximum(c)
    m = iszero(m) ? one(m) : m
    colorview(RGBA{Float64}, 100 * c ./ m)
end


## COPIED FROM Images.jl/src/algorithms.jl

"""
`m = maxabsfinite(A)` calculates the maximum absolute value in `A`, ignoring any values that are not finite (Inf or NaN).
"""
function maxabsfinite(A::AbstractArray{T}) where T
    ret = sentinel_min(typeof(abs(A[1])))
    for a in A
        ret = maxfinite_scalar(abs(a), ret)
    end
    ret
end

minfinite_scalar(a::T, b::T) where {T} = isfinite(a) ? (b < a ? b : a) : b
maxfinite_scalar(a::T, b::T) where {T} = isfinite(a) ? (b > a ? b : a) : b
minfinite_scalar(a::T, b::T) where {T<:Union{Integer,FixedPoint}} = b < a ? b : a
maxfinite_scalar(a::T, b::T) where {T<:Union{Integer,FixedPoint}} = b > a ? b : a
minfinite_scalar(a, b) = minfinite_scalar(promote(a, b)...)
maxfinite_scalar(a, b) = maxfinite_scalar(promote(a, b)...)

function minfinite_scalar(c1::C, c2::C) where C<:AbstractRGB
    C(minfinite_scalar(c1.r, c2.r),
      minfinite_scalar(c1.g, c2.g),
      minfinite_scalar(c1.b, c2.b))
end
function maxfinite_scalar(c1::C, c2::C) where C<:AbstractRGB
    C(maxfinite_scalar(c1.r, c2.r),
      maxfinite_scalar(c1.g, c2.g),
      maxfinite_scalar(c1.b, c2.b))
end

sentinel_min(::Type{T}) where {T<:Union{Integer,FixedPoint}} = typemax(T)
sentinel_max(::Type{T}) where {T<:Union{Integer,FixedPoint}} = typemin(T)
sentinel_min(::Type{T}) where {T<:AbstractFloat} = convert(T, NaN)
sentinel_max(::Type{T}) where {T<:AbstractFloat} = convert(T, NaN)
sentinel_min(::Type{C}) where {C<:AbstractRGB} = _sentinel_min(C, eltype(C))
_sentinel_min(::Type{C},::Type{T}) where {C<:AbstractRGB,T} = (s = sentinel_min(T); C(s,s,s))
sentinel_max(::Type{C}) where {C<:AbstractRGB} = _sentinel_max(C, eltype(C))
_sentinel_max(::Type{C},::Type{T}) where {C<:AbstractRGB,T} = (s = sentinel_max(T); C(s,s,s))
sentinel_min(::Type{C}) where {C<:AbstractGray} = _sentinel_min(C, eltype(C))
_sentinel_min(::Type{C},::Type{T}) where {C<:AbstractGray,T} = C(sentinel_min(T))
sentinel_max(::Type{C}) where {C<:AbstractGray} = _sentinel_max(C, eltype(C))
_sentinel_max(::Type{C},::Type{T}) where {C<:AbstractGray,T} = C(sentinel_max(T))
