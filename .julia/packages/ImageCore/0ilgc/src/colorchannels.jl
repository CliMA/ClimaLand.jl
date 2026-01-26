# Create a special type for permutations. The real point here is to be able to
# unambiguously identify an RRPermArray (see below) so that we may "unwrap" in
# expressions like `channelview(colorview(C, A))`.

"""
ColorChanPerm(perm)

Construct a reordering permutation for the color channel.
This handles swaps between memory layout and constructor argument order for `AbstractRGB` and
various `AlphaChannel` and `ChannelAlpha` color types.
"""
struct ColorChanPerm{N} <: AbstractVector{Int}
    perm::NTuple{N,Int}
end
Base.IndexStyle(::Type{<:ColorChanPerm}) = IndexLinear()
Base.size(v::ColorChanPerm{N}) where N = (N,)
Base.getindex(v::ColorChanPerm, i::Int) = v.perm[i]

dimorder(::Type{<:RGB}) = ColorChanPerm((1, 2, 3))
dimorder(::Type{<:BGR}) = ColorChanPerm((3, 2, 1))
dimorder(::Type{<:XRGB}) = ColorChanPerm((1, 1, 2, 3))
dimorder(::Type{<:RGBX}) = ColorChanPerm((1, 2, 3, 3)) # this causes problems for setindex!, fixed below
dimorder(::Type{<:BGRA}) = ColorChanPerm((3, 2, 1, 4))
dimorder(::Type{<:ABGR}) = ColorChanPerm((4, 3, 2, 1))
dimorder(::Type{<:AlphaColor{<:Color1,T,N}}) where {T,N} = ColorChanPerm((2, 1))
dimorder(::Type{<:AlphaColor{<:Color3,T,N}}) where {T,N} = ColorChanPerm((4, 1, 2, 3))

const ColorChanPermSubArray{T,N,P,I<:Tuple,L} =
    SubArray{T,N,P,I,L}
const RRPermArray{To,From,M,P<:ColorChanPermSubArray} =
    RRArray{To,From,M,P}

# This type exists solely to set multiple values in the color channel axis
struct NVector{T,N} <: AbstractVector{T}
    v::NTuple{N,T}

    NVector{T,N}(x::NTuple{N}) where {T,N} = new{T,N}(x)
end
Base.IndexStyle(::Type{<:NVector}) = IndexLinear()
Base.size(v::NVector{T,N}) where {T,N} = (N,)
Base.getindex(v::NVector, i::Int) = v.v[i]
NVector(x1::T, x::Vararg{T,N}) where {T,N} = NVector{T,N+1}((x1, x...))

@inline Base.setindex!(A::RRPermArray{<:RGBX,<:Number,N}, val::AbstractRGB, i::Vararg{Int,N}) where N =
    setindex!(parent(parent(parent(A))), NVector(red(val), green(val), blue(val)), :, i...)

"""
    channelview(A)

returns a view of `A`, splitting out (if necessary) the color channels
of `A` into a new first dimension.

Of relevance for types like RGB and BGR, the channels of the returned
array will be in constructor-argument order, not memory order (see
`reinterpretc` if you want to use memory order).

# Example
```julia
img = rand(RGB{N0f8}, 10, 10)
A = channelview(img)   # a 3×10×10 array
```

See also: [`colorview`](@ref)
"""
channelview(A::AbstractArray{T}) where {T<:Number} = A
channelview(A::RRArray{<:Colorant,<:Number}) = parent(A)
channelview(A::RRPermArray{<:Colorant,<:Number}) = parent(parent(A))
channelview(A::Base.ReinterpretArray{<:AbstractGray,M,<:Number}) where M = parent(A)
channelview(A::AbstractArray{RGB{T}}) where {T} = reinterpretc(T, A)
function channelview(A::AbstractArray{C}) where {C<:AbstractRGB}
    # BGR, XRGB, etc don't satisfy conditions for reinterpret
    CRGB = RGB{eltype(C)}
    channelview(of_eltype(CRGB, A))
end
channelview(A::AbstractArray{C}) where {C<:Color} = reinterpretc(eltype(C), A)
channelview(A::AbstractArray{C}) where {C<:ColorAlpha} = _channelview(color_type(C), A)
_channelview(::Type{<:RGB}, A) = reinterpretc(eltype(eltype(A)), A)
function _channelview(::Type{C}, A) where {C<:AbstractRGB}
    CRGBA = RGBA{eltype(C)}
    channelview(of_eltype(CRGBA, A))
end
_channelview(::Type{C}, A) where {C<:Color} = reinterpretc(eltype(C), A)
function channelview(A::AbstractArray{AC}) where {AC<:AlphaColor}
    CA = coloralpha(base_color_type(AC)){eltype(AC)}
    channelview(of_eltype(CA, A))
end

"""
    colorview(C, A)

returns a view of the numeric array `A`, interpreting successive
elements of `A` as if they were channels of Colorant `C`.

Of relevance for types like RGB and BGR, the elements of `A` are
interpreted in constructor-argument order, not memory order (see
`reinterpretc` if you want to use memory order).

# Example
```jl
A = rand(3, 10, 10)
img = colorview(RGB, A)
```

See also: [`channelview`](@ref)
"""
colorview(::Type{C}, A::AbstractArray{T}) where {C<:Colorant,T<:Number} =
    _ccolorview(ccolor_number(C, T), A)
_ccolorview(::Type{C}, A::RRPermArray{T,C}) where {C<:Colorant,T<:Number} =
    parent(parent(parent(A)))
_ccolorview(::Type{C}, A::RRArray{T,C}) where {C<:Colorant,T<:Number} =
    parent(parent(A))
_ccolorview(::Type{C}, A::Base.ReinterpretArray{T,M,C,AA,false}) where {C<:RGB,T<:Number,M,AA} =
    reshape(parent(A), Base.tail(axes(parent(A))))
# _ccolorview(::Type{C}, A::Base.ReinterpretArray{T,M,C,AA,true}) where {C<:RGB,T<:Number,M,AA} =
#     parent(A)
_ccolorview(::Type{C}, A::Base.ReinterpretArray{T,M,C,AA,false}) where {C<:Color,T<:Number,M,AA} =
    reshape(parent(A), Base.tail(axes(parent(A))))
# _ccolorview(::Type{C}, A::Base.ReinterpretArray{T,M,C,AA,true}) where {C<:Color,T<:Number,M,AA} =
#     parent(A)
_ccolorview(::Type{C}, A::AbstractArray{T}) where {C<:Colorant,T<:Number} =
    __ccolorview(C, A)  # necessary to avoid ambiguities from dispatch on eltype
__ccolorview(::Type{C}, A::AbstractArray{T}) where {T<:Number,C<:RGB{T}} = reinterpretc(C, A)
__ccolorview(::Type{C}, A::AbstractArray{T}) where {T<:Number,C<:AbstractRGB} =
    _colorview_reorder(C, A)
__ccolorview(::Type{C}, A::AbstractArray{T}) where {T<:Number,C<:Color{T}} = reinterpretc(C, A)
__ccolorview(::Type{C}, A::AbstractArray{T}) where {T<:Number,C<:ColorAlpha} =
    _colorviewalpha(base_color_type(C), C, eltype(C), A)
__ccolorview(::Type{C}, A::AbstractArray{T}) where {T<:Number,C<:AlphaColor} =
    _colorview_reorder(C, A)
_colorviewalpha(::Type{C}, ::Type{CA}, ::Type{T}, A::AbstractArray{T}) where {C<:RGB,CA,T} =
    reinterpretc(CA, A)
_colorviewalpha(::Type{C}, ::Type{CA}, ::Type{T}, A::AbstractArray{T}) where {C<:AbstractRGB,CA,T} =
    _colorview_reorder(CA, A)
_colorviewalpha(::Type{C}, ::Type{CA}, ::Type{T}, A::AbstractArray{T}) where {C<:Color,CA,T} =
    reinterpretc(CA, A)

_colorview_reorder(::Type{C}, A) where C = reinterpretc(C, view(A, dimorder(C), Base.tail(colons(A))...))

colorview(::Type{ARGB32}, A::AbstractArray{BGRA{N0f8}}) = reinterpret(ARGB32, A)

colorview(::Type{C1}, A::AbstractArray{C2}) where {C1<:Colorant,C2<:Colorant} =
    colorview(C1, channelview(A))

colons(A::AbstractArray{T,N}) where {T,N} = ntuple(d->Colon(), Val(N))

"""
    colorview(C, gray1, gray2, ...) -> imgC

Combine numeric/grayscale images `gray1`, `gray2`, etc., into the
separate color channels of an array `imgC` with element type
`C<:Colorant`.

As a convenience, the constant `zeroarray` fills in an array of
matched size with all zeros.

# Example
```julia
imgC = colorview(RGB, r, zeroarray, b)
```

creates an image with `r` in the red chanel, `b` in the blue channel,
and nothing in the green channel.

See also: [`StackedView`](@ref).
"""
function colorview(::Type{C}, gray1, gray2, grays...) where C<:Colorant
    T = _colorview_type(eltype(C), promote_eleltype_all(gray1, gray2, grays...))
    CT = base_colorant_type(C){T}
    axs = firstinds(gray1, gray2, grays...)
    mappedarray(CT, extractchannels, take_zeros(eltype(CT), axs, gray1, gray2, grays...)...)
end

"""
    colorview(C)

Create a function that is equivalent to `(As...) -> colorview(C, Ax...)`.

# Examples

```jldoctest; setup = :(using ImageCore)
julia> ones(Float32, 2, 2) |> colorview(Gray)
2×2 reinterpret(reshape, Gray{Float32}, ::Matrix{Float32}) with eltype Gray{Float32}:
 1.0  1.0
 1.0  1.0
```

This can be slightly convenient when you want to convert a batch of channel data, for example:

```julia
julia> Rs, Gs, Bs = ntuple( i -> [randn(2, 2) for _ in 1:4], 3)

julia> map(colorview(RGB), Rs, Gs, Bs)
```
"""
colorview(::Type{C}) where C<:Colorant = (As...) -> colorview(C, As...)

_colorview_type(::Type{Any}, ::Type{T}) where {T} = T
_colorview_type(::Type{T1}, ::Type{T2}) where {T1,T2} = T1

Base.@pure @inline promote_eleltype_all(gray, grays...) = _promote_eleltype_all(beltype(eltype(gray)), grays...)
@inline function _promote_eleltype_all(::Type{T}, gray, grays...) where T
    _promote_eleltype_all(promote_type(T, beltype(eltype(gray))), grays...)
end
_promote_eleltype_all(::Type{T}) where {T} = T

beltype(::Type{T}) where {T} = eltype(T)
beltype(::Type{Union{}}) = Union{}

extractchannels(c::AbstractGray)    = (gray(c),)
extractchannels(c::TransparentGray) = (gray(c), alpha(c))
extractchannels(c::Color3)          = (comp1(c), comp2(c), comp3(c))
extractchannels(c::Transparent3)    = (comp1(c), comp2(c), comp3(c), alpha(c))

## Tuple & indexing utilities

_size(A::AbstractArray) = map(length, axes(A))

# color->number
@inline channelview_size(parent::AbstractArray{C}) where {C<:Colorant} = (length(C), _size(parent)...)
@inline channelview_axes(parent::AbstractArray{C}) where {C<:Colorant} =
    _cvi(Base.OneTo(length(C)), axes(parent))
_cvi(rc, ::Tuple{}) = (rc,)
_cvi(rc, inds::Tuple{R,Vararg{R}}) where {R<:AbstractUnitRange} = (convert(R, rc), inds...)
@inline channelview_size(parent::AbstractArray{C}) where {C<:Color1} = _size(parent)
@inline channelview_axes(parent::AbstractArray{C}) where {C<:Color1} = axes(parent)

function check_ncolorchan(::AbstractArray{C}, dims) where C<:Colorant
    dims[1] == length(C) || throw(DimensionMismatch("new array has $(dims[1]) color channels, must have $(length(C))"))
end
chanparentsize(::AbstractArray{C}, dims) where {C<:Colorant} = tail(dims)
@inline colparentsize(::Type{C}, dims) where {C<:Colorant} = (length(C), dims...)

channelview_dims_offset(parent::AbstractArray{C}) where {C<:Colorant} = 1

check_ncolorchan(::AbstractArray{C}, dims) where {C<:Color1} = nothing
chanparentsize(::AbstractArray{C}, dims) where {C<:Color1} = dims
colparentsize(::Type{C}, dims) where {C<:Color1} = dims
channelview_dims_offset(parent::AbstractArray{C}) where {C<:Color1} = 0

@inline indexsplit(A::AbstractArray{C}, I) where {C<:Colorant} = I[1], tail(I)
@inline indexsplit(A::AbstractArray{C}, I) where {C<:Color1} = 1, I

# number->color
@inline colorview_size(::Type{C}, parent::AbstractArray) where {C<:Colorant} = tail(_size(parent))
@inline colorview_axes(::Type{C}, parent::AbstractArray) where {C<:Colorant} = tail(axes(parent))
@inline colorview_size(::Type{C}, parent::AbstractArray) where {C<:Color1} = _size(parent)
@inline colorview_axes(::Type{C}, parent::AbstractArray) where {C<:Color1} = axes(parent)

function checkdim1(::Type{C}, inds) where C<:Colorant
    inds[1] == (1:length(C)) || throw(DimensionMismatch("dimension 1 must have indices 1:$(length(C)), got $(inds[1])"))
    nothing
end
checkdim1(::Type{C}, dims) where {C<:Color1} = nothing

parentaxes(::Type, inds) = tail(inds)
parentaxes(::Type{C}, inds) where {C<:Color1} = inds

celtype(::Type{Any}, ::Type{T}) where {T} = T
celtype(::Type{T1}, ::Type{T2}) where {T1,T2} = T1

## Low-level color utilities

tuplify(c::Color1) = (comp1(c),)
tuplify(c::Color3) = (comp1(c), comp2(c), comp3(c))
tuplify(c::Color2) = (comp1(c), alpha(c))
tuplify(c::Color4) = (comp1(c), comp2(c), comp3(c), alpha(c))

"""
    getchannels(P, C::Type, I)

Get a tuple of all channels needed to construct a Colorant of type `C`
from an `P::AbstractArray{<:Number}`.
"""
getchannels
@inline getchannels(P, ::Type{C}, I) where {C<:Color1} = (@inbounds ret = (P[I...],); ret)
@inline getchannels(P, ::Type{C}, I::Real) where {C<:Color1} = (@inbounds ret = (P[I],); ret)
@inline function getchannels(P, ::Type{C}, I) where C<:Color2
    @inbounds ret = (P[1,I...], P[2,I...])
    ret
end
@inline function getchannels(P, ::Type{C}, I) where C<:Color3
    @inbounds ret = (P[1,I...], P[2,I...],P[3,I...])
    ret
end
@inline function getchannels(P, ::Type{C}, I) where C<:Color4
    @inbounds ret = (P[1,I...], P[2,I...], P[3,I...], P[4,I...])
    ret
end

# setchannel (similar to setfield!)
# These don't check bounds since that's already done
"""
    setchannel(c, val, idx)

Equivalent to:

    cc = copy(c)
    cc[idx] = val
    cc

for immutable colors. `idx` is interpreted in the sense of constructor
arguments, so `setchannel(c, 0.5, 1)` would set red color channel for
any `c::AbstractRGB`, even if red isn't the first field in the type.
"""
setchannel(c::Colorant{T,1}, val, Ic::Int) where {T} = typeof(c)(val)

setchannel(c::TransparentColor{C,T,2}, val, Ic::Int) where {C,T} =
    typeof(c)(ifelse(Ic==1,val,comp1(c)),
              ifelse(Ic==2,val,alpha(c)))

setchannel(c::Colorant{T,3}, val, Ic::Int) where {T} = typeof(c)(ifelse(Ic==1,val,comp1(c)),
                                                          ifelse(Ic==2,val,comp2(c)),
                                                          ifelse(Ic==3,val,comp3(c)))
setchannel(c::TransparentColor{C,T,4}, val, Ic::Int) where {C,T} =
    typeof(c)(ifelse(Ic==1,val,comp1(c)),
              ifelse(Ic==2,val,comp2(c)),
              ifelse(Ic==3,val,comp3(c)),
              ifelse(Ic==4,val,alpha(c)))

"""
    setchannels!(P, val, I)

For a color `val`, distribute its channels along `P[:, I...]` for
`P::AbstractArray{<:Number}`.
"""
setchannels!
@inline setchannels!(P, val::Color1, I) = (@inbounds P[I...] = comp1(val); val)
@inline function setchannels!(P, val::Color2, I)
    @inbounds P[1,I...] = comp1(val)
    @inbounds P[2,I...] = alpha(val)
    val
end
@inline function setchannels!(P, val::Color3, I)
    @inbounds P[1,I...] = comp1(val)
    @inbounds P[2,I...] = comp2(val)
    @inbounds P[3,I...] = comp3(val)
    val
end
@inline function setchannels!(P, val::Color4, I)
    @inbounds P[1,I...] = comp1(val)
    @inbounds P[2,I...] = comp2(val)
    @inbounds P[3,I...] = comp3(val)
    @inbounds P[4,I...] = alpha(val)
    val
end
