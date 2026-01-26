"""
    sixel_encode([io], src, [encoder]; kwargs...)

Encode colorant sequence `src` as sixel sequence and write it into a writable io-like object `io`.

# Arguments

- `io`: the output io-like object, which is expected to be writable. Typical choices are: `stdout`,
  `IOBuffer` and `IOStream`. The default value is `stdout`.
- `src`: generic container object(e.g., `AbstractArray`) that contains the colorant sequence.
- `encoder::AbstractSixelEncoder`: the sixel encoder. Currently, only
  [`LibSixelEncoder`](@ref Sixel.LibSixel.LibSixelEncoder) is available.

!!! warning
    `sixel_encode` for tiny images (including vectors) is undefined behavior; sixel format requires
    at least `6` pixels in the row direction. For better visualization result, it's recommended to
    `repeat`(inner) the data so that its size is larger than `(6, 6)` in both dimensions.

!!! note
    For array with `ndims(src) >= 3`, it will be concatenated into a 2d array in row order before
    encoding.

# Parameters

- `transpose::Bool`: whether we need to permute the image's width and height dimension before encoding.
  The default value is `false`.

# Examples

```julia
using Sixel, TestImages

# 2d RGB image
img = testimage("mandril_color")
sixel_encode(img)

# 3d RGB image will be concatenated in row order.
img3d = testimage("mri-stack")
sixel_encode(img3d)

# You can also write the encoded sixel sequences into an `IOBuffer` or `IOStream`(via `open`)
io = IOBuffer()
sixel_encode(io, img)
println(String(take!(io)))

open("out.sixel", "w") do io
    sixel_encode(io, img)
end
```

# References

- [1] VT330/VT340 Programmer Reference Manual, Volume 1: Text Programming
- [2] VT330/VT340 Programmer Reference Manual, Volume 2: Graphics Programming
- [3] https://github.com/saitoha/libsixel
"""
function sixel_encode end

sixel_encode(img::AbstractArray, enc::SEC=default_encoder(img); kwargs...) =
    sixel_encode(stdout, img, enc; kwargs...)

sixel_encode(io::IO, img::AbstractVector, enc::SEC=default_encoder(img); kwargs...) =
    sixel_encode(io, reshape(img, :, 1), enc; kwargs...)

function sixel_encode(io::IO, img::AbstractArray, enc::SEC=default_encoder(img); transpose=false, kwargs...)
    # Conversion from OffsetArrays to Array is not very well supported, so we have to de-offset first.
    img = OffsetArrays.no_offset_view(img)

    # make sure it always tiles along row order
    @assert ndims(img) >= 3
    nrow = transpose ? prod(size(img)[3:end]) : 1
    sixel_encode(io, mosaicview(img; nrow=nrow, rowmajor=true), enc; transpose=transpose, kwargs...)
end

# This is expected to be the only caller of `enc(io, bytes)`; all other methods should eventually
# call this one.
function sixel_encode(
        io::IO,
        img::AbstractMatrix,
        enc::SEC=default_encoder(img);
        transpose=false, kwargs...)
    # Conversion from OffsetArrays to Array is not very well supported, so we have to de-offset first.
    img = OffsetArrays.no_offset_view(img)
    if enc isa LibSixelEncoder
        T = canonical_sixel_eltype(enc, eltype(img))
        AT = Array{T, ndims(img)}
        # libsixel is a C library and assumes row-major memory order, thus `collect` the data into
        # contiguous memeory layout already makes a transpose.
        bytes = transpose ? convert(AT, img) : convert(AT, PermutedDimsArray(img, (2, 1)))
    else
        throw(ArgumentError("Unsupported encoder type $(typeof(enc)). Please open an issue for this."))
    end
    enc(io, bytes; kwargs...)
    return nothing
end
