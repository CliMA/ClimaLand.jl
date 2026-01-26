"""
    jpeg_encode(filename::AbstractString, img; kwargs...) -> Int
    jpeg_encode(io::IO, img; kwargs...) -> Int
    jpeg_encode(img; kwargs...) -> Vector{UInt8}

Encode 2D image `img` as JPEG byte sequences and write to given I/O stream or file. The
return value is number of bytes. If output is not specified, the encoded result is stored
in memory as return value.

# Parameters

- `transpose::Bool`: whether we need to permute the image's width and height dimension
  before encoding. The default value is `false`.
- `quality::Int`: Constructs JPEG quantization tables appropriate for the indicated quality
  setting. The quality value is expressed on the 0..100 scale recommended by IJG. The
  default value is `92`. Pass `quality=nothing` to let libjpeg-turbo dynamicly guess a
  value.

!!! info "Custom compression parameters"
    JPEG has a large number of compression parameters that determine how the image is
    encoded. Most applications don't need or want to know about all these parameters. For
    more detailed information and explaination, please refer to the "Compression parameter
    selection" in [1]. Unsupported custom parameters might cause Julia segmentation fault.

# Examples

```jldoctest
julia> using JpegTurbo, TestImages

julia> img = testimage("cameraman");

julia> jpeg_encode("out.jpg", img) # write to file
51396

julia> buf = jpeg_encode(img); length(buf) # directly write to memory
51396
```

# References

- [1] [libjpeg API Documentation (libjpeg.txt)](https://raw.githubusercontent.com/libjpeg-turbo/libjpeg-turbo/main/libjpeg.txt)
"""
function jpeg_encode(img::AbstractMatrix{T}; transpose=false, kwargs...) where T<:Union{Real, Colorant}
    # quantilized into 8bit sequences first
    CT = T <: Colorant ? n0f8(eltype(img)) : Gray{N0f8}
    AT = Array{CT, ndims(img)}
    clamp01nan!(img)
    # jpegturbo is a C library and assumes row-major memory order, thus `collect` the data into
    # contiguous memeory layout already makes a transpose.
    img = transpose ? convert(AT, img) : convert(AT, PermutedDimsArray(img, (2, 1)))

    return _encode(img; kwargs...)
end

function jpeg_encode(filename::AbstractString, img; kwargs...)
    open(filename, "w") do io
        jpeg_encode(io, img; kwargs...)
    end
end
# TODO(johnnychen94): further improve the performance via asynchronously IO and buffer reuse.
jpeg_encode(io::IO, img; kwargs...) = write(io, jpeg_encode(img; kwargs...))


function _encode(
    img::Matrix{<:Colorant};
    colorspace::Union{Nothing,Type} = nothing,
    # ImageMagick: "the default is to use the estimated quality of your input image if it can be determined, otherwise 92."
    quality::Union{Nothing,Int} = 92,
    arith_code::Union{Nothing,Bool} = nothing,
    optimize_coding::Union{Nothing,Bool} = nothing,
    smoothing_factor::Union{Nothing,Int} = nothing,
    write_JFIF_header::Union{Nothing,Bool} = nothing,
    JFIF_version::Union{Nothing,VersionNumber} = nothing,
    density_unit::Union{Nothing,Int} = nothing,
    X_density::Union{Nothing,Int} = nothing,
    Y_density::Union{Nothing,Int} = nothing,
    write_Adobe_marker::Union{Nothing,Bool} = nothing,
    progressive_mode::Union{Nothing,Bool} = nothing
)
    if prod(size(img)) == 0
        throw(ArgumentError("empty image is not allowed"))
    end

    cinfo = LibJpeg.jpeg_compress_struct()
    cinfo_ref = Ref(cinfo)
    jerr = Ref{LibJpeg.jpeg_error_mgr}()
    cinfo.err = LibJpeg.jpeg_std_error(jerr)
    LibJpeg.jpeg_create_compress(cinfo_ref)

    # set input image information
    cinfo.image_width = size(img, 1)
    cinfo.image_height = size(img, 2)
    cinfo.input_components = jpeg_components(img)
    cinfo.in_color_space = jpeg_color_space(img)

    # set compression keywords
    # it's recommended to call `jpeg_set_defaults` first before setting custom parameters
    # as it's more likely to provide a working parameters and is more likely to be working
    # correctly in the future.
    LibJpeg.jpeg_set_defaults(cinfo_ref)
    isnothing(colorspace) || LibJpeg.jpeg_set_colorspace(cinfo_ref, jpeg_color_space(colorspace))
    isnothing(quality) || LibJpeg.jpeg_set_quality(cinfo_ref, quality, true)
    isnothing(arith_code) || (cinfo.arith_code = arith_code)
    isnothing(optimize_coding) || (cinfo.optimize_coding = optimize_coding)
    isnothing(smoothing_factor) || (cinfo.smoothing_factor = smoothing_factor)
    isnothing(write_JFIF_header) || (cinfo.write_JFIF_header = write_JFIF_header)
    if !isnothing(JFIF_version)
        cinfo.JFIF_major_version = UInt8(JFIF_version.major)
        cinfo.JFIF_minor_version = UInt8(JFIF_version.minor)
    end
    isnothing(density_unit) || (cinfo.density_unit = density_unit)
    isnothing(X_density) || (cinfo.X_density = X_density)
    isnothing(Y_density) || (cinfo.Y_density = Y_density)
    isnothing(write_Adobe_marker) || (cinfo.write_Adobe_marker = write_Adobe_marker)
    if !isnothing(progressive_mode) && progressive_mode
      cinfo.progressive_mode = true
      LibJpeg.jpeg_simple_progression(cinfo_ref)
    end

    # set destination
    # TODO(johnnychen94): allow pre-allocated buffer
    bufsize = Ref{Culong}(0)
    buf_ptr = Ref{Ptr{UInt8}}(C_NULL)
    LibJpeg.jpeg_mem_dest(cinfo_ref, buf_ptr, bufsize)

    # compression stage
    LibJpeg.jpeg_start_compress(cinfo_ref, true)
    row_stride = size(img, 1) * jpeg_components(img)
    row_pointer = Ref{Ptr{UInt8}}(0)
    while (cinfo.next_scanline < cinfo.image_height)
        row_pointer[] = pointer(img) + cinfo.next_scanline * row_stride
        LibJpeg.jpeg_write_scanlines(cinfo_ref, row_pointer, 1);
    end
    LibJpeg.jpeg_finish_compress(cinfo_ref)
    LibJpeg.jpeg_destroy_compress(cinfo_ref)

    return unsafe_wrap(Array, buf_ptr[], bufsize[]; own=true)
end
