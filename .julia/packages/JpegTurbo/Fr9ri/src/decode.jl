"""
    jpeg_decode([T,] filename::AbstractString; kwargs...) -> Matrix{T}
    jpeg_decode([T,] io::IO; kwargs...) -> Matrix{T}
    jpeg_decode([T,] data::Vector{UInt8}; kwargs...) -> Matrix{T}

Decode the JPEG image as colorant matrix. The source data can be either a filename, an IO
, or an in-memory byte sequence.

# parameters

- `transpose::Bool`: whether we need to permute the image's width and height dimension
  before encoding. The default value is `false`.
- `scale_ratio::Real`: scale the image by ratio `scale_ratio` in `M/8` with `M ∈ 1:16`. The
  default value is `1`. For values are not in the range, they will be mapped to the nearest
  value, e.g., `0.3 => 2/8` and `0.35 => 3/8`. `scale_ratio` and `preferred_size` may not be
  used together.
- `preferred_size::Tuple`: infer the minimal `scale_ratio` that `all(size(out) .>=
  preferred_size))` holds. It can optionally be `(op, preferred_size)` format, with compare
  operation `op` be one of `>`, `>=`, `<` or `<=`. If `op in (>=, >)` then it gets the
  minimal `scale_ratio`, otherwise it gets the maximum `scale_ratio` for `op in (<=, <)`.
  `scale_ratio` and `preferred_size` may not be used together. The `preferred_size`
  dimensions are not affected by keyword `transpose`.

# Examples

```jldoctest
julia> using JpegTurbo, TestImages, ImageCore

julia> filename = testimage("earth", download_only=true);

julia> img = jpeg_decode(filename); summary(img)
"3002×3000 Matrix{RGB{N0f8}}"

julia> img = jpeg_decode(Gray, filename; scale_ratio=0.25); summary(img)
"751×750 Matrix{Gray{N0f8}}"
```

For image preview and similar purposes, `T` and `scale_ratio` are useful parameters to
accelerate the JPEG decoding process. For color JPEG image, `jpeg_decode(Gray, filename)` is
faster than `jpeg_decode(filename)` since the color components need not be processed.
Smaller `scale_ratio` permits significantly faster decoding since fewer pixels need be
processed and a simpler IDCT method can be used.

```julia
using BenchmarkTools, TestImages, JpegTurbo
filename = testimage("earth", download_only=true)
# full decompression
@btime jpeg_decode(filename); # 224.760 ms (7 allocations: 51.54 MiB)
# only decompress luminance component
@btime jpeg_decode(Gray, filename); # 91.157 ms (6 allocations: 17.18 MiB)
# process only a few pixels
@btime jpeg_decode(filename; scale_ratio=0.25); # 77.254 ms (8 allocations: 3.23 MiB)
# process only a few pixels for luminance component
@btime jpeg_decode(Gray, filename; scale_ratio=0.25); # 63.119 ms (6 allocations: 1.08 MiB)
```

"""
function jpeg_decode(
        ::Type{CT},
        data::Vector{UInt8};
        transpose=false,
        scale_ratio::Union{Nothing,Real}=nothing,
        preferred_size::Union{Nothing,Tuple}=nothing) where CT<:Colorant
    _jpeg_check_bytes(data)
    out_CT, jpeg_cls = _jpeg_out_color_space(CT)

    cinfo_ref = Ref(LibJpeg.jpeg_decompress_struct())
    jerr = Ref{LibJpeg.jpeg_error_mgr}()
    try
        cinfo = cinfo_ref[]
        cinfo.err = LibJpeg.jpeg_std_error(jerr)
        LibJpeg.jpeg_create_decompress(cinfo_ref)
        LibJpeg.jpeg_mem_src(cinfo_ref, data, length(data))
        LibJpeg.jpeg_read_header(cinfo_ref, true)

        # set decompression parameters, if given
        if !isnothing(preferred_size) && !transpose
            if preferred_size[1] isa Function
                preferred_size = (preferred_size[1], reverse(preferred_size[2]))
            else
                preferred_size = reverse(preferred_size)
            end
        end
        r = _cal_scale_ratio(scale_ratio, preferred_size, cinfo)
        cinfo.scale_num, cinfo.scale_denom = r.num, r.den
        cinfo.out_color_space = jpeg_cls

        # eagerly calculate dimension information so that `output_XXX` fields are valid.
        LibJpeg.jpeg_calc_output_dimensions(cinfo_ref)
        out_size = (Int(cinfo.output_width), Int(cinfo.output_height))
        if !all(x -> x<=65535, out_size)
            error("Suspicious inferred image size $out_size: each dimension is expected to have at most 65535 pixels.")
        end
        out_ndims = Int(cinfo.output_components)
        @assert out_ndims == length(out_CT) "Suspicous output color space: $cinfo.out_color_space"

        out = Matrix{out_CT}(undef, out_size)
        _jpeg_decode!(out, cinfo_ref)

        if out_CT <: CT
            return transpose ? out : permutedims(out, (2, 1))
        else
            return transpose ? CT.(out) : CT.(PermutedDimsArray(out, (2, 1)))
        end
    finally
        LibJpeg.jpeg_destroy_decompress(cinfo_ref)
    end
end
jpeg_decode(data; kwargs...) = jpeg_decode(_default_out_color_space(data), data; kwargs...)

# TODO(johnnychen94): support Progressive JPEG
# TODO(johnnychen94): support partial decoding
function jpeg_decode(::Type{CT}, filename::AbstractString; kwargs...) where CT<:Colorant
    open(filename, "r") do io
        jpeg_decode(CT, io; kwargs...)
    end
end
jpeg_decode(filename::AbstractString; kwargs...) = jpeg_decode(read(filename); kwargs...)

jpeg_decode(io::IO; kwargs...) = jpeg_decode(read(io); kwargs...)
jpeg_decode(::Type{CT}, io::IO; kwargs...) where CT<:Colorant = jpeg_decode(CT, read(io); kwargs...)

function _jpeg_decode!(out::Matrix{<:Colorant}, cinfo_ref::Ref{LibJpeg.jpeg_decompress_struct})
    row_stride = size(out, 1) * length(eltype(out))

    # get a pointer to `out` as if it's a UInt8 array
    out_ptr_ref = Ref(Ptr{UInt8}(pointer(out)))

    cinfo = cinfo_ref[]
    LibJpeg.jpeg_start_decompress(cinfo_ref)
    while cinfo.output_scanline < cinfo.output_height
        GC.@preserve out LibJpeg.jpeg_read_scanlines(cinfo_ref, out_ptr_ref, 1)
        out_ptr_ref[] += row_stride # move the pointer one row
    end
    LibJpeg.jpeg_finish_decompress(cinfo_ref)

    return out
end


_jpeg_decode!(out::Matrix{<:Colorant}, cinfo::LibJpeg.jpeg_decompress_struct) = _jpeg_decode!(out, Ref(cinfo))

# libjpeg-turbo only supports ratio M/8 with M from 1 to 16
const _allowed_scale_ratios = ntuple(i->i//8, 16)
function _cal_scale_ratio(::Real, ::Tuple, cinfo)
    throw(ArgumentError("Keywords `aspect_ratio` and `preferred_size` are exclusive to each other: you should only specify one of them."))
end
_cal_scale_ratio(::Nothing, ::Nothing, cinfo) = 1//1
_cal_scale_ratio(r::Real, ::Nothing, cinfo) = _allowed_scale_ratios[findmin(x->abs(x-r), _allowed_scale_ratios)[2]]
function _cal_scale_ratio(::Nothing, preferred_size::Tuple, cinfo)
    cinfo.scale_num, cinfo.scale_denom = 1, 1
    LibJpeg.jpeg_calc_output_dimensions(Ref(cinfo))
    out_size = (Int(cinfo.output_width), Int(cinfo.output_height))
    op, preferred_size = if preferred_size[1] isa Function
        op = preferred_size[1]
        op in (>, >=, <, <=) || throw(ArgumentError("the compare operation must be one of `>`, `>=`, `<` and `<=`."))
        op, preferred_size[2]
    else
        >=, preferred_size
    end
    if op in (>, >=)
        idx = findfirst(x->all(op.(x .* out_size, preferred_size)), _allowed_scale_ratios)
        if isnothing(idx)
            @warn "Failed to infer appropriate scale ratio, use `scale_ratio=2` instead." actual_size=out_size preferred_size
            idx = length(_allowed_scale_ratios)
        end
    elseif op in (<, <=)
        idx = findlast(x->all(op.(x .* out_size, preferred_size)), _allowed_scale_ratios)
        if isnothing(idx)
            @warn "Failed to infer appropriate scale ratio, use `scale_ratio=1/8` instead." actual_size=out_size preferred_size
            idx = 1
        end
    end
    return _allowed_scale_ratios[idx]
end

function _default_out_color_space(data::Vector{UInt8})
    _jpeg_check_bytes(data)
    cinfo_ref = Ref(LibJpeg.jpeg_decompress_struct())
    try
        jerr = Ref{LibJpeg.jpeg_error_mgr}()
        cinfo_ref[].err = LibJpeg.jpeg_std_error(jerr)
        LibJpeg.jpeg_create_decompress(cinfo_ref)
        LibJpeg.jpeg_mem_src(cinfo_ref, data, length(data))
        LibJpeg.jpeg_read_header(cinfo_ref, true)
        LibJpeg.jpeg_calc_output_dimensions(cinfo_ref)
        return jpeg_color_space(cinfo_ref[].out_color_space)
    finally
        LibJpeg.jpeg_destroy_decompress(cinfo_ref)
    end
end

function _jpeg_out_color_space(::Type{CT}) where CT
    try
        n0f8(CT), jpeg_color_space(n0f8(CT))
    catch e
        @debug "Unsupported libjpeg-turbo color space, fallback to RGB{N0f8}" e
        RGB{N0f8}, jpeg_color_space(RGB{N0f8})
    end
end

# provides some basic integrity check
# TODO(johnnychen94): redirect libjpeg-turbo error to julia
_jpeg_check_bytes(filename::AbstractString) = open(_jpeg_check_bytes, filename, "r")
function _jpeg_check_bytes(io::IO)
    seekend(io)
    nbytes = position(io)
    nbytes > 623 || throw(ArgumentError("Invalid number of bytes."))

    buf = UInt8[]
    seekstart(io)
    readbytes!(io, buf, 623)
    seek(io, nbytes-2)
    append!(buf, read(io, 2))
    return _jpeg_check_bytes(buf)
end
function _jpeg_check_bytes(data::Vector{UInt8})
    length(data) > 623 || throw(ArgumentError("Invalid number of bytes."))
    data[1:2] == [0xff, 0xd8] || throw(ArgumentError("Invalid JPEG byte sequence."))
    data[end-1:end] == [0xff, 0xd9] || @warn "Premature end of JPEG byte sequence."
    return true
end
