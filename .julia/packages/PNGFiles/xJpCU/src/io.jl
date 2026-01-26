"""
    load(fpath::String;
         gamma::Union{Nothing,Float64}=nothing, expand_paletted::Bool=false, background=false)
    load(s::IO;
         gamma::Union{Nothing,Float64}=nothing, expand_paletted::Bool=false, background=false)

Read a PNG image as a julia `Array`.

# Arguments

- `fpath`: the path to the png file to be loaded.
- `s`: an IO stream containing the png file.

# Keywords

- `gamma`: the end-to-end coefficient for gamma correction, can be used to override the
    automatic gamma correction, a value of `1.0` means no gamma correction. By default,
    gamma correction is applied according to a `gAMA` or `sRGB` chunks, if present;
    if neither chunk is present a value of `0.45455` for the file gamma is assumed.
    Screen gamma is currently always assumed to be `2.2`.
- `background`: can be used to set a background color for transparent images.
    Accepted values are `false` for no background, `true` for loading the image
    with a solid background if `bKGD` chunk is set, as well as user-provided backgrounds
    of type `UInt8`, `Gray` or `RGB`, which represent a palette index, gray and true color
    backgrounds respectively.
- `expand_paletted`: when reading in simple paletted images, i.e. having a `PLTE` chunk and
   an 8 bit depth, the image will be represented as an `IndirectArray` with `OffsetArray`
   `values` field. To always get back a plain `Matrix` of colorants, use `expand_paletted=true`.

# Returns

- `Matrix{Gray{N0f8}}`: for grayscale images (without transparency) and bit depth lower or equal to 8.
- `Matrix{Gray{N0f16}}`: for grayscale images (without transparency) and bit depth equal to 16.
- `Matrix{GrayA{N0f8}}`: for grayscale images (with transparency) and bit depth lower or equal to 8.
- `Matrix{GrayA{N0f16}}`: for grayscale images (with transparency) and bit depth equal to 16.
- `Matrix{RGB{N0f8}}`: for true color images (without transparency) and bit depth equal to 8.
    Also for palleted images without transparency when `expand_paletted=true`.
- `Matrix{RGB{N0f16}}`: for true color images (without transparency) and bit depth equal to 16.
- `Matrix{RGBA{N0f8}}`: for true color images (with transparency) and bit depth equal to 8.
    Also for palleted images without transparency when `expand_paletted=true`.
- `Matrix{RGBA{N0f16}}`: for true color images (with transparency) and bit depth equal to 16.
- `IndirectArray{RGB{N0f8}, 2, UInt8, Matrix{UInt8}, OffsetVector{RGB{N0f8}, Vector{RGB{N0f8}}}}`:
    for palleted images with transparency when `expand_paletted=false`.
- `IndirectArray{RGBA{N0f8}, 2, UInt8, Matrix{UInt8}, OffsetVector{RGBA{N0f8}, Vector{RGBA{N0f8}}}}`:
    for palleted images with transparency when `expand_paletted=false`.
"""
function load(fpath::String; gamma::Union{Nothing,Float64}=nothing, expand_paletted::Bool=false, background=false)
    fp = open_png(fpath)
    png_ptr = create_read_struct()
    @debug "Load PNG File:" fpath png_ptr
    info_ptr = create_info_struct(png_ptr)
    png_init_io(png_ptr, fp)
    png_set_sig_bytes(png_ptr, PNG_BYTES_TO_CHECK)
    out = _load(png_ptr, info_ptr, gamma=gamma, expand_paletted=expand_paletted, background=background)
    close_png(fp)
    return out
end

maybe_lock(f, io::IO) = lock(f, io)
# IOStream doesn't support locking...
maybe_lock(f, io::IOStream) = f()
maybe_lock(f, io::IOBuffer) = f()
maybe_lock(f, io::Base64EncodePipe) = f()
maybe_lock(f, io::IOContext) = maybe_lock(f, io.io)

if VERSION < v"1.6.0-DEV.1652"
    _iswritable(x) = Base.iswritable(x)
    _iswritable(pipe::Base64EncodePipe) = Base.iswritable(pipe.io)
else
    _iswritable(x) = iswritable(x)
end

function load(s::IO; gamma::Union{Nothing,Float64}=nothing, expand_paletted::Bool=false, background=false)
    isreadable(s) || throw(ArgumentError("read failed, IOStream is not readable"))
    Base.eof(s) && throw(EOFError())

    png_ptr = create_read_struct()
    info_ptr = create_info_struct(png_ptr)

    maybe_lock(s) do
        if s isa IOBuffer
            png_set_read_fn(png_ptr, pointer_from_objref(s), readcallback_iobuffer_c[])
        else
            png_set_read_fn(png_ptr, s.handle, readcallback_c[])
        end
        # https://stackoverflow.com/questions/22564718/libpng-error-png-unsigned-integer-out-of-range
        png_set_sig_bytes(png_ptr, 0)
        return _load(png_ptr, info_ptr, gamma=gamma, expand_paletted=expand_paletted, background=background)
    end
end

function _readcallback(png_ptr::png_structp, data::png_bytep, length::png_size_t)::Cvoid
    a = png_get_io_ptr(png_ptr)
    ccall(:ios_readall, Csize_t, (Ptr{Cvoid}, Ptr{Cvoid}, Csize_t), a, data, length)
    return
end

function _readcallback_iobuffer(png_ptr::png_structp, data::png_bytep, length::png_size_t)::Cvoid
    a = png_get_io_ptr(png_ptr)
    io = unsafe_pointer_to_objref(a)
    unsafe_read(io, Ptr{UInt8}(data), length)
    return
end

function _load(png_ptr, info_ptr; gamma::Union{Nothing,Float64}=nothing, expand_paletted::Bool=false, background=false)
    png_read_info(png_ptr, info_ptr)
    flags = png_get_valid(
        png_ptr, info_ptr,
        PNG_INFO_sRGB | PNG_INFO_gAMA | PNG_INFO_tRNS | PNG_INFO_cHRM | PNG_INFO_bKGD | PNG_INFO_PLTE
    )
    valid_sRGB = (flags & PNG_INFO_sRGB) != 0
    valid_gAMA = (flags & PNG_INFO_gAMA) != 0
    valid_tRNS = (flags & PNG_INFO_tRNS) != 0
    valid_cHRM = (flags & PNG_INFO_cHRM) != 0
    valid_bKGD = (flags & PNG_INFO_bKGD) != 0
    valid_PLTE = (flags & PNG_INFO_PLTE) != 0

    width = png_get_image_width(png_ptr, info_ptr)
    height = png_get_image_height(png_ptr, info_ptr)
    color_type = color_type_orig = png_get_color_type(png_ptr, info_ptr)
    bit_depth = bit_depth_orig = png_get_bit_depth(png_ptr, info_ptr)
    num_channels = png_get_channels(png_ptr, info_ptr)
    interlace_type = png_get_interlace_type(png_ptr, info_ptr)
    background_color = nothing
    is_transparent = valid_tRNS | (color_type & PNG_COLOR_MASK_ALPHA) != 0
    is_gray = (color_type == PNG_COLOR_TYPE_GRAY_ALPHA) | (color_type == PNG_COLOR_TYPE_GRAY)

    read_as_paletted = !expand_paletted && color_type == PNG_COLOR_TYPE_PALETTE && bit_depth == 8 && valid_PLTE
    read_with_background = _check_background_load(background, is_gray, valid_PLTE, read_as_paletted)

    screen_gamma = PNG_DEFAULT_sRGB
    image_gamma = Ref{Cdouble}(-1.0)
    intent = Ref{Cint}(-1)
    if isnothing(gamma)
        if valid_sRGB
            if png_get_sRGB(png_ptr, info_ptr, intent) != 0
                png_set_gamma(png_ptr, screen_gamma, PNG_DEFAULT_sRGB)
            else
                if png_get_gAMA(png_ptr, info_ptr, image_gamma) == 0
                    image_gamma[] = 0.45455
                end
                png_set_gamma(png_ptr, screen_gamma, image_gamma[])
            end
        elseif valid_gAMA
            if png_get_gAMA(png_ptr, info_ptr, image_gamma) == 0
                image_gamma[] = 0.45455
            end
            png_set_gamma(png_ptr, screen_gamma, image_gamma[])
        end
    elseif gamma != 1
        image_gamma[] = 2.2inv(gamma)
        png_set_gamma(png_ptr, screen_gamma, image_gamma[])
    end

    if !read_as_paletted
        if color_type == PNG_COLOR_TYPE_PALETTE
            png_set_palette_to_rgb(png_ptr)
            color_type = PNG_COLOR_TYPE_RGB
        end

        if color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8
            png_set_expand_gray_1_2_4_to_8(png_ptr)
            png_set_packing(png_ptr)
            bit_depth = 8
        end

        if valid_tRNS
            png_set_tRNS_to_alpha(png_ptr)
            if color_type == PNG_COLOR_TYPE_GRAY || color_type == PNG_COLOR_TYPE_RGB
                color_type |= PNG_COLOR_MASK_ALPHA
            end
        end

        if read_with_background & valid_bKGD & is_transparent
            png_set_filler(png_ptr, 0xff, PNG_FILLER_AFTER)
        end
        buffer_eltype = _buffer_color_type(color_type, bit_depth)

        bit_depth == 16 && png_set_swap(png_ptr)
    else
        buffer_eltype = UInt8
    end

    if read_with_background
        background_color = process_background(png_ptr, info_ptr, _adjust_background_bitdepth(background, bit_depth))
        is_transparent || @warn("Background color for non-transparent image: $(background_color)")
    end

    n_passes = png_set_interlace_handling(png_ptr)
    png_read_update_info(png_ptr, info_ptr)

    # Gamma correction is applied to a palette after `png_read_update_info` is called
    if read_as_paletted
        palette_length = Ref{Cint}(0)
        palette_buffer = Ref{Ptr{PNGFiles.png_color_struct}}(C_NULL)
        png_get_PLTE(png_ptr, info_ptr, palette_buffer, palette_length)
        palette = Vector{RGB{N0f8}}(undef, palette_length[])
        GC.@preserve palette unsafe_copyto!(pointer(palette), Ptr{RGB{N0f8}}(palette_buffer[]), length(palette))
        if valid_tRNS
            alpha_buffer = Ref{Ptr{UInt8}}(C_NULL)
            alphas_cnt = Ref{Cint}(0)
            png_get_tRNS(png_ptr, info_ptr, alpha_buffer, alphas_cnt, C_NULL)
            alpha = fill(one(N0f8), palette_length[]) # palette entries are opaque by default
            @assert palette_length[] >= alphas_cnt[]
            GC.@preserve alpha unsafe_copyto!(pointer(alpha), Ptr{N0f8}(alpha_buffer[]), min(palette_length[], alphas_cnt[]))
            palette = map(RGBA, palette, alpha)
        end
        buffer_eltype = UInt8
    end

    # We transpose to work around libpng expecting row-major arrays
    buffer = Matrix{buffer_eltype}(undef, width, height)

    @debug(
        "Read PNG info:",
        PNG_HEADER_VERSION_STRING,
        png_ptr,
        height,
        width,
        color_type_orig,
        color_type,
        bit_depth_orig,
        bit_depth,
        num_channels,
        interlace_type,
        gamma,
        image_gamma[],
        screen_gamma,
        background,
        intent[],
        read_as_paletted,
        n_passes,
        buffer_eltype,
        valid_sRGB,
        valid_gAMA,
        valid_tRNS,
        valid_cHRM,
        valid_bKGD,
        valid_PLTE,
    )

    GC.@preserve buffer begin
        buffer = Base.invokelatest(_load!, buffer, png_ptr, info_ptr)
    end

    if expand_paletted || color_type != PNG_COLOR_TYPE_PALETTE
        return buffer
    else
        # We got 0-based indices back from libpng and converting to 1-based could overflow UInt8.
        # Using UInt16 for index would cost us large part of the savings provided by IndirectArray.
        return IndirectArray(buffer, OffsetArray(palette, -1))
    end
end

function _load!(buffer::Matrix{T}, png_ptr, info_ptr) where T    # separate to support precompilation of permutedims
    png_read_image(png_ptr, map(a -> Ptr{UInt8}(pointer(a)), eachcol(buffer)))
    png_read_end(png_ptr, info_ptr)
    png_destroy_read_struct(Ref{Ptr{Cvoid}}(png_ptr), Ref{Ptr{Cvoid}}(info_ptr), C_NULL)
    return permutedims(buffer, (2, 1))
end

function _buffer_color_type(color_type, bit_depth)
    bit_depth = Int(bit_depth)
    if color_type == PNG_COLOR_TYPE_GRAY
        colors_type = Gray{bit_depth > 8 ? Normed{UInt16,bit_depth} : Normed{UInt8,bit_depth}}
    elseif color_type == PNG_COLOR_TYPE_PALETTE
        colors_type = RGB{bit_depth == 16 ? N0f16 : N0f8}
    elseif color_type == PNG_COLOR_TYPE_RGB
        colors_type = RGB{bit_depth == 16 ? N0f16 : N0f8}
    elseif color_type == PNG_COLOR_TYPE_RGB_ALPHA
        colors_type = RGBA{bit_depth == 16 ? N0f16 : N0f8}
    elseif color_type == PNG_COLOR_TYPE_GRAY_ALPHA
        colors_type = GrayA{bit_depth > 8 ? Normed{UInt16,bit_depth} : Normed{UInt8,bit_depth}}
    else
        throw(error("Unknown color type: $color_type"))
    end
    return colors_type
end


### Write ##########################################################################################

_dpi_to_ppm(::Nothing) = nothing
function _dpi_to_ppm(dpi::Real)
    @assert dpi > 0
    pixels_per_meter = round(UInt32, dpi / 0.0254)
    return (pixels_per_meter, pixels_per_meter)
end
function _dpi_to_ppm(dpi::Tuple{<:Real,<:Real})
    @assert dpi[1] > 0
    @assert dpi[2] > 0
    return round.(UInt32, dpi ./ 0.0254)
end

const SupportedPaletteColor = Union{
    AbstractRGB{<:Union{N0f8,AbstractFloat}},
    TransparentRGB{T,<:Union{N0f8,AbstractFloat}} where T,
}
"""
    save(fpath::String, image::AbstractArray;
         compression_level::Int=0, compression_strategy::Int=3, filters::Int=4)
    save(s::IO, image::AbstractArray;
         compression_level::Int=0, compression_strategy::Int=3, filters::Int=4)

Write out a julia `Array` as a PNG image.

# Arguments
- `fpath`: the filesystem path where to store the resulting PNG image.
- `s`: the IO stream to which the resulting PNG image should be written to.
- `image`: the julia `Array` which should be encoded as a PNG image. The type, eltype, and
    dimensionality jointly determine the way the PNG file is encoded. The following table
    roughly summarises the encoding rules:
| **dims**    | **type**                                       | **eltype** (`T`)                                     | **PNG color type**          | **PNG bit depth**                                |
|:------------|:-----------------------------------------------|:-----------------------------------------------------|:----------------------------|-------------------------------------------------:|
| `[h, w]`    | `AbstracArray{T,2}`                            | `Union{AbstractFloat, Unsigned, Normed, Gray}`       | `PNG_COLOR_TYPE_GRAY`       | 8 by default, 16 when `N0f16` or `N4f12` is used |
| `[h, w, 1]` | `AbstracArray{T,3}`                            | `Union{AbstractFloat, Unsigned, Normed, Gray}`       | `PNG_COLOR_TYPE_GRAY`       | 8 by default, 16 when `N0f16` or `N4f12` is used |
| `[h, w, 2]` | `AbstracArray{T,3}`                            | `Union{AbstractFloat, Unsigned, Normed, GrayA}`      | `PNG_COLOR_TYPE_GRAY_ALPHA` | 8 by default, 16 when `N0f16` or `N4f12` is used |
| `[h, w, 3]` | `AbstracArray{T,3}`                            | `Union{AbstractFloat, Unsigned, Normed, RGB, BGR}`   | `PNG_COLOR_TYPE_RGB`        | 8 by default, 16 when `N0f16` or `N4f12` is used |
| `[h, w, 4]` | `AbstracArray{T,3}`                            | `Union{AbstractFloat, Unsigned, Normed, RGBA, ARGB}` | `PNG_COLOR_TYPE_RGB_ALPHA`  | 8 by default, 16 when `N0f16` or `N4f12` is used |
| `[h, w]`    | `IndirectArray{T, 2, S, Matrix{S}, Vector{T}}` | `Union{RGB, BGR, RGBA, ARGB}`                        | `PNG_COLOR_TYPE_PALETTE`    |                                                8 |

``Note``: Only `UInt8` and `UInt16` are supported for `Unsigned`, `Float32` and `Float64` for `AbstractFloat`;
    `N0f8`, `N0f16` and `N4f12` for `Normed`. These types are also the only valid types for colorants.

``Note``: `Bool`s and `BitArray` are also supported, but the bit depth of the image will still be 8.

``Note``: When `image` is an `IndirectArray` with up to 256 unique `RGB` colors, the result is encoded as a paletted image.
    Palletes with 16 bit depths are not supported. The palette (`values` field of the `IndirectArray`) could also be represented with an `OffsetArray`.

# Keywords

- `compression_level`: `0` (`Z_NO_COMPRESSION`), `1` (`Z_BEST_SPEED`), ..., `9` (`Z_BEST_COMPRESSION`)
- `compression_strategy`: 0 (`Z_DEFAULT_STRATEGY`), 1 (`Z_FILTERED`), 2 (`Z_HUFFMAN_ONLY`),
    3 (`Z_RLE`), 4 (`Z_FIXED`)
- `filters`: specify a type of preprocessing applied to each row which can increase its compressability.
    Valid values are 0 (`None`), 1 (`Sub`), 2 (`Up`), 3 (`Average`), 4 (`Paeth`).
- `file_gamma`: the value governing the gamma encoding of the image. When `nothing`,
    the image stored as `sRGB`, otherwise the gamma value provided will populate a `gAMA`
    chunk of the image.
- `background`: optional background color to be stored in the `bKGD` chunk. Only meaningful for transparent images.
    Valid values are `nothing` for no background, `UInt8` as a palette index for palleted images,
    `Gray` for grayscale images and `RGB` for true color images.
- `dpi`: stores the pixel density given in dots per inch into the `pHYs` chunk as pixels per meter.
    The density is given as `dpi` for ease of use as pixels per meter is an uncommon format. If set to `nothing`,
    no pixel density is written. If set to a 2-element tuple, a different density is written for x and y, respectively.

# Returns
- `nothing`
"""
function save(
    fpath::String,
    image::S;
    compression_level::Integer = Z_BEST_SPEED,
    compression_strategy::Integer = Z_RLE,
    filters::Integer = Int(PNG_FILTER_PAETH),
    file_gamma::Union{Nothing,Float64} = nothing,
    background::Union{Nothing,UInt8,AbstractGray,AbstractRGB} = nothing,
    dpi::Union{Nothing,Real,Tuple{<:Real,<:Real}} = nothing,
) where {
    T,
    S<:Union{AbstractMatrix{T},AbstractArray{T,3}}
}
    @assert Z_DEFAULT_STRATEGY <= compression_strategy <= Z_FIXED
    @assert Z_NO_COMPRESSION <= compression_level <= Z_BEST_COMPRESSION
    @assert 2 <= ndims(image) <= 3
    @assert size(image, 3) <= 4

    fp = ccall(:fopen, Ptr{Cvoid}, (Cstring, Cstring), fpath, "wb")
    fp == C_NULL && error("Could not open $(fpath) for writing")

    png_ptr = create_write_struct()
    @debug "Save PNG File:" fpath png_ptr
    info_ptr = create_info_struct(png_ptr)

    png_init_io(png_ptr, fp)

    _save(png_ptr, info_ptr, image,
        compression_level=compression_level,
        compression_strategy=compression_strategy,
        filters=filters,
        file_gamma=file_gamma,
        background=background,
        pixels_per_meter = _dpi_to_ppm(dpi),
    )

    close_png(fp)
    return
end
function save(
    s::IO,
    image::S;
    compression_level::Integer = Z_BEST_SPEED,
    compression_strategy::Integer = Z_RLE,
    filters::Integer = Int(PNG_FILTER_PAETH),
    file_gamma::Union{Nothing,Float64} = nothing,
    background::Union{Nothing,UInt8,AbstractGray,AbstractRGB} = nothing,
    dpi::Union{Nothing,Real,Tuple{<:Real,<:Real}} = nothing,
) where {
    S<:Union{AbstractMatrix,AbstractArray{<:Any,3}}
}
    @assert Z_DEFAULT_STRATEGY <= compression_strategy <= Z_FIXED
    @assert Z_NO_COMPRESSION <= compression_level <= Z_BEST_COMPRESSION
    @assert 2 <= ndims(image) <= 3
    @assert size(image, 3) <= 4
    _iswritable(s) || throw(ArgumentError("write failed, IOStream is not writeable"))

    png_ptr = create_write_struct()
    info_ptr = create_info_struct(png_ptr)
    maybe_lock(s) do
        r = Ref{Any}(s)
        GC.@preserve r begin
            png_set_write_fn(png_ptr, r, writecallback_c[], C_NULL)

            _save(png_ptr, info_ptr, image,
                compression_level=compression_level,
                compression_strategy=compression_strategy,
                filters=filters,
                file_gamma=file_gamma,
                background=background,
                pixels_per_meter = _dpi_to_ppm(dpi),
            )
        end
    end
    return
end

function _save(png_ptr, info_ptr, image::S;
    compression_level::Integer = Z_BEST_SPEED,
    compression_strategy::Integer = Z_RLE,
    filters::Integer = Int(PNG_FILTER_PAETH),
    file_gamma::Union{Nothing,Float64} = nothing,
    background::Union{Nothing,UInt8,AbstractGray,AbstractRGB} = nothing,
    pixels_per_meter::Union{Nothing,Tuple{UInt32,UInt32}} = nothing,
) where {
    T,
    S<:Union{AbstractMatrix{T},AbstractArray{T,3}}
}
    image = _enforce_dense_cpu_array(image)
    height, width = size(image)[1:2]
    bit_depth = _get_bit_depth(image)
    color_type = _get_color_type(image)
    is_transparent = (color_type & PNG_COLOR_MASK_ALPHA) != 0
    approx_bytes = round(Int, (height + 1) * width * bit_depth / 8 * (((color_type | PNG_COLOR_MASK_COLOR > 0) ? 3 : 1) + (color_type | PNG_COLOR_MASK_ALPHA > 0)))

    png_set_filter(png_ptr, PNG_FILTER_TYPE_BASE, UInt32(filters))
    png_set_compression_level(png_ptr, compression_level)
    png_set_compression_strategy(png_ptr, compression_strategy)
    png_set_compression_window_bits(png_ptr, min(15, max(8, _nextpow2exp(approx_bytes))))
    if pixels_per_meter !== nothing
        png_set_pHYs(png_ptr, info_ptr, pixels_per_meter..., PNG_RESOLUTION_METER)
    end

    if color_type == PNG_COLOR_TYPE_PALETTE
        # TODO: 1, 2, 4 bit-depth indices for palleted
        _png_check_paletted(image)
        palette = image.values
        is_transparent = eltype(palette) <: TransparentRGB
        color_count = length(palette)
        if is_transparent
            alphas = _palette_alpha(palette)
            alpha_count = color_count
            while (alpha_count > 0) && (alphas[alpha_count] == 1)
                alpha_count -= 1
            end
            png_set_PLTE(png_ptr, info_ptr, _standardize_palette(color.(palette)), color_count)
            png_set_tRNS(png_ptr, info_ptr, alphas, alpha_count, C_NULL)
        else
            png_set_PLTE(png_ptr, info_ptr, _standardize_palette(palette), color_count)
        end
    else
        palette = nothing
        image_eltype = eltype(image)
        if (image_eltype <: BGR || image_eltype <: BGRA || image_eltype <: ABGR || image_eltype <: ARGB32)
            png_set_bgr(png_ptr)
        end

        if (image_eltype <: ABGR || image_eltype <: ARGB)
            png_set_swap_alpha(png_ptr)
        end

        if color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8
            png_set_packing(png_ptr)
            bit_depth = 8  # TODO: support 1, 2, 4 bit-depth gray images
        end
    end
    if !(background === nothing)
        _check_background_save(background, color_type)
        _png_set_bKGD(png_ptr, info_ptr, _adjust_background_bitdepth(background, bit_depth))
    end

    if file_gamma === nothing
        # gAMA and cHRM chunks should be always present for compatibility with older systems
        png_set_sRGB_gAMA_and_cHRM(png_ptr, info_ptr, PNG_sRGB_INTENT_PERCEPTUAL)
    else
        png_set_gAMA(png_ptr, info_ptr, file_gamma)
    end

    @debug(
        "Write PNG info:",
        png_ptr,
        height,
        width,
        bit_depth,
        color_type,
        filters,
        compression_level,
        compression_strategy,
        palette,
        typeof(image)
    )

    # TODO: on error this throws an abort signal because we currently don't handle `longjmp`
    png_set_IHDR(
        png_ptr,
        info_ptr,
        width,
        height,
        bit_depth,
        color_type,
        PNG_INTERLACE_NONE,
        PNG_COMPRESSION_TYPE_BASE,
        PNG_FILTER_TYPE_BASE,
    )

    png_write_info(png_ptr, info_ptr)
    # Handles endianness for 16 bit, must be set after png_write_info
    bit_depth == 16 && png_set_swap(png_ptr)

    # We transpose to work around libpng expecting row-major arrays
    _write_image(permutedims(_prepare_buffer(image), (2, 1)), png_ptr, info_ptr)

    png_destroy_write_struct(Ref(png_ptr), Ref(info_ptr))
end

function _writecallback(png_ptr::png_structp, data::png_bytep, length::png_size_t)::Csize_t
    a = png_get_io_ptr(png_ptr)
    io = unsafe_load(Ptr{Any}(a))
    unsafe_write(io, data, length)
end

function _write_image(buf::AbstractArray{T,2}, png_ptr::Ptr{Cvoid}, info_ptr::Ptr{Cvoid}) where {T}
    GC.@preserve buf begin
        ccall(
            (:png_write_image, libpng),
            Cvoid,
            (Ptr{Cvoid}, Ptr{Ptr{UInt8}}),
            png_ptr,
            map(a -> Ptr{UInt8}(pointer(a)), eachcol(buf)),
        )
    end
    png_write_end(png_ptr, info_ptr)
end

function _png_check_paletted(image)
    palette = image.values
    color_count = length(palette)
    color_type = eltype(palette)

    ndims(image) != 2 && throw(ArgumentError("Only 2D `IndirectArrays` are supported"))
    color_count > 256 && throw(ArgumentError("Maximum size of `image.velues` is 256 colors"))
    if !(color_type <: SupportedPaletteColor)
        throw(ArgumentError(
            "Only 8-bit (transparent) RGB colors are supported for paletted images"
        ))
    end
end

_prepare_buffer(x::IndirectArray) = UInt8.(x.index .- first(axes(x.values, 1)))
_prepare_buffer(x::BitArray) = _prepare_buffer(collect(x))
_prepare_buffer(x::AbstractMatrix{<:T}) where {T<:Colorant{<:Normed}} = x
_prepare_buffer(x::AbstractMatrix{<:T}) where {T<:UInt8} =  reinterpret(Gray{N0f8}, x)
_prepare_buffer(x::AbstractMatrix{<:T}) where {T<:UInt16} = reinterpret(Gray{N0f16}, x)
_prepare_buffer(x::AbstractMatrix{<:T}) where {T<:Normed} = reinterpret(Gray{T}, x)
_prepare_buffer(x::AbstractMatrix{<:T}) where {T<:Gray{Bool}} =  Gray{N0f8}.(x)
_prepare_buffer(x::AbstractMatrix{<:T}) where {T<:Gray{<:AbstractFloat}} =  Gray{N0f8}.(x)
_prepare_buffer(x::AbstractMatrix{<:T}) where {T<:GrayA{<:AbstractFloat}} = GrayA{N0f8}.(x)
_prepare_buffer(x::AbstractMatrix{<:T}) where {T<:RGB{<:AbstractFloat}} =   RGB{N0f8}.(x)
_prepare_buffer(x::AbstractMatrix{<:T}) where {T<:RGBA{<:AbstractFloat}} =  RGBA{N0f8}.(x)
_prepare_buffer(x::AbstractMatrix{<:T}) where {T<:BGR{<:AbstractFloat}} =   BGR{N0f8}.(x)
_prepare_buffer(x::AbstractMatrix{<:T}) where {T<:BGRA{<:AbstractFloat}} =  BGRA{N0f8}.(x)
_prepare_buffer(x::AbstractMatrix{<:T}) where {T<:ARGB{<:AbstractFloat}} =  ARGB{N0f8}.(x)
_prepare_buffer(x::AbstractMatrix{<:T}) where {T<:ABGR{<:AbstractFloat}} =  ABGR{N0f8}.(x)
_prepare_buffer(x::AbstractMatrix{<:T}) where {T<:Union{AbstractFloat,Bool}} = reinterpret(Gray{N0f8}, N0f8.(x))
_prepare_buffer(x::AbstractArray{T,3}) where {T<:Union{AbstractFloat,Bool}} = __prepare_buffer(N0f8.(x))
_prepare_buffer(x::AbstractArray{T,3}) where {T<:Union{UInt8,Int8}} = __prepare_buffer(reinterpret(N0f8, x))
_prepare_buffer(x::AbstractArray{T,3}) where {T<:Union{UInt16,Int16}} = __prepare_buffer(reinterpret(N0f16, x))
_prepare_buffer(x::AbstractArray{T,3}) where {T<:Normed} = __prepare_buffer(x)

function __prepare_buffer(x::AbstractArray{T,3}) where {T}
    nchannels = size(x, 3)
    if nchannels == 1
        ifelse(ndims(x) == 3, _prepare_buffer(dropdims(x, dims=3)), x)
    elseif nchannels == 2
        GrayA.(colorview(GrayA, view(x, :, :, 1), view(x, :, :, 2)))
    elseif nchannels == 3
        RGB.(colorview(RGB, view(x, :, :, 1), view(x, :, :, 2), view(x, :, :, 3)))
    elseif nchannels == 4
        RGBA.(colorview(RGBA, view(x, :, :, 1), view(x, :, :, 2), view(x, :, :, 3), view(x, :, :, 4)))
    else
        error("Unsupported number of channels $(nchannels) in input. Only <= 4 is expected.")
    end
end

# On performance: if the array type has efficient convert method to Array then this is
# almost a no-op
function _enforce_dense_cpu_array(img::AbstractArray)
    if Base.has_offset_axes(img)
        convert(Array, OffsetArrays.no_offset_view(img))
    else
        convert(Array, img)
    end
end
_enforce_dense_cpu_array(img::DenseArray) = img
_enforce_dense_cpu_array(img::OffsetArray) = _enforce_dense_cpu_array(parent(img))
_enforce_dense_cpu_array(img::IndirectArray) = img # PNGFiles has built-in support for this type


_get_color_type(x::AbstractArray{<:Gray}) = PNG_COLOR_TYPE_GRAY
_get_color_type(x::AbstractArray{<:GrayA}) = PNG_COLOR_TYPE_GRAY_ALPHA
_get_color_type(x::AbstractArray{<:AbstractRGB}) = PNG_COLOR_TYPE_RGB
_get_color_type(x::AbstractArray{<:AbstractARGB}) = PNG_COLOR_TYPE_RGBA
_get_color_type(x::AbstractArray{<:AbstractRGBA}) = PNG_COLOR_TYPE_RGBA
_get_color_type(x::IndirectArray) = PNG_COLOR_TYPE_PALETTE
function _get_color_type(
        x::AbstractArray{T, N}
    ) where {
        T<:Union{Normed,Unsigned,Bool,AbstractFloat},
        N
    }
    if N == 2
        return PNG_COLOR_TYPE_GRAY
    elseif N == 3
        d = size(x, 3)
        d == 1 && (return PNG_COLOR_TYPE_GRAY)
        d == 2 && (return PNG_COLOR_TYPE_GRAY_ALPHA)
        d == 3 && (return PNG_COLOR_TYPE_RGB)
        d == 4 && (return PNG_COLOR_TYPE_RGBA)
    end
    error("Number of dimensions $(N) in image not supported (only 2D or 3D Arrays are expected).")
end

_standardize_palette(p::AbstractVector) = _enforce_dense_cpu_array(__standardize_palette(p))
__standardize_palette(p::AbstractVector{<:RGB}) = p
__standardize_palette(p::AbstractVector{<:AbstractRGB}) = RGB.(p)
__standardize_palette(p::AbstractVector{<:RGB{<:AbstractFloat}}) = RGB{N0f8}.(p)
__standardize_palette(p::AbstractVector{<:AbstractRGB{<:AbstractFloat}}) = RGB{N0f8}.(p)
__standardize_palette(p::OffsetArray) = __standardize_palette(parent(p))

_palette_alpha(p::OffsetArray) = _palette_alpha(collect(p))
_palette_alpha(p::AbstractVector{<:TransparentRGB{T,N0f8}}) where {T} = alpha.(p)
_palette_alpha(p::AbstractVector{<:TransparentRGB{T,<:AbstractFloat}}) where {T} = N0f8.(alpha.(p))
