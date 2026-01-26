_nextpow2exp(n) = (8 * sizeof(n) - leading_zeros(n) - (count_ones(n) == 1))

_get_bit_depth(img::BitArray) = 8  # TODO: write 1 bit-depth images
_get_bit_depth(img::AbstractArray{C}) where {C<:Colorant} = __get_bit_depth(eltype(C))
_get_bit_depth(img::AbstractArray{T}) where {T<:Normed} = __get_bit_depth(T)
# __get_bit_depth(::Type{Normed{T,1}}) where T = 1  # TODO: write 1 bit-depth images
# __get_bit_depth(::Type{Normed{T,2}}) where T = 2  # TODO: write 2 bit-depth images
# __get_bit_depth(::Type{Normed{T,4}}) where T = 4  # TODO: write 4 bit-depth images
__get_bit_depth(::Type{Bool}) = 8  # TODO: write 1 bit-depth images
__get_bit_depth(::Type{Normed{T,8}}) where T = 8
__get_bit_depth(::Type{Normed{T,16}}) where T = 16
__get_bit_depth(::Type{Normed{T,N}}) where {T,N} = ifelse(N <= 8, 8, 16)
__get_bit_depth(::Type{<:AbstractFloat}) = 8
_get_bit_depth(img::AbstractArray{T}) where {T<:AbstractFloat} = 8
_get_bit_depth(img::AbstractArray{<:Bool}) = 8  # TODO: write 1 bit-depth images
_get_bit_depth(img::AbstractArray{<:UInt8}) = 8
_get_bit_depth(img::AbstractArray{<:UInt16}) = 16

function _inspect_png_read(fpath, gamma::Union{Nothing,Float64}=nothing)
    fp = open_png(fpath)
    png_ptr = create_read_struct()
    info_ptr = create_info_struct(png_ptr)
    png_init_io(png_ptr, fp)
    png_set_sig_bytes(png_ptr, PNG_BYTES_TO_CHECK)
    png_read_info(png_ptr, info_ptr)

    width = png_get_image_width(png_ptr, info_ptr)
    height = png_get_image_height(png_ptr, info_ptr)
    color_type_orig = png_get_color_type(png_ptr, info_ptr)
    color_type = color_type_orig
    bit_depth_orig = png_get_bit_depth(png_ptr, info_ptr)
    bit_depth = bit_depth_orig
    num_channels = png_get_channels(png_ptr, info_ptr)
    interlace_type = png_get_interlace_type(png_ptr, info_ptr)

    backgroundp = png_color_16p()
    if png_get_bKGD(png_ptr, info_ptr, Ref(backgroundp)[]) != 0
        png_set_background(png_ptr, backgroundp, PNG_BACKGROUND_GAMMA_FILE, 1, 1.0)
    end

    if color_type == PNG_COLOR_TYPE_PALETTE
        png_set_palette_to_rgb(png_ptr)
        color_type = PNG_COLOR_TYPE_RGB
    end

    if color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8
        png_set_expand_gray_1_2_4_to_8(png_ptr)
        png_set_packing(png_ptr)
        bit_depth = UInt8(8)
    end

    if png_get_valid(png_ptr, info_ptr, PNG_INFO_tRNS) != 0
        png_set_tRNS_to_alpha(png_ptr)
        if color_type == PNG_COLOR_TYPE_GRAY || color_type == PNG_COLOR_TYPE_RGB
            color_type |= PNG_COLOR_MASK_ALPHA
        end
    end

    screen_gamma = PNG_DEFAULT_sRGB
    image_gamma = Ref{Cdouble}(0.0)
    intent = Ref{Cint}(-1)
    white_x = Ref{Cdouble}(0.0)
    white_y = Ref{Cdouble}(0.0)
    red_x = Ref{Cdouble}(0.0)
    red_y = Ref{Cdouble}(0.0)
    green_x = Ref{Cdouble}(0.0)
    green_y = Ref{Cdouble}(0.0)
    blue_x = Ref{Cdouble}(0.0)
    blue_y = Ref{Cdouble}(0.0)
    if isnothing(gamma)
        if png_get_valid(png_ptr, info_ptr, PNG_INFO_sRGB) != 0
            if png_get_sRGB(png_ptr, info_ptr, intent) != 0
                png_set_gamma(png_ptr, screen_gamma, PNG_DEFAULT_sRGB);
            else
                if png_get_gAMA(png_ptr, info_ptr, image_gamma) != 0
                    png_set_gamma(png_ptr, screen_gamma, image_gamma[])
                else
                    image_gamma[] = 0.45455
                    png_set_gamma(png_ptr, screen_gamma, image_gamma[])
                end
            end
        elseif png_get_valid(png_ptr, info_ptr, PNG_INFO_gAMA) != 0
            if png_get_gAMA(png_ptr, info_ptr, image_gamma) != 0
                png_set_gamma(png_ptr, screen_gamma, image_gamma[])
            else
                image_gamma[] = 0.45455
                png_set_gamma(png_ptr, screen_gamma, image_gamma[])
            end
        end
    elseif gamma != 1
        image_gamma[] = gamma
        png_set_gamma(png_ptr, screen_gamma, image_gamma[])
    end

    buffer_eltype = _buffer_color_type(color_type, bit_depth)
    bit_depth == 16 && png_set_swap(png_ptr)
    n_passes = png_set_interlace_handling(png_ptr)
    png_read_update_info(png_ptr, info_ptr)
    chunk_info(chunk) = png_get_valid(png_ptr, info_ptr, chunk)

    @show(fpath,
        chunk_info(PNG_INFO_sRGB),
        chunk_info(PNG_INFO_iCCP),
        chunk_info(PNG_INFO_gAMA),
        chunk_info(PNG_INFO_tRNS),
        chunk_info(PNG_INFO_cHRM),
        chunk_info(PNG_INFO_PLTE),
        chunk_info(PNG_INFO_sPLT),
        chunk_info(PNG_INFO_hIST),
        chunk_info(PNG_INFO_sBIT),
        chunk_info(PNG_INFO_bKGD),
        chunk_info(PNG_INFO_pCAL),
        (height, width, num_channels),
        (color_type_orig, color_type),
        (bit_depth_orig, bit_depth),
        (image_gamma[], screen_gamma),
        intent[],
        (white_x[], white_y[], red_x[], red_y[], green_x[], green_y[], blue_x[], blue_y[]),
        (interlace_type, n_passes),
        buffer_eltype,
    )

    png_destroy_read_struct(Ref{Ptr{Cvoid}}(png_ptr), Ref{Ptr{Cvoid}}(info_ptr), C_NULL)
    close_png(fp)
end

mutable struct _png_color_16_struct
    index::png_byte
    red::png_uint_16
    green::png_uint_16
    blue::png_uint_16
    gray::png_uint_16
end
_png_color_16_struct(c::AbstractGray{<:AbstractFloat}) = _png_color_16_struct(0x00, gray(Gray{N0f8}(c)).i, gray(Gray{N0f8}(c)).i, gray(Gray{N0f8}(c)).i, gray(Gray{N0f8}(c)).i)
_png_color_16_struct(c::AbstractGray{<:Normed}) = _png_color_16_struct(0x00, gray(c).i, gray(c).i, gray(c).i, gray(c).i)
_png_color_16_struct(c::UInt8) = _png_color_16_struct(c, 0x0000, 0x0000, 0x0000, 0x0000)
function _png_color_16_struct(c::AbstractRGB{<:AbstractFloat})
    n = RGB{N0f8}(c)
    _png_color_16_struct(0x00, red(n).i, green(n).i, blue(n).i, 0x0000)
end
_png_color_16_struct(c::AbstractRGB{<:Normed}) =  _png_color_16_struct(0x00, red(c).i, green(c).i, blue(c).i, 0x0000)


function _png_get_bKGD(png_ptr, info_ptr, background)
    ccall((:png_get_bKGD, libpng), png_uint_32, (png_const_structrp, png_inforp, Ptr{Ptr{_png_color_16_struct}}), png_ptr, info_ptr, background)
end

function _png_set_background(png_ptr, background_color, background_gamma_code, need_expand, background_gamma)
    ccall((:png_set_background, libpng), Cvoid, (png_structrp, Ptr{_png_color_16_struct}, Cint, Cint, Cdouble), png_ptr, background_color, background_gamma_code, need_expand, background_gamma)
end

function process_background(png_ptr, info_ptr, background::Bool)
    if background
        bg = _png_color_16_struct(0, 0, 0, 0, 0)
        GC.@preserve bg begin
            bgpp = Ref(Ptr{_png_color_16_struct}(pointer_from_objref(bg)))
            if _png_get_bKGD(png_ptr, info_ptr, bgpp) != 0
                _png_set_background(png_ptr, bgpp[], PNG_BACKGROUND_GAMMA_FILE, 1, 1.0)
            end
        end
        return bg
    end
    return
end

function process_background(png_ptr, info_ptr, background::Union{AbstractRGB,AbstractGray})
    bg = _png_color_16_struct(background)
    GC.@preserve bg begin
        bgp = Ptr{_png_color_16_struct}(pointer_from_objref(bg))
        _png_set_background(png_ptr, bgp, PNG_BACKGROUND_GAMMA_FILE, 1, 1.0)
        # _png_set_background(png_ptr, bgp, PNG_BACKGROUND_GAMMA_SCREEN, 0, 1.0)
    end
    return bg
end

function process_background(png_ptr, info_ptr, background::UInt8)
    bg = _png_color_16_struct(background)
    GC.@preserve bg begin
        bgp = Ptr{_png_color_16_struct}(pointer_from_objref(bg))
        _png_set_background(png_ptr, bgp, PNG_BACKGROUND_GAMMA_FILE, 1, 1.0)
        # _png_set_background(png_ptr, bgp, PNG_BACKGROUND_GAMMA_SCREEN, 0, 1.0)
    end
    return bg
end

_check_background_load(background::Bool, is_gray, has_pallette, as_paletted) = background
_check_background_load(background::AbstractGray, is_gray, has_pallette, as_paletted) = true
function _check_background_load(background::UInt8, is_gray, has_pallette, as_paletted)
    has_pallette || throw(ArgumentError("Only paletted images can use a paletted index (background::UInt8)."))
    !as_paletted || throw(ArgumentError("In order to use a (non-default) palette index to specify background, the `expand_paletted` flag must be set."))
    return true
end
function _check_background_load(background::AbstractRGB, is_gray, has_pallette, as_paletted)
    is_gray && throw(ArgumentError("Gray Images cannot use RGB background."))
    return true
end

function _check_background_save(background::AbstractGray, color_type) 
    color_type == PNG_COLOR_TYPE_GRAY || color_type == PNG_COLOR_TYPE_GRAY_ALPHA || 
        throw(ArgumentError("Only gray scale images can use a Gray for a background."))
end
function _check_background_save(background::AbstractRGB, color_type)
   (color_type & PNG_COLOR_MASK_COLOR) != 0 || 
        throw(ArgumentError("Only true color images can use a RGB for background."))
end
function _check_background_save(background::UInt8, color_type) 
    (color_type & PNG_COLOR_MASK_PALETTE) != 0 || 
        throw(ArgumentError("Only paletted images can use a paletted index for background."))
end

function _png_set_bKGD(png_ptr, info_ptr, background)
    bg = _png_color_16_struct(background)
    GC.@preserve bg begin
        bgp = Ptr{_png_color_16_struct}(pointer_from_objref(bg))
        png_set_bKGD(png_ptr, info_ptr, bgp)
    end
    return
end

_adjust_background_bitdepth(x, bit_depth) = x
function _adjust_background_bitdepth(x::AbstractRGB{T}, bit_depth) where {T}
    __get_bit_depth(T) == bit_depth && return x
    return bit_depth == 8 ? RGB{N0f8}(x) : RGB{N0f16}(x)
end
function _adjust_background_bitdepth(x::AbstractGray{T}, bit_depth) where {T}
    __get_bit_depth(T) == bit_depth && return x
    return bit_depth == 8 ? Gray{N0f8}(x) : Gray{N0f16}(x)
end