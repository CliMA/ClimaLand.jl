module Wrapper

using libwebp_jll
export libwebp_jll

using CEnum

struct WebPRGBABuffer
    rgba::Ptr{UInt8}
    stride::Cint
    size::Csize_t
end

struct WebPYUVABuffer
    y::Ptr{UInt8}
    u::Ptr{UInt8}
    v::Ptr{UInt8}
    a::Ptr{UInt8}
    y_stride::Cint
    u_stride::Cint
    v_stride::Cint
    a_stride::Cint
    y_size::Csize_t
    u_size::Csize_t
    v_size::Csize_t
    a_size::Csize_t
end

@cenum WEBP_CSP_MODE::UInt32 begin
    MODE_RGB = 0
    MODE_RGBA = 1
    MODE_BGR = 2
    MODE_BGRA = 3
    MODE_ARGB = 4
    MODE_RGBA_4444 = 5
    MODE_RGB_565 = 6
    MODE_rgbA = 7
    MODE_bgrA = 8
    MODE_Argb = 9
    MODE_rgbA_4444 = 10
    MODE_YUV = 11
    MODE_YUVA = 12
    MODE_LAST = 13
end

struct var"##Ctag#243"
    data::NTuple{80, UInt8}
end

function Base.getproperty(x::Ptr{var"##Ctag#243"}, f::Symbol)
    f === :RGBA && return Ptr{WebPRGBABuffer}(x + 0)
    f === :YUVA && return Ptr{WebPYUVABuffer}(x + 0)
    return getfield(x, f)
end

function Base.getproperty(x::var"##Ctag#243", f::Symbol)
    r = Ref{var"##Ctag#243"}(x)
    ptr = Base.unsafe_convert(Ptr{var"##Ctag#243"}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{var"##Ctag#243"}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end

struct WebPDecBuffer
    data::NTuple{120, UInt8}
end

function Base.getproperty(x::Ptr{WebPDecBuffer}, f::Symbol)
    f === :colorspace && return Ptr{WEBP_CSP_MODE}(x + 0)
    f === :width && return Ptr{Cint}(x + 4)
    f === :height && return Ptr{Cint}(x + 8)
    f === :is_external_memory && return Ptr{Cint}(x + 12)
    f === :u && return Ptr{var"##Ctag#243"}(x + 16)
    f === :pad && return Ptr{NTuple{4, UInt32}}(x + 96)
    f === :private_memory && return Ptr{Ptr{UInt8}}(x + 112)
    return getfield(x, f)
end

function Base.getproperty(x::WebPDecBuffer, f::Symbol)
    r = Ref{WebPDecBuffer}(x)
    ptr = Base.unsafe_convert(Ptr{WebPDecBuffer}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{WebPDecBuffer}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end

mutable struct WebPIDecoder end

struct WebPBitstreamFeatures
    width::Cint
    height::Cint
    has_alpha::Cint
    has_animation::Cint
    format::Cint
    pad::NTuple{5, UInt32}
end

struct WebPDecoderOptions
    bypass_filtering::Cint
    no_fancy_upsampling::Cint
    use_cropping::Cint
    crop_left::Cint
    crop_top::Cint
    crop_width::Cint
    crop_height::Cint
    use_scaling::Cint
    scaled_width::Cint
    scaled_height::Cint
    use_threads::Cint
    dithering_strength::Cint
    flip::Cint
    alpha_dithering_strength::Cint
    pad::NTuple{5, UInt32}
end

struct WebPDecoderConfig
    input::WebPBitstreamFeatures
    output::WebPDecBuffer
    options::WebPDecoderOptions
end

function WebPGetDecoderVersion()
    ccall((:WebPGetDecoderVersion, libwebp), Cint, ())
end

function WebPGetInfo(data, data_size, width, height)
    ccall((:WebPGetInfo, libwebp), Cint, (Ptr{UInt8}, Csize_t, Ptr{Cint}, Ptr{Cint}), data, data_size, width, height)
end

function WebPDecodeRGBA(data, data_size, width, height)
    ccall((:WebPDecodeRGBA, libwebp), Ptr{UInt8}, (Ptr{UInt8}, Csize_t, Ptr{Cint}, Ptr{Cint}), data, data_size, width, height)
end

function WebPDecodeARGB(data, data_size, width, height)
    ccall((:WebPDecodeARGB, libwebp), Ptr{UInt8}, (Ptr{UInt8}, Csize_t, Ptr{Cint}, Ptr{Cint}), data, data_size, width, height)
end

function WebPDecodeBGRA(data, data_size, width, height)
    ccall((:WebPDecodeBGRA, libwebp), Ptr{UInt8}, (Ptr{UInt8}, Csize_t, Ptr{Cint}, Ptr{Cint}), data, data_size, width, height)
end

function WebPDecodeRGB(data, data_size, width, height)
    ccall((:WebPDecodeRGB, libwebp), Ptr{UInt8}, (Ptr{UInt8}, Csize_t, Ptr{Cint}, Ptr{Cint}), data, data_size, width, height)
end

function WebPDecodeBGR(data, data_size, width, height)
    ccall((:WebPDecodeBGR, libwebp), Ptr{UInt8}, (Ptr{UInt8}, Csize_t, Ptr{Cint}, Ptr{Cint}), data, data_size, width, height)
end

function WebPDecodeYUV(data, data_size, width, height, u, v, stride, uv_stride)
    ccall((:WebPDecodeYUV, libwebp), Ptr{UInt8}, (Ptr{UInt8}, Csize_t, Ptr{Cint}, Ptr{Cint}, Ptr{Ptr{UInt8}}, Ptr{Ptr{UInt8}}, Ptr{Cint}, Ptr{Cint}), data, data_size, width, height, u, v, stride, uv_stride)
end

function WebPDecodeRGBAInto(data, data_size, output_buffer, output_buffer_size, output_stride)
    ccall((:WebPDecodeRGBAInto, libwebp), Ptr{UInt8}, (Ptr{UInt8}, Csize_t, Ptr{UInt8}, Csize_t, Cint), data, data_size, output_buffer, output_buffer_size, output_stride)
end

function WebPDecodeARGBInto(data, data_size, output_buffer, output_buffer_size, output_stride)
    ccall((:WebPDecodeARGBInto, libwebp), Ptr{UInt8}, (Ptr{UInt8}, Csize_t, Ptr{UInt8}, Csize_t, Cint), data, data_size, output_buffer, output_buffer_size, output_stride)
end

function WebPDecodeBGRAInto(data, data_size, output_buffer, output_buffer_size, output_stride)
    ccall((:WebPDecodeBGRAInto, libwebp), Ptr{UInt8}, (Ptr{UInt8}, Csize_t, Ptr{UInt8}, Csize_t, Cint), data, data_size, output_buffer, output_buffer_size, output_stride)
end

function WebPDecodeRGBInto(data, data_size, output_buffer, output_buffer_size, output_stride)
    ccall((:WebPDecodeRGBInto, libwebp), Ptr{UInt8}, (Ptr{UInt8}, Csize_t, Ptr{UInt8}, Csize_t, Cint), data, data_size, output_buffer, output_buffer_size, output_stride)
end

function WebPDecodeBGRInto(data, data_size, output_buffer, output_buffer_size, output_stride)
    ccall((:WebPDecodeBGRInto, libwebp), Ptr{UInt8}, (Ptr{UInt8}, Csize_t, Ptr{UInt8}, Csize_t, Cint), data, data_size, output_buffer, output_buffer_size, output_stride)
end

function WebPDecodeYUVInto(data, data_size, luma, luma_size, luma_stride, u, u_size, u_stride, v, v_size, v_stride)
    ccall((:WebPDecodeYUVInto, libwebp), Ptr{UInt8}, (Ptr{UInt8}, Csize_t, Ptr{UInt8}, Csize_t, Cint, Ptr{UInt8}, Csize_t, Cint, Ptr{UInt8}, Csize_t, Cint), data, data_size, luma, luma_size, luma_stride, u, u_size, u_stride, v, v_size, v_stride)
end

function WebPIsPremultipliedMode(mode)
    ccall((:WebPIsPremultipliedMode, libwebp), Cint, (WEBP_CSP_MODE,), mode)
end

function WebPIsAlphaMode(mode)
    ccall((:WebPIsAlphaMode, libwebp), Cint, (WEBP_CSP_MODE,), mode)
end

function WebPIsRGBMode(mode)
    ccall((:WebPIsRGBMode, libwebp), Cint, (WEBP_CSP_MODE,), mode)
end

function WebPInitDecBufferInternal(arg1, arg2)
    ccall((:WebPInitDecBufferInternal, libwebp), Cint, (Ptr{WebPDecBuffer}, Cint), arg1, arg2)
end

function WebPInitDecBuffer(buffer)
    ccall((:WebPInitDecBuffer, libwebp), Cint, (Ptr{WebPDecBuffer},), buffer)
end

function WebPFreeDecBuffer(buffer)
    ccall((:WebPFreeDecBuffer, libwebp), Cvoid, (Ptr{WebPDecBuffer},), buffer)
end

@cenum VP8StatusCode::UInt32 begin
    VP8_STATUS_OK = 0
    VP8_STATUS_OUT_OF_MEMORY = 1
    VP8_STATUS_INVALID_PARAM = 2
    VP8_STATUS_BITSTREAM_ERROR = 3
    VP8_STATUS_UNSUPPORTED_FEATURE = 4
    VP8_STATUS_SUSPENDED = 5
    VP8_STATUS_USER_ABORT = 6
    VP8_STATUS_NOT_ENOUGH_DATA = 7
end

function WebPINewDecoder(output_buffer)
    ccall((:WebPINewDecoder, libwebp), Ptr{WebPIDecoder}, (Ptr{WebPDecBuffer},), output_buffer)
end

function WebPINewRGB(csp, output_buffer, output_buffer_size, output_stride)
    ccall((:WebPINewRGB, libwebp), Ptr{WebPIDecoder}, (WEBP_CSP_MODE, Ptr{UInt8}, Csize_t, Cint), csp, output_buffer, output_buffer_size, output_stride)
end

function WebPINewYUVA(luma, luma_size, luma_stride, u, u_size, u_stride, v, v_size, v_stride, a, a_size, a_stride)
    ccall((:WebPINewYUVA, libwebp), Ptr{WebPIDecoder}, (Ptr{UInt8}, Csize_t, Cint, Ptr{UInt8}, Csize_t, Cint, Ptr{UInt8}, Csize_t, Cint, Ptr{UInt8}, Csize_t, Cint), luma, luma_size, luma_stride, u, u_size, u_stride, v, v_size, v_stride, a, a_size, a_stride)
end

function WebPINewYUV(luma, luma_size, luma_stride, u, u_size, u_stride, v, v_size, v_stride)
    ccall((:WebPINewYUV, libwebp), Ptr{WebPIDecoder}, (Ptr{UInt8}, Csize_t, Cint, Ptr{UInt8}, Csize_t, Cint, Ptr{UInt8}, Csize_t, Cint), luma, luma_size, luma_stride, u, u_size, u_stride, v, v_size, v_stride)
end

function WebPIDelete(idec)
    ccall((:WebPIDelete, libwebp), Cvoid, (Ptr{WebPIDecoder},), idec)
end

function WebPIAppend(idec, data, data_size)
    ccall((:WebPIAppend, libwebp), VP8StatusCode, (Ptr{WebPIDecoder}, Ptr{UInt8}, Csize_t), idec, data, data_size)
end

function WebPIUpdate(idec, data, data_size)
    ccall((:WebPIUpdate, libwebp), VP8StatusCode, (Ptr{WebPIDecoder}, Ptr{UInt8}, Csize_t), idec, data, data_size)
end

function WebPIDecGetRGB(idec, last_y, width, height, stride)
    ccall((:WebPIDecGetRGB, libwebp), Ptr{UInt8}, (Ptr{WebPIDecoder}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}), idec, last_y, width, height, stride)
end

function WebPIDecGetYUVA(idec, last_y, u, v, a, width, height, stride, uv_stride, a_stride)
    ccall((:WebPIDecGetYUVA, libwebp), Ptr{UInt8}, (Ptr{WebPIDecoder}, Ptr{Cint}, Ptr{Ptr{UInt8}}, Ptr{Ptr{UInt8}}, Ptr{Ptr{UInt8}}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}), idec, last_y, u, v, a, width, height, stride, uv_stride, a_stride)
end

function WebPIDecGetYUV(idec, last_y, u, v, width, height, stride, uv_stride)
    ccall((:WebPIDecGetYUV, libwebp), Ptr{UInt8}, (Ptr{WebPIDecoder}, Ptr{Cint}, Ptr{Ptr{UInt8}}, Ptr{Ptr{UInt8}}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}), idec, last_y, u, v, width, height, stride, uv_stride)
end

function WebPIDecodedArea(idec, left, top, width, height)
    ccall((:WebPIDecodedArea, libwebp), Ptr{WebPDecBuffer}, (Ptr{WebPIDecoder}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}), idec, left, top, width, height)
end

function WebPGetFeaturesInternal(arg1, arg2, arg3, arg4)
    ccall((:WebPGetFeaturesInternal, libwebp), VP8StatusCode, (Ptr{UInt8}, Csize_t, Ptr{WebPBitstreamFeatures}, Cint), arg1, arg2, arg3, arg4)
end

function WebPGetFeatures(data, data_size, features)
    ccall((:WebPGetFeatures, libwebp), VP8StatusCode, (Ptr{UInt8}, Csize_t, Ptr{WebPBitstreamFeatures}), data, data_size, features)
end

function WebPInitDecoderConfigInternal(arg1, arg2)
    ccall((:WebPInitDecoderConfigInternal, libwebp), Cint, (Ptr{WebPDecoderConfig}, Cint), arg1, arg2)
end

function WebPInitDecoderConfig(config)
    ccall((:WebPInitDecoderConfig, libwebp), Cint, (Ptr{WebPDecoderConfig},), config)
end

function WebPIDecode(data, data_size, config)
    ccall((:WebPIDecode, libwebp), Ptr{WebPIDecoder}, (Ptr{UInt8}, Csize_t, Ptr{WebPDecoderConfig}), data, data_size, config)
end

function WebPDecode(data, data_size, config)
    ccall((:WebPDecode, libwebp), VP8StatusCode, (Ptr{UInt8}, Csize_t, Ptr{WebPDecoderConfig}), data, data_size, config)
end

mutable struct WebPDemuxer end

@cenum WebPMuxAnimDispose::UInt32 begin
    WEBP_MUX_DISPOSE_NONE = 0
    WEBP_MUX_DISPOSE_BACKGROUND = 1
end

struct WebPData
    bytes::Ptr{UInt8}
    size::Csize_t
end

@cenum WebPMuxAnimBlend::UInt32 begin
    WEBP_MUX_BLEND = 0
    WEBP_MUX_NO_BLEND = 1
end

struct WebPIterator
    frame_num::Cint
    num_frames::Cint
    x_offset::Cint
    y_offset::Cint
    width::Cint
    height::Cint
    duration::Cint
    dispose_method::WebPMuxAnimDispose
    complete::Cint
    fragment::WebPData
    has_alpha::Cint
    blend_method::WebPMuxAnimBlend
    pad::NTuple{2, UInt32}
    private_::Ptr{Cvoid}
end

struct WebPChunkIterator
    chunk_num::Cint
    num_chunks::Cint
    chunk::WebPData
    pad::NTuple{6, UInt32}
    private_::Ptr{Cvoid}
end

struct WebPAnimInfo
    canvas_width::UInt32
    canvas_height::UInt32
    loop_count::UInt32
    bgcolor::UInt32
    frame_count::UInt32
    pad::NTuple{4, UInt32}
end

struct WebPAnimDecoderOptions
    color_mode::WEBP_CSP_MODE
    use_threads::Cint
    padding::NTuple{7, UInt32}
end

function WebPGetDemuxVersion()
    ccall((:WebPGetDemuxVersion, libwebp), Cint, ())
end

@cenum WebPDemuxState::Int32 begin
    WEBP_DEMUX_PARSE_ERROR = -1
    WEBP_DEMUX_PARSING_HEADER = 0
    WEBP_DEMUX_PARSED_HEADER = 1
    WEBP_DEMUX_DONE = 2
end

function WebPDemuxInternal(arg1, arg2, arg3, arg4)
    ccall((:WebPDemuxInternal, libwebp), Ptr{WebPDemuxer}, (Ptr{WebPData}, Cint, Ptr{WebPDemuxState}, Cint), arg1, arg2, arg3, arg4)
end

function WebPDemux(data)
    ccall((:WebPDemux, libwebp), Ptr{WebPDemuxer}, (Ptr{WebPData},), data)
end

function WebPDemuxPartial(data, state)
    ccall((:WebPDemuxPartial, libwebp), Ptr{WebPDemuxer}, (Ptr{WebPData}, Ptr{WebPDemuxState}), data, state)
end

function WebPDemuxDelete(dmux)
    ccall((:WebPDemuxDelete, libwebp), Cvoid, (Ptr{WebPDemuxer},), dmux)
end

@cenum WebPFormatFeature::UInt32 begin
    WEBP_FF_FORMAT_FLAGS = 0
    WEBP_FF_CANVAS_WIDTH = 1
    WEBP_FF_CANVAS_HEIGHT = 2
    WEBP_FF_LOOP_COUNT = 3
    WEBP_FF_BACKGROUND_COLOR = 4
    WEBP_FF_FRAME_COUNT = 5
end

function WebPDemuxGetI(dmux, feature)
    ccall((:WebPDemuxGetI, libwebp), UInt32, (Ptr{WebPDemuxer}, WebPFormatFeature), dmux, feature)
end

function WebPDemuxGetFrame(dmux, frame_number, iter)
    ccall((:WebPDemuxGetFrame, libwebp), Cint, (Ptr{WebPDemuxer}, Cint, Ptr{WebPIterator}), dmux, frame_number, iter)
end

function WebPDemuxNextFrame(iter)
    ccall((:WebPDemuxNextFrame, libwebp), Cint, (Ptr{WebPIterator},), iter)
end

function WebPDemuxPrevFrame(iter)
    ccall((:WebPDemuxPrevFrame, libwebp), Cint, (Ptr{WebPIterator},), iter)
end

function WebPDemuxReleaseIterator(iter)
    ccall((:WebPDemuxReleaseIterator, libwebp), Cvoid, (Ptr{WebPIterator},), iter)
end

function WebPDemuxGetChunk(dmux, fourcc, chunk_number, iter)
    ccall((:WebPDemuxGetChunk, libwebp), Cint, (Ptr{WebPDemuxer}, Ptr{Cchar}, Cint, Ptr{WebPChunkIterator}), dmux, fourcc, chunk_number, iter)
end

function WebPDemuxNextChunk(iter)
    ccall((:WebPDemuxNextChunk, libwebp), Cint, (Ptr{WebPChunkIterator},), iter)
end

function WebPDemuxPrevChunk(iter)
    ccall((:WebPDemuxPrevChunk, libwebp), Cint, (Ptr{WebPChunkIterator},), iter)
end

function WebPDemuxReleaseChunkIterator(iter)
    ccall((:WebPDemuxReleaseChunkIterator, libwebp), Cvoid, (Ptr{WebPChunkIterator},), iter)
end

mutable struct WebPAnimDecoder end

function WebPAnimDecoderOptionsInitInternal(arg1, arg2)
    ccall((:WebPAnimDecoderOptionsInitInternal, libwebp), Cint, (Ptr{WebPAnimDecoderOptions}, Cint), arg1, arg2)
end

function WebPAnimDecoderOptionsInit(dec_options)
    ccall((:WebPAnimDecoderOptionsInit, libwebp), Cint, (Ptr{WebPAnimDecoderOptions},), dec_options)
end

function WebPAnimDecoderNewInternal(arg1, arg2, arg3)
    ccall((:WebPAnimDecoderNewInternal, libwebp), Ptr{WebPAnimDecoder}, (Ptr{WebPData}, Ptr{WebPAnimDecoderOptions}, Cint), arg1, arg2, arg3)
end

function WebPAnimDecoderNew(webp_data, dec_options)
    ccall((:WebPAnimDecoderNew, libwebp), Ptr{WebPAnimDecoder}, (Ptr{WebPData}, Ptr{WebPAnimDecoderOptions}), webp_data, dec_options)
end

function WebPAnimDecoderGetInfo(dec, info)
    ccall((:WebPAnimDecoderGetInfo, libwebp), Cint, (Ptr{WebPAnimDecoder}, Ptr{WebPAnimInfo}), dec, info)
end

function WebPAnimDecoderGetNext(dec, buf, timestamp)
    ccall((:WebPAnimDecoderGetNext, libwebp), Cint, (Ptr{WebPAnimDecoder}, Ptr{Ptr{UInt8}}, Ptr{Cint}), dec, buf, timestamp)
end

function WebPAnimDecoderHasMoreFrames(dec)
    ccall((:WebPAnimDecoderHasMoreFrames, libwebp), Cint, (Ptr{WebPAnimDecoder},), dec)
end

function WebPAnimDecoderReset(dec)
    ccall((:WebPAnimDecoderReset, libwebp), Cvoid, (Ptr{WebPAnimDecoder},), dec)
end

function WebPAnimDecoderGetDemuxer(dec)
    ccall((:WebPAnimDecoderGetDemuxer, libwebp), Ptr{WebPDemuxer}, (Ptr{WebPAnimDecoder},), dec)
end

function WebPAnimDecoderDelete(dec)
    ccall((:WebPAnimDecoderDelete, libwebp), Cvoid, (Ptr{WebPAnimDecoder},), dec)
end

@cenum WebPImageHint::UInt32 begin
    WEBP_HINT_DEFAULT = 0
    WEBP_HINT_PICTURE = 1
    WEBP_HINT_PHOTO = 2
    WEBP_HINT_GRAPH = 3
    WEBP_HINT_LAST = 4
end

struct WebPConfig
    lossless::Cint
    quality::Cfloat
    method::Cint
    image_hint::WebPImageHint
    target_size::Cint
    target_PSNR::Cfloat
    segments::Cint
    sns_strength::Cint
    filter_strength::Cint
    filter_sharpness::Cint
    filter_type::Cint
    autofilter::Cint
    alpha_compression::Cint
    alpha_filtering::Cint
    alpha_quality::Cint
    pass::Cint
    show_compressed::Cint
    preprocessing::Cint
    partitions::Cint
    partition_limit::Cint
    emulate_jpeg_size::Cint
    thread_level::Cint
    low_memory::Cint
    near_lossless::Cint
    exact::Cint
    use_delta_palette::Cint
    use_sharp_yuv::Cint
    qmin::Cint
    qmax::Cint
end

@cenum WebPEncCSP::UInt32 begin
    WEBP_YUV420 = 0
    WEBP_YUV420A = 4
    WEBP_CSP_UV_MASK = 3
    WEBP_CSP_ALPHA_BIT = 4
end

# typedef int ( * WebPWriterFunction ) ( const uint8_t * data , size_t data_size , const WebPPicture * picture )
const WebPWriterFunction = Ptr{Cvoid}

struct WebPAuxStats
    coded_size::Cint
    PSNR::NTuple{5, Cfloat}
    block_count::NTuple{3, Cint}
    header_bytes::NTuple{2, Cint}
    residual_bytes::NTuple{3, NTuple{4, Cint}}
    segment_size::NTuple{4, Cint}
    segment_quant::NTuple{4, Cint}
    segment_level::NTuple{4, Cint}
    alpha_data_size::Cint
    layer_data_size::Cint
    lossless_features::UInt32
    histogram_bits::Cint
    transform_bits::Cint
    cache_bits::Cint
    palette_size::Cint
    lossless_size::Cint
    lossless_hdr_size::Cint
    lossless_data_size::Cint
    pad::NTuple{2, UInt32}
end

@cenum WebPEncodingError::UInt32 begin
    VP8_ENC_OK = 0
    VP8_ENC_ERROR_OUT_OF_MEMORY = 1
    VP8_ENC_ERROR_BITSTREAM_OUT_OF_MEMORY = 2
    VP8_ENC_ERROR_NULL_PARAMETER = 3
    VP8_ENC_ERROR_INVALID_CONFIGURATION = 4
    VP8_ENC_ERROR_BAD_DIMENSION = 5
    VP8_ENC_ERROR_PARTITION0_OVERFLOW = 6
    VP8_ENC_ERROR_PARTITION_OVERFLOW = 7
    VP8_ENC_ERROR_BAD_WRITE = 8
    VP8_ENC_ERROR_FILE_TOO_BIG = 9
    VP8_ENC_ERROR_USER_ABORT = 10
    VP8_ENC_ERROR_LAST = 11
end

# typedef int ( * WebPProgressHook ) ( int percent , const WebPPicture * picture )
const WebPProgressHook = Ptr{Cvoid}

struct WebPPicture
    use_argb::Cint
    colorspace::WebPEncCSP
    width::Cint
    height::Cint
    y::Ptr{UInt8}
    u::Ptr{UInt8}
    v::Ptr{UInt8}
    y_stride::Cint
    uv_stride::Cint
    a::Ptr{UInt8}
    a_stride::Cint
    pad1::NTuple{2, UInt32}
    argb::Ptr{UInt32}
    argb_stride::Cint
    pad2::NTuple{3, UInt32}
    writer::WebPWriterFunction
    custom_ptr::Ptr{Cvoid}
    extra_info_type::Cint
    extra_info::Ptr{UInt8}
    stats::Ptr{WebPAuxStats}
    error_code::WebPEncodingError
    progress_hook::WebPProgressHook
    user_data::Ptr{Cvoid}
    pad3::NTuple{3, UInt32}
    pad4::Ptr{UInt8}
    pad5::Ptr{UInt8}
    pad6::NTuple{8, UInt32}
    memory_::Ptr{Cvoid}
    memory_argb_::Ptr{Cvoid}
    pad7::NTuple{2, Ptr{Cvoid}}
end

struct WebPMemoryWriter
    mem::Ptr{UInt8}
    size::Csize_t
    max_size::Csize_t
    pad::NTuple{1, UInt32}
end

function WebPGetEncoderVersion()
    ccall((:WebPGetEncoderVersion, libwebp), Cint, ())
end

function WebPEncodeRGB(rgb, width, height, stride, quality_factor, output)
    ccall((:WebPEncodeRGB, libwebp), Csize_t, (Ptr{UInt8}, Cint, Cint, Cint, Cfloat, Ptr{Ptr{UInt8}}), rgb, width, height, stride, quality_factor, output)
end

function WebPEncodeBGR(bgr, width, height, stride, quality_factor, output)
    ccall((:WebPEncodeBGR, libwebp), Csize_t, (Ptr{UInt8}, Cint, Cint, Cint, Cfloat, Ptr{Ptr{UInt8}}), bgr, width, height, stride, quality_factor, output)
end

function WebPEncodeRGBA(rgba, width, height, stride, quality_factor, output)
    ccall((:WebPEncodeRGBA, libwebp), Csize_t, (Ptr{UInt8}, Cint, Cint, Cint, Cfloat, Ptr{Ptr{UInt8}}), rgba, width, height, stride, quality_factor, output)
end

function WebPEncodeBGRA(bgra, width, height, stride, quality_factor, output)
    ccall((:WebPEncodeBGRA, libwebp), Csize_t, (Ptr{UInt8}, Cint, Cint, Cint, Cfloat, Ptr{Ptr{UInt8}}), bgra, width, height, stride, quality_factor, output)
end

function WebPEncodeLosslessRGB(rgb, width, height, stride, output)
    ccall((:WebPEncodeLosslessRGB, libwebp), Csize_t, (Ptr{UInt8}, Cint, Cint, Cint, Ptr{Ptr{UInt8}}), rgb, width, height, stride, output)
end

function WebPEncodeLosslessBGR(bgr, width, height, stride, output)
    ccall((:WebPEncodeLosslessBGR, libwebp), Csize_t, (Ptr{UInt8}, Cint, Cint, Cint, Ptr{Ptr{UInt8}}), bgr, width, height, stride, output)
end

function WebPEncodeLosslessRGBA(rgba, width, height, stride, output)
    ccall((:WebPEncodeLosslessRGBA, libwebp), Csize_t, (Ptr{UInt8}, Cint, Cint, Cint, Ptr{Ptr{UInt8}}), rgba, width, height, stride, output)
end

function WebPEncodeLosslessBGRA(bgra, width, height, stride, output)
    ccall((:WebPEncodeLosslessBGRA, libwebp), Csize_t, (Ptr{UInt8}, Cint, Cint, Cint, Ptr{Ptr{UInt8}}), bgra, width, height, stride, output)
end

@cenum WebPPreset::UInt32 begin
    WEBP_PRESET_DEFAULT = 0
    WEBP_PRESET_PICTURE = 1
    WEBP_PRESET_PHOTO = 2
    WEBP_PRESET_DRAWING = 3
    WEBP_PRESET_ICON = 4
    WEBP_PRESET_TEXT = 5
end

function WebPConfigInitInternal(arg1, arg2, arg3, arg4)
    ccall((:WebPConfigInitInternal, libwebp), Cint, (Ptr{WebPConfig}, WebPPreset, Cfloat, Cint), arg1, arg2, arg3, arg4)
end

function WebPConfigInit(config)
    ccall((:WebPConfigInit, libwebp), Cint, (Ptr{WebPConfig},), config)
end

function WebPConfigPreset(config, preset, quality)
    ccall((:WebPConfigPreset, libwebp), Cint, (Ptr{WebPConfig}, WebPPreset, Cfloat), config, preset, quality)
end

function WebPConfigLosslessPreset(config, level)
    ccall((:WebPConfigLosslessPreset, libwebp), Cint, (Ptr{WebPConfig}, Cint), config, level)
end

function WebPValidateConfig(config)
    ccall((:WebPValidateConfig, libwebp), Cint, (Ptr{WebPConfig},), config)
end

function WebPMemoryWriterInit(writer)
    ccall((:WebPMemoryWriterInit, libwebp), Cvoid, (Ptr{WebPMemoryWriter},), writer)
end

function WebPMemoryWriterClear(writer)
    ccall((:WebPMemoryWriterClear, libwebp), Cvoid, (Ptr{WebPMemoryWriter},), writer)
end

function WebPMemoryWrite(data, data_size, picture)
    ccall((:WebPMemoryWrite, libwebp), Cint, (Ptr{UInt8}, Csize_t, Ptr{WebPPicture}), data, data_size, picture)
end

function WebPPictureInitInternal(arg1, arg2)
    ccall((:WebPPictureInitInternal, libwebp), Cint, (Ptr{WebPPicture}, Cint), arg1, arg2)
end

function WebPPictureInit(picture)
    ccall((:WebPPictureInit, libwebp), Cint, (Ptr{WebPPicture},), picture)
end

function WebPPictureAlloc(picture)
    ccall((:WebPPictureAlloc, libwebp), Cint, (Ptr{WebPPicture},), picture)
end

function WebPPictureFree(picture)
    ccall((:WebPPictureFree, libwebp), Cvoid, (Ptr{WebPPicture},), picture)
end

function WebPPictureCopy(src, dst)
    ccall((:WebPPictureCopy, libwebp), Cint, (Ptr{WebPPicture}, Ptr{WebPPicture}), src, dst)
end

function WebPPlaneDistortion(src, src_stride, ref, ref_stride, width, height, x_step, type, distortion, result)
    ccall((:WebPPlaneDistortion, libwebp), Cint, (Ptr{UInt8}, Csize_t, Ptr{UInt8}, Csize_t, Cint, Cint, Csize_t, Cint, Ptr{Cfloat}, Ptr{Cfloat}), src, src_stride, ref, ref_stride, width, height, x_step, type, distortion, result)
end

function WebPPictureDistortion(src, ref, metric_type, result)
    ccall((:WebPPictureDistortion, libwebp), Cint, (Ptr{WebPPicture}, Ptr{WebPPicture}, Cint, Ptr{Cfloat}), src, ref, metric_type, result)
end

function WebPPictureCrop(picture, left, top, width, height)
    ccall((:WebPPictureCrop, libwebp), Cint, (Ptr{WebPPicture}, Cint, Cint, Cint, Cint), picture, left, top, width, height)
end

function WebPPictureView(src, left, top, width, height, dst)
    ccall((:WebPPictureView, libwebp), Cint, (Ptr{WebPPicture}, Cint, Cint, Cint, Cint, Ptr{WebPPicture}), src, left, top, width, height, dst)
end

function WebPPictureIsView(picture)
    ccall((:WebPPictureIsView, libwebp), Cint, (Ptr{WebPPicture},), picture)
end

function WebPPictureRescale(picture, width, height)
    ccall((:WebPPictureRescale, libwebp), Cint, (Ptr{WebPPicture}, Cint, Cint), picture, width, height)
end

function WebPPictureImportRGB(picture, rgb, rgb_stride)
    ccall((:WebPPictureImportRGB, libwebp), Cint, (Ptr{WebPPicture}, Ptr{UInt8}, Cint), picture, rgb, rgb_stride)
end

function WebPPictureImportRGBA(picture, rgba, rgba_stride)
    ccall((:WebPPictureImportRGBA, libwebp), Cint, (Ptr{WebPPicture}, Ptr{UInt8}, Cint), picture, rgba, rgba_stride)
end

function WebPPictureImportRGBX(picture, rgbx, rgbx_stride)
    ccall((:WebPPictureImportRGBX, libwebp), Cint, (Ptr{WebPPicture}, Ptr{UInt8}, Cint), picture, rgbx, rgbx_stride)
end

function WebPPictureImportBGR(picture, bgr, bgr_stride)
    ccall((:WebPPictureImportBGR, libwebp), Cint, (Ptr{WebPPicture}, Ptr{UInt8}, Cint), picture, bgr, bgr_stride)
end

function WebPPictureImportBGRA(picture, bgra, bgra_stride)
    ccall((:WebPPictureImportBGRA, libwebp), Cint, (Ptr{WebPPicture}, Ptr{UInt8}, Cint), picture, bgra, bgra_stride)
end

function WebPPictureImportBGRX(picture, bgrx, bgrx_stride)
    ccall((:WebPPictureImportBGRX, libwebp), Cint, (Ptr{WebPPicture}, Ptr{UInt8}, Cint), picture, bgrx, bgrx_stride)
end

function WebPPictureARGBToYUVA(picture, arg2)
    ccall((:WebPPictureARGBToYUVA, libwebp), Cint, (Ptr{WebPPicture}, WebPEncCSP), picture, arg2)
end

function WebPPictureARGBToYUVADithered(picture, colorspace, dithering)
    ccall((:WebPPictureARGBToYUVADithered, libwebp), Cint, (Ptr{WebPPicture}, WebPEncCSP, Cfloat), picture, colorspace, dithering)
end

function WebPPictureSharpARGBToYUVA(picture)
    ccall((:WebPPictureSharpARGBToYUVA, libwebp), Cint, (Ptr{WebPPicture},), picture)
end

function WebPPictureSmartARGBToYUVA(picture)
    ccall((:WebPPictureSmartARGBToYUVA, libwebp), Cint, (Ptr{WebPPicture},), picture)
end

function WebPPictureYUVAToARGB(picture)
    ccall((:WebPPictureYUVAToARGB, libwebp), Cint, (Ptr{WebPPicture},), picture)
end

function WebPCleanupTransparentArea(picture)
    ccall((:WebPCleanupTransparentArea, libwebp), Cvoid, (Ptr{WebPPicture},), picture)
end

function WebPPictureHasTransparency(picture)
    ccall((:WebPPictureHasTransparency, libwebp), Cint, (Ptr{WebPPicture},), picture)
end

function WebPBlendAlpha(picture, background_rgb)
    ccall((:WebPBlendAlpha, libwebp), Cvoid, (Ptr{WebPPicture}, UInt32), picture, background_rgb)
end

function WebPEncode(config, picture)
    ccall((:WebPEncode, libwebp), Cint, (Ptr{WebPConfig}, Ptr{WebPPicture}), config, picture)
end

mutable struct WebPMux end

@cenum WebPChunkId::UInt32 begin
    WEBP_CHUNK_VP8X = 0
    WEBP_CHUNK_ICCP = 1
    WEBP_CHUNK_ANIM = 2
    WEBP_CHUNK_ANMF = 3
    WEBP_CHUNK_DEPRECATED = 4
    WEBP_CHUNK_ALPHA = 5
    WEBP_CHUNK_IMAGE = 6
    WEBP_CHUNK_EXIF = 7
    WEBP_CHUNK_XMP = 8
    WEBP_CHUNK_UNKNOWN = 9
    WEBP_CHUNK_NIL = 10
end

struct WebPMuxFrameInfo
    bitstream::WebPData
    x_offset::Cint
    y_offset::Cint
    duration::Cint
    id::WebPChunkId
    dispose_method::WebPMuxAnimDispose
    blend_method::WebPMuxAnimBlend
    pad::NTuple{1, UInt32}
end

struct WebPMuxAnimParams
    bgcolor::UInt32
    loop_count::Cint
end

struct WebPAnimEncoderOptions
    anim_params::WebPMuxAnimParams
    minimize_size::Cint
    kmin::Cint
    kmax::Cint
    allow_mixed::Cint
    verbose::Cint
    padding::NTuple{4, UInt32}
end

@cenum WebPMuxError::Int32 begin
    WEBP_MUX_OK = 1
    WEBP_MUX_NOT_FOUND = 0
    WEBP_MUX_INVALID_ARGUMENT = -1
    WEBP_MUX_BAD_DATA = -2
    WEBP_MUX_MEMORY_ERROR = -3
    WEBP_MUX_NOT_ENOUGH_DATA = -4
end

function WebPGetMuxVersion()
    ccall((:WebPGetMuxVersion, libwebp), Cint, ())
end

function WebPNewInternal(arg1)
    ccall((:WebPNewInternal, libwebp), Ptr{WebPMux}, (Cint,), arg1)
end

function WebPMuxNew()
    ccall((:WebPMuxNew, libwebp), Ptr{WebPMux}, ())
end

function WebPMuxDelete(mux)
    ccall((:WebPMuxDelete, libwebp), Cvoid, (Ptr{WebPMux},), mux)
end

function WebPMuxCreateInternal(arg1, arg2, arg3)
    ccall((:WebPMuxCreateInternal, libwebp), Ptr{WebPMux}, (Ptr{WebPData}, Cint, Cint), arg1, arg2, arg3)
end

function WebPMuxCreate(bitstream, copy_data)
    ccall((:WebPMuxCreate, libwebp), Ptr{WebPMux}, (Ptr{WebPData}, Cint), bitstream, copy_data)
end

function WebPMuxSetChunk(mux, fourcc, chunk_data, copy_data)
    ccall((:WebPMuxSetChunk, libwebp), WebPMuxError, (Ptr{WebPMux}, Ptr{Cchar}, Ptr{WebPData}, Cint), mux, fourcc, chunk_data, copy_data)
end

function WebPMuxGetChunk(mux, fourcc, chunk_data)
    ccall((:WebPMuxGetChunk, libwebp), WebPMuxError, (Ptr{WebPMux}, Ptr{Cchar}, Ptr{WebPData}), mux, fourcc, chunk_data)
end

function WebPMuxDeleteChunk(mux, fourcc)
    ccall((:WebPMuxDeleteChunk, libwebp), WebPMuxError, (Ptr{WebPMux}, Ptr{Cchar}), mux, fourcc)
end

function WebPMuxSetImage(mux, bitstream, copy_data)
    ccall((:WebPMuxSetImage, libwebp), WebPMuxError, (Ptr{WebPMux}, Ptr{WebPData}, Cint), mux, bitstream, copy_data)
end

function WebPMuxPushFrame(mux, frame, copy_data)
    ccall((:WebPMuxPushFrame, libwebp), WebPMuxError, (Ptr{WebPMux}, Ptr{WebPMuxFrameInfo}, Cint), mux, frame, copy_data)
end

function WebPMuxGetFrame(mux, nth, frame)
    ccall((:WebPMuxGetFrame, libwebp), WebPMuxError, (Ptr{WebPMux}, UInt32, Ptr{WebPMuxFrameInfo}), mux, nth, frame)
end

function WebPMuxDeleteFrame(mux, nth)
    ccall((:WebPMuxDeleteFrame, libwebp), WebPMuxError, (Ptr{WebPMux}, UInt32), mux, nth)
end

function WebPMuxSetAnimationParams(mux, params)
    ccall((:WebPMuxSetAnimationParams, libwebp), WebPMuxError, (Ptr{WebPMux}, Ptr{WebPMuxAnimParams}), mux, params)
end

function WebPMuxGetAnimationParams(mux, params)
    ccall((:WebPMuxGetAnimationParams, libwebp), WebPMuxError, (Ptr{WebPMux}, Ptr{WebPMuxAnimParams}), mux, params)
end

function WebPMuxSetCanvasSize(mux, width, height)
    ccall((:WebPMuxSetCanvasSize, libwebp), WebPMuxError, (Ptr{WebPMux}, Cint, Cint), mux, width, height)
end

function WebPMuxGetCanvasSize(mux, width, height)
    ccall((:WebPMuxGetCanvasSize, libwebp), WebPMuxError, (Ptr{WebPMux}, Ptr{Cint}, Ptr{Cint}), mux, width, height)
end

function WebPMuxGetFeatures(mux, flags)
    ccall((:WebPMuxGetFeatures, libwebp), WebPMuxError, (Ptr{WebPMux}, Ptr{UInt32}), mux, flags)
end

function WebPMuxNumChunks(mux, id, num_elements)
    ccall((:WebPMuxNumChunks, libwebp), WebPMuxError, (Ptr{WebPMux}, WebPChunkId, Ptr{Cint}), mux, id, num_elements)
end

function WebPMuxAssemble(mux, assembled_data)
    ccall((:WebPMuxAssemble, libwebp), WebPMuxError, (Ptr{WebPMux}, Ptr{WebPData}), mux, assembled_data)
end

mutable struct WebPAnimEncoder end

function WebPAnimEncoderOptionsInitInternal(arg1, arg2)
    ccall((:WebPAnimEncoderOptionsInitInternal, libwebp), Cint, (Ptr{WebPAnimEncoderOptions}, Cint), arg1, arg2)
end

function WebPAnimEncoderOptionsInit(enc_options)
    ccall((:WebPAnimEncoderOptionsInit, libwebp), Cint, (Ptr{WebPAnimEncoderOptions},), enc_options)
end

function WebPAnimEncoderNewInternal(arg1, arg2, arg3, arg4)
    ccall((:WebPAnimEncoderNewInternal, libwebp), Ptr{WebPAnimEncoder}, (Cint, Cint, Ptr{WebPAnimEncoderOptions}, Cint), arg1, arg2, arg3, arg4)
end

function WebPAnimEncoderNew(width, height, enc_options)
    ccall((:WebPAnimEncoderNew, libwebp), Ptr{WebPAnimEncoder}, (Cint, Cint, Ptr{WebPAnimEncoderOptions}), width, height, enc_options)
end

function WebPAnimEncoderAdd(enc, frame, timestamp_ms, config)
    ccall((:WebPAnimEncoderAdd, libwebp), Cint, (Ptr{WebPAnimEncoder}, Ptr{WebPPicture}, Cint, Ptr{WebPConfig}), enc, frame, timestamp_ms, config)
end

function WebPAnimEncoderAssemble(enc, webp_data)
    ccall((:WebPAnimEncoderAssemble, libwebp), Cint, (Ptr{WebPAnimEncoder}, Ptr{WebPData}), enc, webp_data)
end

function WebPAnimEncoderGetError(enc)
    ccall((:WebPAnimEncoderGetError, libwebp), Ptr{Cchar}, (Ptr{WebPAnimEncoder},), enc)
end

function WebPAnimEncoderDelete(enc)
    ccall((:WebPAnimEncoderDelete, libwebp), Cvoid, (Ptr{WebPAnimEncoder},), enc)
end

@cenum WebPFeatureFlags::UInt32 begin
    ANIMATION_FLAG = 2
    XMP_FLAG = 4
    EXIF_FLAG = 8
    ALPHA_FLAG = 16
    ICCP_FLAG = 32
    ALL_VALID_FLAGS = 62
end

function WebPDataInit(webp_data)
    ccall((:WebPDataInit, libwebp), Cvoid, (Ptr{WebPData},), webp_data)
end

function WebPDataClear(webp_data)
    ccall((:WebPDataClear, libwebp), Cvoid, (Ptr{WebPData},), webp_data)
end

function WebPDataCopy(src, dst)
    ccall((:WebPDataCopy, libwebp), Cint, (Ptr{WebPData}, Ptr{WebPData}), src, dst)
end

function WebPMalloc(size)
    ccall((:WebPMalloc, libwebp), Ptr{Cvoid}, (Csize_t,), size)
end

function WebPFree(ptr)
    ccall((:WebPFree, libwebp), Cvoid, (Ptr{Cvoid},), ptr)
end

const WEBP_DECODER_ABI_VERSION = 0x0209

const WEBP_DEMUX_ABI_VERSION = 0x0107

const WEBP_ENCODER_ABI_VERSION = 0x020f

const WEBP_MAX_DIMENSION = 16383

const WEBP_MUX_ABI_VERSION = 0x0108

# Skipping MacroDefinition: WEBP_INLINE inline

# Skipping MacroDefinition: WEBP_EXTERN extern __attribute__ ( ( visibility ( "default" ) ) )

# exports
const PREFIXES = ["WebP"]
for name in names(@__MODULE__; all=true), prefix in PREFIXES
    if startswith(string(name), prefix)
        @eval export $name
    end
end

end # module
