module C

using libsixel_jll
export libsixel_jll

const SIXELSTATUS = Cint

@enum characterSize::UInt32 begin
    CSIZE_7BIT = 0
    CSIZE_8BIT = 1
end

@enum methodForLargest::UInt32 begin
    LARGE_AUTO = 0
    LARGE_NORM = 1
    LARGE_LUM = 2
end

@enum methodForRep::UInt32 begin
    REP_AUTO = 0
    REP_CENTER_BOX = 1
    REP_AVERAGE_COLORS = 2
    REP_AVERAGE_PIXELS = 3
end

@enum methodForDiffuse::UInt32 begin
    DIFFUSE_AUTO = 0
    DIFFUSE_NONE = 1
    DIFFUSE_ATKINSON = 2
    DIFFUSE_FS = 3
    DIFFUSE_JAJUNI = 4
    DIFFUSE_STUCKI = 5
    DIFFUSE_BURKES = 6
    DIFFUSE_A_DITHER = 7
    DIFFUSE_X_DITHER = 8
end

@enum qualityMode::UInt32 begin
    QUALITY_AUTO = 0
    QUALITY_HIGH = 1
    QUALITY_LOW = 2
    QUALITY_FULL = 3
    QUALITY_HIGHCOLOR = 4
end

@enum builtinDither::UInt32 begin
    BUILTIN_MONO_DARK = 0
    BUILTIN_MONO_LIGHT = 1
    BUILTIN_XTERM16 = 2
    BUILTIN_XTERM256 = 3
    BUILTIN_VT340_MONO = 4
    BUILTIN_VT340_COLOR = 5
end

@enum formatType::UInt32 begin
    FORMATTYPE_COLOR = 0
    FORMATTYPE_GRAYSCALE = 64
    FORMATTYPE_PALETTE = 128
end

@enum pixelFormat::UInt32 begin
    PIXELFORMAT_RGB555 = 1
    PIXELFORMAT_RGB565 = 2
    PIXELFORMAT_RGB888 = 3
    PIXELFORMAT_BGR555 = 4
    PIXELFORMAT_BGR565 = 5
    PIXELFORMAT_BGR888 = 6
    PIXELFORMAT_ARGB8888 = 16
    PIXELFORMAT_RGBA8888 = 17
    PIXELFORMAT_G1 = 64
    PIXELFORMAT_G2 = 65
    PIXELFORMAT_G4 = 66
    PIXELFORMAT_G8 = 67
    PIXELFORMAT_AG88 = 83
    PIXELFORMAT_GA88 = 99
    PIXELFORMAT_PAL1 = 128
    PIXELFORMAT_PAL2 = 129
    PIXELFORMAT_PAL4 = 130
    PIXELFORMAT_PAL8 = 131
end

@enum paletteType::UInt32 begin
    PALETTETYPE_AUTO = 0
    PALETTETYPE_HLS = 1
    PALETTETYPE_RGB = 2
end

@enum encodePolicy::UInt32 begin
    ENCODEPOLICY_AUTO = 0
    ENCODEPOLICY_FAST = 1
    ENCODEPOLICY_SIZE = 2
end

@enum methodForResampling::UInt32 begin
    RES_NEAREST = 0
    RES_GAUSSIAN = 1
    RES_HANNING = 2
    RES_HAMMING = 3
    RES_BILINEAR = 4
    RES_WELSH = 5
    RES_BICUBIC = 6
    RES_LANCZOS2 = 7
    RES_LANCZOS3 = 8
    RES_LANCZOS4 = 9
end

# typedef void * ( * sixel_malloc_t ) ( size_t )
const sixel_malloc_t = Ptr{Cvoid}

# typedef void * ( * sixel_calloc_t ) ( size_t , size_t )
const sixel_calloc_t = Ptr{Cvoid}

# typedef void * ( * sixel_realloc_t ) ( void * , size_t )
const sixel_realloc_t = Ptr{Cvoid}

# typedef void ( * sixel_free_t ) ( void * )
const sixel_free_t = Ptr{Cvoid}

mutable struct sixel_allocator end

const sixel_allocator_t = sixel_allocator

function sixel_allocator_new(ppallocator, fn_malloc, fn_calloc, fn_realloc, fn_free)
    ccall((:sixel_allocator_new, libsixel), SIXELSTATUS, (Ptr{Ptr{sixel_allocator_t}}, sixel_malloc_t, sixel_calloc_t, sixel_realloc_t, sixel_free_t), ppallocator, fn_malloc, fn_calloc, fn_realloc, fn_free)
end

function sixel_allocator_ref(allocator)
    ccall((:sixel_allocator_ref, libsixel), Cvoid, (Ptr{sixel_allocator_t},), allocator)
end

function sixel_allocator_unref(allocator)
    ccall((:sixel_allocator_unref, libsixel), Cvoid, (Ptr{sixel_allocator_t},), allocator)
end

function sixel_allocator_malloc(allocator, n)
    ccall((:sixel_allocator_malloc, libsixel), Ptr{Cvoid}, (Ptr{sixel_allocator_t}, Csize_t), allocator, n)
end

function sixel_allocator_calloc(allocator, nelm, elsize)
    ccall((:sixel_allocator_calloc, libsixel), Ptr{Cvoid}, (Ptr{sixel_allocator_t}, Csize_t, Csize_t), allocator, nelm, elsize)
end

function sixel_allocator_realloc(allocator, p, n)
    ccall((:sixel_allocator_realloc, libsixel), Ptr{Cvoid}, (Ptr{sixel_allocator_t}, Ptr{Cvoid}, Csize_t), allocator, p, n)
end

function sixel_allocator_free(allocator, p)
    ccall((:sixel_allocator_free, libsixel), Cvoid, (Ptr{sixel_allocator_t}, Ptr{Cvoid}), allocator, p)
end

mutable struct sixel_output end

const sixel_output_t = sixel_output

# typedef int ( * sixel_write_function ) ( char * data , int size , void * priv )
const sixel_write_function = Ptr{Cvoid}

function sixel_output_new(output, fn_write, priv, allocator)
    ccall((:sixel_output_new, libsixel), SIXELSTATUS, (Ptr{Ptr{sixel_output_t}}, sixel_write_function, Any, Ptr{sixel_allocator_t}), output, fn_write, priv, allocator)
end

function sixel_output_create(fn_write, priv)
    ccall((:sixel_output_create, libsixel), Ptr{sixel_output_t}, (sixel_write_function, Ptr{Cvoid}), fn_write, priv)
end

function sixel_output_destroy(output)
    ccall((:sixel_output_destroy, libsixel), Cvoid, (Ptr{sixel_output_t},), output)
end

function sixel_output_ref(output)
    ccall((:sixel_output_ref, libsixel), Cvoid, (Ptr{sixel_output_t},), output)
end

function sixel_output_unref(output)
    ccall((:sixel_output_unref, libsixel), Cvoid, (Ptr{sixel_output_t},), output)
end

function sixel_output_get_8bit_availability(output)
    ccall((:sixel_output_get_8bit_availability, libsixel), Cint, (Ptr{sixel_output_t},), output)
end

function sixel_output_set_8bit_availability(output, availability)
    ccall((:sixel_output_set_8bit_availability, libsixel), Cvoid, (Ptr{sixel_output_t}, Cint), output, availability)
end

function sixel_output_set_gri_arg_limit(output, value)
    ccall((:sixel_output_set_gri_arg_limit, libsixel), Cvoid, (Ptr{sixel_output_t}, Cint), output, value)
end

function sixel_output_set_penetrate_multiplexer(output, penetrate)
    ccall((:sixel_output_set_penetrate_multiplexer, libsixel), Cvoid, (Ptr{sixel_output_t}, Cint), output, penetrate)
end

function sixel_output_set_skip_dcs_envelope(output, skip)
    ccall((:sixel_output_set_skip_dcs_envelope, libsixel), Cvoid, (Ptr{sixel_output_t}, Cint), output, skip)
end

function sixel_output_set_palette_type(output, palettetype)
    ccall((:sixel_output_set_palette_type, libsixel), Cvoid, (Ptr{sixel_output_t}, Cint), output, palettetype)
end

function sixel_output_set_encode_policy(output, encode_policy)
    ccall((:sixel_output_set_encode_policy, libsixel), Cvoid, (Ptr{sixel_output_t}, Cint), output, encode_policy)
end

mutable struct sixel_dither end

const sixel_dither_t = sixel_dither

function sixel_dither_new(ppdither, ncolors, allocator)
    ccall((:sixel_dither_new, libsixel), SIXELSTATUS, (Ptr{Ptr{sixel_dither_t}}, Cint, Ptr{sixel_allocator_t}), ppdither, ncolors, allocator)
end

function sixel_dither_create(ncolors)
    ccall((:sixel_dither_create, libsixel), Ptr{sixel_dither_t}, (Cint,), ncolors)
end

function sixel_dither_get(builtin_dither)
    ccall((:sixel_dither_get, libsixel), Ptr{sixel_dither_t}, (Cint,), builtin_dither)
end

function sixel_dither_destroy(dither)
    ccall((:sixel_dither_destroy, libsixel), Cvoid, (Ptr{sixel_dither_t},), dither)
end

function sixel_dither_ref(dither)
    ccall((:sixel_dither_ref, libsixel), Cvoid, (Ptr{sixel_dither_t},), dither)
end

function sixel_dither_unref(dither)
    ccall((:sixel_dither_unref, libsixel), Cvoid, (Ptr{sixel_dither_t},), dither)
end

function sixel_dither_initialize(dither, data, width, height, pixelformat, method_for_largest, method_for_rep, quality_mode)
    ccall((:sixel_dither_initialize, libsixel), SIXELSTATUS, (Ptr{sixel_dither_t}, Ptr{Cuchar}, Cint, Cint, Cint, Cint, Cint, Cint), dither, data, width, height, pixelformat, method_for_largest, method_for_rep, quality_mode)
end

function sixel_dither_set_diffusion_type(dither, method_for_diffuse)
    ccall((:sixel_dither_set_diffusion_type, libsixel), Cvoid, (Ptr{sixel_dither_t}, Cint), dither, method_for_diffuse)
end

function sixel_dither_get_num_of_palette_colors(dither)
    ccall((:sixel_dither_get_num_of_palette_colors, libsixel), Cint, (Ptr{sixel_dither_t},), dither)
end

function sixel_dither_get_num_of_histogram_colors(dither)
    ccall((:sixel_dither_get_num_of_histogram_colors, libsixel), Cint, (Ptr{sixel_dither_t},), dither)
end

function sixel_dither_get_num_of_histgram_colors(dither)
    ccall((:sixel_dither_get_num_of_histgram_colors, libsixel), Cint, (Ptr{sixel_dither_t},), dither)
end

function sixel_dither_get_palette(dither)
    ccall((:sixel_dither_get_palette, libsixel), Ptr{Cuchar}, (Ptr{sixel_dither_t},), dither)
end

function sixel_dither_set_palette(dither, palette)
    ccall((:sixel_dither_set_palette, libsixel), Cvoid, (Ptr{sixel_dither_t}, Ptr{Cuchar}), dither, palette)
end

function sixel_dither_set_complexion_score(dither, score)
    ccall((:sixel_dither_set_complexion_score, libsixel), Cvoid, (Ptr{sixel_dither_t}, Cint), dither, score)
end

function sixel_dither_set_body_only(dither, bodyonly)
    ccall((:sixel_dither_set_body_only, libsixel), Cvoid, (Ptr{sixel_dither_t}, Cint), dither, bodyonly)
end

function sixel_dither_set_optimize_palette(dither, do_opt)
    ccall((:sixel_dither_set_optimize_palette, libsixel), Cvoid, (Ptr{sixel_dither_t}, Cint), dither, do_opt)
end

function sixel_dither_set_pixelformat(dither, pixelformat)
    ccall((:sixel_dither_set_pixelformat, libsixel), Cvoid, (Ptr{sixel_dither_t}, Cint), dither, pixelformat)
end

function sixel_dither_set_transparent(dither, transparent)
    ccall((:sixel_dither_set_transparent, libsixel), Cvoid, (Ptr{sixel_dither_t}, Cint), dither, transparent)
end

# typedef void * ( * sixel_allocator_function ) ( size_t size )
const sixel_allocator_function = Ptr{Cvoid}

function sixel_encode(pixels, width, height, depth, dither, context)
    ccall((:sixel_encode, libsixel), SIXELSTATUS, (Ptr{Cuchar}, Cint, Cint, Cint, Ptr{sixel_dither_t}, Ptr{sixel_output_t}), pixels, width, height, depth, dither, context)
end

function sixel_decode_raw(p, len, pixels, pwidth, pheight, palette, ncolors, allocator)
    ccall((:sixel_decode_raw, libsixel), SIXELSTATUS, (Ptr{Cuchar}, Cint, Ptr{Ptr{Cuchar}}, Ptr{Cint}, Ptr{Cint}, Ptr{Ptr{Cuchar}}, Ptr{Cint}, Ptr{sixel_allocator_t}), p, len, pixels, pwidth, pheight, palette, ncolors, allocator)
end

function sixel_decode(sixels, size, pixels, pwidth, pheight, palette, ncolors, fn_malloc)
    ccall((:sixel_decode, libsixel), SIXELSTATUS, (Ptr{Cuchar}, Cint, Ptr{Ptr{Cuchar}}, Ptr{Cint}, Ptr{Cint}, Ptr{Ptr{Cuchar}}, Ptr{Cint}, sixel_allocator_function), sixels, size, pixels, pwidth, pheight, palette, ncolors, fn_malloc)
end

function sixel_helper_set_additional_message(message)
    ccall((:sixel_helper_set_additional_message, libsixel), Cvoid, (Ptr{Cchar},), message)
end

function sixel_helper_get_additional_message()
    ccall((:sixel_helper_get_additional_message, libsixel), Ptr{Cchar}, ())
end

function sixel_helper_format_error(status)
    ccall((:sixel_helper_format_error, libsixel), Ptr{Cchar}, (SIXELSTATUS,), status)
end

function sixel_helper_compute_depth(pixelformat)
    ccall((:sixel_helper_compute_depth, libsixel), Cint, (Cint,), pixelformat)
end

function sixel_helper_normalize_pixelformat(dst, dst_pixelformat, src, src_pixelformat, width, height)
    ccall((:sixel_helper_normalize_pixelformat, libsixel), SIXELSTATUS, (Ptr{Cuchar}, Ptr{Cint}, Ptr{Cuchar}, Cint, Cint, Cint), dst, dst_pixelformat, src, src_pixelformat, width, height)
end

function sixel_helper_scale_image(dst, src, srcw, srch, pixelformat, dstw, dsth, method_for_resampling, allocator)
    ccall((:sixel_helper_scale_image, libsixel), SIXELSTATUS, (Ptr{Cuchar}, Ptr{Cuchar}, Cint, Cint, Cint, Cint, Cint, Cint, Ptr{sixel_allocator_t}), dst, src, srcw, srch, pixelformat, dstw, dsth, method_for_resampling, allocator)
end

@enum imageFormat::UInt32 begin
    FORMAT_GIF = 0
    FORMAT_PNG = 1
    FORMAT_BMP = 2
    FORMAT_JPG = 3
    FORMAT_TGA = 4
    FORMAT_WBMP = 5
    FORMAT_TIFF = 6
    FORMAT_SIXEL = 7
    FORMAT_PNM = 8
    FORMAT_GD2 = 9
    FORMAT_PSD = 10
    FORMAT_HDR = 11
end

@enum loopControl::UInt32 begin
    LOOP_AUTO = 0
    LOOP_FORCE = 1
    LOOP_DISABLE = 2
end

mutable struct sixel_frame end

const sixel_frame_t = sixel_frame

function sixel_frame_new(ppframe, allocator)
    ccall((:sixel_frame_new, libsixel), SIXELSTATUS, (Ptr{Ptr{sixel_frame_t}}, Ptr{sixel_allocator_t}), ppframe, allocator)
end

function sixel_frame_create()
    ccall((:sixel_frame_create, libsixel), Ptr{sixel_frame_t}, ())
end

function sixel_frame_ref(frame)
    ccall((:sixel_frame_ref, libsixel), Cvoid, (Ptr{sixel_frame_t},), frame)
end

function sixel_frame_unref(frame)
    ccall((:sixel_frame_unref, libsixel), Cvoid, (Ptr{sixel_frame_t},), frame)
end

function sixel_frame_init(frame, pixels, width, height, pixelformat, palette, ncolors)
    ccall((:sixel_frame_init, libsixel), SIXELSTATUS, (Ptr{sixel_frame_t}, Ptr{Cuchar}, Cint, Cint, Cint, Ptr{Cuchar}, Cint), frame, pixels, width, height, pixelformat, palette, ncolors)
end

function sixel_frame_get_pixels(frame)
    ccall((:sixel_frame_get_pixels, libsixel), Ptr{Cuchar}, (Ptr{sixel_frame_t},), frame)
end

function sixel_frame_get_palette(frame)
    ccall((:sixel_frame_get_palette, libsixel), Ptr{Cuchar}, (Ptr{sixel_frame_t},), frame)
end

function sixel_frame_get_width(frame)
    ccall((:sixel_frame_get_width, libsixel), Cint, (Ptr{sixel_frame_t},), frame)
end

function sixel_frame_get_height(frame)
    ccall((:sixel_frame_get_height, libsixel), Cint, (Ptr{sixel_frame_t},), frame)
end

function sixel_frame_get_ncolors(frame)
    ccall((:sixel_frame_get_ncolors, libsixel), Cint, (Ptr{sixel_frame_t},), frame)
end

function sixel_frame_get_pixelformat(frame)
    ccall((:sixel_frame_get_pixelformat, libsixel), Cint, (Ptr{sixel_frame_t},), frame)
end

function sixel_frame_get_transparent(frame)
    ccall((:sixel_frame_get_transparent, libsixel), Cint, (Ptr{sixel_frame_t},), frame)
end

function sixel_frame_get_multiframe(frame)
    ccall((:sixel_frame_get_multiframe, libsixel), Cint, (Ptr{sixel_frame_t},), frame)
end

function sixel_frame_get_delay(frame)
    ccall((:sixel_frame_get_delay, libsixel), Cint, (Ptr{sixel_frame_t},), frame)
end

function sixel_frame_get_frame_no(frame)
    ccall((:sixel_frame_get_frame_no, libsixel), Cint, (Ptr{sixel_frame_t},), frame)
end

function sixel_frame_get_loop_no(frame)
    ccall((:sixel_frame_get_loop_no, libsixel), Cint, (Ptr{sixel_frame_t},), frame)
end

function sixel_frame_strip_alpha(frame, bgcolor)
    ccall((:sixel_frame_strip_alpha, libsixel), Cint, (Ptr{sixel_frame_t}, Ptr{Cuchar}), frame, bgcolor)
end

function sixel_frame_resize(frame, width, height, method_for_resampling)
    ccall((:sixel_frame_resize, libsixel), SIXELSTATUS, (Ptr{sixel_frame_t}, Cint, Cint, Cint), frame, width, height, method_for_resampling)
end

function sixel_frame_clip(frame, x, y, width, height)
    ccall((:sixel_frame_clip, libsixel), SIXELSTATUS, (Ptr{sixel_frame_t}, Cint, Cint, Cint, Cint), frame, x, y, width, height)
end

# typedef SIXELSTATUS ( * sixel_load_image_function ) ( sixel_frame_t /* in */ * frame , void /* in/out */ * context )
const sixel_load_image_function = Ptr{Cvoid}

function sixel_helper_load_image_file(filename, fstatic, fuse_palette, reqcolors, bgcolor, loop_control, fn_load, finsecure, cancel_flag, context, allocator)
    ccall((:sixel_helper_load_image_file, libsixel), SIXELSTATUS, (Ptr{Cchar}, Cint, Cint, Cint, Ptr{Cuchar}, Cint, sixel_load_image_function, Cint, Ptr{Cint}, Ptr{Cvoid}, Ptr{sixel_allocator_t}), filename, fstatic, fuse_palette, reqcolors, bgcolor, loop_control, fn_load, finsecure, cancel_flag, context, allocator)
end

function sixel_helper_write_image_file(data, width, height, palette, pixelformat, filename, imageformat, allocator)
    ccall((:sixel_helper_write_image_file, libsixel), SIXELSTATUS, (Ptr{Cuchar}, Cint, Cint, Ptr{Cuchar}, Cint, Ptr{Cchar}, Cint, Ptr{sixel_allocator_t}), data, width, height, palette, pixelformat, filename, imageformat, allocator)
end

mutable struct sixel_encoder end

const sixel_encoder_t = sixel_encoder

function sixel_encoder_new(ppencoder, allocator)
    ccall((:sixel_encoder_new, libsixel), SIXELSTATUS, (Ptr{Ptr{sixel_encoder_t}}, Ptr{sixel_allocator_t}), ppencoder, allocator)
end

function sixel_encoder_create()
    ccall((:sixel_encoder_create, libsixel), Ptr{sixel_encoder_t}, ())
end

function sixel_encoder_ref(encoder)
    ccall((:sixel_encoder_ref, libsixel), Cvoid, (Ptr{sixel_encoder_t},), encoder)
end

function sixel_encoder_unref(encoder)
    ccall((:sixel_encoder_unref, libsixel), Cvoid, (Ptr{sixel_encoder_t},), encoder)
end

function sixel_encoder_set_cancel_flag(encoder, cancel_flag)
    ccall((:sixel_encoder_set_cancel_flag, libsixel), SIXELSTATUS, (Ptr{sixel_encoder_t}, Ptr{Cint}), encoder, cancel_flag)
end

function sixel_encoder_setopt(encoder, arg, optarg)
    ccall((:sixel_encoder_setopt, libsixel), SIXELSTATUS, (Ptr{sixel_encoder_t}, Cint, Ptr{Cchar}), encoder, arg, optarg)
end

function sixel_encoder_encode(encoder, filename)
    ccall((:sixel_encoder_encode, libsixel), SIXELSTATUS, (Ptr{sixel_encoder_t}, Ptr{Cchar}), encoder, filename)
end

function sixel_encoder_encode_bytes(encoder, bytes, width, height, pixelformat, palette, ncolors)
    ccall((:sixel_encoder_encode_bytes, libsixel), SIXELSTATUS, (Ptr{sixel_encoder_t}, Ptr{Cuchar}, Cint, Cint, Cint, Ptr{Cuchar}, Cint), encoder, bytes, width, height, pixelformat, palette, ncolors)
end

mutable struct sixel_decoder end

const sixel_decoder_t = sixel_decoder

function sixel_decoder_new(ppdecoder, allocator)
    ccall((:sixel_decoder_new, libsixel), SIXELSTATUS, (Ptr{Ptr{sixel_decoder_t}}, Ptr{sixel_allocator_t}), ppdecoder, allocator)
end

function sixel_decoder_create()
    ccall((:sixel_decoder_create, libsixel), Ptr{sixel_decoder_t}, ())
end

function sixel_decoder_ref(decoder)
    ccall((:sixel_decoder_ref, libsixel), Cvoid, (Ptr{sixel_decoder_t},), decoder)
end

function sixel_decoder_unref(decoder)
    ccall((:sixel_decoder_unref, libsixel), Cvoid, (Ptr{sixel_decoder_t},), decoder)
end

function sixel_decoder_setopt(decoder, arg, optarg)
    ccall((:sixel_decoder_setopt, libsixel), SIXELSTATUS, (Ptr{sixel_decoder_t}, Cint, Ptr{Cchar}), decoder, arg, optarg)
end

function sixel_decoder_decode(decoder)
    ccall((:sixel_decoder_decode, libsixel), SIXELSTATUS, (Ptr{sixel_decoder_t},), decoder)
end

const LIBSIXEL_VERSION = "1.8.1"

const LIBSIXEL_ABI_VERSION = "1:6:0"

const SIXEL_OUTPUT_PACKET_SIZE = 16384

const SIXEL_PALETTE_MIN = 2

const SIXEL_PALETTE_MAX = 256

const SIXEL_USE_DEPRECATED_SYMBOLS = 1

const SIXEL_OK = 0x0000

const SIXEL_FALSE = 0x1000

const SIXEL_RUNTIME_ERROR = SIXEL_FALSE | 0x0100

const SIXEL_LOGIC_ERROR = SIXEL_FALSE | 0x0200

const SIXEL_FEATURE_ERROR = SIXEL_FALSE | 0x0300

const SIXEL_LIBC_ERROR = SIXEL_FALSE | 0x0400

const SIXEL_CURL_ERROR = SIXEL_FALSE | 0x0500

const SIXEL_JPEG_ERROR = SIXEL_FALSE | 0x0600

const SIXEL_PNG_ERROR = SIXEL_FALSE | 0x0700

const SIXEL_GDK_ERROR = SIXEL_FALSE | 0x0800

const SIXEL_GD_ERROR = SIXEL_FALSE | 0x0900

const SIXEL_STBI_ERROR = SIXEL_FALSE | 0x0a00

const SIXEL_STBIW_ERROR = SIXEL_FALSE | 0x0b00

const SIXEL_INTERRUPTED = SIXEL_OK | 0x0001

const SIXEL_BAD_ALLOCATION = SIXEL_RUNTIME_ERROR | 0x0001

const SIXEL_BAD_ARGUMENT = SIXEL_RUNTIME_ERROR | 0x0002

const SIXEL_BAD_INPUT = SIXEL_RUNTIME_ERROR | 0x0003

const SIXEL_NOT_IMPLEMENTED = SIXEL_FEATURE_ERROR | 0x0001

const SIXEL_LARGE_AUTO = 0x00

const SIXEL_LARGE_NORM = 0x01

const SIXEL_LARGE_LUM = 0x02

const SIXEL_REP_AUTO = 0x00

const SIXEL_REP_CENTER_BOX = 0x01

const SIXEL_REP_AVERAGE_COLORS = 0x02

const SIXEL_REP_AVERAGE_PIXELS = 0x03

const SIXEL_DIFFUSE_AUTO = 0x00

const SIXEL_DIFFUSE_NONE = 0x01

const SIXEL_DIFFUSE_ATKINSON = 0x02

const SIXEL_DIFFUSE_FS = 0x03

const SIXEL_DIFFUSE_JAJUNI = 0x04

const SIXEL_DIFFUSE_STUCKI = 0x05

const SIXEL_DIFFUSE_BURKES = 0x06

const SIXEL_DIFFUSE_A_DITHER = 0x07

const SIXEL_DIFFUSE_X_DITHER = 0x08

const SIXEL_QUALITY_AUTO = 0x00

const SIXEL_QUALITY_HIGH = 0x01

const SIXEL_QUALITY_LOW = 0x02

const SIXEL_QUALITY_FULL = 0x03

const SIXEL_QUALITY_HIGHCOLOR = 0x04

const SIXEL_BUILTIN_MONO_DARK = 0x00

const SIXEL_BUILTIN_MONO_LIGHT = 0x01

const SIXEL_BUILTIN_XTERM16 = 0x02

const SIXEL_BUILTIN_XTERM256 = 0x03

const SIXEL_BUILTIN_VT340_MONO = 0x04

const SIXEL_BUILTIN_VT340_COLOR = 0x05

const SIXEL_BUILTIN_G1 = 0x06

const SIXEL_BUILTIN_G2 = 0x07

const SIXEL_BUILTIN_G4 = 0x08

const SIXEL_BUILTIN_G8 = 0x09

const SIXEL_FORMATTYPE_COLOR = 0

const SIXEL_FORMATTYPE_GRAYSCALE = 1 << 6

const SIXEL_FORMATTYPE_PALETTE = 1 << 7

const SIXEL_PIXELFORMAT_RGB555 = SIXEL_FORMATTYPE_COLOR | 0x01

const SIXEL_PIXELFORMAT_RGB565 = SIXEL_FORMATTYPE_COLOR | 0x02

const SIXEL_PIXELFORMAT_RGB888 = SIXEL_FORMATTYPE_COLOR | 0x03

const SIXEL_PIXELFORMAT_BGR555 = SIXEL_FORMATTYPE_COLOR | 0x04

const SIXEL_PIXELFORMAT_BGR565 = SIXEL_FORMATTYPE_COLOR | 0x05

const SIXEL_PIXELFORMAT_BGR888 = SIXEL_FORMATTYPE_COLOR | 0x06

const SIXEL_PIXELFORMAT_ARGB8888 = SIXEL_FORMATTYPE_COLOR | 0x10

const SIXEL_PIXELFORMAT_RGBA8888 = SIXEL_FORMATTYPE_COLOR | 0x11

const SIXEL_PIXELFORMAT_ABGR8888 = SIXEL_FORMATTYPE_COLOR | 0x12

const SIXEL_PIXELFORMAT_BGRA8888 = SIXEL_FORMATTYPE_COLOR | 0x13

const SIXEL_PIXELFORMAT_G1 = SIXEL_FORMATTYPE_GRAYSCALE | 0x00

const SIXEL_PIXELFORMAT_G2 = SIXEL_FORMATTYPE_GRAYSCALE | 0x01

const SIXEL_PIXELFORMAT_G4 = SIXEL_FORMATTYPE_GRAYSCALE | 0x02

const SIXEL_PIXELFORMAT_G8 = SIXEL_FORMATTYPE_GRAYSCALE | 0x03

const SIXEL_PIXELFORMAT_AG88 = SIXEL_FORMATTYPE_GRAYSCALE | 0x13

const SIXEL_PIXELFORMAT_GA88 = SIXEL_FORMATTYPE_GRAYSCALE | 0x23

const SIXEL_PIXELFORMAT_PAL1 = SIXEL_FORMATTYPE_PALETTE | 0x00

const SIXEL_PIXELFORMAT_PAL2 = SIXEL_FORMATTYPE_PALETTE | 0x01

const SIXEL_PIXELFORMAT_PAL4 = SIXEL_FORMATTYPE_PALETTE | 0x02

const SIXEL_PIXELFORMAT_PAL8 = SIXEL_FORMATTYPE_PALETTE | 0x03

const SIXEL_PALETTETYPE_AUTO = 0

const SIXEL_PALETTETYPE_HLS = 1

const SIXEL_PALETTETYPE_RGB = 2

const SIXEL_ENCODEPOLICY_AUTO = 0

const SIXEL_ENCODEPOLICY_FAST = 1

const SIXEL_ENCODEPOLICY_SIZE = 2

const SIXEL_RES_NEAREST = 0

const SIXEL_RES_GAUSSIAN = 1

const SIXEL_RES_HANNING = 2

const SIXEL_RES_HAMMING = 3

const SIXEL_RES_BILINEAR = 4

const SIXEL_RES_WELSH = 5

const SIXEL_RES_BICUBIC = 6

const SIXEL_RES_LANCZOS2 = 7

const SIXEL_RES_LANCZOS3 = 8

const SIXEL_RES_LANCZOS4 = 9

const SIXEL_FORMAT_GIF = 0x00

const SIXEL_FORMAT_PNG = 0x01

const SIXEL_FORMAT_BMP = 0x02

const SIXEL_FORMAT_JPG = 0x03

const SIXEL_FORMAT_TGA = 0x04

const SIXEL_FORMAT_WBMP = 0x05

const SIXEL_FORMAT_TIFF = 0x06

const SIXEL_FORMAT_SIXEL = 0x07

const SIXEL_FORMAT_PNM = 0x08

const SIXEL_FORMAT_GD2 = 0x09

const SIXEL_FORMAT_PSD = 0x0a

const SIXEL_FORMAT_HDR = 0x0b

const SIXEL_LOOP_AUTO = 0

const SIXEL_LOOP_FORCE = 1

const SIXEL_LOOP_DISABLE = 2

const SIXEL_OPTFLAG_INPUT = Cchar('i')

const SIXEL_OPTFLAG_OUTPUT = Cchar('o')

const SIXEL_OPTFLAG_OUTFILE = Cchar('o')

const SIXEL_OPTFLAG_7BIT_MODE = Cchar('7')

const SIXEL_OPTFLAG_8BIT_MODE = Cchar('8')

const SIXEL_OPTFLAG_HAS_GRI_ARG_LIMIT = Cchar('R')

const SIXEL_OPTFLAG_COLORS = Cchar('p')

const SIXEL_OPTFLAG_MAPFILE = Cchar('m')

const SIXEL_OPTFLAG_MONOCHROME = Cchar('e')

const SIXEL_OPTFLAG_INSECURE = Cchar('k')

const SIXEL_OPTFLAG_INVERT = Cchar('i')

const SIXEL_OPTFLAG_HIGH_COLOR = Cchar('I')

const SIXEL_OPTFLAG_USE_MACRO = Cchar('u')

const SIXEL_OPTFLAG_MACRO_NUMBER = Cchar('n')

const SIXEL_OPTFLAG_COMPLEXION_SCORE = Cchar('C')

const SIXEL_OPTFLAG_IGNORE_DELAY = Cchar('g')

const SIXEL_OPTFLAG_STATIC = Cchar('S')

const SIXEL_OPTFLAG_DIFFUSION = Cchar('d')

const SIXEL_OPTFLAG_FIND_LARGEST = Cchar('f')

const SIXEL_OPTFLAG_SELECT_COLOR = Cchar('s')

const SIXEL_OPTFLAG_CROP = Cchar('c')

const SIXEL_OPTFLAG_WIDTH = Cchar('w')

const SIXEL_OPTFLAG_HEIGHT = Cchar('h')

const SIXEL_OPTFLAG_RESAMPLING = Cchar('r')

const SIXEL_OPTFLAG_QUALITY = Cchar('q')

const SIXEL_OPTFLAG_LOOPMODE = Cchar('l')

const SIXEL_OPTFLAG_PALETTE_TYPE = Cchar('t')

const SIXEL_OPTFLAG_BUILTIN_PALETTE = Cchar('b')

const SIXEL_OPTFLAG_ENCODE_POLICY = Cchar('E')

const SIXEL_OPTFLAG_BGCOLOR = Cchar('B')

const SIXEL_OPTFLAG_PENETRATE = Cchar('P')

const SIXEL_OPTFLAG_PIPE_MODE = Cchar('D')

const SIXEL_OPTFLAG_VERBOSE = Cchar('v')

const SIXEL_OPTFLAG_VERSION = Cchar('V')

const SIXEL_OPTFLAG_HELP = Cchar('H')

end # module
