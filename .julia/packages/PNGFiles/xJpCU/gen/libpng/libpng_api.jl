##############################################
# This file has been automatically generated #
# Do not edit !                              #
##############################################

# TODO: pending https://github.com/JuliaLang/julia/issues/29420

const PNGCAPI = nothing
const PNG_BYTES_TO_CHECK = 8
const tm = Libc.TmStruct
const time_t = UInt
const jmp_buf = Nothing

const png_byte = Cuchar

const png_int_16 = Cshort

const png_uint_16 = Cushort

const png_int_32 = Cint

const png_uint_32 = Cuint

const png_size_t = Csize_t

const png_ptrdiff_t = Cptrdiff_t

const png_alloc_size_t = Csize_t

const png_fixed_point = png_int_32

const png_voidp = Ptr{Cvoid}

const png_const_voidp = Ptr{Cvoid}

const png_bytep = Ptr{png_byte}

const png_const_bytep = Ptr{png_byte}

const png_uint_32p = Ptr{png_uint_32}

const png_const_uint_32p = Ptr{png_uint_32}

const png_int_32p = Ptr{png_int_32}

const png_const_int_32p = Ptr{png_int_32}

const png_uint_16p = Ptr{png_uint_16}

const png_const_uint_16p = Ptr{png_uint_16}

const png_int_16p = Ptr{png_int_16}

const png_const_int_16p = Ptr{png_int_16}

const png_charp = Cstring

const png_const_charp = Cstring

const png_fixed_point_p = Ptr{png_fixed_point}

const png_const_fixed_point_p = Ptr{png_fixed_point}

const png_size_tp = Ptr{Csize_t}

const png_const_size_tp = Ptr{Csize_t}

const png_FILE_p = Ptr{Libc.FILE}

const png_doublep = Ptr{Cdouble}

const png_const_doublep = Ptr{Cdouble}

const png_bytepp = Ptr{Ptr{png_byte}}

const png_uint_32pp = Ptr{Ptr{png_uint_32}}

const png_int_32pp = Ptr{Ptr{png_int_32}}

const png_uint_16pp = Ptr{Ptr{png_uint_16}}

const png_int_16pp = Ptr{Ptr{png_int_16}}

const png_const_charpp = Ptr{Cstring}

const png_charpp = Ptr{Cstring}

const png_fixed_point_pp = Ptr{Ptr{png_fixed_point}}

const png_doublepp = Ptr{Ptr{Cdouble}}

const png_charppp = Ptr{Ptr{Cstring}}

const png_libpng_version_1_6_37 = Cstring

const png_struct_def = Cvoid

const png_struct = png_struct_def

const png_const_structp = Ptr{png_struct}

const png_structp = Ptr{png_struct}

const png_structpp = Ptr{Ptr{png_struct}}

const png_info_def = Cvoid

const png_info = png_info_def

const png_infop = Ptr{png_info}

const png_const_infop = Ptr{png_info}

const png_infopp = Ptr{Ptr{png_info}}

const png_structrp = Ptr{png_struct}

const png_const_structrp = Ptr{png_struct}

const png_inforp = Ptr{png_info}

const png_const_inforp = Ptr{png_info}

struct png_color_struct
    red::png_byte
    green::png_byte
    blue::png_byte
end

const png_color = png_color_struct

const png_colorp = Ptr{png_color}

const png_const_colorp = Ptr{png_color}

const png_colorpp = Ptr{Ptr{png_color}}

struct png_color_16_struct
    index::png_byte
    red::png_uint_16
    green::png_uint_16
    blue::png_uint_16
    gray::png_uint_16
end

const png_color_16 = png_color_16_struct

const png_color_16p = Ptr{png_color_16}

const png_const_color_16p = Ptr{png_color_16}

const png_color_16pp = Ptr{Ptr{png_color_16}}

struct png_color_8_struct
    red::png_byte
    green::png_byte
    blue::png_byte
    gray::png_byte
    alpha::png_byte
end

const png_color_8 = png_color_8_struct

const png_color_8p = Ptr{png_color_8}

const png_const_color_8p = Ptr{png_color_8}

const png_color_8pp = Ptr{Ptr{png_color_8}}

struct png_sPLT_entry_struct
    red::png_uint_16
    green::png_uint_16
    blue::png_uint_16
    alpha::png_uint_16
    frequency::png_uint_16
end

const png_sPLT_entry = png_sPLT_entry_struct

const png_sPLT_entryp = Ptr{png_sPLT_entry}

const png_const_sPLT_entryp = Ptr{png_sPLT_entry}

const png_sPLT_entrypp = Ptr{Ptr{png_sPLT_entry}}

struct png_sPLT_struct
    name::png_charp
    depth::png_byte
    entries::png_sPLT_entryp
    nentries::png_int_32
end

const png_sPLT_t = png_sPLT_struct

const png_sPLT_tp = Ptr{png_sPLT_t}

const png_const_sPLT_tp = Ptr{png_sPLT_t}

const png_sPLT_tpp = Ptr{Ptr{png_sPLT_t}}

struct png_text_struct
    compression::Cint
    key::png_charp
    text::png_charp
    text_length::Csize_t
    itxt_length::Csize_t
    lang::png_charp
    lang_key::png_charp
end

const png_text = png_text_struct

const png_textp = Ptr{png_text}

const png_const_textp = Ptr{png_text}

const png_textpp = Ptr{Ptr{png_text}}

struct png_time_struct
    year::png_uint_16
    month::png_byte
    day::png_byte
    hour::png_byte
    minute::png_byte
    second::png_byte
end

const png_time = png_time_struct

const png_timep = Ptr{png_time}

const png_const_timep = Ptr{png_time}

const png_timepp = Ptr{Ptr{png_time}}

struct png_unknown_chunk_t
    name::NTuple{5, png_byte}
    data::Ptr{png_byte}
    size::Csize_t
    location::png_byte
end

const png_unknown_chunk = png_unknown_chunk_t

const png_unknown_chunkp = Ptr{png_unknown_chunk}

const png_const_unknown_chunkp = Ptr{png_unknown_chunk}

const png_unknown_chunkpp = Ptr{Ptr{png_unknown_chunk}}

struct png_row_info_struct
    width::png_uint_32
    rowbytes::Csize_t
    color_type::png_byte
    bit_depth::png_byte
    channels::png_byte
    pixel_depth::png_byte
end

const png_row_info = png_row_info_struct

const png_row_infop = Ptr{png_row_info}

const png_row_infopp = Ptr{Ptr{png_row_info}}

# typedef PNG_CALLBACK
const png_error_ptr = Ptr{Cvoid}

# typedef PNG_CALLBACK
const png_rw_ptr = Ptr{Cvoid}

# typedef PNG_CALLBACK
const png_flush_ptr = Ptr{Cvoid}

# typedef PNG_CALLBACK
const png_read_status_ptr = Ptr{Cvoid}

# typedef PNG_CALLBACK
const png_write_status_ptr = Ptr{Cvoid}

# typedef PNG_CALLBACK
const png_progressive_info_ptr = Ptr{Cvoid}

# typedef PNG_CALLBACK
const png_progressive_end_ptr = Ptr{Cvoid}

# typedef PNG_CALLBACK
const png_progressive_row_ptr = Ptr{Cvoid}

# typedef PNG_CALLBACK
const png_user_transform_ptr = Ptr{Cvoid}

# typedef PNG_CALLBACK
const png_user_chunk_ptr = Ptr{Cvoid}

# PNG_FUNCTION
const png_longjmp_ptr = Ptr{Cvoid}

# typedef PNG_CALLBACK
const png_malloc_ptr = Ptr{Cvoid}

# typedef PNG_CALLBACK
const png_free_ptr = Ptr{Cvoid}

function png_access_version_number()
    ccall((:png_access_version_number, libpng), png_uint_32, ())
end

function png_set_sig_bytes(png_ptr, num_bytes)
    ccall((:png_set_sig_bytes, libpng), Cvoid, (png_structrp, Cint), png_ptr, num_bytes)
end

function png_sig_cmp(sig, start, num_to_check)
    ccall((:png_sig_cmp, libpng), Cint, (png_const_bytep, Csize_t, Csize_t), sig, start, num_to_check)
end

function png_create_read_struct(user_png_ver, error_ptr, error_fn, warn_fn)
    ccall((:png_create_read_struct, libpng), png_structp, (png_const_charp, png_voidp, png_error_ptr, png_error_ptr), user_png_ver, error_ptr, error_fn, warn_fn)
end

function png_create_write_struct(user_png_ver, error_ptr, error_fn, warn_fn)
    ccall((:png_create_write_struct, libpng), png_structp, (png_const_charp, png_voidp, png_error_ptr, png_error_ptr), user_png_ver, error_ptr, error_fn, warn_fn)
end

function png_get_compression_buffer_size(png_ptr)
    ccall((:png_get_compression_buffer_size, libpng), Csize_t, (png_const_structrp,), png_ptr)
end

function png_set_compression_buffer_size(png_ptr, size)
    ccall((:png_set_compression_buffer_size, libpng), Cvoid, (png_structrp, Csize_t), png_ptr, size)
end

function png_set_longjmp_fn(png_ptr, longjmp_fn, jmp_buf_size)
    ccall((:png_set_longjmp_fn, libpng), Ptr{jmp_buf}, (png_structrp, png_longjmp_ptr, Csize_t), png_ptr, longjmp_fn, jmp_buf_size)
end

function png_longjmp(png_ptr, val)
    ccall((:png_longjmp, libpng), Cvoid, (png_const_structrp, Cint), png_ptr, val)
end

function png_reset_zstream(png_ptr)
    ccall((:png_reset_zstream, libpng), Cint, (png_structrp,), png_ptr)
end

function png_create_read_struct_2(user_png_ver, error_ptr, error_fn, warn_fn, mem_ptr, malloc_fn, free_fn)
    ccall((:png_create_read_struct_2, libpng), png_structp, (png_const_charp, png_voidp, png_error_ptr, png_error_ptr, png_voidp, png_malloc_ptr, png_free_ptr), user_png_ver, error_ptr, error_fn, warn_fn, mem_ptr, malloc_fn, free_fn)
end

function png_create_write_struct_2(user_png_ver, error_ptr, error_fn, warn_fn, mem_ptr, malloc_fn, free_fn)
    ccall((:png_create_write_struct_2, libpng), png_structp, (png_const_charp, png_voidp, png_error_ptr, png_error_ptr, png_voidp, png_malloc_ptr, png_free_ptr), user_png_ver, error_ptr, error_fn, warn_fn, mem_ptr, malloc_fn, free_fn)
end

function png_write_sig(png_ptr)
    ccall((:png_write_sig, libpng), Cvoid, (png_structrp,), png_ptr)
end

function png_write_chunk(png_ptr, chunk_name, data, length)
    ccall((:png_write_chunk, libpng), Cvoid, (png_structrp, png_const_bytep, png_const_bytep, Csize_t), png_ptr, chunk_name, data, length)
end

function png_write_chunk_start(png_ptr, chunk_name, length)
    ccall((:png_write_chunk_start, libpng), Cvoid, (png_structrp, png_const_bytep, png_uint_32), png_ptr, chunk_name, length)
end

function png_write_chunk_data(png_ptr, data, length)
    ccall((:png_write_chunk_data, libpng), Cvoid, (png_structrp, png_const_bytep, Csize_t), png_ptr, data, length)
end

function png_write_chunk_end(png_ptr)
    ccall((:png_write_chunk_end, libpng), Cvoid, (png_structrp,), png_ptr)
end

function png_create_info_struct(png_ptr)
    ccall((:png_create_info_struct, libpng), png_infop, (png_const_structrp,), png_ptr)
end

function png_info_init_3(info_ptr, png_info_struct_size)
    ccall((:png_info_init_3, libpng), Cvoid, (png_infopp, Csize_t), info_ptr, png_info_struct_size)
end

function png_write_info_before_PLTE(png_ptr, info_ptr)
    ccall((:png_write_info_before_PLTE, libpng), Cvoid, (png_structrp, png_const_inforp), png_ptr, info_ptr)
end

function png_write_info(png_ptr, info_ptr)
    ccall((:png_write_info, libpng), Cvoid, (png_structrp, png_const_inforp), png_ptr, info_ptr)
end

function png_read_info(png_ptr, info_ptr)
    ccall((:png_read_info, libpng), Cvoid, (png_structrp, png_inforp), png_ptr, info_ptr)
end

function png_convert_to_rfc1123(png_ptr, ptime)
    ccall((:png_convert_to_rfc1123, libpng), png_const_charp, (png_structrp, png_const_timep), png_ptr, ptime)
end

function png_convert_to_rfc1123_buffer(out, ptime)
    ccall((:png_convert_to_rfc1123_buffer, libpng), Cint, (Ptr{Cchar}, png_const_timep), out, ptime)
end

function png_convert_from_struct_tm(ptime, ttime)
    ccall((:png_convert_from_struct_tm, libpng), Cvoid, (png_timep, Ptr{tm}), ptime, ttime)
end

function png_convert_from_time_t(ptime, ttime)
    ccall((:png_convert_from_time_t, libpng), Cvoid, (png_timep, time_t), ptime, ttime)
end

function png_set_expand(png_ptr)
    ccall((:png_set_expand, libpng), Cvoid, (png_structrp,), png_ptr)
end

function png_set_expand_gray_1_2_4_to_8(png_ptr)
    ccall((:png_set_expand_gray_1_2_4_to_8, libpng), Cvoid, (png_structrp,), png_ptr)
end

function png_set_palette_to_rgb(png_ptr)
    ccall((:png_set_palette_to_rgb, libpng), Cvoid, (png_structrp,), png_ptr)
end

function png_set_tRNS_to_alpha(png_ptr)
    ccall((:png_set_tRNS_to_alpha, libpng), Cvoid, (png_structrp,), png_ptr)
end

function png_set_expand_16(png_ptr)
    ccall((:png_set_expand_16, libpng), Cvoid, (png_structrp,), png_ptr)
end

function png_set_bgr(png_ptr)
    ccall((:png_set_bgr, libpng), Cvoid, (png_structrp,), png_ptr)
end

function png_set_gray_to_rgb(png_ptr)
    ccall((:png_set_gray_to_rgb, libpng), Cvoid, (png_structrp,), png_ptr)
end

function png_set_rgb_to_gray(png_ptr, error_action, red, green)
    ccall((:png_set_rgb_to_gray, libpng), Cvoid, (png_structrp, Cint, Cdouble, Cdouble), png_ptr, error_action, red, green)
end

function png_set_rgb_to_gray_fixed(png_ptr, error_action, red, green)
    ccall((:png_set_rgb_to_gray_fixed, libpng), Cvoid, (png_structrp, Cint, png_fixed_point, png_fixed_point), png_ptr, error_action, red, green)
end

function png_get_rgb_to_gray_status(png_ptr)
    ccall((:png_get_rgb_to_gray_status, libpng), png_byte, (png_const_structrp,), png_ptr)
end

function png_build_grayscale_palette(bit_depth, palette)
    ccall((:png_build_grayscale_palette, libpng), Cvoid, (Cint, png_colorp), bit_depth, palette)
end

function png_set_alpha_mode(png_ptr, mode, output_gamma)
    ccall((:png_set_alpha_mode, libpng), Cvoid, (png_structrp, Cint, Cdouble), png_ptr, mode, output_gamma)
end

function png_set_alpha_mode_fixed(png_ptr, mode, output_gamma)
    ccall((:png_set_alpha_mode_fixed, libpng), Cvoid, (png_structrp, Cint, png_fixed_point), png_ptr, mode, output_gamma)
end

function png_set_strip_alpha(png_ptr)
    ccall((:png_set_strip_alpha, libpng), Cvoid, (png_structrp,), png_ptr)
end

function png_set_swap_alpha(png_ptr)
    ccall((:png_set_swap_alpha, libpng), Cvoid, (png_structrp,), png_ptr)
end

function png_set_invert_alpha(png_ptr)
    ccall((:png_set_invert_alpha, libpng), Cvoid, (png_structrp,), png_ptr)
end

function png_set_filler(png_ptr, filler, flags)
    ccall((:png_set_filler, libpng), Cvoid, (png_structrp, png_uint_32, Cint), png_ptr, filler, flags)
end

function png_set_add_alpha(png_ptr, filler, flags)
    ccall((:png_set_add_alpha, libpng), Cvoid, (png_structrp, png_uint_32, Cint), png_ptr, filler, flags)
end

function png_set_swap(png_ptr)
    ccall((:png_set_swap, libpng), Cvoid, (png_structrp,), png_ptr)
end

function png_set_packing(png_ptr)
    ccall((:png_set_packing, libpng), Cvoid, (png_structrp,), png_ptr)
end

function png_set_packswap(png_ptr)
    ccall((:png_set_packswap, libpng), Cvoid, (png_structrp,), png_ptr)
end

function png_set_shift(png_ptr, true_bits)
    ccall((:png_set_shift, libpng), Cvoid, (png_structrp, png_const_color_8p), png_ptr, true_bits)
end

function png_set_interlace_handling(png_ptr)
    ccall((:png_set_interlace_handling, libpng), Cint, (png_structrp,), png_ptr)
end

function png_set_invert_mono(png_ptr)
    ccall((:png_set_invert_mono, libpng), Cvoid, (png_structrp,), png_ptr)
end

function png_set_background(png_ptr, background_color, background_gamma_code, need_expand, background_gamma)
    ccall((:png_set_background, libpng), Cvoid, (png_structrp, png_const_color_16p, Cint, Cint, Cdouble), png_ptr, background_color, background_gamma_code, need_expand, background_gamma)
end

function png_set_background_fixed(png_ptr, background_color, background_gamma_code, need_expand, background_gamma)
    ccall((:png_set_background_fixed, libpng), Cvoid, (png_structrp, png_const_color_16p, Cint, Cint, png_fixed_point), png_ptr, background_color, background_gamma_code, need_expand, background_gamma)
end

function png_set_scale_16(png_ptr)
    ccall((:png_set_scale_16, libpng), Cvoid, (png_structrp,), png_ptr)
end

function png_set_strip_16(png_ptr)
    ccall((:png_set_strip_16, libpng), Cvoid, (png_structrp,), png_ptr)
end

function png_set_quantize(png_ptr, palette, num_palette, maximum_colors, histogram, full_quantize)
    ccall((:png_set_quantize, libpng), Cvoid, (png_structrp, png_colorp, Cint, Cint, png_const_uint_16p, Cint), png_ptr, palette, num_palette, maximum_colors, histogram, full_quantize)
end

function png_set_gamma(png_ptr, screen_gamma, override_file_gamma)
    ccall((:png_set_gamma, libpng), Cvoid, (png_structrp, Cdouble, Cdouble), png_ptr, screen_gamma, override_file_gamma)
end

function png_set_gamma_fixed(png_ptr, screen_gamma, override_file_gamma)
    ccall((:png_set_gamma_fixed, libpng), Cvoid, (png_structrp, png_fixed_point, png_fixed_point), png_ptr, screen_gamma, override_file_gamma)
end

function png_set_flush(png_ptr, nrows)
    ccall((:png_set_flush, libpng), Cvoid, (png_structrp, Cint), png_ptr, nrows)
end

function png_write_flush(png_ptr)
    ccall((:png_write_flush, libpng), Cvoid, (png_structrp,), png_ptr)
end

function png_start_read_image(png_ptr)
    ccall((:png_start_read_image, libpng), Cvoid, (png_structrp,), png_ptr)
end

function png_read_update_info(png_ptr, info_ptr)
    ccall((:png_read_update_info, libpng), Cvoid, (png_structrp, png_inforp), png_ptr, info_ptr)
end

function png_read_rows(png_ptr, row, display_row, num_rows)
    ccall((:png_read_rows, libpng), Cvoid, (png_structrp, png_bytepp, png_bytepp, png_uint_32), png_ptr, row, display_row, num_rows)
end

function png_read_row(png_ptr, row, display_row)
    ccall((:png_read_row, libpng), Cvoid, (png_structrp, png_bytep, png_bytep), png_ptr, row, display_row)
end

function png_read_image(png_ptr, image)
    ccall((:png_read_image, libpng), Cvoid, (png_structrp, png_bytepp), png_ptr, image)
end

function png_write_row(png_ptr, row)
    ccall((:png_write_row, libpng), Cvoid, (png_structrp, png_const_bytep), png_ptr, row)
end

function png_write_rows(png_ptr, row, num_rows)
    ccall((:png_write_rows, libpng), Cvoid, (png_structrp, png_bytepp, png_uint_32), png_ptr, row, num_rows)
end

function png_write_image(png_ptr, image)
    ccall((:png_write_image, libpng), Cvoid, (png_structrp, png_bytepp), png_ptr, image)
end

function png_write_end(png_ptr, info_ptr)
    ccall((:png_write_end, libpng), Cvoid, (png_structrp, png_inforp), png_ptr, info_ptr)
end

function png_read_end(png_ptr, info_ptr)
    ccall((:png_read_end, libpng), Cvoid, (png_structrp, png_inforp), png_ptr, info_ptr)
end

function png_destroy_info_struct(png_ptr, info_ptr_ptr)
    ccall((:png_destroy_info_struct, libpng), Cvoid, (png_const_structrp, png_infopp), png_ptr, info_ptr_ptr)
end

function png_destroy_read_struct(png_ptr_ptr, info_ptr_ptr, end_info_ptr_ptr)
    ccall((:png_destroy_read_struct, libpng), Cvoid, (png_structpp, png_infopp, png_infopp), png_ptr_ptr, info_ptr_ptr, end_info_ptr_ptr)
end

function png_destroy_write_struct(png_ptr_ptr, info_ptr_ptr)
    ccall((:png_destroy_write_struct, libpng), Cvoid, (png_structpp, png_infopp), png_ptr_ptr, info_ptr_ptr)
end

function png_set_crc_action(png_ptr, crit_action, ancil_action)
    ccall((:png_set_crc_action, libpng), Cvoid, (png_structrp, Cint, Cint), png_ptr, crit_action, ancil_action)
end

function png_set_filter(png_ptr, method, filters)
    ccall((:png_set_filter, libpng), Cvoid, (png_structrp, Cint, Cint), png_ptr, method, filters)
end

function png_set_filter_heuristics(png_ptr, heuristic_method, num_weights, filter_weights, filter_costs)
    ccall((:png_set_filter_heuristics, libpng), Cvoid, (png_structrp, Cint, Cint, png_const_doublep, png_const_doublep), png_ptr, heuristic_method, num_weights, filter_weights, filter_costs)
end

function png_set_filter_heuristics_fixed(png_ptr, heuristic_method, num_weights, filter_weights, filter_costs)
    ccall((:png_set_filter_heuristics_fixed, libpng), Cvoid, (png_structrp, Cint, Cint, png_const_fixed_point_p, png_const_fixed_point_p), png_ptr, heuristic_method, num_weights, filter_weights, filter_costs)
end

function png_set_compression_level(png_ptr, level)
    ccall((:png_set_compression_level, libpng), Cvoid, (png_structrp, Cint), png_ptr, level)
end

function png_set_compression_mem_level(png_ptr, mem_level)
    ccall((:png_set_compression_mem_level, libpng), Cvoid, (png_structrp, Cint), png_ptr, mem_level)
end

function png_set_compression_strategy(png_ptr, strategy)
    ccall((:png_set_compression_strategy, libpng), Cvoid, (png_structrp, Cint), png_ptr, strategy)
end

function png_set_compression_window_bits(png_ptr, window_bits)
    ccall((:png_set_compression_window_bits, libpng), Cvoid, (png_structrp, Cint), png_ptr, window_bits)
end

function png_set_compression_method(png_ptr, method)
    ccall((:png_set_compression_method, libpng), Cvoid, (png_structrp, Cint), png_ptr, method)
end

function png_set_text_compression_level(png_ptr, level)
    ccall((:png_set_text_compression_level, libpng), Cvoid, (png_structrp, Cint), png_ptr, level)
end

function png_set_text_compression_mem_level(png_ptr, mem_level)
    ccall((:png_set_text_compression_mem_level, libpng), Cvoid, (png_structrp, Cint), png_ptr, mem_level)
end

function png_set_text_compression_strategy(png_ptr, strategy)
    ccall((:png_set_text_compression_strategy, libpng), Cvoid, (png_structrp, Cint), png_ptr, strategy)
end

function png_set_text_compression_window_bits(png_ptr, window_bits)
    ccall((:png_set_text_compression_window_bits, libpng), Cvoid, (png_structrp, Cint), png_ptr, window_bits)
end

function png_set_text_compression_method(png_ptr, method)
    ccall((:png_set_text_compression_method, libpng), Cvoid, (png_structrp, Cint), png_ptr, method)
end

function png_init_io(png_ptr, fp)
    ccall((:png_init_io, libpng), Cvoid, (png_structrp, png_FILE_p), png_ptr, fp)
end

function png_set_error_fn(png_ptr, error_ptr, error_fn, warning_fn)
    ccall((:png_set_error_fn, libpng), Cvoid, (png_structrp, png_voidp, png_error_ptr, png_error_ptr), png_ptr, error_ptr, error_fn, warning_fn)
end

function png_get_error_ptr(png_ptr)
    ccall((:png_get_error_ptr, libpng), png_voidp, (png_const_structrp,), png_ptr)
end

function png_set_write_fn(png_ptr, io_ptr, write_data_fn, output_flush_fn)
    ccall((:png_set_write_fn, libpng), Cvoid, (png_structrp, Any, png_rw_ptr, png_flush_ptr), png_ptr, io_ptr, write_data_fn, output_flush_fn)
end

function png_set_read_fn(png_ptr, io_ptr, read_data_fn)
    ccall((:png_set_read_fn, libpng), Cvoid, (png_structrp, png_voidp, png_rw_ptr), png_ptr, io_ptr, read_data_fn)
end

function png_get_io_ptr(png_ptr)
    ccall((:png_get_io_ptr, libpng), png_voidp, (png_const_structrp,), png_ptr)
end

function png_set_read_status_fn(png_ptr, read_row_fn)
    ccall((:png_set_read_status_fn, libpng), Cvoid, (png_structrp, png_read_status_ptr), png_ptr, read_row_fn)
end

function png_set_write_status_fn(png_ptr, write_row_fn)
    ccall((:png_set_write_status_fn, libpng), Cvoid, (png_structrp, png_write_status_ptr), png_ptr, write_row_fn)
end

function png_set_mem_fn(png_ptr, mem_ptr, malloc_fn, free_fn)
    ccall((:png_set_mem_fn, libpng), Cvoid, (png_structrp, png_voidp, png_malloc_ptr, png_free_ptr), png_ptr, mem_ptr, malloc_fn, free_fn)
end

function png_get_mem_ptr(png_ptr)
    ccall((:png_get_mem_ptr, libpng), png_voidp, (png_const_structrp,), png_ptr)
end

function png_set_read_user_transform_fn(png_ptr, read_user_transform_fn)
    ccall((:png_set_read_user_transform_fn, libpng), Cvoid, (png_structrp, png_user_transform_ptr), png_ptr, read_user_transform_fn)
end

function png_set_write_user_transform_fn(png_ptr, write_user_transform_fn)
    ccall((:png_set_write_user_transform_fn, libpng), Cvoid, (png_structrp, png_user_transform_ptr), png_ptr, write_user_transform_fn)
end

function png_set_user_transform_info(png_ptr, user_transform_ptr, user_transform_depth, user_transform_channels)
    ccall((:png_set_user_transform_info, libpng), Cvoid, (png_structrp, png_voidp, Cint, Cint), png_ptr, user_transform_ptr, user_transform_depth, user_transform_channels)
end

function png_get_user_transform_ptr(png_ptr)
    ccall((:png_get_user_transform_ptr, libpng), png_voidp, (png_const_structrp,), png_ptr)
end

function png_get_current_row_number(arg1)
    ccall((:png_get_current_row_number, libpng), png_uint_32, (png_const_structrp,), arg1)
end

function png_get_current_pass_number(arg1)
    ccall((:png_get_current_pass_number, libpng), png_byte, (png_const_structrp,), arg1)
end

function png_set_read_user_chunk_fn(png_ptr, user_chunk_ptr, read_user_chunk_fn)
    ccall((:png_set_read_user_chunk_fn, libpng), Cvoid, (png_structrp, png_voidp, png_user_chunk_ptr), png_ptr, user_chunk_ptr, read_user_chunk_fn)
end

function png_get_user_chunk_ptr(png_ptr)
    ccall((:png_get_user_chunk_ptr, libpng), png_voidp, (png_const_structrp,), png_ptr)
end

function png_set_progressive_read_fn(png_ptr, progressive_ptr, info_fn, row_fn, end_fn)
    ccall((:png_set_progressive_read_fn, libpng), Cvoid, (png_structrp, png_voidp, png_progressive_info_ptr, png_progressive_row_ptr, png_progressive_end_ptr), png_ptr, progressive_ptr, info_fn, row_fn, end_fn)
end

function png_get_progressive_ptr(png_ptr)
    ccall((:png_get_progressive_ptr, libpng), png_voidp, (png_const_structrp,), png_ptr)
end

function png_process_data(png_ptr, info_ptr, buffer, buffer_size)
    ccall((:png_process_data, libpng), Cvoid, (png_structrp, png_inforp, png_bytep, Csize_t), png_ptr, info_ptr, buffer, buffer_size)
end

function png_process_data_pause(arg1, save)
    ccall((:png_process_data_pause, libpng), Csize_t, (png_structrp, Cint), arg1, save)
end

function png_process_data_skip(arg1)
    ccall((:png_process_data_skip, libpng), png_uint_32, (png_structrp,), arg1)
end

function png_progressive_combine_row(png_ptr, old_row, new_row)
    ccall((:png_progressive_combine_row, libpng), Cvoid, (png_const_structrp, png_bytep, png_const_bytep), png_ptr, old_row, new_row)
end

function png_malloc(png_ptr, size)
    ccall((:png_malloc, libpng), png_voidp, (png_const_structrp, png_alloc_size_t), png_ptr, size)
end

function png_calloc(png_ptr, size)
    ccall((:png_calloc, libpng), png_voidp, (png_const_structrp, png_alloc_size_t), png_ptr, size)
end

function png_malloc_warn(png_ptr, size)
    ccall((:png_malloc_warn, libpng), png_voidp, (png_const_structrp, png_alloc_size_t), png_ptr, size)
end

function png_free(png_ptr, ptr)
    ccall((:png_free, libpng), Cvoid, (png_const_structrp, png_voidp), png_ptr, ptr)
end

function png_free_data(png_ptr, info_ptr, free_me, num)
    ccall((:png_free_data, libpng), Cvoid, (png_const_structrp, png_inforp, png_uint_32, Cint), png_ptr, info_ptr, free_me, num)
end

function png_data_freer(png_ptr, info_ptr, freer, mask)
    ccall((:png_data_freer, libpng), Cvoid, (png_const_structrp, png_inforp, Cint, png_uint_32), png_ptr, info_ptr, freer, mask)
end

function png_malloc_default(png_ptr, size)
    ccall((:png_malloc_default, libpng), png_voidp, (png_const_structrp, png_alloc_size_t), png_ptr, size)
end

function png_free_default(png_ptr, ptr)
    ccall((:png_free_default, libpng), Cvoid, (png_const_structrp, png_voidp), png_ptr, ptr)
end

function png_error(png_ptr, error_message)
    ccall((:png_error, libpng), Cvoid, (png_const_structrp, png_const_charp), png_ptr, error_message)
end

function png_chunk_error(png_ptr, error_message)
    ccall((:png_chunk_error, libpng), Cvoid, (png_const_structrp, png_const_charp), png_ptr, error_message)
end

function png_warning(png_ptr, warning_message)
    ccall((:png_warning, libpng), Cvoid, (png_const_structrp, png_const_charp), png_ptr, warning_message)
end

function png_chunk_warning(png_ptr, warning_message)
    ccall((:png_chunk_warning, libpng), Cvoid, (png_const_structrp, png_const_charp), png_ptr, warning_message)
end

function png_benign_error(png_ptr, warning_message)
    ccall((:png_benign_error, libpng), Cvoid, (png_const_structrp, png_const_charp), png_ptr, warning_message)
end

function png_chunk_benign_error(png_ptr, warning_message)
    ccall((:png_chunk_benign_error, libpng), Cvoid, (png_const_structrp, png_const_charp), png_ptr, warning_message)
end

function png_set_benign_errors(png_ptr, allowed)
    ccall((:png_set_benign_errors, libpng), Cvoid, (png_structrp, Cint), png_ptr, allowed)
end

function png_get_valid(png_ptr, info_ptr, flag)
    ccall((:png_get_valid, libpng), png_uint_32, (png_const_structrp, png_const_inforp, png_uint_32), png_ptr, info_ptr, flag)
end

function png_get_rowbytes(png_ptr, info_ptr)
    ccall((:png_get_rowbytes, libpng), Csize_t, (png_const_structrp, png_const_inforp), png_ptr, info_ptr)
end

function png_get_rows(png_ptr, info_ptr)
    ccall((:png_get_rows, libpng), png_bytepp, (png_const_structrp, png_const_inforp), png_ptr, info_ptr)
end

function png_set_rows(png_ptr, info_ptr, row_pointers)
    ccall((:png_set_rows, libpng), Cvoid, (png_const_structrp, png_inforp, png_bytepp), png_ptr, info_ptr, row_pointers)
end

function png_get_channels(png_ptr, info_ptr)
    ccall((:png_get_channels, libpng), png_byte, (png_const_structrp, png_const_inforp), png_ptr, info_ptr)
end

function png_get_image_width(png_ptr, info_ptr)
    ccall((:png_get_image_width, libpng), png_uint_32, (png_const_structrp, png_const_inforp), png_ptr, info_ptr)
end

function png_get_image_height(png_ptr, info_ptr)
    ccall((:png_get_image_height, libpng), png_uint_32, (png_const_structrp, png_const_inforp), png_ptr, info_ptr)
end

function png_get_bit_depth(png_ptr, info_ptr)
    ccall((:png_get_bit_depth, libpng), png_byte, (png_const_structrp, png_const_inforp), png_ptr, info_ptr)
end

function png_get_color_type(png_ptr, info_ptr)
    ccall((:png_get_color_type, libpng), png_byte, (png_const_structrp, png_const_inforp), png_ptr, info_ptr)
end

function png_get_filter_type(png_ptr, info_ptr)
    ccall((:png_get_filter_type, libpng), png_byte, (png_const_structrp, png_const_inforp), png_ptr, info_ptr)
end

function png_get_interlace_type(png_ptr, info_ptr)
    ccall((:png_get_interlace_type, libpng), png_byte, (png_const_structrp, png_const_inforp), png_ptr, info_ptr)
end

function png_get_compression_type(png_ptr, info_ptr)
    ccall((:png_get_compression_type, libpng), png_byte, (png_const_structrp, png_const_inforp), png_ptr, info_ptr)
end

function png_get_pixels_per_meter(png_ptr, info_ptr)
    ccall((:png_get_pixels_per_meter, libpng), png_uint_32, (png_const_structrp, png_const_inforp), png_ptr, info_ptr)
end

function png_get_x_pixels_per_meter(png_ptr, info_ptr)
    ccall((:png_get_x_pixels_per_meter, libpng), png_uint_32, (png_const_structrp, png_const_inforp), png_ptr, info_ptr)
end

function png_get_y_pixels_per_meter(png_ptr, info_ptr)
    ccall((:png_get_y_pixels_per_meter, libpng), png_uint_32, (png_const_structrp, png_const_inforp), png_ptr, info_ptr)
end

function png_get_pixel_aspect_ratio(png_ptr, info_ptr)
    ccall((:png_get_pixel_aspect_ratio, libpng), Cfloat, (png_const_structrp, png_const_inforp), png_ptr, info_ptr)
end

function png_get_pixel_aspect_ratio_fixed(png_ptr, info_ptr)
    ccall((:png_get_pixel_aspect_ratio_fixed, libpng), png_fixed_point, (png_const_structrp, png_const_inforp), png_ptr, info_ptr)
end

function png_get_x_offset_pixels(png_ptr, info_ptr)
    ccall((:png_get_x_offset_pixels, libpng), png_int_32, (png_const_structrp, png_const_inforp), png_ptr, info_ptr)
end

function png_get_y_offset_pixels(png_ptr, info_ptr)
    ccall((:png_get_y_offset_pixels, libpng), png_int_32, (png_const_structrp, png_const_inforp), png_ptr, info_ptr)
end

function png_get_x_offset_microns(png_ptr, info_ptr)
    ccall((:png_get_x_offset_microns, libpng), png_int_32, (png_const_structrp, png_const_inforp), png_ptr, info_ptr)
end

function png_get_y_offset_microns(png_ptr, info_ptr)
    ccall((:png_get_y_offset_microns, libpng), png_int_32, (png_const_structrp, png_const_inforp), png_ptr, info_ptr)
end

function png_get_signature(png_ptr, info_ptr)
    ccall((:png_get_signature, libpng), png_const_bytep, (png_const_structrp, png_const_inforp), png_ptr, info_ptr)
end

function png_get_bKGD(png_ptr, info_ptr, background)
    ccall((:png_get_bKGD, libpng), png_uint_32, (png_const_structrp, png_inforp, Ptr{png_color_16p}), png_ptr, info_ptr, background)
end

function png_set_bKGD(png_ptr, info_ptr, background)
    ccall((:png_set_bKGD, libpng), Cvoid, (png_const_structrp, png_inforp, png_const_color_16p), png_ptr, info_ptr, background)
end

function png_get_cHRM(png_ptr, info_ptr, white_x, white_y, red_x, red_y, green_x, green_y, blue_x, blue_y)
    ccall((:png_get_cHRM, libpng), png_uint_32, (png_const_structrp, png_const_inforp, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), png_ptr, info_ptr, white_x, white_y, red_x, red_y, green_x, green_y, blue_x, blue_y)
end

function png_get_cHRM_XYZ(png_ptr, info_ptr, red_X, red_Y, red_Z, green_X, green_Y, green_Z, blue_X, blue_Y, blue_Z)
    ccall((:png_get_cHRM_XYZ, libpng), png_uint_32, (png_const_structrp, png_const_inforp, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), png_ptr, info_ptr, red_X, red_Y, red_Z, green_X, green_Y, green_Z, blue_X, blue_Y, blue_Z)
end

function png_get_cHRM_fixed(png_ptr, info_ptr, int_white_x, int_white_y, int_red_x, int_red_y, int_green_x, int_green_y, int_blue_x, int_blue_y)
    ccall((:png_get_cHRM_fixed, libpng), png_uint_32, (png_const_structrp, png_const_inforp, Ptr{png_fixed_point}, Ptr{png_fixed_point}, Ptr{png_fixed_point}, Ptr{png_fixed_point}, Ptr{png_fixed_point}, Ptr{png_fixed_point}, Ptr{png_fixed_point}, Ptr{png_fixed_point}), png_ptr, info_ptr, int_white_x, int_white_y, int_red_x, int_red_y, int_green_x, int_green_y, int_blue_x, int_blue_y)
end

function png_get_cHRM_XYZ_fixed(png_ptr, info_ptr, int_red_X, int_red_Y, int_red_Z, int_green_X, int_green_Y, int_green_Z, int_blue_X, int_blue_Y, int_blue_Z)
    ccall((:png_get_cHRM_XYZ_fixed, libpng), png_uint_32, (png_const_structrp, png_const_inforp, Ptr{png_fixed_point}, Ptr{png_fixed_point}, Ptr{png_fixed_point}, Ptr{png_fixed_point}, Ptr{png_fixed_point}, Ptr{png_fixed_point}, Ptr{png_fixed_point}, Ptr{png_fixed_point}, Ptr{png_fixed_point}), png_ptr, info_ptr, int_red_X, int_red_Y, int_red_Z, int_green_X, int_green_Y, int_green_Z, int_blue_X, int_blue_Y, int_blue_Z)
end

function png_set_cHRM(png_ptr, info_ptr, white_x, white_y, red_x, red_y, green_x, green_y, blue_x, blue_y)
    ccall((:png_set_cHRM, libpng), Cvoid, (png_const_structrp, png_inforp, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble), png_ptr, info_ptr, white_x, white_y, red_x, red_y, green_x, green_y, blue_x, blue_y)
end

function png_set_cHRM_XYZ(png_ptr, info_ptr, red_X, red_Y, red_Z, green_X, green_Y, green_Z, blue_X, blue_Y, blue_Z)
    ccall((:png_set_cHRM_XYZ, libpng), Cvoid, (png_const_structrp, png_inforp, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble), png_ptr, info_ptr, red_X, red_Y, red_Z, green_X, green_Y, green_Z, blue_X, blue_Y, blue_Z)
end

function png_set_cHRM_fixed(png_ptr, info_ptr, int_white_x, int_white_y, int_red_x, int_red_y, int_green_x, int_green_y, int_blue_x, int_blue_y)
    ccall((:png_set_cHRM_fixed, libpng), Cvoid, (png_const_structrp, png_inforp, png_fixed_point, png_fixed_point, png_fixed_point, png_fixed_point, png_fixed_point, png_fixed_point, png_fixed_point, png_fixed_point), png_ptr, info_ptr, int_white_x, int_white_y, int_red_x, int_red_y, int_green_x, int_green_y, int_blue_x, int_blue_y)
end

function png_set_cHRM_XYZ_fixed(png_ptr, info_ptr, int_red_X, int_red_Y, int_red_Z, int_green_X, int_green_Y, int_green_Z, int_blue_X, int_blue_Y, int_blue_Z)
    ccall((:png_set_cHRM_XYZ_fixed, libpng), Cvoid, (png_const_structrp, png_inforp, png_fixed_point, png_fixed_point, png_fixed_point, png_fixed_point, png_fixed_point, png_fixed_point, png_fixed_point, png_fixed_point, png_fixed_point), png_ptr, info_ptr, int_red_X, int_red_Y, int_red_Z, int_green_X, int_green_Y, int_green_Z, int_blue_X, int_blue_Y, int_blue_Z)
end

function png_get_eXIf(png_ptr, info_ptr, exif)
    ccall((:png_get_eXIf, libpng), png_uint_32, (png_const_structrp, png_inforp, Ptr{png_bytep}), png_ptr, info_ptr, exif)
end

function png_set_eXIf(png_ptr, info_ptr, exif)
    ccall((:png_set_eXIf, libpng), Cvoid, (png_const_structrp, png_inforp, png_bytep), png_ptr, info_ptr, exif)
end

function png_get_eXIf_1(png_ptr, info_ptr, num_exif, exif)
    ccall((:png_get_eXIf_1, libpng), png_uint_32, (png_const_structrp, png_const_inforp, Ptr{png_uint_32}, Ptr{png_bytep}), png_ptr, info_ptr, num_exif, exif)
end

function png_set_eXIf_1(png_ptr, info_ptr, num_exif, exif)
    ccall((:png_set_eXIf_1, libpng), Cvoid, (png_const_structrp, png_inforp, png_uint_32, png_bytep), png_ptr, info_ptr, num_exif, exif)
end

function png_get_gAMA(png_ptr, info_ptr, file_gamma)
    ccall((:png_get_gAMA, libpng), png_uint_32, (png_const_structrp, png_const_inforp, Ptr{Cdouble}), png_ptr, info_ptr, file_gamma)
end

function png_get_gAMA_fixed(png_ptr, info_ptr, int_file_gamma)
    ccall((:png_get_gAMA_fixed, libpng), png_uint_32, (png_const_structrp, png_const_inforp, Ptr{png_fixed_point}), png_ptr, info_ptr, int_file_gamma)
end

function png_set_gAMA(png_ptr, info_ptr, file_gamma)
    ccall((:png_set_gAMA, libpng), Cvoid, (png_const_structrp, png_inforp, Cdouble), png_ptr, info_ptr, file_gamma)
end

function png_set_gAMA_fixed(png_ptr, info_ptr, int_file_gamma)
    ccall((:png_set_gAMA_fixed, libpng), Cvoid, (png_const_structrp, png_inforp, png_fixed_point), png_ptr, info_ptr, int_file_gamma)
end

function png_get_hIST(png_ptr, info_ptr, hist)
    ccall((:png_get_hIST, libpng), png_uint_32, (png_const_structrp, png_inforp, Ptr{png_uint_16p}), png_ptr, info_ptr, hist)
end

function png_set_hIST(png_ptr, info_ptr, hist)
    ccall((:png_set_hIST, libpng), Cvoid, (png_const_structrp, png_inforp, png_const_uint_16p), png_ptr, info_ptr, hist)
end

function png_get_IHDR(png_ptr, info_ptr, width, height, bit_depth, color_type, interlace_method, compression_method, filter_method)
    ccall((:png_get_IHDR, libpng), png_uint_32, (png_const_structrp, png_const_inforp, Ptr{png_uint_32}, Ptr{png_uint_32}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}), png_ptr, info_ptr, width, height, bit_depth, color_type, interlace_method, compression_method, filter_method)
end

function png_set_IHDR(png_ptr, info_ptr, width, height, bit_depth, color_type, interlace_method, compression_method, filter_method)
    ccall((:png_set_IHDR, libpng), Cvoid, (png_const_structrp, png_inforp, png_uint_32, png_uint_32, Cint, Cint, Cint, Cint, Cint), png_ptr, info_ptr, width, height, bit_depth, color_type, interlace_method, compression_method, filter_method)
end

function png_get_oFFs(png_ptr, info_ptr, offset_x, offset_y, unit_type)
    ccall((:png_get_oFFs, libpng), png_uint_32, (png_const_structrp, png_const_inforp, Ptr{png_int_32}, Ptr{png_int_32}, Ptr{Cint}), png_ptr, info_ptr, offset_x, offset_y, unit_type)
end

function png_set_oFFs(png_ptr, info_ptr, offset_x, offset_y, unit_type)
    ccall((:png_set_oFFs, libpng), Cvoid, (png_const_structrp, png_inforp, png_int_32, png_int_32, Cint), png_ptr, info_ptr, offset_x, offset_y, unit_type)
end

function png_get_pCAL(png_ptr, info_ptr, purpose, X0, X1, type, nparams, units, params)
    ccall((:png_get_pCAL, libpng), png_uint_32, (png_const_structrp, png_inforp, Ptr{png_charp}, Ptr{png_int_32}, Ptr{png_int_32}, Ptr{Cint}, Ptr{Cint}, Ptr{png_charp}, Ptr{png_charpp}), png_ptr, info_ptr, purpose, X0, X1, type, nparams, units, params)
end

function png_set_pCAL(png_ptr, info_ptr, purpose, X0, X1, type, nparams, units, params)
    ccall((:png_set_pCAL, libpng), Cvoid, (png_const_structrp, png_inforp, png_const_charp, png_int_32, png_int_32, Cint, Cint, png_const_charp, png_charpp), png_ptr, info_ptr, purpose, X0, X1, type, nparams, units, params)
end

function png_get_pHYs(png_ptr, info_ptr, res_x, res_y, unit_type)
    ccall((:png_get_pHYs, libpng), png_uint_32, (png_const_structrp, png_const_inforp, Ptr{png_uint_32}, Ptr{png_uint_32}, Ptr{Cint}), png_ptr, info_ptr, res_x, res_y, unit_type)
end

function png_set_pHYs(png_ptr, info_ptr, res_x, res_y, unit_type)
    ccall((:png_set_pHYs, libpng), Cvoid, (png_const_structrp, png_inforp, png_uint_32, png_uint_32, Cint), png_ptr, info_ptr, res_x, res_y, unit_type)
end

function png_get_PLTE(png_ptr, info_ptr, palette, num_palette)
    ccall((:png_get_PLTE, libpng), png_uint_32, (png_const_structrp, png_inforp, Ptr{png_colorp}, Ptr{Cint}), png_ptr, info_ptr, palette, num_palette)
end

function png_set_PLTE(png_ptr, info_ptr, palette, num_palette)
    ccall((:png_set_PLTE, libpng), Cvoid, (png_structrp, png_inforp, png_const_colorp, Cint), png_ptr, info_ptr, palette, num_palette)
end

function png_get_sBIT(png_ptr, info_ptr, sig_bit)
    ccall((:png_get_sBIT, libpng), png_uint_32, (png_const_structrp, png_inforp, Ptr{png_color_8p}), png_ptr, info_ptr, sig_bit)
end

function png_set_sBIT(png_ptr, info_ptr, sig_bit)
    ccall((:png_set_sBIT, libpng), Cvoid, (png_const_structrp, png_inforp, png_const_color_8p), png_ptr, info_ptr, sig_bit)
end

function png_get_sRGB(png_ptr, info_ptr, file_srgb_intent)
    ccall((:png_get_sRGB, libpng), png_uint_32, (png_const_structrp, png_const_inforp, Ptr{Cint}), png_ptr, info_ptr, file_srgb_intent)
end

function png_set_sRGB(png_ptr, info_ptr, srgb_intent)
    ccall((:png_set_sRGB, libpng), Cvoid, (png_const_structrp, png_inforp, Cint), png_ptr, info_ptr, srgb_intent)
end

function png_set_sRGB_gAMA_and_cHRM(png_ptr, info_ptr, srgb_intent)
    ccall((:png_set_sRGB_gAMA_and_cHRM, libpng), Cvoid, (png_const_structrp, png_inforp, Cint), png_ptr, info_ptr, srgb_intent)
end

function png_get_iCCP(png_ptr, info_ptr, name, compression_type, profile, proflen)
    ccall((:png_get_iCCP, libpng), png_uint_32, (png_const_structrp, png_inforp, png_charpp, Ptr{Cint}, png_bytepp, Ptr{png_uint_32}), png_ptr, info_ptr, name, compression_type, profile, proflen)
end

function png_set_iCCP(png_ptr, info_ptr, name, compression_type, profile, proflen)
    ccall((:png_set_iCCP, libpng), Cvoid, (png_const_structrp, png_inforp, png_const_charp, Cint, png_const_bytep, png_uint_32), png_ptr, info_ptr, name, compression_type, profile, proflen)
end

function png_get_sPLT(png_ptr, info_ptr, entries)
    ccall((:png_get_sPLT, libpng), Cint, (png_const_structrp, png_inforp, png_sPLT_tpp), png_ptr, info_ptr, entries)
end

function png_set_sPLT(png_ptr, info_ptr, entries, nentries)
    ccall((:png_set_sPLT, libpng), Cvoid, (png_const_structrp, png_inforp, png_const_sPLT_tp, Cint), png_ptr, info_ptr, entries, nentries)
end

function png_get_text(png_ptr, info_ptr, text_ptr, num_text)
    ccall((:png_get_text, libpng), Cint, (png_const_structrp, png_inforp, Ptr{png_textp}, Ptr{Cint}), png_ptr, info_ptr, text_ptr, num_text)
end

function png_set_text(png_ptr, info_ptr, text_ptr, num_text)
    ccall((:png_set_text, libpng), Cvoid, (png_const_structrp, png_inforp, png_const_textp, Cint), png_ptr, info_ptr, text_ptr, num_text)
end

function png_get_tIME(png_ptr, info_ptr, mod_time)
    ccall((:png_get_tIME, libpng), png_uint_32, (png_const_structrp, png_inforp, Ptr{png_timep}), png_ptr, info_ptr, mod_time)
end

function png_set_tIME(png_ptr, info_ptr, mod_time)
    ccall((:png_set_tIME, libpng), Cvoid, (png_const_structrp, png_inforp, png_const_timep), png_ptr, info_ptr, mod_time)
end

function png_get_tRNS(png_ptr, info_ptr, trans_alpha, num_trans, trans_color)
    ccall((:png_get_tRNS, libpng), png_uint_32, (png_const_structrp, png_inforp, Ptr{png_bytep}, Ptr{Cint}, Ptr{png_color_16p}), png_ptr, info_ptr, trans_alpha, num_trans, trans_color)
end

function png_set_tRNS(png_ptr, info_ptr, trans_alpha, num_trans, trans_color)
    ccall((:png_set_tRNS, libpng), Cvoid, (png_structrp, png_inforp, png_const_bytep, Cint, png_const_color_16p), png_ptr, info_ptr, trans_alpha, num_trans, trans_color)
end

function png_get_sCAL(png_ptr, info_ptr, unit, width, height)
    ccall((:png_get_sCAL, libpng), png_uint_32, (png_const_structrp, png_const_inforp, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}), png_ptr, info_ptr, unit, width, height)
end

function png_get_sCAL_fixed(png_ptr, info_ptr, unit, width, height)
    ccall((:png_get_sCAL_fixed, libpng), png_uint_32, (png_const_structrp, png_const_inforp, Ptr{Cint}, Ptr{png_fixed_point}, Ptr{png_fixed_point}), png_ptr, info_ptr, unit, width, height)
end

function png_get_sCAL_s(png_ptr, info_ptr, unit, swidth, sheight)
    ccall((:png_get_sCAL_s, libpng), png_uint_32, (png_const_structrp, png_const_inforp, Ptr{Cint}, png_charpp, png_charpp), png_ptr, info_ptr, unit, swidth, sheight)
end

function png_set_sCAL(png_ptr, info_ptr, unit, width, height)
    ccall((:png_set_sCAL, libpng), Cvoid, (png_const_structrp, png_inforp, Cint, Cdouble, Cdouble), png_ptr, info_ptr, unit, width, height)
end

function png_set_sCAL_fixed(png_ptr, info_ptr, unit, width, height)
    ccall((:png_set_sCAL_fixed, libpng), Cvoid, (png_const_structrp, png_inforp, Cint, png_fixed_point, png_fixed_point), png_ptr, info_ptr, unit, width, height)
end

function png_set_sCAL_s(png_ptr, info_ptr, unit, swidth, sheight)
    ccall((:png_set_sCAL_s, libpng), Cvoid, (png_const_structrp, png_inforp, Cint, png_const_charp, png_const_charp), png_ptr, info_ptr, unit, swidth, sheight)
end

function png_set_keep_unknown_chunks(png_ptr, keep, chunk_list, num_chunks)
    ccall((:png_set_keep_unknown_chunks, libpng), Cvoid, (png_structrp, Cint, png_const_bytep, Cint), png_ptr, keep, chunk_list, num_chunks)
end

function png_handle_as_unknown(png_ptr, chunk_name)
    ccall((:png_handle_as_unknown, libpng), Cint, (png_const_structrp, png_const_bytep), png_ptr, chunk_name)
end

function png_set_unknown_chunks(png_ptr, info_ptr, unknowns, num_unknowns)
    ccall((:png_set_unknown_chunks, libpng), Cvoid, (png_const_structrp, png_inforp, png_const_unknown_chunkp, Cint), png_ptr, info_ptr, unknowns, num_unknowns)
end

function png_set_unknown_chunk_location(png_ptr, info_ptr, chunk, location)
    ccall((:png_set_unknown_chunk_location, libpng), Cvoid, (png_const_structrp, png_inforp, Cint, Cint), png_ptr, info_ptr, chunk, location)
end

function png_get_unknown_chunks(png_ptr, info_ptr, entries)
    ccall((:png_get_unknown_chunks, libpng), Cint, (png_const_structrp, png_inforp, png_unknown_chunkpp), png_ptr, info_ptr, entries)
end

function png_set_invalid(png_ptr, info_ptr, mask)
    ccall((:png_set_invalid, libpng), Cvoid, (png_const_structrp, png_inforp, Cint), png_ptr, info_ptr, mask)
end

function png_read_png(png_ptr, info_ptr, transforms, params)
    ccall((:png_read_png, libpng), Cvoid, (png_structrp, png_inforp, Cint, png_voidp), png_ptr, info_ptr, transforms, params)
end

function png_write_png(png_ptr, info_ptr, transforms, params)
    ccall((:png_write_png, libpng), Cvoid, (png_structrp, png_inforp, Cint, png_voidp), png_ptr, info_ptr, transforms, params)
end

function png_get_copyright(png_ptr)
    ccall((:png_get_copyright, libpng), png_const_charp, (png_const_structrp,), png_ptr)
end

function png_get_header_ver(png_ptr)
    ccall((:png_get_header_ver, libpng), png_const_charp, (png_const_structrp,), png_ptr)
end

function png_get_header_version(png_ptr)
    ccall((:png_get_header_version, libpng), png_const_charp, (png_const_structrp,), png_ptr)
end

function png_get_libpng_ver(png_ptr)
    ccall((:png_get_libpng_ver, libpng), png_const_charp, (png_const_structrp,), png_ptr)
end

function png_permit_mng_features(png_ptr, mng_features_permitted)
    ccall((:png_permit_mng_features, libpng), png_uint_32, (png_structrp, png_uint_32), png_ptr, mng_features_permitted)
end

function png_set_user_limits(png_ptr, user_width_max, user_height_max)
    ccall((:png_set_user_limits, libpng), Cvoid, (png_structrp, png_uint_32, png_uint_32), png_ptr, user_width_max, user_height_max)
end

function png_get_user_width_max(png_ptr)
    ccall((:png_get_user_width_max, libpng), png_uint_32, (png_const_structrp,), png_ptr)
end

function png_get_user_height_max(png_ptr)
    ccall((:png_get_user_height_max, libpng), png_uint_32, (png_const_structrp,), png_ptr)
end

function png_set_chunk_cache_max(png_ptr, user_chunk_cache_max)
    ccall((:png_set_chunk_cache_max, libpng), Cvoid, (png_structrp, png_uint_32), png_ptr, user_chunk_cache_max)
end

function png_get_chunk_cache_max(png_ptr)
    ccall((:png_get_chunk_cache_max, libpng), png_uint_32, (png_const_structrp,), png_ptr)
end

function png_set_chunk_malloc_max(png_ptr, user_chunk_cache_max)
    ccall((:png_set_chunk_malloc_max, libpng), Cvoid, (png_structrp, png_alloc_size_t), png_ptr, user_chunk_cache_max)
end

function png_get_chunk_malloc_max(png_ptr)
    ccall((:png_get_chunk_malloc_max, libpng), png_alloc_size_t, (png_const_structrp,), png_ptr)
end

function png_get_pixels_per_inch(png_ptr, info_ptr)
    ccall((:png_get_pixels_per_inch, libpng), png_uint_32, (png_const_structrp, png_const_inforp), png_ptr, info_ptr)
end

function png_get_x_pixels_per_inch(png_ptr, info_ptr)
    ccall((:png_get_x_pixels_per_inch, libpng), png_uint_32, (png_const_structrp, png_const_inforp), png_ptr, info_ptr)
end

function png_get_y_pixels_per_inch(png_ptr, info_ptr)
    ccall((:png_get_y_pixels_per_inch, libpng), png_uint_32, (png_const_structrp, png_const_inforp), png_ptr, info_ptr)
end

function png_get_x_offset_inches(png_ptr, info_ptr)
    ccall((:png_get_x_offset_inches, libpng), Cfloat, (png_const_structrp, png_const_inforp), png_ptr, info_ptr)
end

function png_get_x_offset_inches_fixed(png_ptr, info_ptr)
    ccall((:png_get_x_offset_inches_fixed, libpng), png_fixed_point, (png_const_structrp, png_const_inforp), png_ptr, info_ptr)
end

function png_get_y_offset_inches(png_ptr, info_ptr)
    ccall((:png_get_y_offset_inches, libpng), Cfloat, (png_const_structrp, png_const_inforp), png_ptr, info_ptr)
end

function png_get_y_offset_inches_fixed(png_ptr, info_ptr)
    ccall((:png_get_y_offset_inches_fixed, libpng), png_fixed_point, (png_const_structrp, png_const_inforp), png_ptr, info_ptr)
end

function png_get_pHYs_dpi(png_ptr, info_ptr, res_x, res_y, unit_type)
    ccall((:png_get_pHYs_dpi, libpng), png_uint_32, (png_const_structrp, png_const_inforp, Ptr{png_uint_32}, Ptr{png_uint_32}, Ptr{Cint}), png_ptr, info_ptr, res_x, res_y, unit_type)
end

function png_get_io_state(png_ptr)
    ccall((:png_get_io_state, libpng), png_uint_32, (png_const_structrp,), png_ptr)
end

function png_get_io_chunk_type(png_ptr)
    ccall((:png_get_io_chunk_type, libpng), png_uint_32, (png_const_structrp,), png_ptr)
end

function png_get_uint_31(png_ptr, buf)
    ccall((:png_get_uint_31, libpng), png_uint_32, (png_const_structrp, png_const_bytep), png_ptr, buf)
end

function png_save_uint_32(buf, i)
    ccall((:png_save_uint_32, libpng), Cvoid, (png_bytep, png_uint_32), buf, i)
end

function png_save_int_32(buf, i)
    ccall((:png_save_int_32, libpng), Cvoid, (png_bytep, png_int_32), buf, i)
end

function png_save_uint_16(buf, i)
    ccall((:png_save_uint_16, libpng), Cvoid, (png_bytep, Cuint), buf, i)
end

function png_set_check_for_invalid_index(png_ptr, allowed)
    ccall((:png_set_check_for_invalid_index, libpng), Cvoid, (png_structrp, Cint), png_ptr, allowed)
end

function png_get_palette_max(png_ptr, info_ptr)
    ccall((:png_get_palette_max, libpng), Cint, (png_const_structp, png_const_infop), png_ptr, info_ptr)
end

const png_control = Cvoid

const png_controlp = Ptr{png_control}

struct __JL_Ctag_7
    opaque::png_controlp
    version::png_uint_32
    width::png_uint_32
    height::png_uint_32
    format::png_uint_32
    flags::png_uint_32
    colormap_entries::png_uint_32
    warning_or_error::png_uint_32
    message::NTuple{64, Cchar}
end
function Base.getproperty(x::Ptr{__JL_Ctag_7}, f::Symbol)
    f === :opaque && return Ptr{png_controlp}(x + 0)
    f === :version && return Ptr{png_uint_32}(x + 8)
    f === :width && return Ptr{png_uint_32}(x + 12)
    f === :height && return Ptr{png_uint_32}(x + 16)
    f === :format && return Ptr{png_uint_32}(x + 20)
    f === :flags && return Ptr{png_uint_32}(x + 24)
    f === :colormap_entries && return Ptr{png_uint_32}(x + 28)
    f === :warning_or_error && return Ptr{png_uint_32}(x + 32)
    f === :message && return Ptr{NTuple{64, Cchar}}(x + 36)
    return getfield(x, f)
end

function Base.getproperty(x::__JL_Ctag_7, f::Symbol)
    r = Ref{__JL_Ctag_7}(x)
    ptr = Base.unsafe_convert(Ptr{__JL_Ctag_7}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{__JL_Ctag_7}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end


const png_image = __JL_Ctag_7

const png_imagep = Ptr{__JL_Ctag_7}

function png_image_begin_read_from_file(image, file_name)
    ccall((:png_image_begin_read_from_file, libpng), Cint, (png_imagep, Cstring), image, file_name)
end

function png_image_begin_read_from_stdio(image, file)
    ccall((:png_image_begin_read_from_stdio, libpng), Cint, (png_imagep, Ptr{Libc.FILE}), image, file)
end

function png_image_begin_read_from_memory(image, memory, size)
    ccall((:png_image_begin_read_from_memory, libpng), Cint, (png_imagep, png_const_voidp, Csize_t), image, memory, size)
end

function png_image_finish_read(image, background, buffer, row_stride, colormap)
    ccall((:png_image_finish_read, libpng), Cint, (png_imagep, png_const_colorp, Ptr{Cvoid}, png_int_32, Ptr{Cvoid}), image, background, buffer, row_stride, colormap)
end

function png_image_free(image)
    ccall((:png_image_free, libpng), Cvoid, (png_imagep,), image)
end

function png_image_write_to_file(image, file, convert_to_8bit, buffer, row_stride, colormap)
    ccall((:png_image_write_to_file, libpng), Cint, (png_imagep, Cstring, Cint, Ptr{Cvoid}, png_int_32, Ptr{Cvoid}), image, file, convert_to_8bit, buffer, row_stride, colormap)
end

function png_image_write_to_stdio(image, file, convert_to_8_bit, buffer, row_stride, colormap)
    ccall((:png_image_write_to_stdio, libpng), Cint, (png_imagep, Ptr{Libc.FILE}, Cint, Ptr{Cvoid}, png_int_32, Ptr{Cvoid}), image, file, convert_to_8_bit, buffer, row_stride, colormap)
end

function png_image_write_to_memory(image, memory, memory_bytes, convert_to_8_bit, buffer, row_stride, colormap)
    ccall((:png_image_write_to_memory, libpng), Cint, (png_imagep, Ptr{Cvoid}, Ptr{png_alloc_size_t}, Cint, Ptr{Cvoid}, png_int_32, Ptr{Cvoid}), image, memory, memory_bytes, convert_to_8_bit, buffer, row_stride, colormap)
end

function png_set_option(png_ptr, option, onoff)
    ccall((:png_set_option, libpng), Cint, (png_structrp, Cint, Cint), png_ptr, option, onoff)
end

const PNG_LIBPNG_VER_STRING = "1.6.37"

const PNG_HEADER_VERSION_STRING = " libpng version 1.6.37 - April 14, 2019\n"

const PNG_LIBPNG_VER_SONUM = 16

const PNG_LIBPNG_VER_DLLNUM = 16

const PNG_LIBPNG_VER_MAJOR = 1

const PNG_LIBPNG_VER_MINOR = 6

const PNG_LIBPNG_VER_RELEASE = 37

const PNG_LIBPNG_VER_BUILD = 0

const PNG_LIBPNG_BUILD_ALPHA = 1

const PNG_LIBPNG_BUILD_BETA = 2

const PNG_LIBPNG_BUILD_RC = 3

const PNG_LIBPNG_BUILD_STABLE = 4

const PNG_LIBPNG_BUILD_RELEASE_STATUS_MASK = 7

const PNG_LIBPNG_BUILD_PATCH = 8

const PNG_LIBPNG_BUILD_PRIVATE = 16

const PNG_LIBPNG_BUILD_SPECIAL = 32

const PNG_LIBPNG_BUILD_BASE_TYPE = PNG_LIBPNG_BUILD_STABLE

const PNG_LIBPNG_VER = 10637

const PNG_API_RULE = 0

const PNG_DEFAULT_READ_MACROS = 1

const PNG_GAMMA_THRESHOLD_FIXED = 5000

const PNG_INFLATE_BUF_SIZE = 1024

# Skipping MacroDefinition: PNG_LINKAGE_API extern

# Skipping MacroDefinition: PNG_LINKAGE_CALLBACK extern

# Skipping MacroDefinition: PNG_LINKAGE_DATA extern

# Skipping MacroDefinition: PNG_LINKAGE_FUNCTION extern

const PNG_MAX_GAMMA_8 = 11

const PNG_QUANTIZE_BLUE_BITS = 5

const PNG_QUANTIZE_GREEN_BITS = 5

const PNG_QUANTIZE_RED_BITS = 5

const PNG_TEXT_Z_DEFAULT_COMPRESSION = -1

const PNG_TEXT_Z_DEFAULT_STRATEGY = 0

const PNG_USER_CHUNK_CACHE_MAX = 1000

const PNG_USER_CHUNK_MALLOC_MAX = 8000000

const PNG_USER_HEIGHT_MAX = 1000000

const PNG_USER_WIDTH_MAX = 1000000

const PNG_ZBUF_SIZE = 8192

const PNG_ZLIB_VERNUM = 0x12b0

const PNG_Z_DEFAULT_COMPRESSION = -1

const PNG_Z_DEFAULT_NOFILTER_STRATEGY = 0

const PNG_Z_DEFAULT_STRATEGY = 1

const PNG_sCAL_PRECISION = 5

const PNG_sRGB_PROFILE_CHECKS = 2

# Skipping MacroDefinition: PNG_CONST const

const PNGCBAPI = PNGCAPI

const PNGAPI = PNGCAPI

# Skipping MacroDefinition: PNG_USE_RESULT __attribute__ ( ( __warn_unused_result__ ) )

# Skipping MacroDefinition: PNG_NORETURN __attribute__ ( ( __noreturn__ ) )

# Skipping MacroDefinition: PNG_ALLOCATED __attribute__ ( ( __malloc__ ) )

# Skipping MacroDefinition: PNG_DEPRECATED __attribute__ ( ( __deprecated__ ) )

# Skipping MacroDefinition: PNG_PRIVATE __attribute__ ( ( __unavailable__ ( "This function is not exported by libpng." ) ) )

const PNG_LIBPNG_BUILD_TYPE = PNG_LIBPNG_BUILD_BASE_TYPE

const PNG_TEXT_COMPRESSION_NONE_WR = -3

const PNG_TEXT_COMPRESSION_zTXt_WR = -2

const PNG_TEXT_COMPRESSION_NONE = -1

const PNG_TEXT_COMPRESSION_zTXt = 0

const PNG_ITXT_COMPRESSION_NONE = 1

const PNG_ITXT_COMPRESSION_zTXt = 2

const PNG_TEXT_COMPRESSION_LAST = 3

const PNG_HAVE_IHDR = 0x01

const PNG_HAVE_PLTE = 0x02

const PNG_AFTER_IDAT = 0x08

const PNG_UINT_31_MAX = png_uint_32(Clong(0x7fffffff))

const PNG_FP_1 = 100000

const PNG_FP_HALF = 50000

const PNG_FP_MAX = png_fixed_point(Clong(0x7fffffff))

const PNG_FP_MIN = -PNG_FP_MAX

const PNG_COLOR_MASK_PALETTE = 1

const PNG_COLOR_MASK_COLOR = 2

const PNG_COLOR_MASK_ALPHA = 4

const PNG_COLOR_TYPE_GRAY = 0

const PNG_COLOR_TYPE_PALETTE = PNG_COLOR_MASK_COLOR | PNG_COLOR_MASK_PALETTE

const PNG_COLOR_TYPE_RGB = PNG_COLOR_MASK_COLOR

const PNG_COLOR_TYPE_RGB_ALPHA = PNG_COLOR_MASK_COLOR | PNG_COLOR_MASK_ALPHA

const PNG_COLOR_TYPE_GRAY_ALPHA = PNG_COLOR_MASK_ALPHA

const PNG_COLOR_TYPE_RGBA = PNG_COLOR_TYPE_RGB_ALPHA

const PNG_COLOR_TYPE_GA = PNG_COLOR_TYPE_GRAY_ALPHA

const PNG_COMPRESSION_TYPE_BASE = 0

const PNG_COMPRESSION_TYPE_DEFAULT = PNG_COMPRESSION_TYPE_BASE

const PNG_FILTER_TYPE_BASE = 0

const PNG_INTRAPIXEL_DIFFERENCING = 64

const PNG_FILTER_TYPE_DEFAULT = PNG_FILTER_TYPE_BASE

const PNG_INTERLACE_NONE = 0

const PNG_INTERLACE_ADAM7 = 1

const PNG_INTERLACE_LAST = 2

const PNG_OFFSET_PIXEL = 0

const PNG_OFFSET_MICROMETER = 1

const PNG_OFFSET_LAST = 2

const PNG_EQUATION_LINEAR = 0

const PNG_EQUATION_BASE_E = 1

const PNG_EQUATION_ARBITRARY = 2

const PNG_EQUATION_HYPERBOLIC = 3

const PNG_EQUATION_LAST = 4

const PNG_SCALE_UNKNOWN = 0

const PNG_SCALE_METER = 1

const PNG_SCALE_RADIAN = 2

const PNG_SCALE_LAST = 3

const PNG_RESOLUTION_UNKNOWN = 0

const PNG_RESOLUTION_METER = 1

const PNG_RESOLUTION_LAST = 2

const PNG_sRGB_INTENT_PERCEPTUAL = 0

const PNG_sRGB_INTENT_RELATIVE = 1

const PNG_sRGB_INTENT_SATURATION = 2

const PNG_sRGB_INTENT_ABSOLUTE = 3

const PNG_sRGB_INTENT_LAST = 4

const PNG_KEYWORD_MAX_LENGTH = 79

const PNG_MAX_PALETTE_LENGTH = 256

const PNG_INFO_gAMA = Cuint(0x0001)

const PNG_INFO_sBIT = Cuint(0x0002)

const PNG_INFO_cHRM = Cuint(0x0004)

const PNG_INFO_PLTE = Cuint(0x0008)

const PNG_INFO_tRNS = Cuint(0x0010)

const PNG_INFO_bKGD = Cuint(0x0020)

const PNG_INFO_hIST = Cuint(0x0040)

const PNG_INFO_pHYs = Cuint(0x0080)

const PNG_INFO_oFFs = Cuint(0x0100)

const PNG_INFO_tIME = Cuint(0x0200)

const PNG_INFO_pCAL = Cuint(0x0400)

const PNG_INFO_sRGB = Cuint(0x0800)

const PNG_INFO_iCCP = Cuint(0x1000)

const PNG_INFO_sPLT = Cuint(0x2000)

const PNG_INFO_sCAL = Cuint(0x4000)

const PNG_INFO_IDAT = Cuint(0x8000)

const PNG_INFO_eXIf = Cuint(0x00010000)

const PNG_TRANSFORM_IDENTITY = 0x0000

const PNG_TRANSFORM_STRIP_16 = 0x0001

const PNG_TRANSFORM_STRIP_ALPHA = 0x0002

const PNG_TRANSFORM_PACKING = 0x0004

const PNG_TRANSFORM_PACKSWAP = 0x0008

const PNG_TRANSFORM_EXPAND = 0x0010

const PNG_TRANSFORM_INVERT_MONO = 0x0020

const PNG_TRANSFORM_SHIFT = 0x0040

const PNG_TRANSFORM_BGR = 0x0080

const PNG_TRANSFORM_SWAP_ALPHA = 0x0100

const PNG_TRANSFORM_SWAP_ENDIAN = 0x0200

const PNG_TRANSFORM_INVERT_ALPHA = 0x0400

const PNG_TRANSFORM_STRIP_FILLER = 0x0800

const PNG_TRANSFORM_STRIP_FILLER_BEFORE = PNG_TRANSFORM_STRIP_FILLER

const PNG_TRANSFORM_STRIP_FILLER_AFTER = 0x1000

const PNG_TRANSFORM_GRAY_TO_RGB = 0x2000

const PNG_TRANSFORM_EXPAND_16 = 0x4000

const PNG_TRANSFORM_SCALE_16 = 0x8000

const PNG_FLAG_MNG_EMPTY_PLTE = 0x01

const PNG_FLAG_MNG_FILTER_64 = 0x04

const PNG_ALL_MNG_FEATURES = 0x05

const PNG_ERROR_ACTION_NONE = 1

const PNG_ERROR_ACTION_WARN = 2

const PNG_ERROR_ACTION_ERROR = 3

const PNG_RGB_TO_GRAY_DEFAULT = -1

const PNG_ALPHA_PNG = 0

const PNG_ALPHA_STANDARD = 1

const PNG_ALPHA_ASSOCIATED = 1

const PNG_ALPHA_PREMULTIPLIED = 1

const PNG_ALPHA_OPTIMIZED = 2

const PNG_ALPHA_BROKEN = 3

const PNG_DEFAULT_sRGB = -1

const PNG_GAMMA_MAC_18 = -2

const PNG_GAMMA_sRGB = 220000

const PNG_GAMMA_LINEAR = PNG_FP_1

const PNG_FILLER_BEFORE = 0

const PNG_FILLER_AFTER = 1

const PNG_BACKGROUND_GAMMA_UNKNOWN = 0

const PNG_BACKGROUND_GAMMA_SCREEN = 1

const PNG_BACKGROUND_GAMMA_FILE = 2

const PNG_BACKGROUND_GAMMA_UNIQUE = 3

const PNG_GAMMA_THRESHOLD = PNG_GAMMA_THRESHOLD_FIXED * 1.0e-5

const PNG_CRC_DEFAULT = 0

const PNG_CRC_ERROR_QUIT = 1

const PNG_CRC_WARN_DISCARD = 2

const PNG_CRC_WARN_USE = 3

const PNG_CRC_QUIET_USE = 4

const PNG_CRC_NO_CHANGE = 5

const PNG_NO_FILTERS = 0x00

const PNG_FILTER_NONE = 0x08

const PNG_FILTER_SUB = 0x10

const PNG_FILTER_UP = 0x20

const PNG_FILTER_AVG = 0x40

const PNG_FILTER_PAETH = 0x80

const PNG_FAST_FILTERS = (PNG_FILTER_NONE | PNG_FILTER_SUB) | PNG_FILTER_UP

const PNG_ALL_FILTERS = (PNG_FAST_FILTERS | PNG_FILTER_AVG) | PNG_FILTER_PAETH

const PNG_FILTER_VALUE_NONE = 0

const PNG_FILTER_VALUE_SUB = 1

const PNG_FILTER_VALUE_UP = 2

const PNG_FILTER_VALUE_AVG = 3

const PNG_FILTER_VALUE_PAETH = 4

const PNG_FILTER_VALUE_LAST = 5

const PNG_FILTER_HEURISTIC_DEFAULT = 0

const PNG_FILTER_HEURISTIC_UNWEIGHTED = 1

const PNG_FILTER_HEURISTIC_WEIGHTED = 2

const PNG_FILTER_HEURISTIC_LAST = 3

const PNG_DESTROY_WILL_FREE_DATA = 1

const PNG_SET_WILL_FREE_DATA = 1

const PNG_USER_WILL_FREE_DATA = 2

const PNG_FREE_HIST = Cuint(0x0008)

const PNG_FREE_ICCP = Cuint(0x0010)

const PNG_FREE_SPLT = Cuint(0x0020)

const PNG_FREE_ROWS = Cuint(0x0040)

const PNG_FREE_PCAL = Cuint(0x0080)

const PNG_FREE_SCAL = Cuint(0x0100)

const PNG_FREE_UNKN = Cuint(0x0200)

const PNG_FREE_PLTE = Cuint(0x1000)

const PNG_FREE_TRNS = Cuint(0x2000)

const PNG_FREE_TEXT = Cuint(0x4000)

const PNG_FREE_EXIF = Cuint(0x8000)

const PNG_FREE_ALL = Cuint(0xffff)

const PNG_FREE_MUL = Cuint(0x4220)

const PNG_HANDLE_CHUNK_AS_DEFAULT = 0

const PNG_HANDLE_CHUNK_NEVER = 1

const PNG_HANDLE_CHUNK_IF_SAFE = 2

const PNG_HANDLE_CHUNK_ALWAYS = 3

const PNG_HANDLE_CHUNK_LAST = 4

const PNG_IO_NONE = 0x0000

const PNG_IO_READING = 0x0001

const PNG_IO_WRITING = 0x0002

const PNG_IO_SIGNATURE = 0x0010

const PNG_IO_CHUNK_HDR = 0x0020

const PNG_IO_CHUNK_DATA = 0x0040

const PNG_IO_CHUNK_CRC = 0x0080

const PNG_IO_MASK_OP = 0x000f

const PNG_IO_MASK_LOC = 0x00f0

const PNG_INTERLACE_ADAM7_PASSES = 7

const PNG_IMAGE_VERSION = 1

const PNG_IMAGE_WARNING = 1

const PNG_IMAGE_ERROR = 2

const PNG_FORMAT_FLAG_ALPHA = Cuint(0x01)

const PNG_FORMAT_FLAG_COLOR = Cuint(0x02)

const PNG_FORMAT_FLAG_LINEAR = Cuint(0x04)

const PNG_FORMAT_FLAG_COLORMAP = Cuint(0x08)

const PNG_FORMAT_FLAG_BGR = Cuint(0x10)

const PNG_FORMAT_FLAG_AFIRST = Cuint(0x20)

const PNG_FORMAT_FLAG_ASSOCIATED_ALPHA = Cuint(0x40)

const PNG_FORMAT_GRAY = 0

const PNG_FORMAT_GA = PNG_FORMAT_FLAG_ALPHA

const PNG_FORMAT_AG = PNG_FORMAT_GA | PNG_FORMAT_FLAG_AFIRST

const PNG_FORMAT_RGB = PNG_FORMAT_FLAG_COLOR

const PNG_FORMAT_BGR = PNG_FORMAT_FLAG_COLOR | PNG_FORMAT_FLAG_BGR

const PNG_FORMAT_RGBA = PNG_FORMAT_RGB | PNG_FORMAT_FLAG_ALPHA

const PNG_FORMAT_ARGB = PNG_FORMAT_RGBA | PNG_FORMAT_FLAG_AFIRST

const PNG_FORMAT_BGRA = PNG_FORMAT_BGR | PNG_FORMAT_FLAG_ALPHA

const PNG_FORMAT_ABGR = PNG_FORMAT_BGRA | PNG_FORMAT_FLAG_AFIRST

const PNG_FORMAT_LINEAR_Y = PNG_FORMAT_FLAG_LINEAR

const PNG_FORMAT_LINEAR_Y_ALPHA = PNG_FORMAT_FLAG_LINEAR | PNG_FORMAT_FLAG_ALPHA

const PNG_FORMAT_LINEAR_RGB = PNG_FORMAT_FLAG_LINEAR | PNG_FORMAT_FLAG_COLOR

const PNG_FORMAT_LINEAR_RGB_ALPHA = (PNG_FORMAT_FLAG_LINEAR | PNG_FORMAT_FLAG_COLOR) | PNG_FORMAT_FLAG_ALPHA

const PNG_FORMAT_RGB_COLORMAP = PNG_FORMAT_RGB | PNG_FORMAT_FLAG_COLORMAP

const PNG_FORMAT_BGR_COLORMAP = PNG_FORMAT_BGR | PNG_FORMAT_FLAG_COLORMAP

const PNG_FORMAT_RGBA_COLORMAP = PNG_FORMAT_RGBA | PNG_FORMAT_FLAG_COLORMAP

const PNG_FORMAT_ARGB_COLORMAP = PNG_FORMAT_ARGB | PNG_FORMAT_FLAG_COLORMAP

const PNG_FORMAT_BGRA_COLORMAP = PNG_FORMAT_BGRA | PNG_FORMAT_FLAG_COLORMAP

const PNG_FORMAT_ABGR_COLORMAP = PNG_FORMAT_ABGR | PNG_FORMAT_FLAG_COLORMAP

const PNG_IMAGE_FLAG_COLORSPACE_NOT_sRGB = 0x01

const PNG_IMAGE_FLAG_FAST = 0x02

const PNG_IMAGE_FLAG_16BIT_sRGB = 0x04

const PNG_MAXIMUM_INFLATE_WINDOW = 2

const PNG_SKIP_sRGB_CHECK_PROFILE = 4

const PNG_IGNORE_ADLER32 = 8

const PNG_OPTION_NEXT = 12

const PNG_OPTION_UNSET = 0

const PNG_OPTION_INVALID = 1

const PNG_OPTION_OFF = 2

const PNG_OPTION_ON = 3

