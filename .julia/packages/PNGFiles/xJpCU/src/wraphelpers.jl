png_error_handler(::Ptr{Cvoid}, msg::Cstring) = error("Png error: $(unsafe_string(msg))")
png_warn_handler(::Ptr{Cvoid}, msg::Cstring) = @warn("Png warn: $(unsafe_string(msg))")

# Compression strategy
const Z_DEFAULT_STRATEGY = 0
const Z_FILTERED = 1
const Z_HUFFMAN_ONLY = 2
const Z_RLE = 3
const Z_FIXED = 4

# Compression level
const Z_NO_COMPRESSION = 0
const Z_BEST_SPEED = 1
const Z_BEST_COMPRESSION = 9

# Returns the libpng version string
function get_libpng_version()
    ver = png_access_version_number()
    # Version in the format of xxyyzz, where x=major, yy=minor, z=release
    # But on the major version the first x is excluded if 0.
    dig = digits(ver)[end:-1:1]
    length(dig) == 5 ? prepend!(dig, 0) : error("Unknown libpng version: $ver")

    @inbounds major = dig[1] == 0 ? string(dig[2]) : string(dig[1], dig[2])
    @inbounds minor = dig[3] == 0 ? string(dig[4]) : string(dig[3], dig[4])
    @inbounds release = dig[5] == 0 ? string(dig[6]) : string(dig[5], dig[6])

    return "$major.$minor.$release"
end

function open_png(filename::String)
    fp = ccall(:fopen, Ptr{Cvoid}, (Cstring, Cstring), filename, "rb")
    fp == C_NULL && error("Failed to open $filename")

    header = zeros(UInt8, PNG_BYTES_TO_CHECK)
    header_size = ccall(:fread, Csize_t, (Ptr{UInt8}, Cint, Cint, Ptr{Cvoid}), header, 1, PNG_BYTES_TO_CHECK, fp)
    header_size != 8 && error("Failed to read header from $filename")

    is_png = png_sig_cmp(header, 0, PNG_BYTES_TO_CHECK)
    is_png != 0 && error("File $filename is not a png file")

    return fp
end

function create_read_struct()
    png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, C_NULL, png_error_fn_c[], png_warn_fn_c[])
    png_ptr == C_NULL && error("Failed to create png read struct")
    return png_ptr
end

function create_info_struct(png_ptr)
    info_ptr = png_create_info_struct(png_ptr)
    info_ptr == C_NULL && error("Failed to create png info struct")
    return info_ptr
end

close_png(fp::Ptr{Cvoid}) = ccall(:fclose, Cint, (Ptr{Cvoid},), fp)

# Write functions
function create_write_struct()
    png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, C_NULL, png_error_fn_c[], png_warn_fn_c[])
    png_ptr == C_NULL && error("Failed to create png write struct")
    return png_ptr
end
