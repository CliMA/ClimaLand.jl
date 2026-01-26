######################################################################################################
# macros from jpeglib.h

jpeg_create_compress(cinfo) =
    jpeg_CreateCompress(cinfo, JPEG_LIB_VERSION, sizeof(jpeg_compress_struct))
jpeg_create_decompress(cinfo) =
    jpeg_CreateDecompress(cinfo, JPEG_LIB_VERSION, sizeof(jpeg_decompress_struct))
