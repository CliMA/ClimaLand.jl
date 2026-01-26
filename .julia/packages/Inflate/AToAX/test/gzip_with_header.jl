# Use CodecZlib plus one extra ccall into zlib to generate a gzip with
# header information. By default zlib leaves the header empty.
#
# Note: this is only intended for testing of the Inflate
# implementation of Gzip decompression.
#
# See zlib.h in the zlib libary for documentation of the fields in
# this struct.
mutable struct GzHeader
    text::Cint
    time::Culong
    xflags::Cint
    os::Cint
    extra::Ptr{Cuchar}
    extra_len::Cuint
    extra_max::Cuint
    name::Ptr{Cuchar}
    name_max::Cuint
    comment::Ptr{Cuchar}
    comm_max::Cuint
    hcrc::Cint
    done::Cint
end

function GzHeader(mtime, os, extra, name, comment, include_header_crc)
    return GzHeader(true, mtime, 0, os,
                    isempty(extra) ? 0 : pointer(extra), length(extra), 0,
                    isempty(name) ? 0 : pointer(name), 0,
                    isempty(comment) ? 0 : pointer(comment), 0,
                    include_header_crc, 0)
end

function gzip_with_header(text, mtime, os, extra, name, comment,
                          include_header_crc)
    zstream = CodecZlib.ZStream()
    CodecZlib.deflate_init!(zstream, 9, 15 + 16)
    data = Vector{UInt8}(text)
    zstream.next_in = pointer(data)
    zstream.avail_in = length(data)
    out = zeros(UInt8, 1000)
    zstream.next_out = pointer(out)
    zstream.avail_out = length(out)
    gzheader = GzHeader(mtime, os, extra, name, comment, include_header_crc)
    ccall((:deflateSetHeader, CodecZlib.libz), Cint,
          (Ref{CodecZlib.ZStream}, Ref{GzHeader}),
          zstream, gzheader)
    CodecZlib.deflate!(zstream, CodecZlib.Z_FINISH)
    CodecZlib.deflate_end!(zstream)
    return out[1:(length(out) - zstream.avail_out)]
end
