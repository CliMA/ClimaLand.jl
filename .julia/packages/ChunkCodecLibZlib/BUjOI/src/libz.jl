# Constants and c wrapper functions ported to Julia from zlib.h https://github.com/madler/zlib/blob/v1.3.1/zlib.h

const ZLIB_VERSION = "1.3.1"

# constants
const Z_NO_FLUSH      = 0
const Z_PARTIAL_FLUSH = 1
const Z_SYNC_FLUSH    = 2
const Z_FULL_FLUSH    = 3
const Z_FINISH        = 4
const Z_BLOCK         = 5
const Z_TREES         = 6
# Allowed flush values; see deflate() and inflate() below for details 

const Z_OK            = 0
const Z_STREAM_END    = 1
const Z_NEED_DICT     = 2
const Z_ERRNO         = (-1)
const Z_STREAM_ERROR  = (-2)
const Z_DATA_ERROR    = (-3)
const Z_MEM_ERROR     = (-4)
const Z_BUF_ERROR     = (-5)
const Z_VERSION_ERROR = (-6)
# Return codes for the compression/decompression functions. Negative values
# are errors, positive values are used for special but normal events.

const Z_NO_COMPRESSION       = 0
const Z_BEST_SPEED           = 1
const Z_BEST_COMPRESSION     = 9
const Z_DEFAULT_COMPRESSION  = (-1)
# compression levels

const Z_FILTERED            = 1
const Z_HUFFMAN_ONLY        = 2
const Z_RLE                 = 3
const Z_FIXED               = 4
const Z_DEFAULT_STRATEGY    = 0
# compression strategy; see deflateInit2() below for details

const Z_DEFLATED   = 8
# The deflate compression method (the only one supported in this version)

@assert typemax(Csize_t) ≥ typemax(Cuint)

function zalloc(::Ptr{Cvoid}, items::Cuint, size::Cuint)::Ptr{Cvoid}
    s, f = Base.Checked.mul_with_overflow(items, size)
    if f
        C_NULL
    else
        ccall(:jl_malloc, Ptr{Cvoid}, (Csize_t,), s%Csize_t)
    end
end
zfree(::Ptr{Cvoid}, address::Ptr{Cvoid}) = ccall(:jl_free, Cvoid, (Ptr{Cvoid},), address)

mutable struct ZStream
    next_in::Ptr{Cuchar}  # next input byte
    avail_in::Cuint       # number of bytes available at next_in
    total_in::Culong      # total number of input bytes read so far

    next_out::Ptr{Cuchar} # next output byte will go here
    avail_out::Cuint      # remaining free space at next_out
    total_out::Culong     # total number of bytes output so far

    msg::Ptr{Cchar}       # last error message, NULL if no error
    state::Ptr{Cvoid}     # not visible by applications

    zalloc::Ptr{Cvoid}    # used to allocate the internal state
    zfree::Ptr{Cvoid}     # used to free the internal state
    opaque::Ptr{Cvoid}    # private data object passed to zalloc and zfree

    data_type::Cint       # best guess about the data type: binary or text
                          # for deflate, or the decoding state for inflate
    adler::Culong         # Adler-32 or CRC-32 value of the uncompressed data
    reserved::Culong      # reserved for future use

    function ZStream()
        new(
            C_NULL, 0, 0,
            C_NULL, 0, 0,
            C_NULL, C_NULL,
            @cfunction(zalloc, Ptr{Cvoid}, (Ptr{Cvoid}, Cuint, Cuint)),
            @cfunction(zfree, Cvoid, (Ptr{Cvoid}, Ptr{Cvoid})),
            C_NULL,
            0,
            0, 0,
        )
    end
end

function deflateInit2(stream::ZStream, level::Cint, windowBits::Cint)
    memLevel = Cint(8) # default
    ret = ccall(
        (:deflateInit2_, libz),
        Cint,
        (Ref{ZStream}, Cint, Cint, Cint, Cint, Cint, Cstring, Cint),
        stream, level, Z_DEFLATED, windowBits, memLevel, Z_DEFAULT_STRATEGY, ZLIB_VERSION, sizeof(ZStream),
    )
    if ret != Z_OK
        if ret == Z_MEM_ERROR
            throw(OutOfMemoryError())
        elseif ret == Z_STREAM_ERROR
            error("Z_STREAM_ERROR: invalid parameter") # this should be validated before this
        elseif ret == Z_VERSION_ERROR
            error("Z_VERSION_ERROR: zlib version is incompatible")
        else
            error("Unknown zlib error code: $(ret)")
        end
    end
    nothing
end

function deflateEnd(stream::ZStream)
    # free libz stream state if needed
    if stream.state != C_NULL
        ccall(
            (:deflateEnd, libz),
            Cint,
            (Ref{ZStream},),
            stream,
        )
    end
    nothing
end

function inflateInit2(stream::ZStream, windowBits::Cint)
    # indicate that no input data is being provided for future zlib compat
    stream.next_in = C_NULL
    stream.avail_in = 0
    ret = ccall(
        (:inflateInit2_, libz),
        Cint,
        (Ref{ZStream}, Cint, Cstring, Cint),
        stream, windowBits, ZLIB_VERSION, sizeof(ZStream),
    )
    if ret != Z_OK
        if ret == Z_MEM_ERROR
            throw(OutOfMemoryError())
        elseif ret == Z_STREAM_ERROR
            error("Z_STREAM_ERROR: invalid parameter")
        elseif ret == Z_VERSION_ERROR
            error("Z_VERSION_ERROR: zlib version is incompatible")
        else
            error("Unknown zlib error code: $(ret)")
        end
    end
    nothing
end

function inflateEnd(stream::ZStream)
    # free stream state if needed
    if stream.state != C_NULL
        ccall(
            (:inflateEnd, libz),
            Cint,
            (Ref{ZStream},),
            stream,
        )
    end
    nothing
end


# The following is the original license info from zlib.h

#= zlib.h -- interface of the 'zlib' general purpose compression library
  version 1.3.1, January 22nd, 2024

  Copyright (C) 1995-2024 Jean-loup Gailly and Mark Adler

  This software is provided 'as-is', without any express or implied
  warranty.  In no event will the authors be held liable for any damages
  arising from the use of this software.

  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:

  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software. If you use this software
     in a product, an acknowledgment in the product documentation would be
     appreciated but is not required.
  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.
  3. This notice may not be removed or altered from any source distribution.

  Jean-loup Gailly        Mark Adler
  jloup@gzip.org          madler@alumni.caltech.edu


  The data format used by the zlib library is described by RFCs (Request for
  Comments) 1950 to 1952 in the files http://tools.ietf.org/html/rfc1950
  (zlib format), rfc1951 (deflate format) and rfc1952 (gzip format).
=#
