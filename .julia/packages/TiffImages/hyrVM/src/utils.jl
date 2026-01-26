"""
    RemoteData

A placeholder type to describe the location and properties of remote data that
is too large to fit directly in a tag's spot in the IFD. Calling [`TiffImages.load!`](@ref)
on an IFD object replaces all `RemoteData`s with the respective data.

$(FIELDS)
"""
struct RemoteData{O <: Unsigned, T}
    """Position of this data in the stream"""
    position::O

    """The length of the data"""
    count::O

    function RemoteData(position::O, datatype::DataType, count::O) where {O}
        new{O, datatype}(position, count)
    end
end

Base.position(o::R) where {R <: RemoteData} = Int(o.pos)

const filepattern = r"<file (.*)>"
const fdpattern = r"<fd (.*)>"

"""
    extract_filename(io) -> String

Extract the name of the file backing a stream
"""
function extract_filename(io::IOStream)
    name = String(io.name)
    filename = match(filepattern, name)
    if filename !== nothing
        return String(filename[1])
    elseif match(fdpattern, name) !== nothing
        return ""
    else
        error("Can't extract filename from the given stream")
    end
end
extract_filename(io::Stream) = io.filename

"""
    check_bswap(io::Stream)

Check endianness of TIFF file to see if we need to swap bytes
"""
function need_bswap(io::Stream)
    seekstart(io)
    endianness = String(read(io, 2))
    # check if we need to swap byte order
    bswap = endianness != (ENDIAN_BOM == 0x04030201 ? "II" : "MM")
    return bswap
end

function offsetsize(io::Stream)
    seek(io, 2)
    tiff_version = Array{UInt8}(undef, 2)
    read!(io, tiff_version)
    if tiff_version == [0x2a, 0x00] || tiff_version == [0x00, 0x2a]
        return UInt32
    elseif tiff_version == [0x2b, 0x00] || tiff_version == [0x00, 0x2b]
        offset_bytesize = read(io, UInt16)
        constant = read(io, UInt16)
        return UInt64
    else
        error("Unknown TIFF version")
    end
end

# using rephorm's table from
# https://github.com/rephorm/TIFF.jl/blob/master/tiff.jl#L30
const tiff_to_julian = Dict(
  0x0001 => UInt8,            # CHAR
  0x0002 => String,           # ASCII
  0x0003 => UInt16,           # SHORT
  0x0004 => UInt32,           # LONG
  0x0005 => Rational{UInt32}, # RATIONAL
  0x0006 => Int8,             # SBYTE
  0x0007 => Any,              # UNDEFINED
  0x0008 => Int16,            # SSHORT
  0x0009 => Int32,            # SLONG
  0x000a => Rational{Int32},  # SRATIONAL
  0x000b => Float32,          # FLOAT
  0x000c => Float64,          # DOUBLE
  0x000d => UInt32,           # IFD
  0x0010 => UInt64,           # LONG8
  0x0011 => Int64,            # SLONG8
  0x0012 => UInt64,           # IFD8
)

const julian_to_tiff = Dict(
    UInt8               => 0x0001,   # CHAR
    String              => 0x0002,   # ASCII
    UInt16              => 0x0003,   # SHORT
    UInt32              => 0x0004,   # LONG
    Rational{UInt32}    => 0x0005,   # RATIONAL
    Int8                => 0x0006,   # SBYTE
    Any                 => 0x0007,   # UNDEFINED
    Int16               => 0x0008,   # SSHORT
    Int32               => 0x0009,   # SLONG
    Rational{Int32}     => 0x000a,   # SRATIONAL
    Float32             => 0x000b,   # FLOAT
    Float64             => 0x000c,   # DOUBLE
    UInt64              => 0x0010,   # LONG8
    Int64               => 0x0011,   # SLONG8
)

_bswap(a) = bswap(a)
_bswap(c::Colorant{T, N}) where {T, N} = mapc(bswap, c)
_bswap(p::WidePixel{C, X}) where {C, X} = WidePixel(_bswap(p.color), _bswap(p.extra))
_bswap(t::Tuple) = map(bswap, t)

function getstream(fmt, io, name)
    # adapted from https://github.com/JuliaStats/RDatasets.jl/pull/119/
    if isdefined(FileIO, :action)
        # FileIO >= 1.6
        return Stream{fmt}(io, name)
    else
        # FileIO < 1.6
        return Stream(fmt, io, name)
    end
end

getstream(fmt, io::IOBuffer) = getstream(fmt, io, "")
getstream(fmt, io::IOStream) = getstream(fmt, io, extract_filename(io))
# assume OMETIFF if no format given
getstream(io) = getstream(format"TIFF", io)

@static if Sys.iswindows()
    # Be permissive on windows with eager GC to work around
    # https://github.com/tlnagy/TiffImages.jl/pull/79#discussion_r880478304
    function _safe_open(f, filepath::String, mode="r", args...; kwargs...)
        GC.gc()
        try
            open(f, filepath, mode, args...; kwargs...)
        catch err
            if err isa SystemError
                @warn "failed to open file \"$filepath\" in \"$mode\" mode, this may be caused by overwriting a file previously opened with mmap"
            end
            rethrow()
        end
    end
else
    const _safe_open = open
end

function memcpy(dest::Ptr{T}, src::Ptr{T}, n::Int) where T
    ccall(:memcpy, Ptr{T}, (Ptr{T}, Ptr{T}, Int), dest, src, n)
end
