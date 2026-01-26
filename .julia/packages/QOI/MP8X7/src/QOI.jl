module QOI

using FixedPointNumbers
using ColorTypes
using FileIO


#############
# Constants #
#############

const QOI_OP_INDEX = 0x00 # 00xxxxxx
const QOI_OP_DIFF  = 0x40 # 01xxxxxx
const QOI_OP_LUMA  = 0x80 # 10xxxxxx
const QOI_OP_RUN   = 0xc0 # 11xxxxxx
const QOI_OP_RGB   = 0xfe # 11111110
const QOI_OP_RGBA  = 0xff # 11111111
const QOI_MASK_2   = 0xc0 # 11000000
const QOI_MAGIC = UInt32('q') << 24 | UInt32('o') << 16 | UInt32('i') << 8 | UInt32('f')
const QOI_PIXELS_MAX = 400000000
const QOI_PADDING = (0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01)


#########
# Pixel #
#########

struct Pixel
    r::UInt8
    g::UInt8
    b::UInt8
    a::UInt8
    function Pixel(r::UInt8, g::UInt8, b::UInt8, a::UInt8)
        new(r, g, b, a)
    end
end

@inline _qoi_color_hash(p::Pixel) = p.r*3 + p.g*5 + p.b*7 + p.a*11


##############  
# Exceptions #
##############

struct QOIException <: Exception
    msg::String
end
Base.showerror(io::IO, qoi::QOIException) = print(io, qoi.msg)

@noinline throw_magic_bytes_error(magic::UInt32)             = throw(QOIException("invalid magic bytes, got $(repr(magic)), expected $(repr(QOI_MAGIC))"))
@noinline throw_invalid_header_width(width::UInt32)          = throw(QOIException("invalid width in header, got $width"))
@noinline throw_invalid_header_height(height::UInt32)        = throw(QOIException("invalid height in header, got $height"))
@noinline throw_invalid_header_channels(channels::UInt8)     = throw(QOIException("invalid channels in header, got $channels"))
@noinline throw_invalid_header_colorspace(colorspace::UInt8) = throw(QOIException("invalid colorspace in header, got $colorspace"))
@noinline throw_unexpected_eof()                             = throw(QOIException("unexpected end of file"))


#############
# QOIHeader #
#############

@enum QOIChannel::UInt8 begin
    QOI_RGB  = 0x03
    QOI_RGBA = 0x04
end

@enum QOIColorSpace::UInt8 begin
    QOI_SRGB   = 0x00
    QOI_LINEAR = 0x01
end

struct QOIHeader
    width::UInt32
    height::UInt32
    channels::QOIChannel
    colorspace::QOIColorSpace
end
function QOIHeader(width::UInt32, height::UInt32, channels::UInt8, colorspace::UInt8)
    width == 0                   && throw_invalid_header_width(width)
    height == 0                  && throw_invalid_header_height(height)
    channels < 3 || channels > 4 && throw_invalid_header_channels(channels)
    colorspace > 1               && throw_invalid_header_colorspace(colorspace)
    return QOIHeader(width, height, QOIChannel(channels), QOIColorSpace(colorspace))
end


############
# Encoding #
############

mutable struct QOIWriter{V <: AbstractVecOrMat{UInt8}}
    v::V
    pos::Int
end
QOIWriter(v::AbstractVecOrMat{UInt8}) = QOIWriter(v, 0)

@inline function _qoi_write!(qoiw::QOIWriter, v::UInt8)
    qoiw.pos += 1 
    qoiw.pos > length(qoiw.v) && resize!(qoiw.v, max(1, length(qoiw.v) * 2))
    @inbounds qoiw.v[qoiw.pos] = v
end

function _qoi_write_32!(qoiw::QOIWriter, v::UInt32)
    _qoi_write!(qoiw, ((0xff000000 & v) >> 24) % UInt8)
    _qoi_write!(qoiw, ((0x00ff0000 & v) >> 16) % UInt8)
    _qoi_write!(qoiw, ((0x0000ff00 & v) >> 8)  % UInt8)
    _qoi_write!(qoiw, ((0x000000ff & v))       % UInt8)
    return 
end

qoi_encode_raw(image::AbstractVecOrMat{UInt8}, header::QOIHeader) =
    qoi_encode_raw!(Vector{UInt8}(undef, 256), image, header)

function qoi_encode_raw(io::IO, image::AbstractVecOrMat{UInt8}, header::QOIHeader)
    data = qoi_encode_raw(image, header)
    write(io, data)
end

function qoi_encode_raw!(data::AbstractVector{UInt8}, image::AbstractVecOrMat{UInt8}, header::QOIHeader)
    # TODO: Check size of data against QOI_PIXELS_MAX?

    qoiw = QOIWriter(data)

    # Header
    _qoi_write_32!(qoiw, QOI_MAGIC)
    _qoi_write_32!(qoiw, header.width)
    _qoi_write_32!(qoiw, header.height)
    _qoi_write!(qoiw, header.channels |> Integer)
    _qoi_write!(qoiw, header.colorspace |> Integer)

    index = fill(Pixel(0x00, 0x00, 0x00, 0x00), 64)
    run = 0x00
    px_prev = Pixel(0x00, 0x00, 0x00, 0xff)
    px = px_prev
    
    channels = Integer(header.channels)
    px_len = header.width * header.height * channels
    px_end = px_len - channels + 1

    # Data
    for px_pos in 1:channels:px_len
        r = image[px_pos + 0]
        g = image[px_pos + 1]
        b = image[px_pos + 2]
        a = header.channels == QOI_RGBA ? image[px_pos + 3] : 0xff
        px = Pixel(r, g, b, a)

        if px == px_prev
            run += 0x01
            if run == 62 || px_pos == px_end
                _qoi_write!(qoiw, QOI_OP_RUN |  (run-0x01))
                run = 0x00
            end
        else
            if run > 0
                _qoi_write!(qoiw, QOI_OP_RUN | (run-0x01))
                run = 0x00
            end

            index_pos = mod1(_qoi_color_hash(px)+1, 64) % UInt8
            if (@inbounds index[index_pos]) == px
                _qoi_write!(qoiw, QOI_OP_INDEX | (index_pos - 0x01))
            else
                @inbounds index[index_pos] = px
                if px.a == px_prev.a
                    vr = ((px.r) - (px_prev.r)) % Int8
                    vg = ((px.g) - (px_prev.g)) % Int8
                    vb = ((px.b) - (px_prev.b)) % Int8

                    vg_r = vr - vg
                    vg_b = vb - vg
                    if      vr > -3 && vr < 2 &&
                            vg > -3 && vg < 2 &&
                            vb > -3 && vb < 2
                        _qoi_write!(qoiw, QOI_OP_DIFF | ((vr + 0x02) % UInt8) << 4 | ((vg + 0x02) % UInt8) << 2 | (vb + 0x02)  % UInt8)
                    elseif  vg_r > -9 && vg_r < 8  &&
                            vg > -33  && vg   < 32 &&
                            vg_b > -9 && vg_b < 8
                        _qoi_write!(qoiw, QOI_OP_LUMA   | (vg + UInt8(32)) % UInt8)
                        _qoi_write!(qoiw, ((vg_r + 0x08) % UInt8) << 4 | (vg_b + 0x08) % UInt8)
                    else
                        _qoi_write!(qoiw, QOI_OP_RGB)
                        _qoi_write!(qoiw, px.r)
                        _qoi_write!(qoiw, px.g)
                        _qoi_write!(qoiw, px.b)
                    end
                else
                    _qoi_write!(qoiw, QOI_OP_RGBA)
                    _qoi_write!(qoiw, px.r)
                    _qoi_write!(qoiw, px.g)
                    _qoi_write!(qoiw, px.b)
                    _qoi_write!(qoiw, px.a)
                end
            end
        end
        px_prev = px
    end

    # Padding
    for x in QOI_PADDING
        _qoi_write!(qoiw, x)
    end

    sizehint!(data, qoiw.pos)
    resize!(data, qoiw.pos)
    return data
end

function qoi_encode(file::String, image::AbstractMatrix{T}) where T <: Colorant
    if T <: TransparentColor
        if T != RGBA{N0f8}
            image = RGBA{N0f8}.(image)
        end
        channel = QOI_RGBA
    else
        if T != RGB{N0f8}
            image = RGB{N0f8}.(image)
        end 
        channel = QOI_RGB
    end
    header = QOIHeader(size(image, 2), size(image, 1), channel, QOI_SRGB)
    image = permutedims(image)
    open(file, "w") do io
        image_raw = reinterpret(UInt8, image)
        qoi_encode_raw(io, image_raw, header)
    end
end


############
# Decoding #
############

mutable struct QOIReader{V <: AbstractVecOrMat{UInt8}}
    v::V
    pos::Int
end
QOIReader(v::AbstractVecOrMat{UInt8}) = QOIReader(v, 0)

@inline _qoi_read!(qoir::QOIReader) = (qoir.pos+=1; @inbounds qoir.v[qoir.pos])

function _qoi_read_32!(qoir::QOIReader)
    a = UInt32(_qoi_read!(qoir))
    b = UInt32(_qoi_read!(qoir))
    c = UInt32(_qoi_read!(qoir))
    d = UInt32(_qoi_read!(qoir))
    return a << 24 | b << 16 | c << 8 | d
end

@inline _qoi_read_rgb!(qoir::QOIReader) = (_qoi_read!(qoir), _qoi_read!(qoir), _qoi_read!(qoir))
@inline _qoi_read_rgba!(qoir::QOIReader) = (_qoi_read!(qoir), _qoi_read!(qoir), _qoi_read!(qoir), _qoi_read!(qoir))

function qoi_decode_raw(v::AbstractVector{UInt8})
    qoir = QOIReader(v)

    if length(v) < 14 # 4 magic + 10 header
        throw_unexpected_eof()
    end

    # Magic
    magic = _qoi_read_32!(qoir)
    magic == QOI_MAGIC || throw_magic_bytes_error(magic)

    # Header
    width = _qoi_read_32!(qoir)
    height = _qoi_read_32!(qoir)
    channels = _qoi_read!(qoir)
    colorspace = _qoi_read!(qoir)
    header = QOIHeader(width, height, channels, colorspace)

    # Data
    n_pixels = header.width * header.height
    n_values = n_pixels * channels
    data = Vector{UInt8}(undef, n_values)
    index = fill(Pixel(0x00, 0x00, 0x00, 0x00), 64)
    px = Pixel(0x00, 0x00, 0x00, 0xFF)
    px_idx = 1
    run = 0x00
    lenv = length(v)

    for px_idx in 1:channels:n_values
        if qoir.pos > lenv - length(QOI_PADDING)
            throw_unexpected_eof()
        end
        if run > 0
            run -= 0x01
        else
            b1 = _qoi_read!(qoir)
            if b1 == QOI_OP_RGB
                r, g, b = _qoi_read_rgb!(qoir)     
                px = Pixel(r, g, b, px.a)
            elseif b1 == QOI_OP_RGBA
                r, g, b, a = _qoi_read_rgba!(qoir)
                px = Pixel(r, g, b, a)
            elseif b1 & QOI_MASK_2 == QOI_OP_INDEX
                px = index[b1+0x01]
            elseif (b1 & QOI_MASK_2) == QOI_OP_DIFF
                r = px.r + ((b1 >> 0x04) & 0x03) - 0x02
                g = px.g + ((b1 >> 0x02) & 0x03) - 0x02
                b = px.b + ( b1          & 0x03) - 0x02
                px = Pixel(r, g, b, px.a)
            elseif ((b1 & QOI_MASK_2) == QOI_OP_LUMA)
                b2 = _qoi_read!(qoir)
                vg = (b1 & 0x3f) - UInt8(32)
                r = px.r + vg - 0x08 + ((b2 >> 4) & 0x0f)
                g = px.g + vg
                b = px.b + vg - 0x08 +  (b2       & 0x0f)
                px = Pixel(r, g, b, px.a)
            elseif (b1 & QOI_MASK_2) == QOI_OP_RUN
                run = (b1 & 0x3f)
            else
                error("unreachable")
            end
            @inbounds index[mod1(_qoi_color_hash(px)+1, 64)] = px
        end

        @inbounds data[px_idx+0] = px.r
        @inbounds data[px_idx+1] = px.g
        @inbounds data[px_idx+2] = px.b
        if header.channels == QOI_RGBA
            @inbounds data[px_idx+3] = px.a
        end
    end

    if qoir.pos != length(v) - length(QOI_PADDING)
        throw_unexpected_eof()
    end

    # Read padding
    for _ in 1:7
        x = _qoi_read!(qoir) 
        x == 0 || throw_unexpected_eof()
    end  
    x = _qoi_read!(qoir)
    x == 1 || throw_unexpected_eof()

    image = reshape(data, Int(header.width) * channels, Int(header.height))
    return header, image
end

function _to_colortype(header, raw_image)
    T = header.channels == QOI_RGBA ? RGBA{N0f8} : RGB{N0f8} 
    return permutedims(reinterpret(T, raw_image))
end

function qoi_decode(v::AbstractVector{UInt8})
    header, raw_image = qoi_decode_raw(v)
    return _to_colortype(header, raw_image)
end

qoi_decode_raw(f::Union{String, IO}) = qoi_decode_raw(Base.read(f))
qoi_decode(f::Union{String, IO}) = qoi_decode(Base.read(f))


##########
# FileIO #
##########

load(f::File{format"QOI"}) = qoi_decode(f.filename)
save(f::File{format"QOI"}, image::AbstractMatrix{<:Colorant}) = qoi_encode(f.filename, image)

end
