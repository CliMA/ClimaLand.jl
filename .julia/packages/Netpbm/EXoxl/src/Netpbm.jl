module Netpbm

using FileIO, ImageCore, ImageMetadata

# Note: there is no endian standard, but netpbm is big-endian
const is_little_endian = ENDIAN_BOM == 0x04030201
const ufixedtype = Dict(10=>N6f10, 12=>N4f12, 14=>N2f14, 16=>N0f16)

if VERSION < v"1.1"
    isnothing(x) = x === nothing
end

function load(f::Union{File{format"PBMBinary"},File{format"PGMBinary"},File{format"PPMBinary"}})
    open(f) do s
        skipmagic(s)
        load(s)
    end
end

function load(f::Union{File{format"PBMText"},File{format"PGMText"},File{format"PPMText"}})
    open(f) do s
        skipmagic(s)
        load(s)
    end
end

function load(s::Stream{format"PBMBinary"})
    io = stream(s)
    w, h = parse_netpbm_size(io)
    dat = BitArray(undef,w, h)
    nbytes_per_row = ceil(Int, w/8)
    for irow = 1:h, j = 1:nbytes_per_row
        tmp = read(io, UInt8)
        offset = (j-1)*8
        for k = 1:min(8, w-offset)
            dat[offset+k, irow] = (tmp>>>(8-k))&0x01
        end
    end
    PermutedDimsArray(dat, (2,1))
end

function load(s::Stream{format"PGMBinary"})
    io = stream(s)
    w, h = parse_netpbm_size(io)
    maxval = parse_netpbm_maxval(io)
    if maxval <= 255
        dat8 = Array{UInt8}(undef, h, w)
        readio!(io, dat8)
        return reinterpret(Gray{N0f8}, dat8)
    elseif maxval <= typemax(UInt16)
        datraw = Array{UInt16}(undef, h, w)
        readio!(io, datraw)
        # Determine the appropriate Normed type
        T = ufixedtype[ceil(Int, log2(maxval)/2)<<1]
        return reinterpret(Gray{T}, datraw)
    else
        error("Image file may be corrupt. Are there really more than 16 bits in this image?")
    end
end

function load(s::Stream{format"PPMBinary"})
    io = stream(s)
    w, h = parse_netpbm_size(io)
    maxval = parse_netpbm_maxval(io)
    local dat
    if maxval <= 255
        dat8 = Array{UInt8}(undef, 3, h, w)
        readio!(io, dat8)
        return reshape(reinterpret(RGB{N0f8}, dat8), (h, w))
    elseif maxval <= typemax(UInt16)
        datraw = Array{UInt16}(undef, 3, h, w)
        readio!(io, datraw)
        # Determine the appropriate Normed type
        T = ufixedtype[ceil(Int, log2(maxval)/2)<<1]
        return reshape(reinterpret(RGB{T}, datraw), (h, w))
    else
        error("Image file may be corrupt. Are there really more than 16 bits in this image?")
    end
end

function load(s::Stream{format"PBMText"})
    io = stream(s)
    w, h = parse_netpbm_size(io)
    dat = BitArray(undef, h, w)
    nbytes_per_row = ceil(Int, w/8)
    readtextio!(io, dat)
    return(dat)
end

function load(s::Stream{format"PGMText"})
    io = stream(s)
    w, h = parse_netpbm_size(io)
    maxval = parse_netpbm_maxval(io)
    if maxval <= 255
        dat8 = Array{UInt8}(undef, h, w)
        readtextio!(io, dat8)
        return reinterpret(Gray{N0f8}, dat8)
    elseif maxval <= typemax(UInt16)
        datraw = Array{UInt16}(undef, h, w)
        readtextio!(io, datraw)
        # Determine the appropriate Normed type
        T = ufixedtype[ceil(Int, log2(maxval)/2)<<1]
        return reinterpret(Gray{T}, datraw)
    else
        error("Image file may be corrupt. Are there really more than 16 bits in this image?")
    end
end

function load(s::Stream{format"PPMText"})
    io = stream(s)
    w, h = parse_netpbm_size(io)
    maxval = parse_netpbm_maxval(io)
    local dat
    if maxval <= 255
        dat8 = Array{UInt8}(undef, 3, h, w)
        readtextio!(io, dat8)
        return reshape(reinterpret(RGB{N0f8}, dat8), (h, w))
    elseif maxval <= typemax(UInt16)
        datraw = Array{UInt16}(undef, 3, h, w)
        readtextio!(io, datraw)
        # Determine the appropriate Normed type
        T = ufixedtype[ceil(Int, log2(maxval)/2)<<1]
        return reshape(reinterpret(RGB{T}, datraw), (h, w))
    else
        error("Image file may be corrupt. Are there really more than 16 bits in this image?")
    end
end

@noinline function readio!(io, dat::AbstractMatrix{T}) where {T}
    axh, axw = axes(dat)
    for i = axh, j = axw  # io is stored in row-major format
        dat[i,j] = default_swap(read(io, T))
    end
    dat
end

@noinline function readio!(io, dat::AbstractArray{T,3}) where {T}
    size(dat, 1) == 3 || throw(DimensionMismatch("must be of size 3 in first dimension, got $(size(dat, 1))"))
    axh, axw = axes(dat, 2), axes(dat, 3)
    for i = axh, j = axw, k = 1:3  # io is stored row-major, color-first
        dat[k,i,j] = default_swap(read(io, T))
    end
    dat
end

@noinline function readtextio!(io, dat::AbstractMatrix{T}) where {T}
    axh, axw = axes(dat)
    for i = axh, j = axw  # io is stored in row-major format
        dat[i,j] = parsenextint(io)
    end
    dat
end

@noinline function readtextio!(io, dat::AbstractArray{T,3}) where {T}
    size(dat, 1) == 3 || throw(DimensionMismatch("must be of size 3 in first dimension, got $(size(dat, 1))"))
    axh, axw = axes(dat, 2), axes(dat, 3)
    for i = axh, j = axw, k = 1:3  # io is stored row-major, color-first
        dat[k,i,j] = parsenextint(io)
    end
    dat
end

function save(filename::File{format"PBMBinary"}, img; mapf=identity, mapi=nothing)
    mapf = kwrename(:mapf, mapf, :mapi, mapi, :save)
    comment = _read_metadata(img)
    open(filename, "w") do s
        io = stream(s)
        write(io, "P4\n")
        write(io, isnothing(comment) ? "# pbm file written by Julia" : comment, "\n")
        save(s, img, mapf=mapf)
    end
end

function save(filename::File{format"PGMBinary"}, img; mapf=identity, mapi=nothing)
    mapf = kwrename(:mapf, mapf, :mapi, mapi, :save)
    comment = _read_metadata(img)
    open(filename, "w") do s
        io = stream(s)
        write(io, "P5\n")
        write(io, isnothing(comment) ? "# pgm file written by Julia" : comment, "\n")
        save(s, img, mapf=mapf)
    end
end

function save(filename::File{format"PPMBinary"}, img; mapf=identity, mapi=nothing)
    mapf = kwrename(:mapf, mapf, :mapi, mapi, :save)
    comment = _read_metadata(img)
    open(filename, "w") do s
        io = stream(s)
        write(io, "P6\n")
        write(io, isnothing(comment) ? "# ppm file written by Julia" : comment, "\n")
        save(s, img, mapf=mapf)
    end
end

function save(filename::File{format"PBMText"}, img; mapf=identity, mapi=nothing)
    mapf = kwrename(:mapf, mapf, :mapi, mapi, :save)
    comment = _read_metadata(img)
    open(filename, "w") do s
        io = stream(s)
        write(io, "P1\n")
        write(io, isnothing(comment) ? "# pbm file written by Julia" : comment, "\n")
        save(s, img, mapf=mapf)
    end
end

function save(filename::File{format"PGMText"}, img; mapf=identity, mapi=nothing)
    mapf = kwrename(:mapf, mapf, :mapi, mapi, :save)
    comment = _read_metadata(img)
    open(filename, "w") do s
        io = stream(s)
        write(io, "P2\n")
        write(io, isnothing(comment) ? "# pbm file written by Julia" : comment, "\n")
        save(s, img, mapf=mapf)
    end
end

function save(filename::File{format"PPMText"}, img; mapf=identity, mapi=nothing)
    mapf = kwrename(:mapf, mapf, :mapi, mapi, :save)
    comment = _read_metadata(img)
    open(filename, "w") do s
        io = stream(s)
        write(io, "P3\n")
        write(io, isnothing(comment) ? "# ppm file written by Julia" : comment, "\n")
        save(s, img, mapf=mapf)
    end
end

function save(s::Stream, img::AbstractMatrix; mapf=identity, mapi=nothing)
    mapf = kwrename(:mapf, mapf, :mapi, mapi, :save)
    save(s, img, mapf)
end

@noinline function save(s::Stream{format"PBMBinary"}, img::AbstractMatrix{T}, mapf) where {T<:Number}
    axh, axw = axes(img)
    w = length(axw)
    write(s, "$(w) $(length(axh))\n")
    dat = PermutedDimsArray(img, (2,1))
    nbytes_per_row = ceil(Int, w/8)
    for irow = axh, j = 1:nbytes_per_row
        tmp = UInt8(0)
        offset = (j-1)*8
        for k = 1:min(8, w-offset)
            tmp |= round(UInt8, mapf(dat[offset+k, irow]) << (8-k))
        end
        write(s, UInt8(tmp))
    end
    nothing
end

@noinline function save(s::Stream{format"PGMBinary"}, img::AbstractMatrix{T}, mapf) where {T<:Union{Gray,Number}}
    axh, axw = axes(img)
    Tout, mx = pnmmax(img)
    if sizeof(Tout) > 2
        error("element type $Tout (from $T) not supported")
    end
    write(s, "$(length(axw)) $(length(axh))\n$mx\n")
    for i = axh, j = axw  # s is stored in row-major format
        write(s, default_swap(round(Tout, mx*gray(mapf(img[i,j])))))
    end
    nothing
end

@noinline function save(s::Stream{format"PPMBinary"}, img::AbstractMatrix{T}, mapf) where {T<:Color}
    axh, axw = axes(img)
    Tout, mx = pnmmax(img)
    if sizeof(Tout) > 2
        error("element type $Tout (from $T) not supported")
    end
    write(s, "$(length(axw)) $(length(axh))\n$mx\n")
    for i = axh, j = axw  # io is stored row-major, color-first
        c = RGB(mapf(img[i,j]))
        write(s, default_swap(round(Tout, mx*red(c))))
        write(s, default_swap(round(Tout, mx*green(c))))
        write(s, default_swap(round(Tout, mx*blue(c))))
    end
    nothing
end

@noinline function save(s::Stream{format"PBMText"}, img::AbstractMatrix{T}, mapf) where {T<:Number}
    axh, axw = axes(img)
    w = length(axw)
    write(s, "$(w) $(length(axh))\n")
    linelen = offset = 0
    maxdigits = 1
    for i = axh, j = axw  # s is stored in row-major format
        if linelen + maxdigits + 1 > 70
            write(s, "\n")
            linelen = 0
        elseif linelen > 0
            write(s, " ")
            linelen += 1
        end
        write(s, Bool(mapf(round(img[i,j]))) ? "1" : "0")
        linelen += maxdigits
        offset += 1
        if offset == w
            write(s, "\n")
            linelen = 0
        end
    end
    nothing
end

@noinline function save(s::Stream{format"PGMText"}, img::AbstractMatrix{T}, mapf) where {T<:Union{Gray,Number}}
    axh, axw = axes(img)
    Tout, mx = pnmmax(img)
    if sizeof(Tout) > 2
        error("element type $Tout (from $T) not supported")
    end
    w = length(axw)
    write(s, "$(w) $(length(axh))\n$mx\n")
    linelen = offset = 0
    maxdigits = round(Int, log10(mx + 1), RoundUp)
    for i = axh, j = axw  # s is stored in row-major format
        value = round(Tout, mx*gray(mapf(img[i,j])))
        if linelen + maxdigits + 1 > 70
            write(s, "\n")
            linelen = 0
        elseif linelen > 0
            write(s, " ")
            linelen += 1
        end
        write(s, lpad("$(Int(value))", maxdigits))
        linelen += maxdigits
        offset += 1
        if offset == w
            write(s, "\n")
            linelen = 0
        end
    end
    nothing
end

@noinline function save(s::Stream{format"PPMText"}, img::AbstractMatrix{T}, mapf) where {T<:Color}
    axh, axw = axes(img)
    Tout, mx = pnmmax(img)
    if sizeof(Tout) > 2
        error("element type $Tout (from $T) not supported")
    end
    w = length(axw)
    write(s, "$(w) $(length(axh))\n$mx\n")
    linelen = offset = 0
    maxdigits = round(Int, log10(mx + 1), RoundUp)
    for i = axh, j = axw  # io is stored row-major, color-first
        c = RGB(mapf(img[i,j]))
        r,g,b = round.(Tout, mx.*(red(c),green(c),blue(c)))
        if linelen + 3*maxdigits + 1 > 70
            write(s, "\n")
            linelen = 0
        elseif linelen > 0
            write(s, "  ")
            linelen += 1
        end
        write(s, string(lpad("$(Int(r))", maxdigits), ' ', lpad("$(Int(g))", maxdigits), ' ', lpad("$(Int(b))", maxdigits)))
        linelen += 2+3*maxdigits
        offset += 1
        if offset == w
            write(s, "\n")
            linelen = 0
        end
    end
    nothing
end

function parse_netpbm_size(stream::IO)
    (parsenextint(stream), parsenextint(stream))
end

function parse_netpbm_maxval(stream::IO)
    parsenextint(stream)
end

function parsenextint(stream::IO)
    # ikirill: ugly, but I can't figure out a better way
    skipchars(isspace, stream, linecomment='#')
    from = position(stream)
    mark(stream)
    skipchars(isdigit, stream)
    to = position(stream)
    reset(stream)
    parse(Int, String(read(stream, to-from+1)))
end

function pnmmax(img::AbstractArray{T}) where {T}
    if isconcretetype(T)
        return pnmmax(eltype(T))
    end
    # Determine the concrete type that can hold all the elements
    S = typeof(first(img))
    for val in img
        S = promote_type(S, typeof(val))
    end
    pnmmax(eltype(S))
end

pnmmax(::Type{T}) where {T<:AbstractFloat} = UInt8, 255
function pnmmax(::Type{U}) where {U<:Normed}
    FixedPointNumbers.rawtype(U), reinterpret(one(U))
end
pnmmax(::Type{T}) where {T<:Unsigned} = T, typemax(T)

mybswap(i::Integer)  = bswap(i)
mybswap(i::Normed)   = bswap(i)
mybswap(c::Colorant) = mapc(bswap, c)
mybswap(c::RGB24) = c

const default_swap = is_little_endian ? mybswap : identity

function kwrename(newname, newval, oldname, oldval, caller::Symbol)
    if oldval !== nothing
        Base.depwarn("keyword $oldname has been renamed $newname", caller)
        return oldval
    end
    newval
end

_read_metadata(img::AbstractArray) = nothing

function _read_metadata(img::ImageMeta)
    isempty(properties(img)) && return nothing
    buf = IOBuffer()
    for (i, (k, v)) in enumerate(properties(img))
        i > 1 && write(buf, '\n')
        print(buf, "# ", k, ": ", v)
    end
    return String(take!(buf))
end

end # module
