"""
    $(TYPEDEF)

In-memory representation of Tiff Tags, which are essentially key value pairs.
The `data` field can either be a `String`, a `Number`, an Array of bitstypes, or
a [`RemoteData`](@ref) type.

$FIELDS
"""
struct Tag{T}
    tag::UInt16
    data::T
end

Tag(tag::TiffTag, data) = Tag(UInt16(tag), data)
Tag(tag::UInt16, data::String) = Tag{String}(tag, data)
Tag(tag::UInt16, data::T) where {T <: Enum} = Tag{UInt16}(tag, UInt16(data))
function Tag(tag::UInt16, data::AbstractVector{T}) where {T}
    if length(data) == 1
        Tag{T}(tag, first(data))
    else
        Tag{Vector{T}}(tag, data)
    end
end

function Base.getproperty(t::Tag{T}, f::Symbol) where {T}
    if f === :data
        isa(t, Tag{Vector{UInt32}}) ? getfield(t, :data)::Vector{UInt32} :
        isa(t, Tag{Vector{UInt16}}) ? getfield(t, :data)::Vector{UInt16} :
        isa(t, Tag{Vector{UInt8}})  ? getfield(t, :data)::Vector{UInt8} :
        isa(t, Tag{UInt32})         ? getfield(t, :data)::UInt32 :
        isa(t, Tag{UInt16})         ? getfield(t, :data)::UInt16 :
        isa(t, Tag{UInt8})          ? getfield(t, :data)::UInt8 :
        isa(t, Tag{Int64})          ? getfield(t, :data)::Int64 :
        getfield(t, :data)
    else
        getfield(t, :tag)::UInt16
    end
end

function Base.getproperty(t::Tag{<: AbstractString}, f::Symbol)
    if f === :data
        # strip NUL terminators from the end of strings
        strip(getfield(t, :data)::String, '\0')
    else
        getfield(t, :tag)::UInt16
    end
end
Base.getproperty(t::Tag{<: RemoteData}, f::Symbol) = (f == :data) ? error("Data hasn't been loaded, use `load!` on this IFD") : getfield(t, f)

isloaded(t::Tag) = true
isloaded(t::Tag{<: RemoteData}) = false

Base.length(t::Tag{<: AbstractVector}) = length(t.data)
Base.length(t::Tag{<: AbstractString}) = (endswith(t.data, '\0') ? ncodeunits(t.data) : ncodeunits(t.data) + 1)
Base.length(t::Tag{<: RemoteData}) = getfield(t, :data).count
Base.length(t::Tag) = 1

Base.eltype(::Tag{T}) where {T} = T
Base.eltype(::Tag{<: AbstractVector{T}}) where {T} = T
Base.eltype(t::Tag{RemoteData{O, T}}) where {O, T} = T

"""
    sizeof(tag::TiffImages.Tag)

Minimum number of bytes that the _data_ in `tag` will use on disk.

!!! note 
    Actual space on disk will be different because the tag's representation depends on
    the file's offset. For example, given a 2 bytes of data in `tag` and a file with
    `UInt32` offsets, the actual usage on disk will be `sizeof(UInt32)=4` for the
    data + tag overhead
"""
Base.sizeof(t::Tag{T}) where {T} = bytes(T) * length(t)

load(::TiffFile, t::Tag) = t

function load(tf::TiffFile{O}, t::Tag{RemoteData{O, T}}) where {O <: Unsigned, T}
    t_data = getfield(t, :data)
    N = t_data.count
    nb = bytes(T)::Int
    rawdata = Vector{UInt8}(undef, nb*N)

    pos = position(tf.io)
    seek(tf, t_data.position)
    read!(tf, rawdata)

    # if this datatype is comprised of multiple bytes and this file needs to be
    # bitswapped then we'll need to reverse the byte order inside each datatype
    # unit
    if tf.need_bswap && nb >= 2
        reverse!(rawdata)
        data = Array{T}(reverse(reinterpret(T, rawdata)))
    elseif T == String
        data = String(rawdata)
    elseif T == Any
        data = Array{Any}(rawdata)
    else
        data = Array{T}(reinterpret(T, rawdata))
    end

    if N == 1
        data = first(data)
    end
    seek(tf, pos)

    Tag(t.tag, data)::Tag
end

bytes(x::Type) = sizeof(x)
bytes(::Type{Any}) = 1
bytes(::Type{String}) = 1
bytes(::Type{RemoteData{O, T}}) where {O, T} = bytes(T)
bytes(::Type{<: AbstractVector{T}}) where {T} = bytes(T)

function Base.read(tf::TiffFile{O}, ::Type{Tag}) where O <: Unsigned
    tag = read(tf, UInt16)
    datatype = read(tf, UInt16)
    count = read(tf, O)
    data = Vector{UInt8}(undef, sizeof(O))
    read!(tf, data)

    T = Any
    if datatype in keys(tiff_to_julian)
        T = tiff_to_julian[datatype]
    end

    nbytes = bytes(T)::Int * count
    if nbytes <= sizeof(O)
        if tf.need_bswap
            reverse!(view(data, 1:nbytes))
        end
        if T === String
            return Tag(tag, String(data))
        elseif T === Any
            return Tag(tag, Vector{Any}(data))
        elseif count == 1
            return Tag(tag, first(reinterpret(T, data)))::Tag
        else
            return Tag(tag, reinterpret(T, data)[1:Int(count)])
        end
    else
        (tf.need_bswap) && reverse!(data)
        return Tag(tag, RemoteData(first(reinterpret(O, data)), T, count))
    end
end

_showdata(io::IO, t::Tag{RemoteData{O, T}}) where {O, T} = print(io, "REMOTE@$(getfield(t, :data).position) $(T)[] len=$(getfield(t, :data).count)")
_showdata(io::IO, t::Tag{<: AbstractString}) = print(io, "\"", first(t.data, 20), (length(t.data) > 20) ? "..." : "", "\"")
_showdata(io::IO, t::Tag{<: AbstractVector}) = print(io, "$(eltype(t.data))[", join(t.data[1:min(5, end)], ", "), (length(t.data) > 5) ? ", ..." : "", "]")
_showdata(io::IO, t::Tag) = print(io, t.data)

function Base.show(io::IO, t::Tag)
    print(io, "Tag(")
    try
        print(io, TiffTag(t.tag), ", ")
    catch
        print(io, "UNKNOWN($(Int(t.tag))), ")
    end
    if t.tag == Int(COMPRESSION)
        print(io, CompressionType(t.data))
    else
        _showdata(io, t)
    end
    print(io, ")")
end

"""
    write(tf, t)

Write tag `t` to the tiff file `tf`. Returns `true` if the tag data fit
entirely in the IFD space and was written to disk. Otherwise it returns false.
"""
function Base.write(tf::TiffFile{O}, t::Tag{T}) where {O <: Unsigned, T}
    # if the data are too large to fit then we'll need to skip writing this tag
    # for now until we know the length of the entire IFD
    if sizeof(t) > sizeof(O)
        _writeblank(tf)
        return false
    end

    data = t.data
    # add NUL terminator to the end of Strings that don't have it already
    if T === String
        data = data::SubString{String}
        if !endswith(data, '\0')
            data *= "\0"
        end
    end

    write(tf, t.tag)
    write(tf, julian_to_tiff[eltype(t)])
    write(tf, O(length(t)))
    nbytes = write(tf.io, data)::Int

    # write padding
    if nbytes < sizeof(O)
        write(tf, fill(0x00, sizeof(O) - nbytes))
    end
    true
end

function Base.write(tf::TiffFile{O}, t::Tag{RemoteData{O, T}}) where {O <: Unsigned, T}
    write(tf, t.tag)
    write(tf, julian_to_tiff[T])
    t_data = getfield(t, :data)
    write(tf, t_data.count)
    write(tf, t_data.position)
    true
end


function Base.write(tf::TiffFile{O}, t::Tag, offset) where {O <: Unsigned}
    write(tf, Tag(t.tag, RemoteData(O(offset), eltype(t), O(length(t)))))
end

# Base.write(tf::TiffFile, t::Tag) = error("Tag offsets must agree with file offsets")

function _writeblank(tf::TiffFile{O}) where O
    write(tf, UInt32(0))
    write(tf, zero(O))
    write(tf, zero(O))
end

const tagfields = fieldnames(Tag)
function Base.:(==)(t1::Tag{O1}, t2::Tag{O2}) where {O1, O2}
    return O1 == O2 && getproperty.(Ref(t1), tagfields) == getproperty.(Ref(t2), tagfields)
end
