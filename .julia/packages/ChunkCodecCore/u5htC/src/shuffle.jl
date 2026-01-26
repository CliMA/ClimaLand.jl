# The Shuffle codec

"""
    struct ShuffleCodec <: Codec
    ShuffleCodec(element_size::Integer)

Byte shuffle. The element size is required
to be able to decode the shuffle.

For example:

```julia
[
0x11, 0x12, 0x13,
0x21, 0x22, 0x23,
0x31, 0x32, 0x33,
0x41, 0x42, 0x43,
]
```

with element size 3 will be encoded as

```julia
[
0x11, 0x21, 0x31, 0x41,
0x12, 0x22, 0x32, 0x42,
0x13, 0x23, 0x33, 0x43,
]
```

If the length of the data is not evenly divisible
by the element size, the remainder data is appended.

For example "12312312312345" with element size 3 will be encoded as
"11112222333345"

A `ShuffleCodec` can be used as an encoder or decoder.
"""
struct ShuffleCodec <: Codec
    element_size::Int64
    function ShuffleCodec(element_size::Integer)
        check_in_range(Int64(1):typemax(Int64); element_size)
        new(Int64(element_size))
    end
end

decode_options(x::ShuffleCodec) = ShuffleDecodeOptions(;codec=x) # default decode options

# Allow ShuffleCodec to be used as an encoder
decoded_size_range(::ShuffleCodec) = Int64(0):Int64(1):typemax(Int64)-Int64(1)

encode_bound(::ShuffleCodec, src_size::Int64)::Int64 = src_size

function try_encode!(e::ShuffleCodec, dst::AbstractVector{UInt8}, src::AbstractVector{UInt8}; kwargs...)::MaybeSize
    dst_size::Int64 = length(dst)
    src_size::Int64 = length(src)
    element_size = e.element_size
    check_in_range(decoded_size_range(e); src_size)
    if dst_size < src_size
        NOT_SIZE
    else
        if src_size>>1 < element_size || element_size == 1
            copyto!(dst, src)
            return src_size
        end
        n_elements, n_remainder = fldmod(src_size, element_size)
        @inbounds for i in 0:(element_size-1)
            for j in 0:(n_elements-1)
                dst[begin + j + n_elements*i] = src[begin + i + element_size*j]
            end
        end
        offset = n_elements*element_size
        for i in 0:(n_remainder-1)
            dst[begin + offset + i] = src[begin + offset + i]
        end
        return src_size
    end
end

"""
    struct ShuffleEncodeOptions <: EncodeOptions
    ShuffleEncodeOptions(; kwargs...)

Byte shuffle encoding.

# Keyword Arguments

- `codec::ShuffleCodec`
"""
struct ShuffleEncodeOptions <: EncodeOptions
    codec::ShuffleCodec
end
function ShuffleEncodeOptions(;
        codec::ShuffleCodec,
        kwargs...
    )
    ShuffleEncodeOptions(codec)
end

is_thread_safe(::ShuffleEncodeOptions) = true

decoded_size_range(x::ShuffleEncodeOptions) = decoded_size_range(x.codec)

encode_bound(x::ShuffleEncodeOptions, src_size::Int64)::Int64 = encode_bound(x.codec, src_size)

function try_encode!(x::ShuffleEncodeOptions, dst::AbstractVector{UInt8}, src::AbstractVector{UInt8}; kwargs...)::MaybeSize
    try_encode!(x.codec, dst, src)
end

"""
    struct ShuffleDecodeOptions <: DecodeOptions
    ShuffleDecodeOptions(; kwargs...)

Byte shuffle decoding.

# Keyword Arguments

- `codec::ShuffleCodec`
"""
struct ShuffleDecodeOptions <: DecodeOptions
    codec::ShuffleCodec
end
function ShuffleDecodeOptions(;
        codec::ShuffleCodec,
        kwargs...
    )
    ShuffleDecodeOptions(codec)
end

is_thread_safe(::ShuffleDecodeOptions) = true

function try_find_decoded_size(::ShuffleDecodeOptions, src::AbstractVector{UInt8})::Int64
    length(src)
end

function try_decode!(d::ShuffleDecodeOptions, dst::AbstractVector{UInt8}, src::AbstractVector{UInt8}; kwargs...)::MaybeSize
    dst_size::Int64 = length(dst)
    src_size::Int64 = length(src)
    element_size = d.codec.element_size
    if dst_size < src_size
        NOT_SIZE
    else
        if src_size>>1 < element_size || element_size == 1
            copyto!(dst, src)
            return MaybeSize(src_size)
        end
        n_elements, n_remainder = fldmod(src_size, element_size)
        @inbounds for i in 0:(element_size-1)
            for j in 0:(n_elements-1)
                dst[begin + i + element_size*j] = src[begin + j + n_elements*i]
            end
        end
        offset = n_elements*element_size
        for i in 0:(n_remainder-1)
            dst[begin + offset + i] = src[begin + offset + i]
        end
        return MaybeSize(src_size)
    end
end
