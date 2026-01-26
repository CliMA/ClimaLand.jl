# The Noop codec

"""
    struct NoopCodec <: Codec
    NoopCodec()

Copies the input.

See also [`NoopEncodeOptions`](@ref) and [`NoopDecodeOptions`](@ref)
"""
struct NoopCodec <: Codec end
decode_options(::NoopCodec) = NoopDecodeOptions() # default decode options

"""
    struct NoopEncodeOptions <: EncodeOptions
    NoopEncodeOptions(; kwargs...)

Copies the input.

# Keyword Arguments

- `codec::NoopCodec=NoopCodec()`
"""
struct NoopEncodeOptions <: EncodeOptions
    codec::NoopCodec
end
function NoopEncodeOptions(;
        codec::NoopCodec=NoopCodec(),
        kwargs...
    )
    NoopEncodeOptions(codec)
end

is_thread_safe(::NoopEncodeOptions) = true

decoded_size_range(::NoopEncodeOptions) = Int64(0):Int64(1):typemax(Int64)-Int64(1)

encode_bound(::NoopEncodeOptions, src_size::Int64)::Int64 = src_size

function try_encode!(e::NoopEncodeOptions, dst::AbstractVector{UInt8}, src::AbstractVector{UInt8}; kwargs...)::MaybeSize
    dst_size::Int64 = length(dst)
    src_size::Int64 = length(src)
    check_in_range(decoded_size_range(e); src_size)
    if dst_size < src_size
        NOT_SIZE
    else
        copyto!(dst, src)
        src_size
    end
end

"""
    struct NoopDecodeOptions <: DecodeOptions
    NoopDecodeOptions(; kwargs...)

Copies the input.

# Keyword Arguments

- `codec::NoopCodec=NoopCodec()`
"""
struct NoopDecodeOptions <: DecodeOptions
    codec::NoopCodec
end
function NoopDecodeOptions(;
        codec::NoopCodec=NoopCodec(),
        kwargs...
    )
    NoopDecodeOptions(codec)
end

is_thread_safe(::NoopDecodeOptions) = true

can_concatenate(::NoopDecodeOptions) = true

function try_find_decoded_size(::NoopDecodeOptions, src::AbstractVector{UInt8})::Int64
    length(src)
end

function try_decode!(::NoopDecodeOptions, dst::AbstractVector{UInt8}, src::AbstractVector{UInt8}; kwargs...)::MaybeSize
    dst_size::Int64 = length(dst)
    src_size::Int64 = length(src)
    if dst_size < src_size
        NOT_SIZE
    else
        copyto!(dst, src)
        MaybeSize(src_size)
    end
end
