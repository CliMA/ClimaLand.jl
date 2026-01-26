"""
    abstract type Codec

Required information to decode encoded data.

Properties are public for reading.

Required methods for a type `T <: Codec` to implement:
- `decode_options(::T)::DecodeOptions`
"""
abstract type Codec end

"""
    abstract type EncodeOptions

Options for encoding data.

Properties are public for reading.
All `EncodeOptions` have a keyword argument constructor
that accept all properties as arguments.
All `EncodeOptions` have a `codec::Codec` property.

Required methods for a type `T <: EncodeOptions` to implement:
- `decoded_size_range(::T)::StepRange{Int64, Int64}`
- `encode_bound(::T, src_size::Int64)::Int64`
- `try_encode!(::T, dst::AbstractVector{UInt8}, src::AbstractVector{UInt8}; kwargs...)::MaybeSize`

Optional methods to implement:
- `is_thread_safe(::T)::Bool`: defaults to `false`.
- `is_lossless(::T)::Bool`: defaults to `true`.
"""
abstract type EncodeOptions end


"""
    abstract type DecodeOptions

Options for decoding data.

Properties are public for reading.
All `DecodeOptions` have a keyword argument constructor
that accept all properties as arguments.
All `DecodeOptions` have a `codec::Codec` property.

Required methods for a type `T <: DecodeOptions` to implement:
- `try_find_decoded_size(::T, src::AbstractVector{UInt8})::Union{Nothing, Int64}`
- `try_decode!(::T, dst::AbstractVector{UInt8}, src::AbstractVector{UInt8}; kwargs...)::MaybeSize`

Optional methods to implement:
- `is_thread_safe(::T)::Bool`: defaults to `false`.
- `can_concatenate(::T)::Bool`: defaults to `false`.
- `try_resize_decode!(::T, dst::AbstractVector{UInt8}, src::AbstractVector{UInt8}, max_size::Int64; kwargs...)::MaybeSize`: defaults to using `try_decode!` and `try_find_decoded_size`
"""
abstract type DecodeOptions end