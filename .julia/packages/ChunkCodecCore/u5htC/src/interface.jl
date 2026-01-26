"""
    encode(e, src::AbstractVector{UInt8})::Vector{UInt8}

Encode the input vector `src` using encoder `e`.
`e` must implement [`decoded_size_range`](@ref),
[`encode_bound`](@ref), and [`try_encode!`](@ref).

Throw an error if `length(src)` is not in `decoded_size_range(e)`.

Otherwise throw an error.

Precondition: `src` has a stride of 1.

See also [`EncodeOptions`](@ref) and [`decode`](@ref).
"""
function encode(e, src::AbstractVector{UInt8})::Vector{UInt8}
    src_size::Int64 = length(src)
    check_in_range(decoded_size_range(e)::StepRange{Int64, Int64}; src_size)
    dst_size_bound = encode_bound(e, src_size)::Int64
    if dst_size_bound == typemax(Int64)
        throw(ArgumentError("`encode_bound(e, $(src_size))` saturated"))
    end
    dst = Vector{UInt8}(undef, dst_size_bound)
    real_dst_size = Int64(try_encode!(e, dst, src)::MaybeSize)
    @assert real_dst_size ∈ 0:dst_size_bound
    if real_dst_size < dst_size_bound
        resize!(dst, real_dst_size)
    end
    dst
end

"""
    decode(d, src::AbstractVector{UInt8}; max_size::Integer=typemax(Int64), size_hint::Integer=Int64(0))::Vector{UInt8}

Decode the input data `src` using decoder `d`.
`d` must implement [`try_find_decoded_size`](@ref), [`try_decode!`](@ref), and optionally [`try_resize_decode!`](@ref).

Throw a [`DecodedSizeError`](@ref) if decoding fails because the output size would be greater than `max_size`.

Throw a [`DecodingError`](@ref) if decoding fails because the input data is not valid.

Otherwise throw an error.

If you have a good idea of what the decoded size is, using the `size_hint` keyword argument
can greatly improve performance.

Precondition: `src` has a stride of 1.

See also [`DecodeOptions`](@ref) and [`encode`](@ref)
"""
function decode(
        d,
        src::AbstractVector{UInt8};
        max_size::Integer=typemax(Int64),
        size_hint::Integer=Int64(0),
    )::Vector{UInt8}
    _clamp_max_size::Int64 = clamp(max_size, Int64)
    if _clamp_max_size < Int64(0)
        throw(DecodedSizeError(_clamp_max_size, NOT_SIZE))
    end
    _clamp_size_hint::Int64 = clamp(size_hint, Int64(0), _clamp_max_size)
    dst = Vector{UInt8}(undef, _clamp_size_hint)
    maybe_dst_size = try_resize_decode!(d, dst, src, _clamp_max_size)::MaybeSize
    if !is_size(maybe_dst_size)
        throw(DecodedSizeError(_clamp_max_size, maybe_dst_size))
    end
    real_dst_size = Int64(maybe_dst_size)
    @assert real_dst_size ∈ 0:_clamp_max_size
    @assert real_dst_size ≤ length(dst)
    resize!(dst, real_dst_size)
    dst
end

"""
    decode!(d, dst::AbstractVector{UInt8}, src::AbstractVector{UInt8}) -> dst

Decode the input data `src` into `dst` using decoder `d`.
`d` must implement [`try_decode!`](@ref).

Throw a [`DecodedSizeError`](@ref) if the decoded output size is not exactly `length(dst)`.

Throw a [`DecodingError`](@ref) if decoding fails because the input data is not valid.

Otherwise throw an error.

Precondition: `dst` and `src` do not overlap in memory.

Precondition: `src` and `dst` have a stride of 1.

See also [`decode`](@ref) and [`encode`](@ref)
"""
function decode!(d, dst::AbstractVector{UInt8}, src::AbstractVector{UInt8})
    expected_size::Int64 = length(dst)
    n = try_decode!(d, dst, src)::MaybeSize
    if n.val != expected_size
        throw(DecodedSizeError(expected_size, n))
    else
        dst
    end
end

"""
    decode_options(::Codec)::DecodeOptions

Return the default decode options for the codec.
"""
function decode_options end

"""
    decoded_size_range(e)::StepRange{Int64, Int64}

Return the range of allowed `src` sizes for encoding.

See also [`encode_bound`](@ref)
"""
function decoded_size_range end

"""
    encode_bound(e, src_size::Int64)::Int64

Return the size of `dst` required to ensure [`try_encode!`](@ref)
succeeds regardless of `src`'s content if `src_size` is also in
[`decoded_size_range(e)`](@ref) and there is enough memory.

On the domain of `0:typemax(Int64)` this function must not error and must be monotonically increasing.
"""
function encode_bound end

"""
    try_encode!(e, dst::AbstractVector{UInt8}, src::AbstractVector{UInt8}; kwargs...)::MaybeSize

Try to encode `src` into `dst` using encoder `e`.

Return the size of the encoded output in `dst` if successful.

If `dst` is too small, return `NOT_SIZE`.

Throw an error if `length(src)` is not in `decoded_size_range(e)`

Otherwise throw an error.

Precondition: `dst` and `src` do not overlap in memory.

Precondition: `src` and `dst` have a stride of 1.

All of `dst` can be written to or used as scratch space by the encoder.
Only the initial returned number of bytes are valid output.

See also [`encode_bound`](@ref) and [`decoded_size_range`](@ref)
"""
function try_encode! end

"""
    is_thread_safe(::Union{Codec, DecodeOptions, EncodeOptions})::Bool

Return `true` if it is safe to use the options to encode or decode concurrently in multiple threads.
"""
is_thread_safe(::EncodeOptions) = false
is_thread_safe(::DecodeOptions) = false
is_thread_safe(c::Codec) = is_thread_safe(decode_options(c))

"""
    is_lossless(e)::Bool

Return `true` if the encoder is lossless.
"""
is_lossless(::Any) = true

"""
    try_find_decoded_size(d, src::AbstractVector{UInt8})::Union{Nothing, Int64}

Try to return the size of the decoded output of `src` using `d`.

If the size cannot be quickly determined, return `nothing`.

If the encoded data is found to be invalid, throw a `DecodingError`.

If an `Int64` is returned, it must be the exact size of the decoded output.
If [`try_decode!`](@ref) is called with a `dst` of this size, it must succeed and return the same size, or throw an error.
"""
function try_find_decoded_size end

"""
    try_decode!(d, dst::AbstractVector{UInt8}, src::AbstractVector{UInt8}; kwargs...)::MaybeSize

Try to decode `src` into `dst` using decoder `d`.

Return the size of the decoded output in `dst` if successful.

If `dst` is too small to fit the decoded output, return `NOT_SIZE`.
If `dst` is too small but a positive hint of the required size can be found, return `MaybeSize(-hint)`.

Throw a [`DecodingError`](@ref) if decoding fails because the input data is not valid.

Otherwise throw an error.

Precondition: `dst` and `src` do not overlap in memory.

Precondition: `src` and `dst` have a stride of 1.

All of `dst` can be written to or used as scratch space by the decoder.
Only the initial returned number of bytes are valid output.
"""
function try_decode! end

"""
    try_resize_decode!(d, dst::AbstractVector{UInt8}, src::AbstractVector{UInt8}, max_size::Int64; kwargs...)::MaybeSize

Try to decode the input `src` into `dst` using decoder `d`.

Return the size of the decoded output in `dst` if successful.

`dst` can be grown using the `resize!` function to any size between one more than the original length of `dst` and `max_size`.

Return `NOT_SIZE` if the size of `dst` is too small to contain the decoded output and cannot be grown due to the `max_size` restriction.
If a positive hint of the size can be found, return `MaybeSize(-hint)`.

Precondition: `dst` and `src` do not overlap in memory.

Precondition: `src` and `dst` have a stride of 1.

All of `dst` can be written to or used as scratch space by the decoder.
Only the initial returned number of bytes are valid output.
"""
function try_resize_decode!(d, dst::AbstractVector{UInt8}, src::AbstractVector{UInt8}, max_size::Int64; kwargs...)::MaybeSize
    init_dst_size::Int64 = length(dst)
    maybe_decoded_size = try_find_decoded_size(d, src)::Union{Nothing, Int64}
    if isnothing(maybe_decoded_size)
        while true
            ds = try_decode!(d, dst, src)::MaybeSize
            if !is_size(ds)
                if length(dst) ≥ max_size
                    return ds
                end
                local hint = -ds.val
                if hint ≤ length(dst)
                    grow_dst!(dst, max_size)
                else
                    resize!(dst, min(hint, max_size))
                end
            else
                @assert Int64(ds) ∈ 0:length(dst)
                return ds
            end
        end
    else
        decoded_size = Int64(maybe_decoded_size)
        @assert !signbit(decoded_size)
        @assert !signbit(init_dst_size)
        if decoded_size > init_dst_size
            if decoded_size > max_size
                return MaybeSize(-decoded_size)
            end
            resize!(dst, decoded_size)
        end
        real_dst_size = try_decode!(d, dst, src)::MaybeSize
        @assert Int64(real_dst_size) == decoded_size
        return real_dst_size
    end
end

"""
    can_concatenate(d)::Bool

Return `true` if the decoder has concatenation transparency.

If `true`, and some encoded data `a` and `b` decode to `x` and `y` respectively, then
the concatenation of `a` and `b` will
decode to the concatenation of `x` and `y`
"""
can_concatenate(::Any) = false

# allow passing codec to decode
try_find_decoded_size(c::Codec, src::AbstractVector{UInt8}) = try_find_decoded_size(decode_options(c), src)
try_decode!(c::Codec, dst::AbstractVector{UInt8}, src::AbstractVector{UInt8}; kwargs...) = try_decode!(decode_options(c), dst, src; kwargs...)
try_resize_decode!(c::Codec, dst::AbstractVector{UInt8}, src::AbstractVector{UInt8}, max_size::Int64; kwargs...) = try_resize_decode!(decode_options(c), dst, src, max_size; kwargs...)
can_concatenate(c::Codec) = can_concatenate(decode_options(c))

"""
    check_contiguous(x::AbstractVector{UInt8})

Check if the given vector is contiguous in memory.
Throw an error if the vector is not contiguous or if its length cannot be represented as Int64.
"""
function check_contiguous(x::AbstractVector{UInt8})
    y = Base.cconvert(Ptr{UInt8}, x)
    GC.@preserve y Base.unsafe_convert(Ptr{UInt8}, y)
    isone(only(strides(x))) || throw(ArgumentError("vector is not contiguous in memory"))
    Int64(length(x))
    @assert !signbit(length(x))
    nothing
end

"""
    check_in_range(range; kwargs...)

Check if all keyword arguments are within the specified range.
Throw an `ArgumentError` if any value is outside the range.

# Arguments
- `range`: The allowed range of values
- `kwargs...`: Keyword arguments to check against the range
"""
function check_in_range(range; kwargs...)
    for (k, v) in kwargs
        if v ∉ range
            throw(ArgumentError("$(k) ∈ $(range) must hold. Got\n$(k) => $(v)"))
        end
    end
end

"""
    grow_dst!(dst::AbstractVector{UInt8}, max_size::Int64)::Union{Nothing, Int64}

Grow the destination vector `dst` to a size between its current size and `max_size`.
Return the new size of `dst` if it was grown, or `nothing` if it could not be grown due to the `max_size` restriction.
"""
function grow_dst!(dst::AbstractVector{UInt8}, max_size::Int64)::Union{Nothing, Int64}
    cur_size::Int64 = length(dst)
    if cur_size ≥ max_size
        return
    end
    # This inequality prevents overflow
    next_size::Int64 = if max_size - cur_size ≤ cur_size
        max_size
    else
        max(2*cur_size, Int64(1))
    end
    resize!(dst, next_size)
    return next_size
end
