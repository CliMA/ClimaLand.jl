"""
    abstract type DecodingError <: Exception

Generic error for data that cannot be decoded.
"""
abstract type DecodingError <: Exception end

"""
    struct MaybeSize
        val::Int64
    end

If `val ≥ 0` it is a size, and can be converted back and forth with `Int64`.
If `val < 0` converting to and from `Int64` will error.
If `val == typemin(Int64)` it is an unknown size.
Otherwise it is a size hint of `-val`.

"""
struct MaybeSize
    val::Int64
end

"""
    const NOT_SIZE = MaybeSize(typemin(Int64))
"""
const NOT_SIZE = MaybeSize(typemin(Int64))
function is_size(x::MaybeSize)::Bool
    !signbit(x.val)
end
function Base.Int64(x::MaybeSize)
    if !is_size(x)
        throw(InexactError(:Int64, Int64, x))
    else
        x.val
    end
end
function Base.convert(::Type{Int64}, x::MaybeSize)
    Int64(x)
end
function Base.convert(::Type{MaybeSize}, x::Int64)::MaybeSize
    if signbit(x)
        throw(InexactError(:convert, MaybeSize, x))
    else
        MaybeSize(x)
    end
end

"""
    struct DecodedSizeError <: Exception
    DecodedSizeError(max_size, decoded_size)

Exception thrown when the decoded data size doesn't match expectations or exceeds limits.

# Fields
- `max_size::Int64`: The maximum allowed or expected size in bytes
- `decoded_size::MaybeSize`: The actual decoded size, size hint, or `NOT_SIZE` if unknown

This error can occur in several scenarios:
1. Decoded size exceeds the maximum allowed size
2. Decoded size is less than expected
3. Decoder provides a size hint when the decoded size exceeds limits
4. Decoded size is completely unknown but exceeds limits
"""
struct DecodedSizeError <: Exception
    max_size::Int64
    decoded_size::MaybeSize
end

function Base.showerror(io::IO, err::DecodedSizeError)
    print(io, "DecodedSizeError: ")
    if err.decoded_size === NOT_SIZE
        print(io, "decoded size > ")
        print(io, err.max_size)
    elseif !is_size(err.decoded_size)
        suggested_size = -err.decoded_size.val
        print(io, "decoded size > ")
        print(io, err.max_size)
        print(io, ", try max_size = ")
        print(io, suggested_size)
    else
        decoded_size = err.decoded_size.val
        if decoded_size < err.max_size
            print(io, "decoded size ")
            print(io, decoded_size)
            print(io, " < expected ")
            print(io, err.max_size)
        else
            print(io, "decoded size ")
            print(io, decoded_size)
            print(io, " > ")
            print(io, err.max_size)
        end
    end
    nothing
end
