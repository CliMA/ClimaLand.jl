"""
    sixel_decode([T=RGB{N0f8}], src, [decoder]; kwargs...) -> img::IndirectArray

Decode the sixel format sequence provided by `src` and output as an indexed image.

# Arguments

- `T`: output eltype. By default it is `RGB{N0f8}`.
- `src`: the input sixel data source. It can be either an `IO`, `AbstractVector{UInt8}`, or `AbstractString`.
- `decoder::AbstractSixelDecoder`: the sixel decoder. Currently only `LibSixelDecoder` is available.

!!! warning
    `sixel_decode` for tiny images (including vectors) is undefined behavior; sixel format requires
    at least `6` pixels in the row direction. For example, you should not expect it returning a
    `Vector` data even if it is.

# Parameters

- `transpose::Bool`: whether we need to permute the image's width and height dimension before encoding.
  The default value is `false`.

# References

- [1] VT330/VT340 Programmer Reference Manual, Volume 1: Text Programming
- [2] VT330/VT340 Programmer Reference Manual, Volume 2: Graphics Programming
- [3] https://github.com/saitoha/libsixel
"""
function sixel_decode end

sixel_decode(src, dec=default_decoder(); kwargs...) = sixel_decode(RGB{N0f8}, src, dec; kwargs...)

sixel_decode(::Type{T}, data::AbstractString, dec=default_decoder(); kwargs...) where T =
    sixel_decode(T, collect(UInt8, convert(String, data)), dec; kwargs...)

function sixel_decode(::Type{T}, bytes::AbstractArray, dec=default_decoder(); transpose=false) where {T}
    bytes = convert(Vector{UInt8}, bytes)

    expected_size = read_sixel_size(bytes)
    index, values = dec(bytes)
    values = eltype(values) == T ? values : map(T, values)

    if dec isa LibSixelDecoder
        # Julia uses column-major order while libsixel uses row-major order,
        # thus transpose=true means no permutation.
        index = transpose ? index : PermutedDimsArray(index, (2, 1))
        expected_size = transpose ? (expected_size[2], expected_size[1]) : expected_size
    else
        throw(ArgumentError("Unsupported decoder type $(typeof(enc)). Please open an issue for this."))
    end

    check_size(size(index), expected_size)
    # We use IndirectArray to mark it as an indexed image so as to avoid unnecessary memory allocation
    # Users that expect a dense Array can always call `collect(rst)` or `convert(Array, rst)` on this.
    return IndirectArray(index, values)
end

function sixel_decode(::Type{T}, io::IO, dec=default_decoder(); kwargs...) where T
    # TODO: This is actually a duplicated check since the bytes method also checks it
    #       but this is a quite fast operation...
    expected_size = read_sixel_size(io)

    bytes = read(io)
    img = sixel_decode(T, bytes, dec; kwargs...)

    check_size(size(img), expected_size)
    return img
end


"""
    read_sixel_size(io::IO)
    read_sixel_size(bytes::Vector{UInt8})

Read the header from a sixel sequence `io`/`bytes` and return the size of the sixel image.
"""
function read_sixel_size(bytes::Vector{UInt8})
    # There's absolutely something wrong if the header content exceeds 50 bytes
    max_header_length = 50

    p_end = findfirst(isequal(UInt8('#')), bytes)
    p_end = isnothing(p_end) ? length(bytes) : p_end - 1
    buffer = view(bytes, 1:min(p_end, max_header_length))
    seps = findall(isequal(UInt8(';')), buffer)
    length(seps) == 3 || throw(ArgumentError("The input data is not recognizable as sixel sequence."))
    w = parse(Int, String(buffer[seps[2]+1:seps[3]-1]))
    h = parse(Int, String(buffer[seps[3]+1:end]))
    # (h, w) in column-major order convention
    return h, w
end

function read_sixel_size(io::IO)
    # Sixel sequence format always start with this
    #     \ePq"1;1;w;h#
    # where `#` indicates the first palette value
    buffer = UInt8[]
    p = position(io)
    try
        # There's absolutely something wrong if the header content exceeds 50 bytes
        i, max_header_length = 1, 50
        ch = read(io, Cuchar)
        # Cuchar('#') == 0x23
        while ch != 0x23 && i < max_header_length
            push!(buffer, ch)
            ch = read(io, Cuchar)
            i += 1
        end
    catch e
        e isa EOFError && throw(ArgumentError("The input data is not recognizable as sixel sequence."))
        rethrow(e)
    finally
        seek(io, p)
    end
    return read_sixel_size(buffer)
end

function check_size(actual_size, expected_size)
    valid = length(actual_size)==length(expected_size)
    valid = valid && all(zip(actual_size, expected_size)) do x
        # libsixel encode/decode for small images is undefined behavior to us
        minimum(x) <= 6 && return true
        x[1] == x[2]
    end
    if !valid
        @warn "Output size mismatch during decoding sixel sequences" actual_size expected_size
    end
end
