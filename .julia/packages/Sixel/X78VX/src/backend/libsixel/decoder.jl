struct LibSixelDecoder <: AbstractSixelDecoder
    # (Experimental) internal fields

    # We need `allocator` to constructor libsixel objects, however, users of this
    # Julia package is not expected to use this field as it really should just live
    # in the C world.
    allocator::SixelAllocator

    function LibSixelDecoder(allocator=SixelAllocator())
        new(allocator)
    end
end

function (dec::LibSixelDecoder)(bytes::Vector{UInt8}; transpose=false)
    pixels = Ref{Ptr{Cuchar}}()
    palette = Ref{Ptr{Cuchar}}()
    pwidth = Ref(Cint(0))
    pheight = Ref(Cint(0))
    ncolors = Ref(Cint(0))

    status = C.sixel_decode_raw(
        bytes, length(bytes),
        pixels, pwidth, pheight,
        palette, ncolors,
        dec.allocator.ptr
    )
    check_status(status)

    palette_length = 256
    ncolors[] <= palette_length || throw(ArgumentError("palette `ncolors` should be less than $(palette_length), instead it is $(ncolors[])"))

    index = unsafe_wrap(Matrix{Cuchar}, pixels[], (pwidth[], pheight[]); own=false)

    # libsixel declares it to be ARGB but it's actually RGB{N0f8}
    PT = RGB{N0f8}
    pvalues = convert(Ptr{PT}, palette[])
    # TODO: we can acutally further compress the data with `values` of length `ncolors`
    #       Does it worth it? Maybe, but not very much, though.
    values = unsafe_wrap(Vector{PT}, pvalues, (palette_length, ); own=false)
    values = OffsetArray(values, OffsetArrays.Origin(0)) # again, libsixel assumes 0-based indexing

    return index, values
end
