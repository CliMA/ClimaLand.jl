# Provide the interface as a module so that different sub-modules can reuse the symbols for dispatching
module SixelInterface

export AbstractSixelDecoder,
       AbstractSixelEncoder,
       canonical_sixel_eltype

abstract type AbstractSixelEncoder end
(enc::AbstractSixelEncoder)(io, src) = error("The encoder functor method for inputs (::$(typeof(io)), ::$(typeof(src)) is not implemented.")

abstract type AbstractSixelDecoder end
(enc::AbstractSixelDecoder)(io, src) = error("The decoder functor method for inputs (::$(typeof(io)), ::$(typeof(src)) is not implemented.")

"""
    canonical_sixel_eltype(enc, CT1) -> CT2

Given input type `CT1`, infer the expected colorant type `CT2` for encoder `enc`. This function is
used to early convert the image data into formats that backend encoder (e.g., the C libsixel) understands.
"""
function canonical_sixel_eltype(::AbstractSixelEncoder, ::DataType) end

end # module
