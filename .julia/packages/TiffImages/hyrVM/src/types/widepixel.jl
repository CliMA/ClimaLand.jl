"""
    $(TYPEDEF)

Pixel type used for images with extra (non-alpha) channels

See also [`nchannels`](@ref), [`channel`](@ref), [`color`](@ref)
"""
struct WidePixel{C <: Colorant, X <: Tuple}
    color::C
    extra::X
end

Base.zero(::WidePixel{C, X}) where {C, X} = WidePixel(zero(C), zero.(fieldtypes(X)))

function Base.show(io::IO, x::WidePixel{C,X}) where {C, X}
    if get(io, :compact, false)
        print(io, "WidePixel($(x.color))")
    else
        len = length(x.extra)
        print(io, "WidePixel($(x.color)) + $len extra channel$(len > 1 ? "s" : "")")
    end
end