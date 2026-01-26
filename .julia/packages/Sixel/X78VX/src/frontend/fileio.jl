using FileIO

function fileio_load(f::File{format"SIXEL"}; kwargs...)
    open(f, "r") do io
        sixel_decode(io)
    end
end

function fileio_load(io::Stream{format"SIXEL"}; kwargs...)
    sixel_decode(RGB{N0f8}, io; kwargs...)
end

sixel_decode(::Type{T}, io::Stream{format"SIXEL"}, dec::SDC=default_decoder(); kwargs...) where T =
    sixel_decode(T, io.io, dec; kwargs...)

function fileio_save(f::File{format"SIXEL"}, img::AbstractArray; kwargs...)
    open(f, "w") do io
        sixel_encode(io, img; kwargs...)
    end
end

function fileio_save(io::Stream{format"SIXEL"}, img::AbstractArray; kwargs...)
    sixel_encode(io, img; kwargs...)
end

sixel_encode(io::Stream{format"SIXEL"}, img::AbstractArray, enc::SEC=default_encoder(img); kwargs...) =
    sixel_encode(io.io, img, enc; kwargs...)
