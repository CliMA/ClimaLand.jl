using FileIO

function fileio_load(f::File{format"JPEG"}; kwargs...)
    open(f.filename, "r") do io
        jpeg_decode(io; kwargs...)
    end
end
fileio_load(io::Stream{format"JPEG"}; kwargs...) = jpeg_decode(read(io); kwargs...)

function fileio_save(f::File{format"JPEG"}, img::AbstractArray; kwargs...)
    open(f.filename, "w") do io
        jpeg_encode(io, img; kwargs...)
    end
end
function fileio_save(io::Stream{format"JPEG"}, img::AbstractArray; kwargs...)
    jpeg_encode(io.io, img; kwargs...)
end
