module ImageIO

using UUIDs
using FileIO: File, DataFormat, Stream, stream, Formatted

import IndirectArrays: IndirectArray

using LazyModules # @lazy macro is used to delay the package loading to its first usage

@lazy import Sixel = "45858cf5-a6b0-47a3-bbea-62219f50df47"
@lazy import Netpbm = "f09324ee-3d7c-5217-9330-fc30815ba969"
@lazy import PNGFiles = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
@lazy import TiffImages = "731e570b-9d59-4bfa-96dc-6df516fadf69"
@lazy import OpenEXR = "52e1d378-f018-4a11-a4be-720524705ac7"
@lazy import QOI = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
@lazy import JpegTurbo = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
@lazy import WebP = "e3aaa7dc-3e4b-44e0-be63-ffb868ccd7c1"

# Enforce a type conversion to be backend independent (issue #25)
# Note: If the backend does not provide efficient `convert` implementation,
#       there will be an extra memory allocation and thus hurt the performance.
for FMT in (
    :PBMBinary, :PGMBinary, :PPMBinary, :PBMText, :PGMText, :PPMText,
    :TIFF,
    :EXR,
    :QOI,
    :SIXEL,
    :JPEG,
    :WebP,
)
    @eval canonical_type(::DataFormat{$(Expr(:quote, FMT))}, ::AbstractArray{T, N}) where {T,N} =
        Array{T,N}
end
@inline canonical_type(::DataFormat{:PNG}, ::AbstractArray{T, N}) where {T,N} = Union{Array{T,N}, IndirectArray{T,N}}
@inline canonical_type(::Formatted{T}, data) where T = canonical_type(T(), data)

function enforce_canonical_type(f, data)
    AT = canonical_type(f, data)
    # This may not be type stable if `AT` is not a concrete type,
    # but it's not an issue for `load`; it can never be type stable.

    # work around the invokelatest overhead with an eager type check
    if data isa AT
        return data
    else
        # the backend might provide its own convert method
        # use invokelatest to avoid world age issues
        # See issue #34
        return Base.invokelatest(convert, AT, data)
    end
end

# For lazy loads by either of
#   load(io; mmap=true)
#   load(io; lazyio=true)
# do not canonicalize by default. A package that supports both is TiffImages v0.6+.
uses_lazy(kwargs) = get(kwargs, :mmap, false) || get(kwargs, :lazyio, false)

## PNGs


function load(f::File{DataFormat{:PNG}}; kwargs...)
    data = PNGFiles.load(f.filename; kwargs...)
    return enforce_canonical_type(f, data)
end
function load(s::Stream{DataFormat{:PNG}}; kwargs...)
    data = PNGFiles.load(stream(s); kwargs...)
    return enforce_canonical_type(s, data)
end

function save(f::File{DataFormat{:PNG}}, image::S; kwargs...) where {T, S<:Union{AbstractMatrix, AbstractArray{T,3}}}
    PNGFiles.save(f.filename, image; kwargs...)
end

function save(s::Stream{DataFormat{:PNG}}, image::S; permute_horizontal=false, mapi=identity, kwargs...) where {T, S<:Union{AbstractMatrix, AbstractArray{T,3}}}
    imgout = map(mapi, image)
    if permute_horizontal
        perm = ndims(imgout) == 2 ? (2, 1) : ndims(imgout) == 3 ? (2, 1, 3) : error("$(ndims(imgout)) dims array is not supported")
        PNGFiles.save(stream(s), PermutedDimsArray(imgout, perm); kwargs...)
    else
        PNGFiles.save(stream(s), imgout; kwargs...)
    end
end

# Netpbm types

for NETPBMFORMAT in (:PBMBinary, :PGMBinary, :PPMBinary, :PBMText, :PGMText, :PPMText)
    @eval begin
        function load(f::File{DataFormat{$(Expr(:quote,NETPBMFORMAT))}})
            data = Netpbm.load(f)
            return enforce_canonical_type(f, data)
        end

        function load(s::Stream{DataFormat{$(Expr(:quote,NETPBMFORMAT))}})
            data = Netpbm.load(s)
            return enforce_canonical_type(s, data)
        end

        function save(f::File{DataFormat{$(Expr(:quote,NETPBMFORMAT))}}, image::S; kwargs...) where {S<:AbstractMatrix}
            return Netpbm.save(f, image; kwargs...)
        end

        function save(s::Stream{DataFormat{$(Expr(:quote,NETPBMFORMAT))}}, image::S; kwargs...) where {S<:AbstractMatrix}
            return Netpbm.save(s, image; kwargs...)
        end
    end
end

## TIFFs

function load(f::File{DataFormat{:TIFF}}; canonicalize::Union{Bool,Nothing}=nothing, kwargs...)
    canonicalize = something(canonicalize, !uses_lazy(kwargs))
    data = TiffImages.load(f.filename; kwargs...)
    return canonicalize ? enforce_canonical_type(f, data) : data
end
function load(s::Stream{DataFormat{:TIFF}}; canonicalize::Union{Bool,Nothing}=nothing, kwargs...)
    canonicalize = something(canonicalize, !uses_lazy(kwargs))
    data = TiffImages.load(stream(s); kwargs...)
    return canonicalize ? enforce_canonical_type(s, data) : data
end

function save(f::File{DataFormat{:TIFF}}, image::S) where {T, S<:Union{AbstractMatrix, AbstractArray{T,3}}}
    TiffImages.save(f.filename, image)
end

function save(s::Stream{DataFormat{:TIFF}}, image::S; permute_horizontal=false, mapi=identity) where {T, S<:Union{AbstractMatrix, AbstractArray{T,3}}}
    imgout = map(mapi, image)
    if permute_horizontal
        perm = ndims(imgout) == 2 ? (2, 1) : ndims(imgout) == 3 ? (2, 1, 3) : error("$(ndims(imgout)) dims array is not supported")
        TiffImages.save(stream(s), PermutedDimsArray(imgout, perm))
    else
        TiffImages.save(stream(s), imgout)
    end
end

## OpenEXR

function load(f::File{DataFormat{:EXR}}; kwargs...)
    data = OpenEXR.load(f; kwargs...)
    return enforce_canonical_type(f, data)
end

function save(f::File{DataFormat{:EXR}}, args...; kwargs...)
    OpenEXR.save(f, args...; kwargs...)
end

## QOI

function load(f::File{DataFormat{:QOI}}; kwargs...)
    data = QOI.qoi_decode(f.filename; kwargs...)
    return enforce_canonical_type(f, data)
end

function save(f::File{DataFormat{:QOI}}, args...; kwargs...)
    QOI.qoi_encode(f.filename, args...; kwargs...)
end

## Sixel
# Sixel.jl itself provides `fileio_load`/`fileio_save` so we simply delegate everything to it
function load(f::File{DataFormat{:SIXEL}}; kwargs...)
    data = Sixel.fileio_load(f, kwargs...)
    return  enforce_canonical_type(f, data)
end
function load(s::Stream{DataFormat{:SIXEL}}; kwargs...)
    data = Sixel.fileio_load(s, kwargs...)
    return  enforce_canonical_type(s, data)
end
function save(f::File{DataFormat{:SIXEL}}, image::AbstractArray; kwargs...)
    Sixel.fileio_save(f, image; kwargs...)
end
function save(s::Stream{DataFormat{:SIXEL}}, image::AbstractArray; kwargs...)
    Sixel.fileio_save(s, image; kwargs...)
end

## JPEG
function load(f::File{DataFormat{:JPEG}}; kwargs...)
    data = JpegTurbo.fileio_load(f, kwargs...)
    return enforce_canonical_type(f, data)
end
function load(s::Stream{DataFormat{:JPEG}}; kwargs...)
    data = JpegTurbo.fileio_load(s, kwargs...)
    return enforce_canonical_type(s, data)
end
function save(f::File{DataFormat{:JPEG}}, image::AbstractArray; kwargs...)
    JpegTurbo.fileio_save(f, image; kwargs...)
end
function save(s::Stream{DataFormat{:JPEG}}, image::AbstractArray; kwargs...)
    JpegTurbo.fileio_save(s, image; kwargs...)
end

## WebP
function load(f::File{DataFormat{:WebP}}; kwargs...)
    data = WebP.fileio_load(f, kwargs...)
    return enforce_canonical_type(f, data)
end
function load(s::Stream{DataFormat{:WebP}}; kwargs...)
    data = WebP.fileio_load(s, kwargs...)
    return enforce_canonical_type(s, data)
end
function save(f::File{DataFormat{:WebP}}, image::AbstractArray; kwargs...)
    WebP.fileio_save(f, image; kwargs...)
end
function save(s::Stream{DataFormat{:WebP}}, image::AbstractArray; kwargs...)
    WebP.fileio_save(s, image; kwargs...)
end

## Function names labelled for FileIO. Makes FileIO lookup quicker
const fileio_save = save
const fileio_load = load

end # module
