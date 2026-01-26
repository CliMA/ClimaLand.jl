module OpenEXR

export load_exr, save_exr

using FileIO, Colors

module C
using OpenEXR_jll
using Colors
const ImfHalf = Float16
const ImfRgba = RGBA{ImfHalf}
include("OpenEXR_common.jl")
include("OpenEXR_api.jl")
const IMF_WRITE_RGBA = IMF_WRITE_RGB + IMF_WRITE_A
end  # module C

const MAGIC = Cint(C.IMF_MAGIC)

@enum Compression::Cint begin
    NO_COMPRESSION = C.IMF_NO_COMPRESSION
    RLE_COMPRESSION = C.IMF_RLE_COMPRESSION
    ZIPS_COMPRESSION = C.IMF_ZIPS_COMPRESSION
    ZIP_COMPRESSION = C.IMF_ZIP_COMPRESSION
    PIZ_COMPRESSION = C.IMF_PIZ_COMPRESSION
    PXR24_COMPRESSION = C.IMF_PXR24_COMPRESSION
    B44_COMPRESSION = C.IMF_B44_COMPRESSION
    B44A_COMPRESSION = C.IMF_B44A_COMPRESSION
    DWAA_COMPRESSION = C.IMF_DWAA_COMPRESSION
    DWAB_COMPRESSION = C.IMF_DWAB_COMPRESSION
end

@enum RgbaChannels::Cint begin
    WRITE_R = C.IMF_WRITE_R
    WRITE_G = C.IMF_WRITE_G
    WRITE_B = C.IMF_WRITE_B
    WRITE_A = C.IMF_WRITE_A
    WRITE_Y = C.IMF_WRITE_Y
    WRITE_C = C.IMF_WRITE_C
    WRITE_RGB = C.IMF_WRITE_RGB
    WRITE_RGBA = C.IMF_WRITE_RGBA
    WRITE_YC = C.IMF_WRITE_YC
    WRITE_YA = C.IMF_WRITE_YA
    WRITE_YCA = C.IMF_WRITE_YCA
end

function check(ret)
    ret == typeof(ret)(0) && error(unsafe_string(C.ImfErrorMessage()))
end

"""
    load_exr(filename)::(Array{RGBA{Float16},2}, RgbaChannels)

Returns the image data contained in `filename` along with flags representing the on-disk
storage format.
"""
function load_exr(filename)
    infile = C.ImfOpenInputFile(filename)  # open the file
    check(infile)
    chans = RgbaChannels(C.ImfInputChannels(infile))
    try
        # get the header
        hdr = C.ImfInputHeader(infile)

        # read its data window
        xmin = Ref{Cint}()
        ymin = Ref{Cint}()
        xmax = Ref{Cint}()
        ymax = Ref{Cint}()
        C.ImfHeaderDataWindow(hdr, xmin, ymin, xmax, ymax)

        # compute the window size
        width = xmax[] - xmin[] + 1
        height = ymax[] - ymin[] + 1

        # allocate space for the result and get its strides
        data = Array{C.ImfRgba,2}(undef, height, width)
        (xstride, ystride) = strides(data)

        # get the pointer to the data, shifting it according to the expected window
        dataptr =
            Base.unsafe_convert(Ptr{C.ImfRgba}, data) - xmin[] * xstride - ymin[] * ystride

        # copy the data
        check(C.ImfInputSetFrameBuffer(infile, dataptr, ystride, xstride))
        check(C.ImfInputReadPixels(infile, ymin[], ymax[]))

        # return the loaded raster along with the channels
        return (data, chans)
    finally
        check(C.ImfCloseInputFile(infile))
    end
end

"""
    save_exr(filename, image[, compression[, channels]])

Save `image` as an OpenEXR file in `filename`, storing the data in a format
indicated by `channels` using the `compression` algorithm.
"""
function save_exr(
    filename,
    image::AbstractArray{C.ImfRgba,2},
    compression::Compression = ZIP_COMPRESSION,
    channels::RgbaChannels = WRITE_RGBA,
)
    # get the size of the data
    (height, width) = size(image)

    # create a new header
    hdr = C.ImfNewHeader()
    check(hdr)
    try
        # set the compression
        C.ImfHeaderSetCompression(hdr, compression)

        # set the correct window sizes
        C.ImfHeaderSetDataWindow(hdr, 0, 0, width - 1, height - 1)
        C.ImfHeaderSetDisplayWindow(hdr, 0, 0, width - 1, height - 1)

        # open the output file
        outfile = C.ImfOpenOutputFile(filename, hdr, channels)
        check(outfile)
        try
            # get the strides and a pointer to the raster
            (xstride, ystride) = strides(image)
            dataptr = Base.unsafe_convert(Ptr{C.ImfRgba}, image)

            # copy the data
            check(C.ImfOutputSetFrameBuffer(outfile, dataptr, ystride, xstride))
            check(C.ImfOutputWritePixels(outfile, height))
        finally
            check(C.ImfCloseOutputFile(outfile))
        end
    finally
        C.ImfDeleteHeader(hdr)
    end
    nothing
end

function save_exr(
    filename,
    image::AbstractArray{T,2},
    compression::Compression = ZIP_COMPRESSION,
    channels::RgbaChannels = WRITE_RGBA,
) where {T}
    save_exr(filename, (c -> convert(C.ImfRgba, c)).(image), compression, channels)
end

"""
    load(filename)::Array{[RGB|RGBA|Gray|GrayA]{Float16},2}

Returns the image data contained in `filename`.
"""
function load(filename::AbstractString)
    (rgba, chans) = load_exr(filename)
    if chans == WRITE_YA
        return (c -> convert(GrayA{Float16}, c)).(rgba)
    elseif chans == WRITE_Y
        return (c -> convert(Gray{Float16}, c)).(rgba)
    elseif UInt32(chans) & UInt32(WRITE_A) == 0x00
        return (c -> convert(RGB{Float16}, c)).(rgba)
    else
        return rgba
    end
end

"""
    save(filename, image[, compression])

Save `image` as an OpenEXR file in `filename` using the `compression` algorithm.
"""
function save(
    filename::AbstractString,
    image::AbstractArray{T,2},
    compression::Compression = ZIP_COMPRESSION,
) where {T<:Transparent3}
    save_exr(filename, image, compression, WRITE_RGBA)
end

function save(
    filename::AbstractString,
    image::AbstractArray{T,2},
    compression::Compression = ZIP_COMPRESSION,
) where {T<:Color3}
    save_exr(filename, image, compression, WRITE_RGB)
end

function save(
    filename::AbstractString,
    image::AbstractArray{T,2},
    compression::Compression = ZIP_COMPRESSION,
) where {T<:TransparentGray}
    save_exr(filename, image, compression, WRITE_YA)
end

function save(
    filename::AbstractString,
    image::AbstractArray{T,2},
    compression::Compression = ZIP_COMPRESSION,
) where {T<:AbstractGray}
    save_exr(filename, image, compression, WRITE_Y)
end

# FileIO interface
load(f::File{DataFormat{:EXR}}, args...) = load(f.filename, args...)
save(f::File{DataFormat{:EXR}}, args...) = save(f.filename, args...)

end  # module OpenEXR
