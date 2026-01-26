# Julia wrapper for header: ImfCRgbaFile.h
# Automatically generated using Clang.jl


function ImfFloatToHalf(f, h)
    ccall((:ImfFloatToHalf, libOpenEXR), Cvoid, (Cfloat, Ptr{ImfHalf}), f, h)
end

function ImfFloatToHalfArray(n, f, h)
    ccall((:ImfFloatToHalfArray, libOpenEXR), Cvoid, (Cint, Ptr{Cfloat}, Ptr{ImfHalf}), n, f, h)
end

function ImfHalfToFloat(h)
    ccall((:ImfHalfToFloat, libOpenEXR), Cfloat, (ImfHalf,), h)
end

function ImfHalfToFloatArray(n, h, f)
    ccall((:ImfHalfToFloatArray, libOpenEXR), Cvoid, (Cint, Ptr{ImfHalf}, Ptr{Cfloat}), n, h, f)
end

function ImfNewHeader()
    ccall((:ImfNewHeader, libOpenEXR), Ptr{ImfHeader}, ())
end

function ImfDeleteHeader(hdr)
    ccall((:ImfDeleteHeader, libOpenEXR), Cvoid, (Ptr{ImfHeader},), hdr)
end

function ImfCopyHeader(hdr)
    ccall((:ImfCopyHeader, libOpenEXR), Ptr{ImfHeader}, (Ptr{ImfHeader},), hdr)
end

function ImfHeaderSetDisplayWindow(hdr, xMin, yMin, xMax, yMax)
    ccall((:ImfHeaderSetDisplayWindow, libOpenEXR), Cvoid, (Ptr{ImfHeader}, Cint, Cint, Cint, Cint), hdr, xMin, yMin, xMax, yMax)
end

function ImfHeaderDisplayWindow(hdr, xMin, yMin, xMax, yMax)
    ccall((:ImfHeaderDisplayWindow, libOpenEXR), Cvoid, (Ptr{ImfHeader}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}), hdr, xMin, yMin, xMax, yMax)
end

function ImfHeaderSetDataWindow(hdr, xMin, yMin, xMax, yMax)
    ccall((:ImfHeaderSetDataWindow, libOpenEXR), Cvoid, (Ptr{ImfHeader}, Cint, Cint, Cint, Cint), hdr, xMin, yMin, xMax, yMax)
end

function ImfHeaderDataWindow(hdr, xMin, yMin, xMax, yMax)
    ccall((:ImfHeaderDataWindow, libOpenEXR), Cvoid, (Ptr{ImfHeader}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}), hdr, xMin, yMin, xMax, yMax)
end

function ImfHeaderSetPixelAspectRatio(hdr, pixelAspectRatio)
    ccall((:ImfHeaderSetPixelAspectRatio, libOpenEXR), Cvoid, (Ptr{ImfHeader}, Cfloat), hdr, pixelAspectRatio)
end

function ImfHeaderPixelAspectRatio(hdr)
    ccall((:ImfHeaderPixelAspectRatio, libOpenEXR), Cfloat, (Ptr{ImfHeader},), hdr)
end

function ImfHeaderSetScreenWindowCenter(hdr, x, y)
    ccall((:ImfHeaderSetScreenWindowCenter, libOpenEXR), Cvoid, (Ptr{ImfHeader}, Cfloat, Cfloat), hdr, x, y)
end

function ImfHeaderScreenWindowCenter(hdr, x, y)
    ccall((:ImfHeaderScreenWindowCenter, libOpenEXR), Cvoid, (Ptr{ImfHeader}, Ptr{Cfloat}, Ptr{Cfloat}), hdr, x, y)
end

function ImfHeaderSetScreenWindowWidth(hdr, width)
    ccall((:ImfHeaderSetScreenWindowWidth, libOpenEXR), Cvoid, (Ptr{ImfHeader}, Cfloat), hdr, width)
end

function ImfHeaderScreenWindowWidth(hdr)
    ccall((:ImfHeaderScreenWindowWidth, libOpenEXR), Cfloat, (Ptr{ImfHeader},), hdr)
end

function ImfHeaderSetLineOrder(hdr, lineOrder)
    ccall((:ImfHeaderSetLineOrder, libOpenEXR), Cvoid, (Ptr{ImfHeader}, Cint), hdr, lineOrder)
end

function ImfHeaderLineOrder(hdr)
    ccall((:ImfHeaderLineOrder, libOpenEXR), Cint, (Ptr{ImfHeader},), hdr)
end

function ImfHeaderSetCompression(hdr, compression)
    ccall((:ImfHeaderSetCompression, libOpenEXR), Cvoid, (Ptr{ImfHeader}, Cint), hdr, compression)
end

function ImfHeaderCompression(hdr)
    ccall((:ImfHeaderCompression, libOpenEXR), Cint, (Ptr{ImfHeader},), hdr)
end

function ImfHeaderSetIntAttribute(hdr, name, value)
    ccall((:ImfHeaderSetIntAttribute, libOpenEXR), Cint, (Ptr{ImfHeader}, Ptr{UInt8}, Cint), hdr, name, value)
end

function ImfHeaderIntAttribute(hdr, name, value)
    ccall((:ImfHeaderIntAttribute, libOpenEXR), Cint, (Ptr{ImfHeader}, Ptr{UInt8}, Ptr{Cint}), hdr, name, value)
end

function ImfHeaderSetFloatAttribute(hdr, name, value)
    ccall((:ImfHeaderSetFloatAttribute, libOpenEXR), Cint, (Ptr{ImfHeader}, Ptr{UInt8}, Cfloat), hdr, name, value)
end

function ImfHeaderSetDoubleAttribute(hdr, name, value)
    ccall((:ImfHeaderSetDoubleAttribute, libOpenEXR), Cint, (Ptr{ImfHeader}, Ptr{UInt8}, Cdouble), hdr, name, value)
end

function ImfHeaderFloatAttribute(hdr, name, value)
    ccall((:ImfHeaderFloatAttribute, libOpenEXR), Cint, (Ptr{ImfHeader}, Ptr{UInt8}, Ptr{Cfloat}), hdr, name, value)
end

function ImfHeaderDoubleAttribute(hdr, name, value)
    ccall((:ImfHeaderDoubleAttribute, libOpenEXR), Cint, (Ptr{ImfHeader}, Ptr{UInt8}, Ptr{Cdouble}), hdr, name, value)
end

function ImfHeaderSetStringAttribute(hdr, name, value)
    ccall((:ImfHeaderSetStringAttribute, libOpenEXR), Cint, (Ptr{ImfHeader}, Ptr{UInt8}, Ptr{UInt8}), hdr, name, value)
end

function ImfHeaderStringAttribute(hdr, name, value)
    ccall((:ImfHeaderStringAttribute, libOpenEXR), Cint, (Ptr{ImfHeader}, Ptr{UInt8}, Ptr{Cstring}), hdr, name, value)
end

function ImfHeaderSetBox2iAttribute(hdr, name, xMin, yMin, xMax, yMax)
    ccall((:ImfHeaderSetBox2iAttribute, libOpenEXR), Cint, (Ptr{ImfHeader}, Ptr{UInt8}, Cint, Cint, Cint, Cint), hdr, name, xMin, yMin, xMax, yMax)
end

function ImfHeaderBox2iAttribute(hdr, name, xMin, yMin, xMax, yMax)
    ccall((:ImfHeaderBox2iAttribute, libOpenEXR), Cint, (Ptr{ImfHeader}, Ptr{UInt8}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}), hdr, name, xMin, yMin, xMax, yMax)
end

function ImfHeaderSetBox2fAttribute(hdr, name, xMin, yMin, xMax, yMax)
    ccall((:ImfHeaderSetBox2fAttribute, libOpenEXR), Cint, (Ptr{ImfHeader}, Ptr{UInt8}, Cfloat, Cfloat, Cfloat, Cfloat), hdr, name, xMin, yMin, xMax, yMax)
end

function ImfHeaderBox2fAttribute(hdr, name, xMin, yMin, xMax, yMax)
    ccall((:ImfHeaderBox2fAttribute, libOpenEXR), Cint, (Ptr{ImfHeader}, Ptr{UInt8}, Ptr{Cfloat}, Ptr{Cfloat}, Ptr{Cfloat}, Ptr{Cfloat}), hdr, name, xMin, yMin, xMax, yMax)
end

function ImfHeaderSetV2iAttribute(hdr, name, x, y)
    ccall((:ImfHeaderSetV2iAttribute, libOpenEXR), Cint, (Ptr{ImfHeader}, Ptr{UInt8}, Cint, Cint), hdr, name, x, y)
end

function ImfHeaderV2iAttribute(hdr, name, x, y)
    ccall((:ImfHeaderV2iAttribute, libOpenEXR), Cint, (Ptr{ImfHeader}, Ptr{UInt8}, Ptr{Cint}, Ptr{Cint}), hdr, name, x, y)
end

function ImfHeaderSetV2fAttribute(hdr, name, x, y)
    ccall((:ImfHeaderSetV2fAttribute, libOpenEXR), Cint, (Ptr{ImfHeader}, Ptr{UInt8}, Cfloat, Cfloat), hdr, name, x, y)
end

function ImfHeaderV2fAttribute(hdr, name, x, y)
    ccall((:ImfHeaderV2fAttribute, libOpenEXR), Cint, (Ptr{ImfHeader}, Ptr{UInt8}, Ptr{Cfloat}, Ptr{Cfloat}), hdr, name, x, y)
end

function ImfHeaderSetV3iAttribute(hdr, name, x, y, z)
    ccall((:ImfHeaderSetV3iAttribute, libOpenEXR), Cint, (Ptr{ImfHeader}, Ptr{UInt8}, Cint, Cint, Cint), hdr, name, x, y, z)
end

function ImfHeaderV3iAttribute(hdr, name, x, y, z)
    ccall((:ImfHeaderV3iAttribute, libOpenEXR), Cint, (Ptr{ImfHeader}, Ptr{UInt8}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}), hdr, name, x, y, z)
end

function ImfHeaderSetV3fAttribute(hdr, name, x, y, z)
    ccall((:ImfHeaderSetV3fAttribute, libOpenEXR), Cint, (Ptr{ImfHeader}, Ptr{UInt8}, Cfloat, Cfloat, Cfloat), hdr, name, x, y, z)
end

function ImfHeaderV3fAttribute(hdr, name, x, y, z)
    ccall((:ImfHeaderV3fAttribute, libOpenEXR), Cint, (Ptr{ImfHeader}, Ptr{UInt8}, Ptr{Cfloat}, Ptr{Cfloat}, Ptr{Cfloat}), hdr, name, x, y, z)
end

function ImfHeaderSetM33fAttribute(hdr, name, m)
    ccall((:ImfHeaderSetM33fAttribute, libOpenEXR), Cint, (Ptr{ImfHeader}, Ptr{UInt8}, Ptr{NTuple{3, Cfloat}}), hdr, name, m)
end

function ImfHeaderM33fAttribute(hdr, name, m)
    ccall((:ImfHeaderM33fAttribute, libOpenEXR), Cint, (Ptr{ImfHeader}, Ptr{UInt8}, Ptr{NTuple{3, Cfloat}}), hdr, name, m)
end

function ImfHeaderSetM44fAttribute(hdr, name, m)
    ccall((:ImfHeaderSetM44fAttribute, libOpenEXR), Cint, (Ptr{ImfHeader}, Ptr{UInt8}, Ptr{NTuple{4, Cfloat}}), hdr, name, m)
end

function ImfHeaderM44fAttribute(hdr, name, m)
    ccall((:ImfHeaderM44fAttribute, libOpenEXR), Cint, (Ptr{ImfHeader}, Ptr{UInt8}, Ptr{NTuple{4, Cfloat}}), hdr, name, m)
end

function ImfOpenOutputFile(name, hdr, channels)
    ccall((:ImfOpenOutputFile, libOpenEXR), Ptr{ImfOutputFile}, (Ptr{UInt8}, Ptr{ImfHeader}, Cint), name, hdr, channels)
end

function ImfCloseOutputFile(out)
    ccall((:ImfCloseOutputFile, libOpenEXR), Cint, (Ptr{ImfOutputFile},), out)
end

function ImfOutputSetFrameBuffer(out, base, xStride, yStride)
    ccall((:ImfOutputSetFrameBuffer, libOpenEXR), Cint, (Ptr{ImfOutputFile}, Ptr{ImfRgba}, Csize_t, Csize_t), out, base, xStride, yStride)
end

function ImfOutputWritePixels(out, numScanLines)
    ccall((:ImfOutputWritePixels, libOpenEXR), Cint, (Ptr{ImfOutputFile}, Cint), out, numScanLines)
end

function ImfOutputCurrentScanLine(out)
    ccall((:ImfOutputCurrentScanLine, libOpenEXR), Cint, (Ptr{ImfOutputFile},), out)
end

function ImfOutputHeader(out)
    ccall((:ImfOutputHeader, libOpenEXR), Ptr{ImfHeader}, (Ptr{ImfOutputFile},), out)
end

function ImfOutputChannels(out)
    ccall((:ImfOutputChannels, libOpenEXR), Cint, (Ptr{ImfOutputFile},), out)
end

function ImfOpenTiledOutputFile(name, hdr, channels, xSize, ySize, mode, rmode)
    ccall((:ImfOpenTiledOutputFile, libOpenEXR), Ptr{ImfTiledOutputFile}, (Ptr{UInt8}, Ptr{ImfHeader}, Cint, Cint, Cint, Cint, Cint), name, hdr, channels, xSize, ySize, mode, rmode)
end

function ImfCloseTiledOutputFile(out)
    ccall((:ImfCloseTiledOutputFile, libOpenEXR), Cint, (Ptr{ImfTiledOutputFile},), out)
end

function ImfTiledOutputSetFrameBuffer(out, base, xStride, yStride)
    ccall((:ImfTiledOutputSetFrameBuffer, libOpenEXR), Cint, (Ptr{ImfTiledOutputFile}, Ptr{ImfRgba}, Csize_t, Csize_t), out, base, xStride, yStride)
end

function ImfTiledOutputWriteTile(out, dx, dy, lx, ly)
    ccall((:ImfTiledOutputWriteTile, libOpenEXR), Cint, (Ptr{ImfTiledOutputFile}, Cint, Cint, Cint, Cint), out, dx, dy, lx, ly)
end

function ImfTiledOutputWriteTiles(out, dxMin, dxMax, dyMin, dyMax, lx, ly)
    ccall((:ImfTiledOutputWriteTiles, libOpenEXR), Cint, (Ptr{ImfTiledOutputFile}, Cint, Cint, Cint, Cint, Cint, Cint), out, dxMin, dxMax, dyMin, dyMax, lx, ly)
end

function ImfTiledOutputHeader(out)
    ccall((:ImfTiledOutputHeader, libOpenEXR), Ptr{ImfHeader}, (Ptr{ImfTiledOutputFile},), out)
end

function ImfTiledOutputChannels(out)
    ccall((:ImfTiledOutputChannels, libOpenEXR), Cint, (Ptr{ImfTiledOutputFile},), out)
end

function ImfTiledOutputTileXSize(out)
    ccall((:ImfTiledOutputTileXSize, libOpenEXR), Cint, (Ptr{ImfTiledOutputFile},), out)
end

function ImfTiledOutputTileYSize(out)
    ccall((:ImfTiledOutputTileYSize, libOpenEXR), Cint, (Ptr{ImfTiledOutputFile},), out)
end

function ImfTiledOutputLevelMode(out)
    ccall((:ImfTiledOutputLevelMode, libOpenEXR), Cint, (Ptr{ImfTiledOutputFile},), out)
end

function ImfTiledOutputLevelRoundingMode(out)
    ccall((:ImfTiledOutputLevelRoundingMode, libOpenEXR), Cint, (Ptr{ImfTiledOutputFile},), out)
end

function ImfOpenInputFile(name)
    ccall((:ImfOpenInputFile, libOpenEXR), Ptr{ImfInputFile}, (Ptr{UInt8},), name)
end

function ImfCloseInputFile(in)
    ccall((:ImfCloseInputFile, libOpenEXR), Cint, (Ptr{ImfInputFile},), in)
end

function ImfInputSetFrameBuffer(in, base, xStride, yStride)
    ccall((:ImfInputSetFrameBuffer, libOpenEXR), Cint, (Ptr{ImfInputFile}, Ptr{ImfRgba}, Csize_t, Csize_t), in, base, xStride, yStride)
end

function ImfInputReadPixels(in, scanLine1, scanLine2)
    ccall((:ImfInputReadPixels, libOpenEXR), Cint, (Ptr{ImfInputFile}, Cint, Cint), in, scanLine1, scanLine2)
end

function ImfInputHeader(in)
    ccall((:ImfInputHeader, libOpenEXR), Ptr{ImfHeader}, (Ptr{ImfInputFile},), in)
end

function ImfInputChannels(in)
    ccall((:ImfInputChannels, libOpenEXR), Cint, (Ptr{ImfInputFile},), in)
end

function ImfInputFileName(in)
    ccall((:ImfInputFileName, libOpenEXR), Cstring, (Ptr{ImfInputFile},), in)
end

function ImfOpenTiledInputFile(name)
    ccall((:ImfOpenTiledInputFile, libOpenEXR), Ptr{ImfTiledInputFile}, (Ptr{UInt8},), name)
end

function ImfCloseTiledInputFile(in)
    ccall((:ImfCloseTiledInputFile, libOpenEXR), Cint, (Ptr{ImfTiledInputFile},), in)
end

function ImfTiledInputSetFrameBuffer(in, base, xStride, yStride)
    ccall((:ImfTiledInputSetFrameBuffer, libOpenEXR), Cint, (Ptr{ImfTiledInputFile}, Ptr{ImfRgba}, Csize_t, Csize_t), in, base, xStride, yStride)
end

function ImfTiledInputReadTile(in, dx, dy, lx, ly)
    ccall((:ImfTiledInputReadTile, libOpenEXR), Cint, (Ptr{ImfTiledInputFile}, Cint, Cint, Cint, Cint), in, dx, dy, lx, ly)
end

function ImfTiledInputReadTiles(in, dxMin, dxMax, dyMin, dyMax, lx, ly)
    ccall((:ImfTiledInputReadTiles, libOpenEXR), Cint, (Ptr{ImfTiledInputFile}, Cint, Cint, Cint, Cint, Cint, Cint), in, dxMin, dxMax, dyMin, dyMax, lx, ly)
end

function ImfTiledInputHeader(in)
    ccall((:ImfTiledInputHeader, libOpenEXR), Ptr{ImfHeader}, (Ptr{ImfTiledInputFile},), in)
end

function ImfTiledInputChannels(in)
    ccall((:ImfTiledInputChannels, libOpenEXR), Cint, (Ptr{ImfTiledInputFile},), in)
end

function ImfTiledInputFileName(in)
    ccall((:ImfTiledInputFileName, libOpenEXR), Cstring, (Ptr{ImfTiledInputFile},), in)
end

function ImfTiledInputTileXSize(in)
    ccall((:ImfTiledInputTileXSize, libOpenEXR), Cint, (Ptr{ImfTiledInputFile},), in)
end

function ImfTiledInputTileYSize(in)
    ccall((:ImfTiledInputTileYSize, libOpenEXR), Cint, (Ptr{ImfTiledInputFile},), in)
end

function ImfTiledInputLevelMode(in)
    ccall((:ImfTiledInputLevelMode, libOpenEXR), Cint, (Ptr{ImfTiledInputFile},), in)
end

function ImfTiledInputLevelRoundingMode(in)
    ccall((:ImfTiledInputLevelRoundingMode, libOpenEXR), Cint, (Ptr{ImfTiledInputFile},), in)
end

function ImfNewRound12logLut(channels)
    ccall((:ImfNewRound12logLut, libOpenEXR), Ptr{ImfLut}, (Cint,), channels)
end

function ImfNewRoundNBitLut(n, channels)
    ccall((:ImfNewRoundNBitLut, libOpenEXR), Ptr{ImfLut}, (UInt32, Cint), n, channels)
end

function ImfDeleteLut(lut)
    ccall((:ImfDeleteLut, libOpenEXR), Cvoid, (Ptr{ImfLut},), lut)
end

function ImfApplyLut(lut, data, nData, stride)
    ccall((:ImfApplyLut, libOpenEXR), Cvoid, (Ptr{ImfLut}, Ptr{ImfRgba}, Cint, Cint), lut, data, nData, stride)
end

function ImfErrorMessage()
    ccall((:ImfErrorMessage, libOpenEXR), Cstring, ())
end
