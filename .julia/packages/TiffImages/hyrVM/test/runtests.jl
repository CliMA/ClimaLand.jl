using ColorTypes
using ColorVectorSpace
using Documenter
using FixedPointNumbers
using OffsetArrays
import AxisArrays: AxisArray # No, we don't want to import `AxisArrays.axes`
using Statistics
using Test
using TiffImages

include("Aqua.jl")

DocMeta.setdocmeta!(TiffImages, :DocTestSetup, :(using TiffImages); recursive=true)
doctest(TiffImages)

_wrap(name) = "https://github.com/tlnagy/exampletiffs/blob/master/$name?raw=true"

using Downloads
get_example(x) = Downloads.download(_wrap(x))

@testset "Tag loading" begin
    include("tags.jl")
end

@testset "Gray with ExtraSample" begin
    filepath = get_example("house.tif")
    img = TiffImages.load(filepath)
    @test size(img) == (512, 512)
    @test eltype(img) == GrayA{N0f8}
    @test img[50,50] == GrayA{N0f8}(0.804, 1.0) # value from ImageMagick.jl
    img[50:300, 50:150] .= 0.0
    @test img[50, 50] == GrayA{N0f8}(0.0, 1.0)

    # Efficient convert method
    img_cvt = convert(Array{eltype(img), ndims(img)}, img)
    @test img_cvt === img.data
end

@testset "MRI stack" begin
    filepath = get_example("mri.tif")
    img = TiffImages.load(filepath)
    @test size(img) == (128, 128, 27)
    @test eltype(img) == RGB{N0f16}

    # TODO: inefficient convert method
    img_cvt = convert(Array{eltype(img), ndims(img)}, img)
    @test img_cvt == img.data
    @test img_cvt !== img.data

    @testset "Multiple planes" begin
        # test slice based access
        @test red(img[92, 21, 2]) == 0.74118N0f16

        # make sure IFDs are included if whole slice is grabbed
        @test typeof(ifds(img[:, :, 2])) <: TiffImages.IFD
        @test length(ifds(img[:, :, 2:5])) == 4

        # make sure that subslicing in XY drops IFD info
        @test typeof(img[:, 50:60, 27]) == Array{RGB{Normed{UInt16,16}},2}
    end
end

@testset "Floating point RGB image" begin
    filepath = get_example("spring.tif")
    img = TiffImages.load(filepath)
    @test size(img) == (619, 858)
    @test eltype(img) == RGB{Float16}
    
    # Efficient convert method
    img_cvt = convert(Array{eltype(img), ndims(img)}, img)
    @test img_cvt === img.data
end

@testset "Packbits image" begin
    filepath = get_example("coffee.tif")
    img = TiffImages.load(filepath)
    @test size(img) == (378, 504)
    @test eltype(img) == Gray{N0f8}
    
    # Efficient convert method
    img_cvt = convert(Array{eltype(img), ndims(img)}, img)
    @test img_cvt === img.data
end

@testset "Adobe Deflate image" begin
    filepath = get_example("underwater_bmx.tif")
    img = TiffImages.load(filepath)
    @test size(img) == (773, 1076)
    @test eltype(img) == RGB{N0f8}

    # Efficient convert method
    img_cvt = convert(Array{eltype(img), ndims(img)}, img)
    @test img_cvt === img.data

    # Dirty patch the TIFF_COMPRESSION tag from
    # COMPRESSION_ADOBE_DEFLATE (8) to COMPRESSION_DEFLATE (32946).
    # Those are in fact exactly the same so nothing else needs to be
    # adjusted.
    data = read(filepath)
    @assert data[1063947:1063948] == [0x08, 0x00]
    data[1063947:1063948] .= reinterpret(UInt8, [UInt16(32946)])
    _, io = mktemp()
    write(io, data)
    seekstart(io)
    img = TiffImages.load(io)
    @test size(img) == (773, 1076)
    @test eltype(img) == RGB{N0f8}
end

@testset "Bilevel image" begin
    filepath = get_example("capitol.tif")
    img = TiffImages.load(filepath)
    @test size(img) == (378, 504)
    @test eltype(img) == Gray{N7f1}
    
    # Efficient convert method
    img_cvt = convert(Array{eltype(img), ndims(img)}, img)
    @test img_cvt === img.data
end

@testset "Striped bilevel image" begin
    filepath = get_example("capitol2.tif")
    img = TiffImages.load(filepath)
    @test size(img) == (378, 504)
    @test eltype(img) == Gray{N7f1}

    # Efficient convert method
    img_cvt = convert(Array{eltype(img), ndims(img)}, img)
    @test img_cvt === img.data
end

@testset "Signed integer type" begin
    filepath = get_example("4D-series.ome.tif")
    img = TiffImages.load(filepath)

    @test size(img) == (167, 439, 35)
    expected_rng = reinterpret.(Q0f7, Int8.((-1, 96)))
    @test extrema(img[:, :, 1]) == expected_rng

    # Efficient convert method
    img_cvt = convert(Array{eltype(img), ndims(img)}, img)
    @test img_cvt === img.data
end

@testset "Discontiguous striped image, Issue #38" begin
    filepath = get_example("julia.tif")
    img = TiffImages.load(filepath)

    @test size(img) == (300, 500)
    # if discontiguous striping is broken then the garbage padding data will
    # leak into the actual image
    @test all(img[3, 1:50] .== RGB{N0f8}(1, 1, 1))

    # Efficient convert method
    img_cvt = convert(Array{eltype(img), ndims(img)}, img)
    @test img_cvt === img.data
end

@testset "Issue #12" begin
    @testset "Big endian striped file" begin
        filepath = get_example("flagler.tif")
        img = TiffImages.load(filepath)

        # verify that strip offsets are correctly bswapped for endianness
        img_stripoffsets = Int.(ifds(img)[TiffImages.STRIPOFFSETS].data)
        @test img_stripoffsets == [8, 129848, 259688, 389528]

        # Efficient convert method
        img_cvt = convert(Array{eltype(img), ndims(img)}, img)
        @test img_cvt === img.data
    end

    @testset "Little endian striped file" begin
        filepath = get_example("house.tif")
        img = TiffImages.load(filepath)

        # verify that strip offsets are correctly bswapped for endianness
        img_stripoffsets = Int.(ifds(img)[TiffImages.STRIPOFFSETS].data)
        @test issorted(img_stripoffsets)

        # Efficient convert method
        img_cvt = convert(Array{eltype(img), ndims(img)}, img)
        @test img_cvt === img.data
    end
end

@testset "Mmap & lazyio" begin
    include("mmap_lazyio.jl")
end

@testset "Writing" begin
    include("writer.jl")
end

@testset "Interpreting IFD layouts" begin
    include("layouts.jl")
end

@testset "Issue #69" begin
    rawarray = Gray.(zeros(10, 10, 2))
    ifds = [TiffImages.IFD(UInt64), TiffImages.IFD(UInt64)]
    # test that the constructor can handle small images using 64bit offsets 
    @test size(TiffImages.DenseTaggedImage(rawarray, ifds)) == size(rawarray)
end

@testset "Issue #96" begin
    # Handle TIFFs where multiple slices are stored in the same strip. These are
    # still contiguous since the entire slice is within one strip.
    ifd = TiffImages.IFD(UInt32)
    ifd[TiffImages.IMAGELENGTH] = 16
    ifd[TiffImages.ROWSPERSTRIP] = 256

    @test TiffImages.iscontiguous(ifd)
end

@testset "LZW" begin
    uncompressed = get_example("shapes_uncompressed.tif")
    compressed = get_example("shapes_lzw.tif")
    @test TiffImages.load(uncompressed) == TiffImages.load(compressed)
end

@testset "Tiled" begin
    uncompressed = get_example("shapes_uncompressed.tif")
    compressed_tiled = get_example("shapes_lzw_tiled.tif")
    compressed_tiled_multi = get_example("shapes_tiled_multi.tif")
    uncompressed_loaded = TiffImages.load(uncompressed)
    compressed_tiled_multi_loaded = TiffImages.load(compressed_tiled_multi; lazyio=true)
    @test uncompressed_loaded == TiffImages.load(compressed_tiled)
    @test uncompressed_loaded[:] == TiffImages.load(compressed_tiled; lazyio=true)[:]
    @test compressed_tiled_multi_loaded[:,:,1] == uncompressed_loaded
    @test compressed_tiled_multi_loaded[:,:,2] == uncompressed_loaded
    @test compressed_tiled_multi_loaded[:,:,3] == uncompressed_loaded
end

@testset "Planar" begin
    ref = get_example("shapes_uncompressed.tif")
    planar = get_example("shapes_lzw_planar.tif")
    tiled_planar = get_example("shapes_lzw_tiled_planar.tif")
    uncompressed_tiled_planar = get_example("shapes_uncompressed_tiled_planar.tif")
    @test TiffImages.load(ref) == TiffImages.load(planar)
    @test TiffImages.load(ref) == TiffImages.load(tiled_planar)
    @test TiffImages.load(ref) == TiffImages.load(uncompressed_tiled_planar)

    for typ in [Int8,UInt16,Float32]
        for planes in 1:33
            for size in 64:164
                out = Vector{typ}(undef, size * planes)
                a=reduce(vcat,[fill(typ(x),size) for x in 1:planes])
                b=copy(a)
                TiffImages.deplane!(out, a, Val(planes))
                @test a == TiffImages.deplane_slow(b, planes)
            end
        end
    end

    for typ in [Int8,UInt16,Float32]
        for planes in 1:33
            for size in 1:164
                a=reduce(vcat,[fill(typ(x),size) for x in 1:planes])
                b=copy(a)
                TiffImages.deplane!(a, planes)
                @test a == TiffImages.deplane_slow(b, planes)
            end
        end
    end
end

if VERSION < v"1.7"
    partition(s::String, n) = map(xs -> reduce(*, xs), Iterators.partition(s, n))
else
    partition = Iterators.partition
end

@testset "Arbitrary Bit Depth" begin
    ref = TiffImages.load(get_example("shapes_uncompressed.tif"))
    other = TiffImages.load(get_example("shapes_lzw_12bps.tif"))
    m = sum(ref.data .- other.data) ./ length(ref)
    @test m.r < 0.001 && m.g < 0.001 && m.b < 0.001

    @test ifds(other)[TiffImages.BITSPERSAMPLE].data == UInt16[12,12,12]

    other = TiffImages.load(get_example("shapes_lzw_14bps.tif"))
    m = sum(ref.data .- other.data) ./ length(ref)
    @test m.r < 0.001 && m.g < 0.001 && m.b < 0.001

    xs = [rand(UInt8) & 0x7f for _ in 1:160];
    bytes = parse.(UInt8, partition(reduce(*,string.(xs; base=2, pad=7)),8); base=2)
    resize!(bytes, fld(length(bytes) * 8, 7))
    unpacked = TiffImages.unpack_integers_simd(bytes, Val(7))

    @test xs == unpacked

    xs = [rand(UInt32) & 0x1fffffff for _ in 1:256];
    bytes = parse.(UInt8, partition(reduce(*,string.(xs; base=2, pad=29)),8); base=2)
    resize!(bytes, fld(length(bytes) * 32, 29))
    unpacked = TiffImages.unpack_integers_simd(reinterpret(UInt32, bytes), Val(29))

    @test xs == unpacked

    xs = [rand(UInt8) & 0x0f for _ in 1:128];
    bytes = parse.(UInt8, partition(reduce(*,string.(xs; base=2, pad=4)),8); base=2)
    resize!(bytes, fld(length(bytes) * 8, 4))
    unpacked = TiffImages.unpack_integers_simd(bytes, Val(4))

    @test xs == unpacked

    xs = [rand(UInt16) & 0x1ff for _ in 1:124];
    bytes = parse.(UInt8, partition(rpad(reduce(*,string.(xs; base=2, pad=9)), cld(124 * 9, 8) * 8, '0'),8); base=2)
    resize!(bytes, fld(length(bytes) * 16, 9))
    unpacked = TiffImages.unpack_integers_slow(reinterpret(UInt16, bytes), 1, 124, 9)

    @test xs == unpacked

    xs = [rand(UInt32) & 0x7ffffff for _ in 1:124];
    bytes = parse.(UInt8, partition(rpad(reduce(*,string.(xs; base=2, pad=27)), cld(124 * 27, 8) * 8, '0'),8); base=2)
    resize!(bytes, fld(length(bytes) * 32, 27))
    unpacked = TiffImages.unpack_integers(reinterpret(UInt32, bytes), 27)

    @test xs == unpacked
end

@testset "predictor == 3" begin
    original = get_example("shapes_uncompressed.tif")
    encoded = get_example("shapes_lzw_predictor3.tif")
    @test TiffImages.load(original) == TiffImages.load(encoded)
end

@testset "Issue #148" begin
    im = get_example("earthlab.tif")
    @test size(TiffImages.load(im)) == (2400, 2400) # no error
end

@testset "Ragged" begin
    original = TiffImages.load(get_example("shapes_uncompressed.tif"))
    half = TiffImages.load(get_example("shapes_uncompressed_half.tif"))
    palette = TiffImages.load(get_example("shapes_lzw_palette.tif"))
    multisize = TiffImages.load(get_example("shapes_multi_size.tif"))
    multicolor = TiffImages.load(get_example("shapes_multi_color.tif"))

    @test original == multisize[1]
    @test half == multisize[2]

    @test original == multicolor[1]
    @test palette == multicolor[2]
    @test original == multicolor[3]
    @test sum(convert.(Float64, reinterpret(N0f8, original)) .- convert.(Float64, reinterpret(N4f12, multicolor[4]))) ./ (216*128) < 0.0001
    @test original == multicolor[5]
end

@testset "WidePixel" begin
    original = TiffImages.load(get_example("shapes_uncompressed.tif"))
    hyper = TiffImages.load(get_example("shapes_hyper.tif"))

    @test TiffImages.color(hyper) == original
    @test TiffImages.color(original) == original

    @test TiffImages.nchannels(original) == 3
    @test TiffImages.nchannels(hyper) == 7

    @test TiffImages.channel(hyper[100], 4) == 0.1f0
    @test TiffImages.channel(hyper[200], 5) == 0.2f0
    @test TiffImages.channel(hyper[300], 6) == 0.3f0
    @test TiffImages.channel(hyper[400], 7) == 0.5f0

    @test isapprox(sum(sum.(TiffImages.channel.(Ref(original), 1:3))), 23145.05882)
    @test isapprox(sum(sum.(TiffImages.channel.(Ref(hyper), 1:7))), 33282.65886)
end

@testset "Palette color images" begin
    img = TiffImages.load(get_example("poppies.tif"))
    
    @test ifds(img)[TiffImages.BITSPERSAMPLE].data == UInt16(8)
    @test ifds(img)[TiffImages.SAMPLESPERPIXEL].data == UInt16(1)
    @test ifds(img)[TiffImages.PHOTOMETRIC].data == UInt16(TiffImages.PHOTOMETRIC_PALETTE)
end