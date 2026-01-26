using QOI
using Test
using Downloads
using p7zip_jll
using ImageIO
using FileIO
using Scratch
using ColorTypes
using FixedPointNumbers

const DOWNLOAD_SHA = "869a6433a3af7ce84fc55fda6a5387d6c2113c3e8231153549a6407ed1e71696"

function get_testset()
    scratch = get_scratch!(Scratch, DOWNLOAD_SHA)
    if isempty(readdir(scratch))
        p7z = p7zip_jll.p7zip()
        file = Downloads.download("https://github.com/KristofferC/QOI.jl/releases/download/v0.0.0/qoi_test_images.zip")
        run(`$p7z e $file -o$scratch`)
    end
    return scratch
end

function check_roundtrip_qoi(file)
    header, image = QOI.qoi_decode_raw(file)
    data = QOI.qoi_encode_raw(image, header)
    return data == read(file)
end

testset = get_testset()    

# Test roundtrip of raw decode/encode
for file in readdir(testset; join=true)
    if endswith(file, ".qoi")
        @info "testing roundtrip of $file"
        @test check_roundtrip_qoi(file)
    end
end

# Test correctness vs PNG
for file_name in Set(getindex.(splitext.(readdir(testset; join=true)), 1))
    qoi = file_name * ".qoi"
    isfile(qoi) || continue
    @info "comparing to png $qoi"
    png = file_name * ".png"
    img_qoi = QOI.qoi_decode(qoi)
    img_png = load(png)
    @test img_qoi == img_png

    t = tempname()
    QOI.qoi_encode(t, img_qoi)
    @test read(t) == read(qoi)
end

# Test saving imges not in N0f8 format
img_gray = rand(Gray{N0f8}, 5, 5)
t = tempname()
QOI.qoi_encode(t, img_gray)
img_qoi = QOI.qoi_decode(t)
@test convert(Matrix{RGB{N0f8}}, img_gray) == img_qoi

img_24bit = rand(RGBA{Float64}, 5, 5)
t = tempname()
QOI.qoi_encode(t, img_24bit)
img_qoi = QOI.qoi_decode(t)
@test convert(Matrix{RGBA{N0f8}}, img_24bit) == img_qoi


# Invalid images

# Unexpected end
io = IOBuffer()
write(io, "qoif")
write(io, hton(UInt32(10)))
write(io, hton(UInt32(0))) 
@test_throws QOI.QOIException("unexpected end of file") QOI.qoi_decode_raw(take!(io))


# Invalid width
io = IOBuffer()
write(io, "qoif")
write(io, hton.(UInt32.([0, 10])))
write(io, hton.(UInt8.([3, 1])))
@test_throws QOI.QOIException("invalid width in header, got 0") QOI.qoi_decode_raw(take!(io))

# Invalid height
io = IOBuffer()
write(io, "qoif")
write(io, hton.(UInt32.([10, 0])))
write(io, hton.(UInt8.([3, 1])))
@test_throws QOI.QOIException("invalid height in header, got 0") QOI.qoi_decode_raw(take!(io))

# Invalid channels
io = IOBuffer()
write(io, "qoif")
write(io, hton.(UInt32.([10, 5])))
write(io, hton.(UInt8.([5, 1])))
@test_throws QOI.QOIException("invalid channels in header, got 5") QOI.qoi_decode_raw(take!(io))

# Invalid colorspace
io = IOBuffer()
write(io, "qoif")
write(io, hton.(UInt32.([10, 5])))
write(io, hton.(UInt8.([3, 2])))
@test_throws QOI.QOIException("invalid colorspace in header, got 2") QOI.qoi_decode_raw(take!(io))

# Too little data after header
io = IOBuffer()
write(io, "qoif")
write(io, hton.(UInt32.([10, 5])))
write(io, hton.(UInt8.([3, 1])))
@test_throws QOI.QOIException("unexpected end of file") QOI.qoi_decode_raw(take!(io))

# Invalid magic bytes
io = IOBuffer()
write(io, "qoiz")
write(io, hton.(UInt32.([10, 5])))
write(io, hton.(UInt8.([3, 1])))
@test_throws QOI.QOIException("invalid magic bytes, got 0x716f697a, expected 0x716f6966") QOI.qoi_decode_raw(take!(io))


const f_logo = joinpath(@__DIR__, "qoi_logo.qoi")

# Test IO API
data_io = open(QOI.qoi_decode, f_logo)
data_file = QOI.qoi_decode(f_logo)
@test data_io == data_file

# Test FileIO integration (enable when https://github.com/JuliaIO/FileIO.jl/pull/354 is merged)
@test data_file == FileIO.load(f_logo)
t = tempname() * ".qoi"
FileIO.save(t, data_file)
@test read(f_logo) == read(t)
