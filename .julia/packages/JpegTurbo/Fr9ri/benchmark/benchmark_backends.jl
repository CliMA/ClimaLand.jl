# compare JpegTurbo with other IO backends that supports JPEG format
using BenchmarkTools
using TestImages
using ImageMagick
using QuartzImageIO
using JpegTurbo
using ImageCore
using FileIO
using ProgressMeter
using ImageQualityIndexes
using Printf
using PyCall
using InteractiveUtils

const tmpname = "test.jpg"
use_quartz_imageio = Sys.isapple()

const cv = pyimport("cv2")
const skimageio = pyimport("skimage.io")
const skimage = pyimport("skimage")

@pyinclude joinpath(@__DIR__, "pyhelpers.py")

function print_cv_version(io)
    cv_jpeg_ver = match(r"\d+\.\d+\.\d+-\d+", filter(split(cv.getBuildInformation(), "\n")) do line
        occursin("build-libjpeg-turbo", line)
    end[1]).match

    println(io, "OpenCV version: $(cv.__version__)")
    println(io, "OpenCV libjpeg-turbo version: $(cv_jpeg_ver)")
end

img_list = Any[
    # Gray
    "moonsurface" => testimage("moonsurface"),
    "cameraman" => testimage("cameraman"),
    "pirate" => testimage("pirate"),
    "house" => Gray.(testimage("house")),
    "rand" => rand(Gray, 512, 512),
    "rand" => rand(Gray, 4096, 4096),
    # RGB
    "fabio" => testimage("fabio_color_512"),
    "barbara" => testimage("barbara_color"),
    "mandril" => testimage("mandril_color"),
    "coffee" => testimage("coffee"),
    "lighthouse" => testimage("lighthouse"),
    "earth_apollo" => RGB.(testimage("earth_apollo")),
    "rand" => rand(RGB, 512, 512),
    "rand" => rand(RGB, 4096, 4096)
]

open("report.md", "w") do out_io
    println(out_io, "# JPEG backends comparison")
    println(out_io)

    println(out_io, "```")
    println(out_io, "Julia versioninfo:")
    versioninfo(out_io)
    println(out_io)

    println(out_io, "JpegTurbo.jl versioninfo:")
    JpegTurbo.versioninfo(out_io)
    println(out_io)

    print_cv_version(out_io)
    println(out_io)

    println(out_io, "Scikit-image version: $(skimage.__version__)")
    println(out_io, "```")
    println(out_io)

    @showprogress for (name, img) in img_list
        io = IOBuffer()
        println(io, "## $name $(eltype(img)) $(size(img))")
        println(io)
        println(io, "| Backend | encode time(ms) | decode time(ms) | encoded size(KB) | PSNR(dB) |")
        println(io, "| ------- | --------------- | --------------- | ---------------- | -------- |")

        # JpegTurbo
        @info "JpegTurbo: $name"
        t_encode = 1000 * @belapsed jpeg_encode(tmpname, $img) seconds=2
        sz = filesize(tmpname) / 1024
        t_decode = 1000 * @belapsed jpeg_decode(tmpname) seconds=2
        psnr = assess_psnr(img, eltype(img).(jpeg_decode(tmpname)))
        row_msg = @sprintf "| JpegTurbo.jl | %.2f | %.2f | %.2f | %.4f |" t_encode t_decode sz psnr
        println(io, row_msg)

        # ImageMagick
        @info "ImageMagick: $name"
        t_encode = 1000 * @belapsed ImageMagick.save(tmpname, $img) seconds=2
        sz = filesize(tmpname) / 1024
        t_decode = 1000 * @belapsed ImageMagick.load(tmpname) seconds=2
        psnr = assess_psnr(img, eltype(img).(ImageMagick.load(tmpname)))
        row_msg = @sprintf "| ImageMagick.jl | %.2f | %.2f | %.2f | %.4f |" t_encode t_decode sz psnr
        println(io, row_msg)

        # QuartzImageIO
        if use_quartz_imageio
            @info "QuartzImageIO: $name"
            t_encode = 1000 * @belapsed QuartzImageIO.save(File{format"JPEG"}(tmpname), $img) seconds=2
            sz = filesize(tmpname) / 1024
            t_decode = 1000 * @belapsed QuartzImageIO.load(File{format"JPEG"}(tmpname)) seconds=2
            psnr = assess_psnr(img, eltype(img).(QuartzImageIO.load(File{format"JPEG"}(tmpname))))
            row_msg = @sprintf "| QuartzImageIO.jl | %.2f | %.2f | %.2f | %.4f |" t_encode t_decode sz psnr
            println(io, row_msg)
        end

        # OpenCV (Python)
        @info "OpenCV (Python): $name"
        save("tmp.png", img)
        flag = eltype(img) <: Gray ? cv.IMREAD_GRAYSCALE : cv.IMREAD_COLOR
        t_encode, t_decode = py"belapsed_cv2"(tmpname, flag, number=100)
        img_cv = cv.imread("tmp.png", flag)
        cv.imwrite(tmpname, img_cv)
        sz = filesize(tmpname) / 1024
        img_cv_new = cv.imread(tmpname, flag)
        if eltype(img) <: RGB
            img_cv_new = collect(reinterpret(reshape, RGB{N0f8}, permutedims(img_cv_new, (3, 1, 2))))
            img_cv = collect(reinterpret(reshape, RGB{N0f8}, permutedims(img_cv, (3, 1, 2))))
        elseif eltype(img) <: Gray
            img_cv_new = collect(reinterpret(reshape, Gray{N0f8}, img_cv_new))
            img_cv = collect(reinterpret(reshape, Gray{N0f8}, img_cv))
        end
        psnr = assess_psnr(img_cv, img_cv_new)
        row_msg = @sprintf "| OpenCV (Python) | %.2f | %.2f | %.2f | %.4f |" t_encode t_decode sz psnr
        println(io, row_msg)

        # Scikit-image
        @info "Scikit-image: $name"
        save("tmp.png", img)
        t_encode, t_decode = py"belapsed_skimage"(tmpname, number=100)
        img_ski = skimageio.imread("tmp.png")
        skimageio.imsave(tmpname, img_ski)
        sz = filesize(tmpname) / 1024
        img_ski_new = skimageio.imread(tmpname)
        if eltype(img) <: RGB
            img_ski_new = collect(reinterpret(reshape, RGB{N0f8}, permutedims(img_ski_new, (3, 1, 2))))
            img_ski = collect(reinterpret(reshape, RGB{N0f8}, permutedims(img_ski, (3, 1, 2))))
        elseif eltype(img) <: Gray
            img_ski_new = collect(reinterpret(reshape, Gray{N0f8}, img_ski_new))
            img_ski = collect(reinterpret(reshape, Gray{N0f8}, img_ski))
        end
        psnr = assess_psnr(img_ski, img_ski_new)
        row_msg = @sprintf "| Scikit-image | %.2f | %.2f | %.2f | %.4f |" t_encode t_decode sz psnr
        println(io, row_msg)

        println(io)
        write(out_io, take!(io))
    end
    isfile(tmpname) && rm(tmpname)
    isfile("tmp.png") && rm("tmp.png")
end
