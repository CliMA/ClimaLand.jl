using Logging
using Glob
using ImageCore
using ImageMagick
using IndirectArrays
using Test
using PNGFiles
import PNGFiles._inspect_png_read

const DEBUG_FILE = joinpath(@__DIR__, "_debugfile.jmd")

logger = ConsoleLogger(stdout, Logging.Info)
global_logger(logger)


PNG_TEST_PATH = joinpath(@__DIR__, "temp")
isdir(PNG_TEST_PATH) && rm(PNG_TEST_PATH, recursive = true)
mkdir(PNG_TEST_PATH)

include("utils.jl")

open(DEBUG_FILE, "w") do f
    write(f, "```julia\n")
    write(f, """include(\"$(joinpath(@__DIR__, "utils.jl"))\")\n """)
    write(f, "```\n")
end
function _add_debugging_entry(fpath, case, imdiff_val=missing)
    apath = joinpath(@__DIR__, abspath(fpath))
    open(DEBUG_FILE, "a") do f
        write(f, "\n\n")
        write(f, "# $(case) ############################################################\n")
        write(f, "```julia\n")
        write(f, "p = _standardize_grayness(load(\"$(apath)\"))\n")
        write(f, "i = _standardize_grayness(ImageMagick.load(\"$(apath)\"))\n")
        write(f, "plotdiffs(p, i)\n")
        write(f, "```\n")
        write(f, "```julia\n")
        write(f, "imdiff(i, p) # $(imdiff_val)\n")
        write(f, "_inspect_png_read(\"$(apath)\");\n")
        write(f, "```\n")
    end
end

function ensure_imagemagick()
    # ImageMagick backend itself might be broken
    # https://github.com/JuliaIO/ImageMagick.jl/issues/206
    tmpfile = joinpath(PNG_TEST_PATH, "tmp.png")
    img = Gray{N0f8}.(reshape(collect(0.0:0.01:0.49), 5, 10))
    ImageMagick.save(tmpfile, img)
    @test ImageMagick.load(tmpfile) == img
end

ensure_imagemagick()
@testset "PNGFiles" begin
    include("test_pngsuite.jl")
    include("test_paletted_images.jl")
    include("test_testimages.jl")
    include("test_synthetic_images.jl")
    include("test_invalid_inputs.jl")
    include("test_images_with_background.jl")
    include("test_io.jl")
    include("test_various_array_types.jl")
    include("test_dpi.jl")
end

# Cleanup
if !is_ci()
    isdir(PNG_TEST_PATH) && rm(PNG_TEST_PATH, recursive = true)
end
