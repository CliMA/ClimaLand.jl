module JpegTurbo

using ImageCore
using TOML
const project_info = TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))

include("../lib/LibJpeg.jl")
using .LibJpeg
include("libjpeg_utils.jl")

include("common.jl")
include("encode.jl")
include("decode.jl")

include("fileio.jl")

export jpeg_encode, jpeg_decode

end
