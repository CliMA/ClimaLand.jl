using Clang.JLLEnvs
using Clang.Generators
using JpegTurbo_jll

cd(@__DIR__)

include_dir = normpath(JpegTurbo_jll.artifact_dir, "include")


function generate_prologue(options)
    # Macro "#define LIBJPEG_TURBO_VERSION 2.1.0" is not well-defined values
    # for Clang to parse. Thus we manually parse and generate it as prologue.
    function get_jpegturbo_version()
        lines = readlines(joinpath(include_dir, "jconfig.h"))
        prefix = "#define LIBJPEG_TURBO_VERSION "
        version_line = filter(lines) do line
            startswith(line, prefix)
        end[1]
        version_str = strip(split(version_line, prefix)[2])
        return VersionNumber(version_str)
    end

    outfile = options["general"]["prologue_file_path"]
    open(outfile, "w") do io
        println(io, "const LIBJPEG_TURBO_VERSION = v\"$(get_jpegturbo_version())\"")
        println(io)

        # https://github.com/libjpeg-turbo/libjpeg-turbo/blob/14ce28a92d45e4e22b643bd845ba6c543ebcd388/win/jconfig.h.in#L14-L18
        println(io, "const boolean = @static Sys.iswindows() ? Cuchar : Cint")
        println(io)
    end
end

options = load_options(joinpath(@__DIR__, "generator.toml"))

args = get_default_args()
push!(args, "-I$include_dir")
push!(args, "--include=stdio.h")

generate_prologue(options)

# https://github.com/JuliaInterop/Clang.jl/discussions/373
headers = [
    joinpath(include_dir, "jpeglib.h"),
    joinpath(include_dir, "turbojpeg.h")
]

ctx = create_context(headers, args, options)

build!(ctx)
