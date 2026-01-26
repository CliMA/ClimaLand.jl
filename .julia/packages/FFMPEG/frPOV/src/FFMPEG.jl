module FFMPEG

using FFMPEG_jll

av_version(v) = VersionNumber(v >> 16, (v >> 8) & 0xff, v & 0xff)

_avcodec_version()    = av_version(ccall((:avcodec_version,    libavcodec),    UInt32, ()))
_avformat_version()   = av_version(ccall((:avformat_version,   libavformat),   UInt32, ()))
_avutil_version()     = av_version(ccall((:avutil_version,     libavutil),     UInt32, ()))
_swscale_version()    = av_version(ccall((:swscale_version,    libswscale),    UInt32, ()))
_avdevice_version()   = av_version(ccall((:avdevice_version,   libavdevice),   UInt32, ()))
_avfilter_version()   = av_version(ccall((:avfilter_version,   libavfilter),   UInt32, ()))
_swresample_version() = av_version(ccall((:swresample_version, libswresample), UInt32, ()))

function versioninfo()
    println("Using ffmpeg")
    println("AVCodecs version $(_avcodec_version())")
    println("AVFormat version $(_avformat_version())")
    println("AVUtil version $(_avutil_version())")
    println("SWScale version $(_swscale_version())")
    println("AVDevice version $(_avdevice_version())")
    println("AVFilters version $(_avfilter_version())")
    println("SWResample version $(_swresample_version())")
end

"""
    @ffmpeg_env arg

Runs `arg` within the build environment of FFMPEG.

## Examples

```jldoctest
julia> @ffmpeg_env run(`\$ffmpeg -version`)
ffmpeg version 4.1 Copyright (c) 2000-2018 the FFmpeg developers
built with clang version 6.0.1 (tags/RELEASE_601/final)
[...]
```
"""
macro ffmpeg_env(arg)
    return esc(quote
        ffmpeg() do ffmpeg_path
            $(arg)
        end
    end)
end

"""
    exe(args...)

Execute the given commands as arguments to the given executable.

## Examples

```jldoctest
julia> FFMPEG.exe("-version")
ffmpeg version 4.1 Copyright (c) 2000-2018 the FFmpeg developers
built with clang version 6.0.1 (tags/RELEASE_601/final)
[...]
```
"""
exe(args::AbstractString...; command = FFMPEG.ffmpeg, collect = false) = exe(Cmd([args...]), command=command, collect=collect)

"""
    collectexecoutput(exec::Cmd) -> Array of output lines

Takes the dominant output std from ffmpeg.
"""
function collectexecoutput(exec::Cmd)
    out_s, err_s = readexecoutput(exec)
    return (length(out_s) > length(err_s)) ? out_s : err_s
end

"""
    readexecoutput(exec::Cmd) -> (out, err)

Takes the output stdout and stderr from the input command.  

Returns a Tuple of String vectors.
"""
function readexecoutput(exec::Cmd)
    out = Pipe(); err = Pipe()
    p = Base.open(pipeline(ignorestatus(exec), stdout=out, stderr=err))
    close(out.in); close(err.in)
    err_s = readlines(err); out_s = readlines(out)
    return out_s, err_s
end

"""
    exe(arg)

Execute the given command literal as an argument to the given executable.

## Examples

```jldoctest
julia> FFMPEG.exe(`-version`)
ffmpeg version 4.1 Copyright (c) 2000-2018 the FFmpeg developers
built with clang version 6.0.1 (tags/RELEASE_601/final)
[...]
```
"""
function exe(arg::Cmd; command = ffmpeg, collect = false)
    f = collect ? collectexecoutput : Base.run
    @static if VERSION ≥ v"1.6"
        f(`$(command()) $arg`)
    else
        command() do command_path
            f(`$command_path $arg`)
        end
    end
end

"""
    ffmpeg_exe(arg::Cmd)
    ffmpeg_exe(args::String...)

Execute the given arguments as arguments to the `ffmpeg` executable.
"""
ffmpeg_exe(args...) = exe(args...; command = ffmpeg)

"""
    ffprobe_exe(arg::Cmd)
    ffprobe_exe(args::String...)

Execute the given arguments as arguments to the `ffprobe` executable.
"""
ffprobe_exe(args...) = exe(args...; command = ffprobe)

"""
    ffmpeg\`<ARGS>\`

Execute the given arguments as arguments to the `ffmpeg` executable.
"""
macro ffmpeg_cmd(arg)
    esc(:(ffmpeg_exe($arg)))
end

"""
    ffprobe\`<ARGS>\`

Execute the given arguments as arguments to the `ffprobe` executable.
"""
macro ffprobe_cmd(arg)
    esc(:(ffprobe_exe($arg)))
end

export ffmpeg_exe, @ffmpeg_env, ffprobe_exe, ffmpeg, ffprobe, @ffmpeg_cmd, @ffprobe_cmd, libavcodec, libavformat, libavutil, libswscale, libavfilter, libavdevice

end # module
