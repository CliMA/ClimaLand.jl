# FFMPEG

[![Codecov](https://codecov.io/gh/JuliaIO/FFMPEG.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaIO/FFMPEG.jl)

This package simply offers:

```julia
# a simple way to invoke ffmpeg:
FFMPEG.exe("-version")
FFMPEG.exe("-version", collect=true) #collect output lines into an array of strings
FFMPEG.exe("-version", command=FFMPEG.ffprobe, collect=true) #collect ffprobe output lines into an array of strings (defaults to ffmpeg)

@ffmpeg_env run(`$(FFMPEG.ffmpeg) -version`) #Manually sets up the shared lib environment location. Note the $(FFMPEG.ffmpeg)

ffmpeg_exe("-version") #takes strings
ffmpeg_exe(`-version`) #or command strings

ffprobe_exe("-version") # we wrap FFPROBE too!

ffmpeg`-version` # Cmd string macros too
ffprobe`-version`

# the AV libraries (exported too):
FFMPEG.libavcodec
FFMPEG.libavformat
FFMPEG.libavutil
FFMPEG.libswscale
FFMPEG.libavdevice
FFMPEG.libavfilter

# and for good measure:
FFMPEG.versioninfo()
```

For a high level API to the AV libraries in `libav`, have a look at [VideoIO.jl](https://github.com/JuliaIO/VideoIO.jl/).

This package is made to be included into packages that just need the ffmpeg binaries + executables, and don't want to take on the 3.6 second load time of VideoIO.
