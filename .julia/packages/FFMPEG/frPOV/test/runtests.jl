using FFMPEG
using FFMPEG_jll
using Test

text_execute(f) = try
    f()
    return true
catch e
    @warn "can't execute" exception=e
    return false
end

@testset "FFMPEG.jl" begin
    FFMPEG.versioninfo()
    
    # Test run and parse output
    out = FFMPEG.exe(`-version`, collect=true)
    @test occursin("ffmpeg version ",out[1])
    
    out = FFMPEG.exe(`-version`, command=FFMPEG.ffprobe, collect=true)
    @test occursin("ffprobe version ",out[1])
    
    # Test different invokation methods
    @test text_execute(() -> FFMPEG.exe("-version"))
    @test text_execute(() -> FFMPEG.exe(`-version`))
    @test text_execute(() -> FFMPEG.exe(`-version`, collect=true))
    @test text_execute(() -> FFMPEG.ffmpeg_exe(`-version`))
    @test text_execute(() -> FFMPEG.ffprobe_exe(`-version`))
    @test text_execute(() -> ffmpeg`-version`)
    @test text_execute(() -> ffprobe`-version`)
    @test text_execute(() -> @ffmpeg_env run(`$(FFMPEG_jll.ffmpeg_path) -version`))
end
