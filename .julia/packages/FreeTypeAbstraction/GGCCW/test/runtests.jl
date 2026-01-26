ENV["FREETYPE_ABSTRACTION_FONT_PATH"] = @__DIR__  # coverage

using FreeTypeAbstraction, Colors, ColorVectorSpace, GeometryBasics
using GeometryBasics: Vec2f
import FreeTypeAbstraction
const FA = FreeTypeAbstraction
using FreeType
using Test

@testset "init and done" begin
    @test_throws ErrorException FA.ft_init()
    @test FA.ft_done()
    @test_throws ErrorException FA.ft_done()
    @test FA.ft_init()
end

face = FA.findfont("hack")

@testset "basics" begin
    @test :size in propertynames(face)
    @test repr(face) == "FTFont (family = Hack, style = Regular)"
    @test Broadcast.broadcastable(face)[] == face

    @test FA.ascender(face) isa Real
    @test FA.descender(face) isa Real

    FA.set_pixelsize(face, 64) # should be the default
    img, extent = FA.renderface(face, 'C', 64)
    @test size(img) == (30, 49)
    @test typeof(img) == Array{UInt8,2}

    @test FA.hadvance(extent) == 39
    @test FA.vadvance(extent) == 62
    @test FA.inkheight(extent) == 49
    @test FA.hbearing_ori_to_left(extent) == 4
    @test FA.hbearing_ori_to_top(extent) == 48
    @test FA.leftinkbound(extent) == 4
    @test FA.rightinkbound(extent) == 34
    @test FA.bottominkbound(extent) == -1
    @test FA.topinkbound(extent) == 48
    @test FA.inkboundingbox(extent) == HyperRectangle{2, Float32}(Float32[4.0, -1.0], Float32[30.0, 49.0])
    @test_broken FA.height_insensitive_boundingbox(extent, face) == HyperRectangle{2, Float32}(Float32[4.0, 64 * -0.23583984], Float32[30.0, 64 * 1.2006836])

    a = renderstring!(zeros(UInt8, 20, 100), "helgo", face, 10, 10, 10)

    @test any(a[3:12, :] .!= 0)
    @test all(a[vcat(1:2, 13:20), :] .== 0)
    @test any(a[:, 11:40] .!= 0)
    @test all(a[:, vcat(1:10, 41:100)] .== 0)
    a = renderstring!(zeros(UInt8, 20, 100), "helgo", face, 10, 15, 70)
    @test any(a[8:17, :] .!= 0)
    @test all(a[vcat(1:7, 18:20), :] .== 0)
    @test any(a[:, 71:100] .!= 0)
    @test all(a[:, 1:70] .== 0)

    a = renderstring!(zeros(Float32, 20, 100), "helgo", face, 10, 10, 50)
    @test maximum(a) <= 1.0
    a = renderstring!(zeros(Float64, 20, 100), "helgo", face, 10, 10, 50)
    @test maximum(a) <= 1.0

    renderstring!(zeros(UInt8, 20, 100), "helgo", face, 10, 25, 80)
    @test_logs (:warn, "using tuple for pixelsize is deprecated, please use one integer") renderstring!(zeros(UInt8, 20, 100), "helgo", face, (10, 10), 1, 1)
end

@testset "ways to access glyphs" begin
    i = FA.glyph_index(face, 'A')
    @test FA.glyph_index(face, i) == i
    @test FA.glyph_index(face, "A") == i
end

@testset "alignements" begin
    a = renderstring!(
        zeros(UInt8, 20, 100),
        "helgo",
        face,
        10,
        10,
        50,
        valign = :vtop,
    )

    @test all(a[1:10, :] .== 0)
    @test any(a[11:20, :] .!= 0)
    a = renderstring!(
        zeros(UInt8, 20, 100),
        "helgo",
        face,
        10,
        10,
        50,
        valign = :vcenter,
    )
    @test all(a[vcat(1:5, 16:end), :] .== 0)
    @test any(a[6:15, :] .!= 0)
    a = renderstring!(
        zeros(UInt8, 20, 100),
        "helgo",
        face,
        10,
        10,
        50,
        valign = :vbaseline,
    )
    @test all(a[vcat(1:2, 13:end), :] .== 0)
    @test any(a[3:12, :] .!= 0)
    a = renderstring!(
        zeros(UInt8, 20, 100),
        "helgo",
        face,
        10,
        10,
        50,
        valign = :vbottom,
    )
    @test any(a[1:10, :] .!= 0)
    @test all(a[11:20, :] .== 0)
    a = renderstring!(
        zeros(UInt8, 20, 100),
        "helgo",
        face,
        10,
        10,
        50,
        halign = :hleft,
    )
    @test all(a[:, 1:50] .== 0)
    @test any(a[:, 51:100] .!= 0)
    a = renderstring!(
        zeros(UInt8, 20, 100),
        "helgo",
        face,
        10,
        10,
        50,
        halign = :hcenter,
    )
    @test all(a[:, vcat(1:35, 66:end)] .== 0)
    @test any(a[:, 36:65] .!= 0)
    a = renderstring!(
        zeros(UInt8, 20, 100),
        "helgo",
        face,
        10,
        10,
        50,
        halign = :hright,
    )
    @test any(a[:, 1:50] .!= 0)
    @test all(a[:, 51:100] .== 0)

    renderstring!(
        zeros(UInt8, 20, 100),
        "helgo",
        face,
        10,
        1,
        1,
        halign = :hcenter,
        valign = :vcenter,
    )
end

@testset "foreground / background colors" begin
    a = renderstring!(
        zeros(UInt8, 20, 100),
        "helgo",
        face,
        10,
        10,
        50,
        fcolor = 0x80,
    )
    @test maximum(a) <= 0x80
    a = renderstring!(
        zeros(UInt8, 20, 100),
        "helgo",
        face,
        10,
        10,
        50,
        fcolor = 0x80,
        bcolor = 0x40,
    )
    @test any(a .== 0x40)
    a = renderstring!(
        fill(0x01, 20, 100),
        "helgo",
        face,
        10,
        10,
        50,
        bcolor = nothing,
    )
    @test !any(a .== 0x00)
end

@testset "array of grays" begin
    renderstring!(zeros(Gray, 20, 100), "helgo", face, 10, 10, 50)
    renderstring!(
        zeros(Gray{Float64}, 20, 100),
        "helgo",
        face,
        10,
        10,
        50,
        fcolor = Gray(0.5),
    )
    @test true
end

@testset "per char background / colors" begin
    for str in ("helgo", collect("helgo"))
        fcolor = [RGB{Float32}(rand(3)...) for _ ∈ 1:length(str)]
        renderstring!(zeros(RGB{Float32}, 20, 100), str, face, 10, 1, 1; fcolor = fcolor)
        gcolor = [RGB{Float32}(rand(3)...) for _ ∈ 1:length(str)]
        gstr = fill('█', length(str))
        renderstring!(
            zeros(RGB{Float32}, 20, 100), str, face, 10, 1, 1;
            fcolor = fcolor, gcolor = gcolor, gstr = gstr
        )
    end
    @test true
end

@testset "draw bounding boxes" begin
    renderstring!(
        zeros(RGB{Float32}, 20, 100), "helgo", face, 10, 1, 1;
        gcolor = RGB{Float32}(0., 1., 0.), bbox = RGB{Float32}(1., 0., 0.), bbox_glyph = RGB{Float32}(0., 0., 1.)
    )
    @test true
end

# Find fonts
# these fonts should be available on all platforms:

# debug travis... does it even have fonts?
fontpaths = FA.fontpaths()

isempty(fontpaths) && println("OS doesn't have any font folder")

if Sys.islinux()
    fonts = ["dejavu sans"]
    # apple on gh-actions doesn't seem to have any fonts...
elseif Sys.isapple()
    fonts = []
else # windows have some more fonts installed per default
    fonts = [
        "Times New Roman",
        "Arial",
        "Comic Sans MS",
        "Impact",
        "Tahoma",
        "Trebuchet MS",
        "Verdana",
        "Courier New",
    ]
end

@testset "finding fonts" begin
    valid_fontpaths = filter(x -> x != @__DIR__, FA.valid_fontpaths)
    empty!(FA.valid_fontpaths)
    append!(FA.valid_fontpaths, valid_fontpaths)
    for font in fonts
        @testset "finding $font" begin
            @test findfont(font) !== nothing
        end
    end
    @testset "find in additional dir" begin
        @test findfont("Hack") === nothing
        @test findfont("Hack", additional_fonts = @__DIR__) !== nothing
    end
end

@testset "loading lots of fonts" begin
    for i = 1:10, font in fonts
        @time findfont(font)
    end
    @test true
end

@testset "Font extent" begin
    f1 = FontExtent(Vec2f(1, 2), Vec2f(3, 4), Vec2f(5, 6), Vec2f(7, 8))
    f2 = FA.broadcasted(x -> 2 * x, f1)
    @test f2 == FontExtent(Vec2f(2, 4), Vec2f(6, 8), Vec2f(10, 12), Vec2f(14, 16))
    f3 = FA.broadcasted(*, f1, Vec2f(2, 3))
    @test f3 == FontExtent(Vec2f(2, 4), Vec2f(9, 12), Vec2f(10, 18), Vec2f(14, 24))
end

@testset "Boundingbox" begin
    for glyph in ('a', FA.glyph_index(face, 'a'), "a")
        bb, extent = FA.metrics_bb(glyph, face, 64)
        bb2 = FA.boundingbox(glyph, face, 64)
        @test bb == bb2
        w = GeometryBasics.widths(bb2)
        @test w == FA.glyph_ink_size(glyph, face, 64)
    end
end

@testset "Unicode in path" begin
    mktempdir() do dir
        path = joinpath(dir, "𝓜𝔂 𝔽𝕒𝕟𝕔𝕪 𝐅𝐨𝐧𝐭 🅵🅸🅻🅴.ttf")
        cp(joinpath(@__DIR__, "hack_regular.ttf"), path)
        f = FreeTypeAbstraction.FTFont(path)
        @test repr(f) == "FTFont (family = Hack, style = Regular)"
    end
end

@testset "Thread safety" begin
    mktempdir() do dir
        n = 100
        fontfiles = map(1:n) do i
            p = joinpath(dir, "hack_regular_$i.ttf")
            cp(joinpath(@__DIR__, "hack_regular.ttf"), p)
            p
        end
        Threads.@threads for f in fontfiles
            fo = FreeTypeAbstraction.FTFont(f)
            Threads.@threads for i in 1:100
                g = FreeTypeAbstraction.load_glyph(fo, i)
                g = FreeTypeAbstraction.loadglyph(fo, i, 64)
                FreeTypeAbstraction.renderface(fo, i, 16)
            end
        end
    end
end
