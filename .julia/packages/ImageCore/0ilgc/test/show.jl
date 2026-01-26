using ImageCore, Colors, FixedPointNumbers, OffsetArrays, Test

sumsz(img) = Base.dims2string(size(img)) * ' '

const rrstr = "reshape, "
rrdim(n) = n-1

# N0f8 is shown as either N0f8 or Normed{UInt8, 8}
# RGB is shown as ColorTypes.RGB or RGB
function typestring(::Type{T}) where T
    buf = IOBuffer()
    show(buf, T)
    String(take!(buf))
end
N0f8_str = typestring(N0f8)
N0f16_str = typestring(N0f16)
RGB_str = typestring(RGB)

@testset "show" begin
    rgb32 = rand(RGB{Float32}, 3, 5)
    v = view(rgb32, 2:3, :)
    @test summary(v) == "2×5 view(::Matrix{ColorTypes.RGB{Float32}}, 2:3, :) with eltype ColorTypes.RGB{Float32}"
    a = channelview(rgb32)
    @test summary(a) == "3×3×5 reinterpret(reshape, Float32, ::Matrix{ColorTypes.RGB{Float32}}) with eltype Float32"
    num64 = rand(3,5)
    b = colorview(RGB, num64)
    str = summary(b)
    @test occursin("5-element", str) && occursin("reinterpret", str) && occursin("reshape", str) && occursin("RGB{Float64}", str) &&
          occursin("::$(typeof(num64))", str) && occursin("with eltype", str)
    rgb8 = rand(RGB{N0f8}, 3, 5)
    c = rawview(channelview(rgb8))
    str = summary(c)
    @test_broken occursin("3×3×5 rawview(reinterpret($(rrstr)", str) && occursin("N0f8", str) &&
          occursin("reshape, N0f8, ::Matrix{RGB{N0f8}}", str) && occursin("with eltype UInt8", str)
    @test_broken summary(rgb8) == "3×5 Matrix{RGB{N0f8}}"
    rand8 = rand(UInt8, 3, 5)
    d = normedview(PermutedDimsArray(rand8, (2,1)))
    @test summary(d) == "5×3 normedview(N0f8, PermutedDimsArray(::$(typeof(rand8)), (2, 1))) with eltype $(N0f8_str)"
    e = PermutedDimsArray(normedview(rand8), (2,1))
    str = summary(e)
    @test occursin("5×3 PermutedDimsArray(reinterpret", str) && occursin("N0f8", str) &&
          occursin("::$(typeof(rand8))), (2, 1)", str) && occursin("with eltype $(N0f8_str)", str)
    rand16 = rand(UInt16, 3, 5)
    f = PermutedDimsArray(normedview(N0f16, rand16), (2,1))
    str = summary(f)
    @test occursin("5×3 PermutedDimsArray(reinterpret", str) && occursin("N0f16", str) &&
          occursin("::$(typeof(rand16))), (2, 1)", str) && occursin("with eltype $(N0f16_str)", str)
    g = channelview(rgb8)
    str = summary(g)
    @test_broken occursin("3×3×5 reinterpret($(rrstr)", str) && occursin("N0f8", str) &&
          occursin("::Array{RGB{N0f8},$(rrdim(3))}", str)
    @test occursin("with eltype", str)
    h = OffsetArray(rgb8, -1:1, -2:2)
    @test_broken summary(h) == "$(sumsz(h))OffsetArray(::Array{RGB{N0f8},2}, -1:1, -2:2) with eltype $(RGB_str){$(N0f8_str)} with indices -1:1×-2:2"
    i = channelview(h)
    str = summary(i)
    @test_broken occursin("$(sumsz(i))reinterpret($(rrstr)", str) && occursin("N0f8", str) &&
          occursin("OffsetArray(::Array{RGB{N0f8},$(rrdim(3))}", str) && occursin("-1:1, -2:2", str) && occursin("with indices", str)
    c = channelview(rand(RGB{N0f8}, 2))
    o = OffsetArray(c, -1:1, 0:1)
    str = summary(o)
    @test_broken occursin("$(sumsz(o))OffsetArray(reinterpret($(rrstr)", str) && occursin("N0f8,", str) &&
          occursin("::Array{RGB{N0f8},$(rrdim(2))}", str) && occursin("-1:1, 0:1", str) && occursin("with eltype $(N0f8_str)", str) &&
          occursin("with indices", str)
    # Issue #45
    a = collect(tuple())
    @test summary(a) == "0-element $(typeof(a))"
    b = view(a, :)
    @test summary(b) == "0-element view(::$(typeof(a)), :) with eltype Union{}"
end

nothing
