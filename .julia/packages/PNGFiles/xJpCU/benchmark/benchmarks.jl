# Designed to be used with PkgBenchmark
# i.e.
# using PkgBenchmark
# import PNGFiles
# benchmarkpkg(PNGFiles)

using BenchmarkTools
using PNGFiles
using PNGFiles.ImageCore

include("../test/test_images/synth_images.jl")

SUITE = BenchmarkGroup()

SUITE["synth_imgs"] = BenchmarkGroup()
for (name, img) in synth_imgs
    SUITE["synth_imgs"][name] = BenchmarkGroup()
    SUITE["synth_imgs"][name]["save"] = @benchmarkable PNGFiles.save("tst.png", $img) teardown=(rm("tst.png"))
    SUITE["synth_imgs"][name]["load"] = @benchmarkable PNGFiles.load("tst.png") setup=(PNGFiles.save("tst.png", $img)) teardown=(rm("tst.png"))
end