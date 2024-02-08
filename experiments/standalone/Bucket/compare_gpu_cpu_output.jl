using DelimitedFiles
using Statistics
import ClimaLand
outdir = joinpath(pkgdir(ClimaLand), "experiments/standalone/Bucket/artifacts")
cpu_state = readdlm(joinpath(outdir, "tf_state_cpu.txt"), ',')
gpu_state = readdlm(joinpath(outdir, "tf_state_gpu.txt"), ',')
@show mean(cpu_state .- gpu_state)
@show eps(Float64)
@assert mean(cpu_state .- gpu_state) < eps(Float64) * 5e2
