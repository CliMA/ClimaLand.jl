using DelimitedFiles
using Statistics
import ClimaLand
outdir = joinpath(pkgdir(ClimaLand), "experiments/standalone/Bucket/artifacts")
cpu_state = readdlm(joinpath(outdir, "tf_state_cpu.txt"), ',')
gpu_state = readdlm(joinpath(outdir, "tf_state_gpu.txt"), ',')
@show abs(maximum(cpu_state .- gpu_state))
@show abs(median(cpu_state .- gpu_state))
@show abs(mean(cpu_state .- gpu_state))
@assert isapprox(cpu_state, gpu_state)
