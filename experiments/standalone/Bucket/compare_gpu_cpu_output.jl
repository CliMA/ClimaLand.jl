using DelimitedFiles
using Statistics
import ClimaLand
function check(job)
    outdir = joinpath(
        pkgdir(ClimaLand),
        "experiments/standalone/Bucket/artifacts_$job",
    )
    cpu_state = readdlm(joinpath(outdir, "tf_state_cpu_$job.txt"), ',')
    gpu_state = readdlm(joinpath(outdir, "tf_state_gpu_$job.txt"), ',')
    @show abs(maximum(cpu_state .- gpu_state))
    @show abs(median(cpu_state .- gpu_state))
    @show abs(mean(cpu_state .- gpu_state))
    @assert isapprox(cpu_state, gpu_state)
end

check("function")
check("staticmap")
check("temporalmap")
