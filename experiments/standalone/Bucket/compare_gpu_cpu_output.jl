using DelimitedFiles
using Statistics
import ClimaLand

function check(job)
    cpu_state = readdlm(
        joinpath(
            pkgdir(ClimaLand),
            "experiments/standalone/Bucket/artifacts_$(job)_cpu/output_0000/tf_state_cpu_$job.txt",
        ),
        ',',
    )
    gpu_state = readdlm(
        joinpath(
            pkgdir(ClimaLand),
            "experiments/standalone/Bucket/artifacts_$(job)_gpu/output_0000/tf_state_gpu_$job.txt",
        ),
        ',',
    )

    # gpu_state = readdlm(joinpath(outdir, "tf_state_gpu_$job.txt"), ',')
    @show abs(maximum(cpu_state .- gpu_state))
    @show abs(median(cpu_state .- gpu_state))
    @show abs(mean(cpu_state .- gpu_state))
    @assert isapprox(cpu_state, gpu_state)
end

check("function")
check("era5")
check("temporalmap")
