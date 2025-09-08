# This script contains functions to assess performance, either via Nsight or by running the
# model multiple times withj the integrated profiler and collecting statistics

# You can choose between the profiling type
# by requestion --profiler nsight or --profiler integrated. If not provided,
# the integrated profiler is used. If benchmarking on cpu with the integrated profiler,
# flamegraphs are created and saved. If benchmarking on gpu with the integrated profiler, the results
# are saved to csv files. If nsight is used, the simulation is profiled for two steps.



import CUDA
using ArgParse
using Test
using Statistics
import Profile, ProfileCanvas
using Dates
import ClimaLand
using ClimaComms
using CSV

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--profiler"
        help = "Profiler option: nsight or flamegraph"
        arg_type = String
        default = "integrated"
    end
    return parse_args(s)
end

function run_benchmarks(
    device,
    setup_simulation,
    profiler,
    PREVIOUS_GPU_TIME,
    outdir,
)
    # run once just to ensure everything is compiled
    simulation = setup_simulation()
    ClimaLand.Simulations.solve!(simulation)
    if profiler == "integrated"
        @info "Running timing benchmarks"
        (average_timing_s, std_timing_s) =
            run_timing_benchmarks(device, setup_simulation)
        @info "Done with timing benchmarks. Starting integrated profiling"
        run_integrated_profiler(device, setup_simulation, outdir)
        if get(ENV, "BUILDKITE_PIPELINE_SLUG", nothing) ==
           "climaland-benchmark" && device isa ClimaComms.CUDADevice
            if average_timing_s > PREVIOUS_GPU_TIME + std_timing_s
                @info "Possible performance regression, previous average time was $(PREVIOUS_GPU_TIME)"
            elseif average_timing_s < PREVIOUS_GPU_TIME - std_timing_s
                @info "Possible significant performance improvement, please update PREVIOUS_GPU_TIME in $(@__DIR__)"
            end
            @testset "Performance" begin
                @test PREVIOUS_GPU_TIME - 2std_timing_s <=
                      average_timing_s <=
                      PREVIOUS_GPU_TIME + 2std_timing_s
            end
        end
    elseif profiler == "nsight"
        run_external_profiler(setup_simulation)
    else
        @error "Unknown profiler: $profiler"
    end
    return

end

# integrated profiling for cpu - flame graphs are hard to interpret when gpu is used
function run_integrated_profiler(device, setup_simulation, outdir)
    simulation = setup_simulation()
    Profile.@profile ClimaLand.Simulations.solve!(simulation)
    results = Profile.fetch()
    flame_file = joinpath(outdir, "flame_$device_suffix.html")
    ProfileCanvas.html_file(flame_file, results)
    @info "Saved compute flame to $flame_file"

    simulation = setup_simulation()
    Profile.Allocs.@profile sample_rate = 0.0025 ClimaLand.Simulations.solve!(
        simulation,
    )
    results = Profile.Allocs.fetch()
    profile = ProfileCanvas.view_allocs(results)
    alloc_flame_file = joinpath(outdir, "alloc_flame_$device_suffix.html")
    ProfileCanvas.html_file(alloc_flame_file, profile)
    @info "Saved allocation flame to $alloc_flame_file"
    return
end

# integrated profiling for CUDA
function run_integrated_profiler(
    device::ClimaComms.CUDADevice,
    setup_simulation,
    outdir,
)
    simulation = setup_simulation()
    p = CUDA.@profile begin
        for i in 1:20
            ClimaLand.Simulations.step!(simulation)
        end
    end
    envs = ("COLUMNS" => 120,)
    withenv(envs...) do
        io = IOContext(
            stdout,
            :crop => :horizontal,
            :limit => true,
            :displaysize => displaysize(),
        )
        show(io, p)
    end
    println()
    host_file = joinpath(outdir, "host.csv")
    device_file = joinpath(outdir, "device.csv")
    nvtx_file = joinpath(outdir, "nvtx.csv")
    CSV.write(host_file, p.host)
    CSV.write(device_file, p.device)
    CSV.write(nvtx_file, p.nvtx)
    @info "CUDA.jl integrated profiling results saved in $outdir"
    return
end

# time `solve!(simulation)` - simulation setup is not timed
function run_timing_benchmarks(
    device,
    setup_simulation;
    MAX_PROFILING_TIME_SECONDS = 400,
)
    MAX_PROFILING_SAMPLES = 100
    time_now = time()
    timings_s = Float64[]
    while (time() - time_now) < MAX_PROFILING_TIME_SECONDS &&
        length(timings_s) < MAX_PROFILING_SAMPLES
        simulation = setup_simulation()
        push!(
            timings_s,
            ClimaComms.@elapsed device ClimaLand.Simulations.solve!(simulation)
        )
    end
    num_samples = length(timings_s)
    average_timing_s = round(sum(timings_s) / num_samples, sigdigits = 3)
    max_timing_s = round(maximum(timings_s), sigdigits = 3)
    min_timing_s = round(minimum(timings_s), sigdigits = 3)
    std_timing_s = round(
        sqrt(sum(((timings_s .- average_timing_s) .^ 2) / num_samples)),
        sigdigits = 3,
    )
    @info "Num samples: $num_samples"
    @info "Average time: $(average_timing_s) s"
    @info "Max time: $(max_timing_s) s"
    @info "Min time: $(min_timing_s) s"
    @info "Standard deviation time: $(std_timing_s) s"
    @info "Done profiling"
    return (average_timing_s, std_timing_s)
end

function run_external_profiler(setup_simulation)
    simulation = setup_simulation()
    # step to ensure compilation
    ClimaLand.Simulations.step!(simulation)
    CUDA.@profile external = true begin
        for i in 1:3
            ClimaLand.Simulations.step!(simulation)
        end
    end
    return
end
