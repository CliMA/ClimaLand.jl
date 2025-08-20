# This script contains functions to assess performance, either via Nsight or by running the
# model multiple times and collecting statistics for execution time and allocations
# to make flame graphs. You can choose between the two
# by requestion --profiler nsight or --profiler integrated. If not provided,
# the integrated profiler is used, and if benchmarking on cpu, flamegraphs are created.

import CUDA
import Profile, ProfileCanvas
using Dates
import ClimaLand
using ClimaComms


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
        @info "Benchmarking with Integrated Profiler"
        run_integrated_profiler(device, setup_simulation, outdir)
        @info "Done with integrated profiler"
        @info "Running timing benchmarks"
        (average_timing_s, std_timing_s) =
            run_timing_benchmarks(device, setup_simulation)
        @show (average_timing_s, std_timing_s)
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
        simulation = setup_simulation()
        ClimaLand.Simulations.step!(simulation)
        CUDA.@profile external = true begin
            for i in 1:5
                ClimaLand.Simulations.step!(simulation)
            end
        end
    else
        @error "Unknown profiler: $profiler"
    end
    return

end

function run_integrated_profiler(device, setup_simulation, outdir)
    simulation = setup_simulation()
    Profile.@profile ClimaLand.Simulations.solve!(simulation)
    results = Profile.fetch()
    flame_file = joinpath(outdir, "flame_$device_suffix.html")
    ProfileCanvas.html_file(flame_file, results)
    @info "Saved compute flame to $flame_file"
    return
end

function run_integrated_profiler(
    device::ClimaComms.CUDADevice,
    setup_simulation,
)
    simulation = setup_simulation()
    p = CUDA.@profile ClimaLand.Simulations.solve!(simulation)
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
    return
end

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
