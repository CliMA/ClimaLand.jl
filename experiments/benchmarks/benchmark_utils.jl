import CUDA
using Test
using Statistics
import Profile, ProfileCanvas
using Dates
import ClimaLand
using ClimaComms
using CSV

"""
    profile_and_benchmark(
        setup_simulation::Function,
        device::ClimaComms.AbstractDevice,
        reference_time::Number,
        outdir::AbstractString,
    )

Profile and benchmark the simulation returned by `setup_simulation()`. Also test that
wall time it takes to run the simulation is significantly greater than `reference_time`.
If profiling produces any outputs, they are saved in `outdir`
"""
function profile_and_benchmark(
    setup_simulation::Function,
    device::ClimaComms.AbstractDevice,
    reference_time::Number,
    outdir::AbstractString,
)
    # setup and run for a bit just to ensure everything is compiled
    setup_time = ClimaComms.@elapsed device setup_simulation()
    simulation = setup_simulation()
    first_step_time =
        ClimaComms.@elapsed device ClimaLand.Simulations.step!(simulation)
    model_name = ClimaLand.name(simulation.model)
    @info "$model_name setup took $setup_time (s)"
    @info "First step time: $first_step_time (s)"
    step_or_solve!(simulation, device)
    simulation = nothing
    GC.gc()
    # only run timing benchmarks if no external CUDA profiler is active
    use_external_profiler = CUDA.Profile.detect_cupti()
    if !use_external_profiler
        (average_timing_s, std_timing_s) =
            run_timing_benchmarks(setup_simulation, device)
        GC.gc()
        # if the profiler is ran before the benchmarks, the timing seems less consistent
        run_profiler(setup_simulation, device, outdir)
        if get(ENV, "BUILDKITE_PIPELINE_SLUG", nothing) == "climaland-benchmark"
            if average_timing_s > reference_time + std_timing_s
                @info "Possible performance regression, previous average time was $(reference_time)"
                extra_note_string = "If this is unexpected, it may be caused by someone "
                extra_note_string *= "using a gpu on the clima cluster without a slurm reservation. "
                extra_note_string *= "Please use nvtop and squeue to check"
                @info extra_note_string
            elseif average_timing_s < reference_time - std_timing_s
                device_suffix =
                    device isa ClimaComms.AbstractCPUDevice ? "CPU" : "GPU"
                @info "Possible significant performance improvement, please update PREVIOUS_$(device_suffix)_TIME in $(@__DIR__)"
            end
            # cpu benchmarks seem to have more variation because less samples are taken
            N_ALLOWED_STDV = device isa ClimaComms.CUDADevice ? 2 : 3
            @testset "Performance" begin
                @test reference_time - N_ALLOWED_STDV * std_timing_s <=
                      average_timing_s <=
                      reference_time + N_ALLOWED_STDV * std_timing_s
            end
        end
    else
        run_profiler(setup_simulation, device, outdir)
    end
    return
end

"""
    run_profiler(
        setup_simulation::Function,
        device::ClimaComms.CPUSingleThreaded,
        outdir::AbstractString,
    )

Use the Profile package to profile 5 steps of the simulation returned by `setup_simulation`.
The profiling results are used to produce compute and allocation flames which
are saved in `outdir` as flame_cpu.html and alloc_flame_cpu.html.
"""
function run_profiler(
    setup_simulation::Function,
    device::ClimaComms.CPUSingleThreaded,
    outdir::AbstractString,
)

    simulation = setup_simulation()
    Profile.@profile step_or_solve!(simulation, device)
    results = Profile.fetch()
    flame_file = joinpath(outdir, "flame_cpu.html")
    ProfileCanvas.html_file(flame_file, results)
    @info "Saved compute flame to $flame_file"

    simulation = setup_simulation()
    GC.gc()
    Profile.Allocs.@profile sample_rate = 0.0025 step_or_solve!(simulation, device)
    results = Profile.Allocs.fetch()
    profile = ProfileCanvas.view_allocs(results)
    alloc_flame_file = joinpath(outdir, "alloc_flame_cpu.html")
    ProfileCanvas.html_file(alloc_flame_file, profile)
    @info "Saved allocation flame to $alloc_flame_file"
    return
end

"""
    run_profiler(
        setup_simulation::Function,
        device::ClimaComms.CUDADevice,
        outdir::AbstractString;
        use_external_profiler = CUDA.Profile.detect_cupti(),
    )

Profile 3 steps of the simulation returned by `setup_simulation`. If `use_external_profiler`
is true, the output should be handeled by the external profiler, Otherwise, the following are
saved in `outdir`:
    host.csv - walltime spent in each CUDA API call
    device.csv - walltime spent in each CUDA kernel
    nvtx.csv - walltime spent in each NVTX range
"""
function run_profiler(
    setup_simulation::Function,
    device::ClimaComms.CUDADevice,
    outdir::AbstractString;
    use_external_profiler = CUDA.Profile.detect_cupti(),
)

    simulation = setup_simulation()
    # step once for closures that may recompile
    ClimaLand.Simulations.step!(simulation)
    if use_external_profiler
        @info "Profiling using external profiler"
        CUDA.@profile external = true step_N_times!(simulation, 3)
        @info "Profiling complete"
    else
        @info "Profiling using internal profiler"
        p = CUDA.@profile step_N_times!(simulation, 3)
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
    end
    # solve to close file readers
    step_or_solve!(simulation, device)
    return
end

"""
    run_timing_benchmarks(
        setup_simulation::Function,
        device::ClimaComms.AbstractDevice;
        MAX_PROFILING_TIME_SECONDS::Number = 400,
        MAX_PROFILING_SAMPLES::Number = 100,
    )

Benchmark solving the simulation returned by `setup_simulation()`, with
`MAX_PROFILING_SAMPLES` samples. If the total benchmark time takes longer than
`MAX_PROFILING_SAMPLES`, the benchmark will terminate early.
"""
function run_timing_benchmarks(
    setup_simulation::Function,
    device::ClimaComms.AbstractDevice;
    MAX_PROFILING_TIME_SECONDS::Number = 400,
    MAX_PROFILING_SAMPLES::Number = 100,
)
    time_now = time()
    timings_s = Float64[]
    @info "Running benchmarks"
    while (time() - time_now) < MAX_PROFILING_TIME_SECONDS &&
        length(timings_s) < MAX_PROFILING_SAMPLES
        GC.gc()
        simulation = setup_simulation()
        # step once for closures that may recompile
        ClimaLand.Simulations.step!(simulation)
        push!(
            timings_s,
            ClimaComms.@elapsed device step_or_solve!(simulation, device)
        )
    end
    length(timings_s) == 1 && error("Only one sample was obtained. Try increasing MAX_PROFILING_TIME_SECONDS")
    popfirst!(timings_s) # the first sample is always slower, but I'm not sure why
    num_samples = length(timings_s)
    average_timing_s = round(sum(timings_s) / num_samples, sigdigits = 3)
    max_timing_s = round(maximum(timings_s), sigdigits = 3)
    min_timing_s = round(minimum(timings_s), sigdigits = 3)
    std_timing_s = round(
        sqrt(sum(((timings_s .- average_timing_s) .^ 2) / num_samples)),
        sigdigits = 3,
    )
    @info "Benchmarks complete"
    @info "Num samples: $num_samples"
    @info "Average time: $(average_timing_s) s"
    @info "Max time: $(max_timing_s) s"
    @info "Min time: $(min_timing_s) s"
    @info "Standard deviation time: $(std_timing_s) s"
    return (average_timing_s, std_timing_s)
end

"""
    step_N_times!(sim::ClimaLand.Simulations.LandSimulation, N::Number)

Steps `sim` `N` times.
"""
function step_N_times!(sim::ClimaLand.Simulations.LandSimulation, N::Number)
    for i in 1:N
        ClimaLand.Simulations.step!(sim)
    end
    return
end

"""
    step_or_solve!(sim::ClimaLand.Simulations.LandSimulation, device)

If device is a ClimaComms.CUDADevice, solve `sim`. Otherwise step `sim` 5 times and close
all output writers and file readers.
"""
step_or_solve!(sim::ClimaLand.Simulations.LandSimulation, device::ClimaComms.CUDADevice) =
    ClimaLand.Simulations.solve!(sim)
function step_or_solve!(sim::ClimaLand.Simulations.LandSimulation, device)
    step_N_times!(sim, 5)
    ClimaLand.Simulations.close_output_writers(sim.diagnostics)
    ClimaLand.Simulations.close_all_ncfiles()
    return
end
