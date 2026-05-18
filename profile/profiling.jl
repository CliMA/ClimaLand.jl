include("sim.jl")

# Profiling
import Profile
import ProfileCanvas

duration_days = 90

sim = setup_simulation(; site_ID = "US-MOz", duration_days)
solve!(sim)
GC.gc()

# No warm up is done
function run_profiler(
    setup_simulation::Function,
    outdir::AbstractString,
)
    simulation = setup_simulation()
    Profile.@profile solve!(simulation)
    results = Profile.fetch()
    flame_file = joinpath(outdir, "flame_cpu.html")
    ProfileCanvas.html_file(flame_file, results)
    @info "Saved compute flame to $flame_file"

    simulation = setup_simulation()
    GC.gc()
    Profile.Allocs.@profile sample_rate = 0.0001 solve!(simulation)
    results = Profile.Allocs.fetch()
    profile = ProfileCanvas.view_allocs(results)
    alloc_flame_file = joinpath(outdir, "alloc_flame_cpu.html")
    ProfileCanvas.html_file(alloc_flame_file, profile)
    @info "Saved allocation flame to $alloc_flame_file"
    return
end

if ClimaComms.device() isa ClimaComms.CPUSingleThreaded
    run_profiler(
        () -> setup_simulation(; site_ID = "US-MOz", duration_days),
        "profile",
    )
end

true

