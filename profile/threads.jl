import ClimaComms
ClimaComms.@import_required_backends
using Dates
using ClimaLand
import ClimaLand.FluxnetSimulations as FluxnetSimulations
using ClimaLand.Simulations: solve!
import CUDA

include(joinpath(@__DIR__, "sim.jl"))

site_IDs = ["US-MOz", "US-NR1", "US-Ha1"]
duration_days = 90

# Create one stream per site
streams = [CUDA.CuStream() for _ in site_IDs]

# Launch all solves concurrently on separate streams
@info "Launching concurrent solves..."
t_start = time()
@sync for (stream, site_ID) in zip(streams, site_IDs)
    Threads.@spawn begin
        @info "START" site_ID thread=Threads.threadid() t=time() - t_start
        sim = setup_simulation(; site_ID, duration_days)
        CUDA.stream!(stream) do
            solve!(sim)
        end
        @info "END" site_ID thread=Threads.threadid() t=time() - t_start
    end
end
t_elapsed = time() - t_start
@info "All done" t_elapsed
