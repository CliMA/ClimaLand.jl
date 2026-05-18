using Distributed

nworkers = 10
addprocs(nworkers)

@everywhere begin
    include(joinpath(@__DIR__, "sim.jl"))

    using CUDA
    function run_on_stream(site_ID, duration_days)
        stream = CUDA.CuStream()
        CUDA.stream!(stream) do
            sim = setup_simulation(; site_ID, duration_days)
            @info "Starting simulation for worker $(Distributed.myid())"
            solve!(sim)
            @info "Finished simulation for worker $(Distributed.myid())"
        end
    end
end

site_IDs = repeat(["US-MOz", "US-NR1", "US-Ha1"], 100)[1:nworkers]
duration_days = 90

# pmap distributes one site per worker and waits for all to finish
pmap(site_IDs) do site_ID
    run_on_stream(site_ID, duration_days)
end

# 3 workers
# About 110 SYPD for each worker (in total 300 sypd vs 200 spyd for an
# individual worker)

# 10 workers
# About 30 SYPD for each workers (in total 300 sypd)


# This might requires nvidia-cuda-mps-control to actually allow parallelism
# between multiple processes
