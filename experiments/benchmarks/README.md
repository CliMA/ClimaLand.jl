# ClimaLand Benchmarks

This folder contains scripts that setup and benchmark/profile certain ClimaLand models.

## `benchmark_utils.jl`

This contains functions that are used to profile individual simulations. The utilities in this
file all require a simulation setup function to be provided by the user. This function is
called with no arguments, and is expected to return a `ClimaLand.Simulations.LandSimulation`,
which will be benchmarked/profiled.

## Adding a Performance Test

The general format of the performance tests in this folder is:

```julia
delete!(ENV, "JULIA_CUDA_MEMORY_POOL")
import ClimaLand

######################################################################
## This result is from a benchmark ran on an A100 on the clima cluster
const PREVIOUS_GPU_TIME_S = 1.0
## This result is from a benchmark ran with a single process on the clima cluster
const PREVIOUS_CPU_TIME_S = 100.0
######################################################################

include("benchmark_utils.jl")

context = ClimaComms.context()
ClimaComms.init(context)
device = ClimaComms.device()
device_suffix = device isa ClimaComms.CPUSingleThreaded ? "cpu" : "gpu"

outdir = "mymodel_benchmark_$(device_suffix)"
!ispath(outdir) && mkpath(outdir)

function setup_mymodel()
    # setup my model
    # return a LandSimulation
end

reference_time = device isa ClimaComms.AbstractCPUDevice ? PREVIOUS_CPU_TIME_S : PREVIOUS_GPU_TIME_S
profile_and_benchmark(setup_bucket, device, reference_time, outdir)
```

### Adding to the buildkite short_perf pipeline

This pipeline is intended to be used as a benchmark that is run for each PR. It only runs
a GPU benchmark for models who's source files have been changed. This uses less resources than
the full benchmark pipeline, but will not necessarily run the benchmark for a model if a method
is added outside of the tracked source files.

When the short_perf pipeline is triggered, it generates a list of files that are different
from the branch point off of main. Each individual benchmark can then be triggered if any file
in its "watched" files are modified. To add a benchmark, add an entry to the "watch".
The "path" key for the entry should be a list of files or paths that can be a glob pattern. The normal
buildkite step entry can then be added within the "config" key.



## `profile_and_benchmark`

`profile_and_benchmark` behaves differently if an external profiler is used.

### Usage with the internal profiler

If the ClimaComms device is a `ClimaComms.CPUSingleThreaded`, steps 2-7 are benchmarked. Then the entire simulation
is profiled, and a flame graph and allocations flame graph are saved.

If the ClimaComms device is a `ClimaComms.CUDADevice`, the entire simulation, excluding the first step, is benchmarked
by repeatedly running it. Then three steps of the simulation
are profiled. The results are printed to stdout and saved as csv files.

If ran on buildkite, the mean simulation
time is compared to the previous best time for that device. If it is significantly slower,
the test fails.

### Usage with an external profiler

If an external CUDA profiler is detected, the simulation isprofiled for three steps,
and no benchmarks are ran. This only supports a `ClimaComms.CUDADevice` device.
