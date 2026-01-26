# ClimaComms

`ClimaComms.jl` is a small package that provides abstractions for different
computing devices and environments. `ClimaComms.jl` is use extensively by
`CliMA` packages to control where and how simulations are run (e.g., one on
core, on multiple GPUs, et cetera).

This page highlights the most important user-facing `ClimaComms` concepts. If
you are using `ClimaComms` as a developer, refer to the [Developing with
`ClimaComms`](@ref) page. For a detailed list of all the functions and objects
implemented, the [APIs](@ref) page collects all of them.

## `Device`s and `Context`s

The two most important objects in `ClimaComms.jl` are the [`Device`](@ref
ClimaComms.AbstractDevice) and the [`Context`](@ref ClimaComms.AbstractCommsContext).

A `Device` identifies a computing device, a piece of hardware that will be
executing some code. The `Device`s currently implemented are
- [`CPUSingleThreaded`](@ref ClimaComms.CPUSingleThreaded) for a CPU core with a single thread,
- [`CPUMultiThreaded`](@ref ClimaComms.CPUMultiThreaded) for a CPU core with multiple threads,
- [`CUDADevice`](@ref ClimaComms.CUDADevice) for a single CUDA-enabled GPU.

`Device`s are part of [`Context`](@ref ClimaComms.AbstractCommsContext)s,
objects that contain information require for multiple `Device`s to communicate.
Implemented `Context`s are
- [`SingletonCommsContext`](@ref ClimaComms.SingletonCommsContext), when there is no parallelism;
- [`MPICommsContext`](@ref ClimaComms.MPICommsContext) , for a MPI-parallelized runs.

To choose a device and a context, most `CliMA` packages use the
[`device()`](@ref ClimaComms.device()) and [`context()`](@ref ClimaComms.context()) functions. These functions look at
specific environment variables and set the `device` and `context` accordingly.
By default, the [`CPUSingleThreaded`](@ref ClimaComms.CPUSingleThreaded) device is chosen and the context is
set to [`SingletonCommsContext`](@ref ClimaComms.SingletonCommsContext) unless `ClimaComms` detects being run in
a standard MPI launcher (as `srun` or `mpiexec`).

For example, to run a simulation on a GPU, run `julia` as
```bash
export CLIMACOMMS_DEVICE="CUDA"
export CLIMACOMMS_CONTEXT="SINGLETON"
# call/open julia as usual
```

!!! note
    There might be other ways to control the device and context. Please,
    refer to the documentation of the specific package to learn more.

## Running with MPI/CUDA

`CliMA` packages do not depend directly on `MPI` or `CUDA`, so, if you want to
run your simulation in parallel mode and/or on GPUs, you will need to install
some packages separately.

For parallel simulations, [`MPI.jl`](https://github.com/JuliaParallel/MPI.jl), and
for GPU runs, [`CUDA.jl`](https://github.com/JuliaGPU/CUDA.jl). You can install
these packages in your base environment
```bash
julia -E "using Pkg; Pkg.add(\"CUDA\"); Pkg.add(\"MPI\")"
```
Some packages come with environments that includes all possible backends
(typically `.buildkite`). You can also consider directly using those
environments.

# Writing generic kernels

To implement kernels that work across all `Device`s, `ClimaComms.jl` provides
the [`@threaded`](@ref ClimaComms.@threaded) macro. Like the standard library's
[`Threads.@threads`](https://docs.julialang.org/en/v1/base/multi-threading/#Base.Threads.@threads)
macro, this can be placed in front of any `for`-loop to parallelize its
iterations across threads. For example, given two vectors `a` and `b` of equal
size, a threaded version of the `copyto!` function can be implemented as
```julia-repl
julia> threaded_copyto!(a, b) = ClimaComms.@threaded for i in axes(a, 1)
           a[i] = b[i]
       end
threaded_copyto! (generic function with 1 method)

julia> threaded_copyto!(a, b)
```
This macro offers several options for fine-tuning performance:
- the `Device` can be specified to reduce compilation time
- the level of thread coarsening (number of loop iterations evaluated in each
  thread) can be specified to reduce the overhead of launching many threads
- on GPUs, the size of each block (a collection of threads launched on a single
  multiprocessor) can be specified to improve GPU
  [utilization](https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#multiprocessor-level)

!!! note
    Executing code on GPUs requires *static compilation*, which means that the
    compiler must be able to infer the types of all variables used in a threaded
    loop. For example, global variables and type variables defined outside a
    loop should be avoided when the loop is parallelized on a GPU. See the
    docstring of [`@threaded`](@ref ClimaComms.@threaded) for more information.
