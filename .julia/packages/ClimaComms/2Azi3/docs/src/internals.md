# Developing with `ClimaComms`

This page goes into more depth about `ClimaComms` in `CliMA` packages.

First, we will describe what `Device`s and `Context`s are.

### `Device`s

`Device`s identify a specific type of computing hardware (e.g., a CPU/a NVidia
GPU, et cetera). The `Device`s implemented are
- [`CPUSingleThreaded`](@ref ClimaComms.CPUSingleThreaded), for a CPU core with a single thread;
- [`CUDADevice`](@ref ClimaComms.CUDADevice), for a single CUDA GPU.

`Device`s in `ClimaComms` are
[singletons](https://docs.julialang.org/en/v1/manual/types/#man-singleton-types),
meaning that they do not contain any data: they are used exclusively to
determine what implementation of a given function or data structure should be
used. For example, let us implement a function to allocate an array `what` on
the CPU or on the GPU depending on the `Device`:

```julia
import ClimaComms: CPUSingleThreaded, CUDADevice
import CUDA

function allocate(device::CPUSingleThreaded, what)
    return Array(what)
end

function allocate(device::CUDADevice, what)
    return CUDA.CuArray(what)
end
```

If we want to allocate memory on the GPU, we would do something like:
```julia-repl
julia> allocate(CUDADevice(), [1, 2, 3])
CUDA.CuArray([1, 2, 3])
```

Low-level `CliMA` code often needs to implement different methods for different
`Device`s (e.g., in [ClimaCore](https://github.com/CliMA/ClimaCore.jl)), but
this level of specialization is often not required at higher levels.

Higher-level code often interacts with `Device`s through `ClimaComms` functions
such [`time`](@ref ClimaComms.@time) or [`sync`](@ref ClimaComms.@sync). These functions implement device-agnostic
operations. For instance, the proper way to compute how long a given expression takes to compute is
```julia
import ClimaComms: @time

device = ClimaComms.device()  # or the device of interest

@time device my_expr
```
This will ensure that the correct time is computed (e.g., the time to run a GPU kernel, and not the time to launch the kernel).

For a complete list of such functions, consult the [APIs](@ref) page.


### `Context`s

A `Context` contains the information needed for multiple devices to communicate.
For simulations with only one device, [`SingletonCommsContext`](@ref ClimaComms.SingletonCommsContext) simply
contains an instance of an [`AbstractDevice`](@ref ClimaComms.AbstractDevice). For `MPI` simulations, the
context contains the MPI communicator as well.

`Context`s specify devices and form of parallelism, so they are often passed
around in both low-level and higher-level code.

`ClimaComms` provide functions that are context-agnostic. For instance,
[`reduce`](@ref ClimaComms.reduce) applies a given function to an array across difference
processes and collects the result. Let us see an example
```julia
import ClimaComms
ClimaComms.@import_required_backends

context = ClimaComms.context()  # Default context (read from environment variables)
device = ClimaComms.device(context)  # Default device

mypid, nprocs = ClimaComms.init(context)

ArrayType = ClimaComms.array_type(device)

my_array = mypid * ArrayType([1, 1, 1])

reduced_array = ClimaComms.reduce(context, my_array, +)
ClimaComms.iamroot(context) && @show reduced_array
```

[`@import_required_backends`](@ref ClimaComms.@import_required_backends) is responsible for loading relevant
libraries, for more information refer to the [Backends and extensions](@ref)
section.

In this snippet, we obtained the default context from environment variables
using the [`context`](@ref ClimaComms.context) function. As developers, we do not know whether this
code is being run on a single process or multiple, so took the more generic
stance that the code _might_ be run on several processes. When several processes
are being used, the same code is being run by parallel Julia instances, each
with a different process id (`pid`). The line ```julia mypid, nprocs =
ClimaComms.init(context) ``` assigns different `mypid` to the various processes
and returns `nprocs` so that we can use this information if needed. This
function is also responsible for distributing GPUs across processes, if
relevant.

In this example, we used `mypid` to set up `my_array` in such a way that it
would be different on different processes. We set up the array with `ArrayType`,
obtained with [`array_type`](@ref ClimaComms.array_type). This function provides the type to allocate
the array on the correct device (CPU or GPU). Then, we applied [`reduce`](@ref ClimaComms.reduce)
to sum them all. [`reduce`](@ref ClimaComms.reduce) collects the result to the _root_ process, the
one with `pid = 1`. For single-process runs, the only process is also a root
process.

The code above works independently on the number of processes and over all the
devices supported by `ClimaComms`.

### Backends and extensions

Except the most basic ones, `ClimaComms` computing devices and contexts are
implemented as independent backends. For instance, `ClimaComms` provides an
`AbstractDevice` interface for which `CUDADevice` is an implementation that
depends on [`CUDA.jl`](https://github.com/JuliaGPU/CUDA.jl). Scripts that use
`ClimaComms` have to load the packages that power the desired backend (e.g.,
`CUDA.jl` has to be explicitly loaded if one wants to use `CUDADevice`s).

When using the default context and device (as read from environment variables),
`ClimaComms` can automatically load the required packages. To instruct
`ClimaComms` to do, add [`ClimaComms.@import_required_backends`](@ref) at the
beginning of your scripts. `ClimaComms` can also identify some cases when a
package is required but not loaded and warn you about it.

Using non-trivial backends might require you to install `CUDA.jl` and/or
`MPI.jl` in your environment.

> Note: When using [`context`](@ref ClimaComms.context) to select the context, it is safe to always
> add [`ClimaComms.@import_required_backends`](@ref) at the top of your scripts.
> *Do not add* [`ClimaComms.@import_required_backends`](@ref) to library code
> (i.e., in `src`) because the macro requires dependencies that should not be
> satisfied in that instance.

Technically, `ClimaComms` backends are implemented in Julia
[extensions](https://pkgdocs.julialang.org/v1/creating-packages/#Conditional-loading-of-code-in-packages-(Extensions)).
There are two main reasons for this:
- optional packages (such as `CUDA` and `MPI`) should not be hard dependencies
  of `ClimaComms`;
- implementing code in extensions reduces the loading time for `ClimaComms` and
  downstream packages because it skips compilation of code that is not needed.

If you are implementing features for `ClimaComms`, make sure that your
backend-specific code is in a Julia extension (in the `ext` folder).
