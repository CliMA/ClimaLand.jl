# NVTX.jl

NVTX.jl provides Julia bindings to the [NVIDIA Tools Extension Library (NVTX)](https://nvidia.github.io/NVTX/doxygen/index.html) for instrumenting Julia code for use with the [Nsight systems profiler](https://developer.nvidia.com/nsight-systems).

## Installation

NVTX.jl includes the necessary NVTX library for supported platforms, however you will need to install the [NVIDIA Nsight Systems](https://developer.nvidia.com/nsight-systems) to actually run the profiler: it is available for Linux (x86\_64, Aarch64, Power) and Windows (x86_64), and the profile viewer is available on Linux, Windows and MacOS.

NVTX.jl can be loaded on any platform, even those without the NVTX library, and so can safely be included as a package dependency.

## Instrumenting Julia code

The simplest way to use this package is through the macros [`NVTX.@mark`](@ref), [`NVTX.@range`](@ref), [`NVTX.@annotate`](@ref).

```julia
NVTX.@mark "my message"

NVTX.@range "my message" begin
    # code to measure
end

NVTX.@annotate function myfunction(args...)
    # function body
end
```

## Instrumenting Julia internals

Additionally, it is possible to annotate the Julia garbage collector and inference by calling [`NVTX.enable_gc_hooks()`](@ref) and [`NVTX.enable_inference_hook()`](@ref), or setting the `JULIA_NVTX_CALLBACKS` environment variable to a comma (`,`) or bar (`|`) separated list of:

- `gc`: instrument the Julia garbage collector
- `alloc`: instrument calls to alloc
- `free`: instrument calls to free
- `inference`: instrument calls to compiler inference.

## Running the profiler

To run the Nsight Systems profiler on a script, use
```
nsys profile julia script.jl
```

See [Nsight Systems User Manual](https://docs.nvidia.com/nsight-systems/UserGuide/index.html) for more information.
