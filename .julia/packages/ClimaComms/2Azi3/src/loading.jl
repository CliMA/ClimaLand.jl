import ..ClimaComms

export @import_required_backends

"""
    mpi_is_required()

Returns a Bool indicating if MPI should be loaded, based on the
`ENV["CLIMACOMMS_CONTEXT"]`. See [`ClimaComms.context`](@ref) for
more information.

```julia
mpi_is_required() && using MPI
```
"""
mpi_is_required() = context_type() == :MPICommsContext

"""
    cuda_is_required()

Returns a Bool indicating if CUDA should be loaded, based on the
`ENV["CLIMACOMMS_DEVICE"]`. See [`ClimaComms.device`](@ref) for
more information.

```julia
cuda_is_required() && using CUDA
```
"""
cuda_is_required() = device_type() == :CUDADevice

cuda_ext_is_loaded() =
    !isnothing(Base.get_extension(ClimaComms, :ClimaCommsCUDAExt))

mpi_ext_is_loaded() =
    !isnothing(Base.get_extension(ClimaComms, :ClimaCommsMPIExt))

"""
    ClimaComms.@import_required_backends

If the desired context is MPI (as determined by `ClimaComms.context()`), try loading MPI.jl.
If the desired device is CUDA (as determined by `ClimaComms.device()`), try loading CUDA.jl.
"""
macro import_required_backends()
    return quote
        @static if $mpi_is_required()
            @debug "Loading MPI via `import MPI`..."
            import MPI
        end
        @static if $cuda_is_required()
            @debug "Loading CUDA via `import CUDA`..."
            import CUDA
        end
    end
end
