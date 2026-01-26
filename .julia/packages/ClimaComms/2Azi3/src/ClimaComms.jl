"""
    ClimaComms

Abstracts the communications interface for the various CliMA components
in order to:
- support different computational backends (CPUs, GPUs)
- enable the use of different backends as transports (MPI, SharedArrays,
  etc.), and
- transparently support single or double buffering for GPUs, depending
  on whether the transport has the ability to access GPU memory.
"""
module ClimaComms

include("devices.jl")
include("context.jl")
include("singleton.jl")
include("mpi.jl")
include("loading.jl")
include("adapt.jl")
include("logging.jl")

end # module
