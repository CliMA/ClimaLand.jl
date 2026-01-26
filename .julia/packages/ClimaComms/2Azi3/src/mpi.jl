"""
    MPICommsContext()
    MPICommsContext(device)
    MPICommsContext(device, comm)

A MPI communications context, used for distributed runs.
[`AbstractCPUDevice`](@ref) and [`CUDADevice`](@ref) device options are currently supported.
"""
struct MPICommsContext{D <: AbstractDevice, C} <: AbstractCommsContext
    device::D
    mpicomm::C
end

function MPICommsContext end

"""
    local_communicator(ctx::MPICommsContext)

Internal function to create a new MPI communicator for processes on the same physical node.

Sample Usage:
```
ClimaComms.local_communicator(ctx) do local_comm
    ClimaComms._assign_device(ClimaComms.device(ctx), MPI.Comm_rank(local_comm))
end
```
"""
function local_communicator end
