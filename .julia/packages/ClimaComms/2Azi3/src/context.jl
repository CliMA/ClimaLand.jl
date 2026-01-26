import ..ClimaComms

function context_type()
    name = get(ENV, "CLIMACOMMS_CONTEXT", nothing)
    if !isnothing(name)
        if name == "MPI"
            return :MPICommsContext
        elseif name == "SINGLETON"
            return :SingletonCommsContext
        else
            error("Invalid context: $name")
        end
    end
    # detect common environment variables used by MPI launchers
    #   PMI_RANK appears to be used by MPICH and srun
    #   OMPI_COMM_WORLD_RANK appears to be used by OpenMPI
    if haskey(ENV, "PMI_RANK") || haskey(ENV, "OMPI_COMM_WORLD_RANK")
        return :MPICommsContext
    else
        return :SingletonCommsContext
    end
end

"""
    ClimaComms.context(device=device())

Construct a default communication context.

By default, it will try to determine if it is running inside an MPI environment
variables are set; if so it will return a [`MPICommsContext`](@ref); otherwise
it will return a [`SingletonCommsContext`](@ref).

Behavior can be overridden by setting the `CLIMACOMMS_CONTEXT` environment variable
to either `MPI` or `SINGLETON`.
"""
function context(device = device(); target_context = context_type())
    if target_context == :MPICommsContext && !mpi_ext_is_loaded()
        error(
            "Loading MPI.jl is required to use MPICommsContext. You might want to call ClimaComms.@import_required_backends",
        )
    end
    ContextConstructor = getproperty(ClimaComms, target_context)
    return ContextConstructor(device)
end

"""
    AbstractCommsContext

The base type for a communications context. Each backend defines a
concrete subtype of this.
"""
abstract type AbstractCommsContext end

"""
    (pid, nprocs) = init(ctx::AbstractCommsContext)

Perform any necessary initialization for the specified backend. Return a
tuple of the processor ID and the number of participating processors.
"""
function init end

"""
    mypid(ctx::AbstractCommsContext)

Return the processor ID.
"""
function mypid end

"""
    iamroot(ctx::AbstractCommsContext)

Return `true` if the calling processor is the root processor.
"""
function iamroot end

"""
    nprocs(ctx::AbstractCommsContext)

Return the number of participating processors.
"""
function nprocs end


"""
    barrier(ctx::CC) where {CC <: AbstractCommsContext}

Perform a global synchronization across all participating processors.
"""
function barrier end
barrier(::Nothing) = nothing

"""
    reduce(ctx::CC, val, op) where {CC <: AbstractCommsContext}

Perform a reduction across all participating processors, using `op` as
the reduction operator and `val` as this rank's reduction value. Return
the result to the first processor only.
"""
function reduce end
reduce(::Nothing, val, op) = val

"""
    reduce!(ctx::CC, sendbuf, recvbuf, op)
    reduce!(ctx::CC, sendrecvbuf, op)

Performs elementwise reduction using the operator `op` on the buffer `sendbuf`, storing the result in the `recvbuf` of the process.
If only one `sendrecvbuf` buffer is provided, then the operation is performed in-place.

"""
function reduce! end

"""
    allreduce(ctx::CC, sendbuf, op)

Performs elementwise reduction using the operator `op` on the buffer `sendbuf`, allocating a new array for the result.
`sendbuf` can also be a scalar, in which case `recvbuf` will be a value of the same type.
"""
function allreduce end

"""
    allreduce!(ctx::CC, sendbuf, recvbuf, op)
    allreduce!(ctx::CC, sendrecvbuf, op)

Performs elementwise reduction using the operator `op` on the buffer `sendbuf`, storing the result in the `recvbuf` of all processes in the group.
`Allreduce!` is equivalent to a `Reduce!` operation followed by a `Bcast!`, but can lead to better performance.
If only one `sendrecvbuf` buffer is provided, then the operation is performed in-place.

"""
function allreduce! end

"""
    gather(ctx::AbstractCommsContext, array)

Gather an array of values from all processors into a single array,
concattenating along the last dimension.
"""
gather(::Nothing, array) = array

"""
    bcast(ctx::AbstractCommsContext, object)

Broadcast `object` from the root process to all other processes.
The value of `object` on non-root processes is ignored.
"""
function bcast end

"""
    abort(ctx::CC, status::Int) where {CC <: AbstractCommsContext}

Terminate the caller and all participating processors with the specified
`status`.
"""
function abort end
abort(::Nothing, status::Int) = exit(status)

"""
    AbstractGraphContext

A context for communicating between processes in a graph.
"""
abstract type AbstractGraphContext end

"""
    graph_context(context::AbstractCommsContext,
                  sendarray, sendlengths, sendpids,
                  recvarray, recvlengths, recvpids)

Construct a communication context for exchanging neighbor data via a graph.

Arguments:
- `context`: the communication context on which to construct the graph context.
- `sendarray`: array containing data to send
- `sendlengths`: list of lengths of data to send to each process ID
- `sendpids`: list of processor IDs to send
- `recvarray`: array to receive data into
- `recvlengths`: list of lengths of data to receive from each process ID
- `recvpids`: list of processor IDs to receive from

This should return an `AbstractGraphContext` object.
"""
function graph_context end


"""
    start(ctx::AbstractGraphContext)

Initiate graph data exchange.
"""
function start end

"""
    progress(ctx::AbstractGraphContext)

Drive communication. Call after `start` to ensure that communication
proceeds asynchronously.
"""
function progress end

"""
    finish(ctx::AbstractGraphContext)

Complete the communications step begun by `start()`. After this returns,
data received from all neighbors will be available in the stage areas of
each neighbor's receive buffer.
"""
function finish end
