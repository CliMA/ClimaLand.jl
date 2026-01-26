"""
    SingletonCommsContext(device=device())

A singleton communications context, used for single-process runs.
[`ClimaComms.AbstractCPUDevice`](@ref) and [`ClimaComms.CUDADevice`](@ref)
device options are currently supported.
"""
struct SingletonCommsContext{D <: AbstractDevice} <: AbstractCommsContext
    device::D
end

SingletonCommsContext() = SingletonCommsContext(device())

device(ctx::SingletonCommsContext) = ctx.device

init(::SingletonCommsContext) = (1, 1)

mypid(::SingletonCommsContext) = 1
iamroot(::SingletonCommsContext) = true
nprocs(::SingletonCommsContext) = 1
barrier(::SingletonCommsContext) = nothing
reduce(::SingletonCommsContext, val, op) = val
gather(::SingletonCommsContext, array) = array
allreduce(::SingletonCommsContext, sendbuf, op) = sendbuf
bcast(::SingletonCommsContext, object) = object

function reduce!(::SingletonCommsContext, sendbuf, recvbuf, op)
    copyto!(recvbuf, sendbuf)
    return nothing
end
function reduce!(::SingletonCommsContext, sendrecvbuf, op)
    return nothing
end

function allreduce!(::SingletonCommsContext, sendbuf, recvbuf, op)
    copyto!(recvbuf, sendbuf)
    return nothing
end
function allreduce!(::SingletonCommsContext, sendrecvbuf, op)
    return nothing
end

struct SingletonGraphContext <: AbstractGraphContext
    context::SingletonCommsContext
end

graph_context(ctx::SingletonCommsContext, args...) = SingletonGraphContext(ctx)

start(gctx::SingletonGraphContext) = nothing
progress(gctx::SingletonGraphContext) = nothing
finish(gctx::SingletonGraphContext) = nothing

function Base.summary(io::IO, ctx::SingletonCommsContext)
    println(io, "Context: $(nameof(typeof(ctx)))")
    println(io, "Device: $(typeof(device(ctx)))")
end
