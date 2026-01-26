module ClimaCommsCUDAExt

import CUDA

import Adapt
import ClimaComms
import ClimaComms: CUDADevice

function ClimaComms._assign_device(::CUDADevice, rank_number)
    CUDA.device!(rank_number % CUDA.ndevices())
    return nothing
end

function Base.summary(io::IO, ::CUDADevice)
    dev = CUDA.device()
    name = CUDA.name(dev)
    uuid = CUDA.uuid(dev)
    return "$name ($uuid)"
end

function ClimaComms.device_functional(::CUDADevice)
    return CUDA.functional()
end

function Adapt.adapt_structure(
    to::Type{<:CUDA.CuArray},
    ctx::ClimaComms.AbstractCommsContext,
)
    return ClimaComms.context(Adapt.adapt(to, ClimaComms.device(ctx)))
end

Adapt.adapt_structure(
    ::Type{<:CUDA.CuArray},
    device::ClimaComms.AbstractDevice,
) = ClimaComms.CUDADevice()

ClimaComms.array_type(::CUDADevice) = CUDA.CuArray
ClimaComms.free_memory(::CUDADevice) = CUDA.free_memory()
ClimaComms.total_memory(::CUDADevice) = CUDA.total_memory()
ClimaComms.allowscalar(f, ::CUDADevice, args...; kwargs...) =
    CUDA.@allowscalar f(args...; kwargs...)

# Extending ClimaComms methods that operate on expressions (cannot use dispatch here)
ClimaComms.sync(f::F, ::CUDADevice, args...; kwargs...) where {F} =
    CUDA.@sync f(args...; kwargs...)
ClimaComms.cuda_sync(f::F, ::CUDADevice, args...; kwargs...) where {F} =
    CUDA.@sync f(args...; kwargs...)
ClimaComms.time(f::F, ::CUDADevice, args...; kwargs...) where {F} =
    CUDA.@time f(args...; kwargs...)
ClimaComms.elapsed(f::F, ::CUDADevice, args...; kwargs...) where {F} =
    CUDA.@elapsed f(args...; kwargs...)
ClimaComms.assert(::CUDADevice, cond::C, text::T) where {C, T} =
    isnothing(text) ? (CUDA.@cuassert cond()) : (CUDA.@cuassert cond() text())

# The number of threads in the kernel being executed by the calling thread.
threads_in_kernel() = CUDA.blockDim().x * CUDA.gridDim().x

# The index of the calling thread, which is between 1 and threads_in_kernel().
thread_index() =
    (CUDA.blockIdx().x - 1) * CUDA.blockDim().x + CUDA.threadIdx().x

# The maximum number of blocks that can fit on the GPU used for this kernel.
grid_size_limit(kernel) = CUDA.attribute(
    CUDA.device(kernel.fun.mod.ctx),
    CUDA.DEVICE_ATTRIBUTE_MAX_GRID_DIM_X,
)

# Either the first value if it is available, or the maximum number of threads
# that can fit in one block of this kernel (cuOccupancyMaxPotentialBlockSize).
# With enough blocks, the latter value will maximize the occupancy of the GPU.
block_size_limit(max_threads_in_block::Int, _) = max_threads_in_block
block_size_limit(::Val{:auto}, kernel) =
    CUDA.launch_configuration(kernel.fun).threads

function ClimaComms.run_threaded(
    f::F,
    ::CUDADevice,
    ::Val,
    itr;
    block_size,
) where {F}
    n_items = length(itr)
    n_items > 0 || return nothing

    function call_f_from_thread()
        item_index = thread_index()
        item_index <= n_items &&
            @inbounds f(itr[firstindex(itr) + item_index - 1])
        return nothing
    end
    kernel = CUDA.@cuda always_inline=true launch=false call_f_from_thread()
    max_blocks = grid_size_limit(kernel)
    max_threads_in_block = block_size_limit(block_size, kernel)

    # If there are too many items, coarsen by the smallest possible amount.
    n_items <= max_blocks * max_threads_in_block ||
        return ClimaComms.run_threaded(f, CUDADevice(), 1, itr; block_size)

    threads_in_block = min(max_threads_in_block, n_items)
    blocks = cld(n_items, threads_in_block)
    kernel(; blocks, threads = threads_in_block)
end

function ClimaComms.run_threaded(
    f::F,
    ::CUDADevice,
    min_items_in_thread::Int,
    itr;
    block_size,
) where {F}
    min_items_in_thread > 0 || throw(ArgumentError("`coarsen` is not positive"))
    n_items = length(itr)
    n_items > 0 || return nothing

    # Maximize memory coalescing with a "grid-stride loop"; for reference, see
    # https://developer.nvidia.com/blog/cuda-pro-tip-write-flexible-kernels-grid-stride-loops
    call_f_from_thread() =
        for item_index in thread_index():threads_in_kernel():n_items
            @inbounds f(itr[firstindex(itr) + item_index - 1])
        end
    kernel = CUDA.@cuda always_inline=true launch=false call_f_from_thread()
    max_blocks = grid_size_limit(kernel)
    max_threads_in_block = block_size_limit(block_size, kernel)

    # If there are too many items to use the specified coarsening, increase it
    # by the smallest possible amount.
    max_required_threads = cld(n_items, min_items_in_thread)
    items_in_thread =
        max_required_threads <= max_blocks * max_threads_in_block ?
        min_items_in_thread : cld(n_items, max_blocks * max_threads_in_block)

    threads_in_block = min(max_threads_in_block, max_required_threads)
    blocks = cld(n_items, items_in_thread * threads_in_block)
    kernel(; blocks, threads = threads_in_block)
end

end
