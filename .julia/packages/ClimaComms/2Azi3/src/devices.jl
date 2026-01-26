import ..ClimaComms
import Adapt

"""
    AbstractDevice

The base type for a device.
"""
abstract type AbstractDevice end

"""
    AbstractCPUDevice()

Abstract device type for single-threaded and multi-threaded CPU runs.
"""
abstract type AbstractCPUDevice <: AbstractDevice end


"""
    CPUSingleThreaded()

Use the CPU with single thread.
"""
struct CPUSingleThreaded <: AbstractCPUDevice end

"""
    CPUMultiThreaded()

Use the CPU with multiple thread.
"""
struct CPUMultiThreaded <: AbstractCPUDevice end

"""
    CUDADevice()

Use NVIDIA GPU accelarator
"""
struct CUDADevice <: AbstractDevice end

"""
    ClimaComms.device_functional(device)

Return true when the `device` is correctly set up.
"""
function device_functional end

device_functional(::CPUSingleThreaded) = true
device_functional(::CPUMultiThreaded) = true

function device_type()
    env_var = get(ENV, "CLIMACOMMS_DEVICE", "CPU")
    if env_var == "CPU"
        return Threads.nthreads() > 1 ? :CPUMultiThreaded : :CPUSingleThreaded
    elseif env_var == "CPUSingleThreaded"
        return :CPUSingleThreaded
    elseif env_var == "CPUMultiThreaded"
        return :CPUMultiThreaded
    elseif env_var == "CUDA"
        return :CUDADevice
    else
        error("Invalid CLIMACOMMS_DEVICE: $env_var")
    end
end

"""
    ClimaComms.device()

Determine the device to use depending on the `CLIMACOMMS_DEVICE` environment variable.

Allowed values:
- `CPU`, single-threaded or multi-threaded depending on the number of threads;
- `CPUSingleThreaded`,
- `CPUMultiThreaded`,
- `CUDA`.

The default is `CPU`.
"""
function device()
    target_device = device_type()
    if target_device == :CUDADevice && !cuda_ext_is_loaded()
        error(
            "Loading CUDA.jl is required to use CUDADevice. You might want to call ClimaComms.@import_required_backends",
        )
    end
    DeviceConstructor = getproperty(ClimaComms, target_device)
    return DeviceConstructor()
end

Base.summary(io::IO, device::AbstractDevice) = string(device_type())

"""
    ClimaComms.array_type(::AbstractDevice)

The base array type used by the specified device (currently `Array` or `CuArray`).
"""
array_type(::AbstractCPUDevice) = Array

"""
Internal function that can be used to assign a device to a process.

Currently used to assign CUDADevices to MPI ranks.
"""
_assign_device(device, id) = nothing

"""
    ClimaComms.free_memory(device)

Bytes of memory that are currently available for allocation on the `device`.
"""
free_memory(::AbstractCPUDevice) = Sys.free_memory()

"""
    ClimaComms.total_memory(device)

Bytes of memory that are theoretically available for allocation on the `device`.
"""
total_memory(::AbstractCPUDevice) = Sys.total_memory()

"""
    @time f(args...; kwargs...)

Device-flexible `@time`:

Calls
```julia
@time f(args...; kwargs...)
```
for CPU devices and
```julia
CUDA.@time f(args...; kwargs...)
```
for CUDA devices.
"""
function time(f::F, device::AbstractCPUDevice, args...; kwargs...) where {F}
    Base.@time begin
        f(args...; kwargs...)
    end
end

"""
    elapsed(f::F, device::AbstractDevice, args...; kwargs...)

Device-flexible `elapsed`.

Calls
```julia
@elapsed f(args...; kwargs...)
```
for CPU devices and
```julia
CUDA.@elapsed f(args...; kwargs...)
```
for CUDA devices.
"""
function elapsed(f::F, device::AbstractCPUDevice, args...; kwargs...) where {F}
    Base.@elapsed begin
        f(args...; kwargs...)
    end
end

"""
    sync(f, ::AbstractDevice, args...; kwargs...)

Device-flexible function that calls `@sync`.

Calls
```julia
@sync f(args...; kwargs...)
```
for CPU devices and
```julia
CUDA.@sync f(args...; kwargs...)
```
for CUDA devices.

An example use-case of this might be:
```julia
BenchmarkTools.@benchmark begin
    if ClimaComms.device() isa ClimaComms.CUDADevice
        CUDA.@sync begin
            launch_cuda_kernels_or_spawn_tasks!(...)
        end
    elseif ClimaComms.device() isa ClimaComms.CPUMultiThreading
        Base.@sync begin
            launch_cuda_kernels_or_spawn_tasks!(...)
        end
    end
end
```

If the CPU version of the above example does not leverage
spawned tasks (which require using `Base.sync` or `Threads.wait`
to synchronize), then you may want to simply use [`cuda_sync`](@ref).
"""
function sync(f::F, ::AbstractCPUDevice, args...; kwargs...) where {F}
    Base.@sync begin
        f(args...; kwargs...)
    end
end

"""
    cuda_sync(f, ::AbstractDevice, args...; kwargs...)

Device-flexible function that (may) call `CUDA.@sync`.

Calls
```julia
f(args...; kwargs...)
```
for CPU devices and
```julia
CUDA.@sync f(args...; kwargs...)
```
for CUDA devices.
"""
function cuda_sync(f::F, ::AbstractCPUDevice, args...; kwargs...) where {F}
    f(args...; kwargs...)
end

"""
    allowscalar(f, ::AbstractDevice, args...; kwargs...)

Device-flexible version of `CUDA.@allowscalar`.

Lowers to
```julia
f(args...)
```
for CPU devices and
```julia
CUDA.@allowscalar f(args...)
```
for CUDA devices.

This is usefully written with closures via
```julia
allowscalar(device) do
    f()
end
```
"""
allowscalar(f, ::AbstractCPUDevice, args...; kwargs...) = f(args...; kwargs...)

"""
    @time device expr

Device-flexible `@time`.

Lowers to
```julia
@time expr
```
for CPU devices and
```julia
CUDA.@time expr
```
for CUDA devices.
"""
macro time(device, expr)
    __CC__ = ClimaComms
    return :($__CC__.time(() -> $(esc(expr)), $(esc(device))))
end

"""
    @sync device expr

Device-flexible `@sync`.

Lowers to
```julia
@sync expr
```
for CPU devices and
```julia
CUDA.@sync expr
```
for CUDA devices.

An example use-case of this might be:
```julia
BenchmarkTools.@benchmark begin
    if ClimaComms.device() isa ClimaComms.CUDADevice
        CUDA.@sync begin
            launch_cuda_kernels_or_spawn_tasks!(...)
        end
    elseif ClimaComms.device() isa ClimaComms.CPUMultiThreading
        Base.@sync begin
            launch_cuda_kernels_or_spawn_tasks!(...)
        end
    end
end
```

If the CPU version of the above example does not leverage
spawned tasks (which require using `Base.sync` or `Threads.wait`
to synchronize), then you may want to simply use [`@cuda_sync`](@ref).
"""
macro sync(device, expr)
    __CC__ = ClimaComms
    return :($__CC__.sync(() -> $(esc(expr)), $(esc(device))))
end

"""
    @cuda_sync device expr

Device-flexible `CUDA.@sync`.

Lowers to
```julia
expr
```
for CPU devices and
```julia
CUDA.@sync expr
```
for CUDA devices.
"""
macro cuda_sync(device, expr)
    __CC__ = ClimaComms
    return :($__CC__.cuda_sync(() -> $(esc(expr)), $(esc(device))))
end

"""
    @elapsed device expr

Device-flexible `@elapsed`.

Lowers to
```julia
@elapsed expr
```
for CPU devices and
```julia
CUDA.@elapsed expr
```
for CUDA devices.
"""
macro elapsed(device, expr)
    __CC__ = ClimaComms
    return :($__CC__.elapsed(() -> $(esc(expr)), $(esc(device))))
end

"""
    @assert device cond [text]

Device-flexible `@assert`.

Lowers to
```julia
@assert cond [text]
```
for CPU devices and
```julia
CUDA.@cuassert cond [text]
```
for CUDA devices.
"""
macro assert(device, cond, text = nothing)
    text_func = isnothing(text) ? nothing : :(() -> $(esc(text)))
    return :($assert($(esc(device)), () -> $(esc(cond)), $text_func))
end
assert(::AbstractCPUDevice, cond::C, text::T) where {C, T} =
    isnothing(text) ? (Base.@assert cond()) : (Base.@assert cond() text())


"""
    @threaded [device] [coarsen=...] [block_size=...] for ... end

Device-flexible generalization of `Threads.@threads`, which distributes the
iterations of a for-loop across multiple threads, with the option to control
thread coarsening and GPU kernel configuration. Coarsening makes each thread
evaluate more than one iteration of the loop, which can improve performance by
reducing the runtime overhead of launching additional threads (though too much
coarsening worsens performance because it reduces parallelization). The `device`
is either inferred by calling `ClimaComms.device()`, or it can be specified
manually, with the following device-dependent behavior:

 - When `device` is a `CPUSingleThreaded()`, the loop is evaluated as-is. This
   avoids the runtime overhead of calling `Threads.@threads` with a single
   thread, and, when the device type is statically inferrable, it also avoids
   compilation overhead.

 - When `device` is a `CPUMultiThreaded()`, the loop is passed to
   `Threads.@threads`. This supports three different kinds of "schedulers" for
   determining how many iterations of the loop to evaluate in each thread:
     1) (default) a "dynamic" scheduler that changes the number of iterations as
        new threads are launched,
     2) a "static" scheduler that evaluates a fixed number of iterations per
        thread, and
     3) a "greedy" scheduler that uses a small number of threads, continuously
        evaluating iterations in each thread until the loop is completed (only
        available as of Julia 1.11).
   Setting `coarsen` to `:dynamic` or `:greedy` launches threads with those
   schedulers. Setting it to `:static` or an integer value launches threads with
   static scheduling (using `:static` is similar to using `1`, but slightly more
   performant). To read more about multi-threading, see the documentation for
   [`Threads.@threads`](https://docs.julialang.org/en/v1/base/multi-threading/#Base.Threads.@threads).

 - When `device` is a `CUDADevice()`, the loop is compiled with `CUDA.@cuda` and
   run with `CUDA.@sync`. Since CUDA launches all threads at the same time, only
   static scheduling can be used. Setting `coarsen` to any symbol causes each
   thread to evaluate a single iteration (default), and setting it to an integer
   value causes each thread to evaluate that number of iterations (the default
   is similar to using `1`, but slightly more performant). If the total number
   of iterations in the loop is extremely large, the specified coarsening may
   require more threads than can be simultaneously launched on the GPU, in which
   case the amount of coarsening is automatically increased.

   The optional argument `block_size` is also available for manually specifying
   the size of each block on a GPU. The default value of `:auto` sets the number
   of threads in each block to the largest possible value that permits a high
   GPU "occupancy" (the number of active thread warps in each multiprocessor
   executing the kernel). An integer can be used instead of `:auto` to override
   this default value. If the specified value exceeds the total number of
   threads, it is automatically decreased to avoid idle threads.

Any iterator with methods for `firstindex`, `length`, and `getindex` can be used
in a `@threaded` loop. All lazy iterators from `Base` and `Base.Iterators`, such
as `zip`, `enumerate`, `Iterators.product`, and generator expressions, are also
compatible with `@threaded`. (Although these iterators do not define methods for
`getindex`, they are automatically modified by `threadable` to support
`getindex`.) Using multiple iterators with `@threaded` is equivalent to looping
over a single `Iterators.product`, with the innermost iterator of the loop
appearing first in the product, and the outermost iterator appearing last.

NOTE: When a value in the body of the loop has a type that cannot be inferred by
the compiler, an `InvalidIRError` will be thrown during compilation for a
`CUDADevice()`. In particular, global variables are not inferrable, so
`@threaded` must be wrapped in a function whenever it is used in the REPL:

```julia-repl
julia> a = CUDA.CuArray{Int}(undef, 100); b = similar(a);

julia> threaded_copyto!(a, b) = ClimaComms.@threaded for i in axes(a, 1)
           a[i] = b[i]
       end
threaded_copyto! (generic function with 1 method)

julia> threaded_copyto!(a, b)

julia> ClimaComms.@threaded for i in axes(a, 1)
           a[i] = b[i]
       end
ERROR: InvalidIRError: ...
```

Moreover, type variables are not inferrable across function boundaries, so types
used in a threaded loop cannot be precomputed before the loop:

```julia-repl
julia> threaded_add_epsilon!(a) = ClimaComms.@threaded for i in axes(a, 1)
           FT = eltype(a)
           a[i] += eps(FT)
       end
threaded_add_epsilon! (generic function with 1 method)

julia> threaded_add_epsilon!(a)

julia> function threaded_add_epsilon!(a)
           FT = eltype(a)
           ClimaComms.@threaded for i in axes(a, 1)
               a[i] += eps(FT)
           end
       end
threaded_add_epsilon! (generic function with 1 method)

julia> threaded_add_epsilon!(a)
ERROR: InvalidIRError: ...
```

To fix other kinds of inference issues on GPUs, especially ones brought about by
indexing into iterators with nonuniform element types, see
[`UnrolledUtilities.jl`](https://github.com/CliMA/UnrolledUtilities.jl).
"""
macro threaded(args...)
    usage_string = "Usage: @threaded [device] [coarsen=...] [block_size=...] for ... end"
    (device_and_kwarg_exprs..., loop_block_expr) = args

    if all(expr -> Meta.isexpr(expr, :(=)), device_and_kwarg_exprs)
        device_expr = :($ClimaComms.device())
        kwarg_exprs = device_and_kwarg_exprs
    elseif all(expr -> Meta.isexpr(expr, :(=)), device_and_kwarg_exprs[2:end])
        (device_expr, kwarg_exprs...) = device_and_kwarg_exprs
    else
        throw(ArgumentError(usage_string))
    end
    for kwarg_expr in kwarg_exprs
        if kwarg_expr.args[2] isa QuoteNode
            kwarg_expr.args[2] = :(Val($(kwarg_expr.args[2])))
        end
    end

    loop_expr = if Meta.isexpr(loop_block_expr, :for)
        loop_block_expr
    elseif (
        Meta.isexpr(loop_block_expr, :block) &&
        length(loop_block_expr.args) == 2 &&
        Meta.isexpr(loop_block_expr.args[2], :for)
    )
        loop_block_expr.args[2]
    else
        throw(ArgumentError(usage_string))
    end
    (var_and_itr_expr, loop_body) = loop_expr.args
    if Meta.isexpr(var_and_itr_expr, :(=))
        var_exprs = Any[var_and_itr_expr.args[1]]
        itr_exprs = Any[var_and_itr_expr.args[2]]
    elseif (
        Meta.isexpr(var_and_itr_expr, :block) &&
        all(expr -> Meta.isexpr(expr, :(=)), var_and_itr_expr.args)
    )
        # Reverse the order of iterators to match standard for-loop behavior
        # (for-loops go from outermost iterator to innermost iterator, whereas
        # Iterators.product goes from innermost iterator to outermost iterator).
        var_exprs =
            Base.mapany(expr -> expr.args[1], reverse(var_and_itr_expr.args))
        itr_exprs =
            Base.mapany(expr -> expr.args[2], reverse(var_and_itr_expr.args))
    else
        throw(ArgumentError("Invalid for-loop expression: $var_and_itr_expr"))
    end

    # Improve latency by inlining the loop on single-threaded CPU devices.
    return quote
        device = $(esc(device_expr))
        device isa CPUSingleThreaded ? $(esc(loop_expr)) :
        threaded(
            ($(Base.mapany(esc, var_exprs)...),) -> $(esc(loop_body)),
            device,
            $(Base.mapany(esc, itr_exprs)...);
            $(Base.mapany(esc, kwarg_exprs)...),
        )
    end
end

"""
    threaded(f, device, itrs...; kwargs...)

Functional form of `@threaded`. If there are `n` iterators and `f` is a function
of `n` arguments, the `threaded` function is similar to

```julia
@threaded device [kwargs...] for xₙ in itrs[n], ..., x₂ in itrs[2], x₁ in itrs[1]
    f(x₁, x₂, ..., xₙ)
end
```

On single-threaded CPU devices, the `@threaded` macro inlines the for-loop
without any intermediate function calls, so that it has a lower latency than the
`threaded` function. On other devices, the only difference between the macro and
the function is that keyword argument symbols like `:dynamic` and `:auto` must
be wrapped in `Val`s for the function.
"""
threaded(f::F, device, itrs...; kwargs...) where {F} =
    threaded(splat(f), device, Iterators.product(itrs...); kwargs...)

threaded(
    f::F,
    device,
    itr;
    coarsen = Val(:dynamic),
    block_size = Val(:auto),
) where {F} =
    run_threaded(f, device, coarsen, threadable(device, itr); block_size)

run_threaded(f::F, device::AbstractCPUDevice, coarsen, itr; _...) where {F} =
    run_threaded(f, device, coarsen, itr)

run_threaded(f::F, ::CPUSingleThreaded, _, itr) where {F} =
    for item in itr
        f(item)
    end

run_threaded(f::F, ::CPUMultiThreaded, ::Val{:dynamic}, itr) where {F} =
    Threads.@threads :dynamic for item in itr
        f(item)
    end

run_threaded(f::F, ::CPUMultiThreaded, ::Val{:static}, itr) where {F} =
    Threads.@threads :static for item in itr
        f(item)
    end

@static if VERSION >= v"1.11"
    run_threaded(f::F, ::CPUMultiThreaded, ::Val{:greedy}, itr) where {F} =
        Threads.@threads :greedy for item in itr
            f(item)
        end
end

function run_threaded(
    f::F,
    ::CPUMultiThreaded,
    items_in_thread::Int,
    itr,
) where {F}
    items_in_thread > 0 || throw(ArgumentError("`coarsen` is not positive"))
    n_items = length(itr)
    Threads.@threads :static for thread_index in 1:cld(n_items, items_in_thread)
        first_item_index = items_in_thread * (thread_index - 1) + 1
        last_item_index = items_in_thread * thread_index
        for item_index in first_item_index:min(last_item_index, n_items)
            @inbounds f(itr[firstindex(itr) + item_index - 1])
        end
    end
end

"""
    threadable(device, itr)

Modifies an iterator to ensure that it can be used in a `@threaded` loop. This
will typically return `itr` or a `ThreadableWrapper` of `itr`.
"""
threadable(_, itr) = itr

"""
    ThreadableWrapper

Wrapper for an iterator from `Base` or `Iterators` that can be used in
`@threaded`, with methods for `firstindex`, `length`, and `getindex`. The
`getindex` method only supports linear indices between `firstindex` and
`firstindex + length - 1`. For the `ThreadableWrapper` of `Iterators.product`,
`getindex` converts each linear index to a Cartesian index using regular integer
division on CPUs and `Base.multiplicativeinverse` on GPUs.
"""
abstract type ThreadableWrapper end
Base.firstindex(::ThreadableWrapper) = 1
Base.iterate(wrapper::ThreadableWrapper, state = 1) =
    state > length(wrapper) ? nothing : (wrapper[state], state + 1)

struct ThreadableGenerator{F, I} <: ThreadableWrapper
    f::F
    itr::I
end
Adapt.@adapt_structure ThreadableGenerator
threadable(device, (; f, iter)::Base.Generator) =
    ThreadableGenerator(f, threadable(device, iter))
Base.length((; itr)::ThreadableGenerator) = length(itr)
Base.@propagate_inbounds Base.getindex(
    (; f, itr)::ThreadableGenerator,
    i::Integer,
) = f(itr[firstindex(itr) + i - 1])

struct ThreadableReverse{I} <: ThreadableWrapper
    itr::I
end
Adapt.@adapt_structure ThreadableReverse
threadable(device, (; itr)::Iterators.Reverse) =
    ThreadableReverse(threadable(device, itr))
Base.length((; itr)::ThreadableReverse) = length(itr)
Base.@propagate_inbounds Base.getindex((; itr)::ThreadableReverse, i::Integer) =
    itr[firstindex(itr) + length(itr) - i]

struct ThreadableEnumerate{I} <: ThreadableWrapper
    itr::I
end
Adapt.@adapt_structure ThreadableEnumerate
threadable(device, (; itr)::Iterators.Enumerate) =
    ThreadableEnumerate(threadable(device, itr))
Base.length((; itr)::ThreadableEnumerate) = length(itr)
Base.@propagate_inbounds Base.getindex(
    (; itr)::ThreadableEnumerate,
    i::Integer,
) = (i, itr[firstindex(itr) + i - 1])

struct ThreadableZip{I} <: ThreadableWrapper
    itrs::I
end
Adapt.@adapt_structure ThreadableZip
threadable(device, (; is)::Iterators.Zip) =
    ThreadableZip(map(itr -> threadable(device, itr), is))
Base.length((; itrs)::ThreadableZip) = minimum(length, itrs)
Base.@propagate_inbounds Base.getindex((; itrs)::ThreadableZip, i::Integer) =
    map(itr -> itr[firstindex(itr) + i - 1], itrs)

struct ThreadableProduct{I, D} <: ThreadableWrapper
    itrs::I
    divisors::D
end
Adapt.@adapt_structure ThreadableProduct
function threadable(device, (; iterators)::Iterators.ProductIterator)
    itrs = map(itr -> threadable(device, itr), iterators)
    divisors = reverse(cumprod((1, map(length, itrs[1:(end - 1)])...)))
    device_optimized_divisors =
        device isa AbstractCPUDevice ? divisors :
        map(Base.multiplicativeinverse, divisors) # Avoid Int division on GPUs.
    return ThreadableProduct(itrs, device_optimized_divisors)
end
Base.length((; itrs)::ThreadableProduct) = prod(length, itrs)
Base.@propagate_inbounds function Base.getindex(
    (; itrs, divisors)::ThreadableProduct,
    i::Integer,
)
    reversed_offset_remainder_pairs =
        accumulate(divisors; init = (0, i - 1)) do (_, remainder), divisor
            divrem(remainder, divisor)
        end
    offsets = map(first, reverse(reversed_offset_remainder_pairs))
    return map((itr, offset) -> itr[firstindex(itr) + offset], itrs, offsets)
end

struct ThreadableFlatten{I, O} <: ThreadableWrapper
    itrs::I
    offsets::O
end
Adapt.@adapt_structure ThreadableFlatten
function threadable(device, (; it)::Iterators.Flatten)
    itrs = map(itr -> threadable(device, itr), it)
    offsets = cumsum((0, map(length, itrs[1:(end - 1)])...))
    return ThreadableFlatten(itrs, offsets)
end
Base.length((; itrs)::ThreadableFlatten) = sum(length, itrs)
Base.@propagate_inbounds function Base.getindex(
    (; itrs, offsets)::ThreadableFlatten,
    i::Integer,
)
    # TODO: Implement binary search (searchsortedfirst doesn't work for tuples).
    n = 0
    while n < length(offsets) && i > offsets[n + 1]
        n += 1
    end
    return itrs[n][firstindex(itrs[n]) - offsets[n] + i - 1]
end

struct ThreadablePartition{I} <: ThreadableWrapper
    itr::I
    partition_size::Int
end
Adapt.@adapt_structure ThreadablePartition
threadable(device, (; c, n)::Iterators.PartitionIterator) =
    ThreadablePartition(threadable(device, c), n)
Base.length((; itr, partition_size)::ThreadablePartition) =
    cld(length(itr), partition_size)
Base.@propagate_inbounds function Base.getindex(
    (; itr, partition_size)::ThreadablePartition,
    i::Integer,
)
    first_index_in_partition = firstindex(itr) + (i - 1) * partition_size
    last_index_in_partition =
        firstindex(itr) + min(i * partition_size, length(itr)) - 1
    return (itr[i] for i in first_index_in_partition:last_index_in_partition)
end

# TODO: Check whether conversion of every Int to Int32 improves GPU performance.
