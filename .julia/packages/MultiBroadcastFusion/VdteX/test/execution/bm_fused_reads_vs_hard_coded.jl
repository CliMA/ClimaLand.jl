#=
using Revise; include(joinpath("test", "execution", "bm_fused_reads_vs_hard_coded.jl"))
=#
include("utils_test.jl")
include("utils_setup.jl")
include("utils_benchmark.jl")

import MultiBroadcastFusion as MBF

# =========================================== hard-coded implementations
perf_kernel_hard_coded!(X, Y) = perf_kernel_hard_coded!(X, Y, MBF.device(X.x1))

function perf_kernel_hard_coded!(X, Y, ::MBF.MBF_CPU)
    (; x1, x2, x3, x4) = X
    (; y1, y2, y3, y4) = Y
    @inbounds for i in eachindex(x1)
        y1[i] = x1[i] + x2[i] + x3[i] + x4[i]
        y2[i] = x1[i] * x2[i] * x3[i] * x4[i]
        y3[i] = x1[i] + x2[i] - x3[i] + x4[i]
        y4[i] = x1[i] * x2[i] + x3[i] * x4[i]
    end
end

@static get(ENV, "USE_CUDA", nothing) == "true" && using CUDA
use_cuda = @isdefined(CUDA) && CUDA.has_cuda() # will be true if you first run `using CUDA`
@static if use_cuda
    function perf_kernel_hard_coded!(X, Y, ::MBF.MBF_CUDA)
        x1 = X.x1
        nitems = length(parent(x1))
        max_threads = 256 # can be higher if conditions permit
        nthreads = min(max_threads, nitems)
        nblocks = cld(nitems, nthreads)
        CUDA.@cuda threads = (nthreads) blocks = (nblocks) knl_multi_copyto_hard_coded!(
            X,
            Y,
            Val(nitems),
        )
    end
    function knl_multi_copyto_hard_coded!(X, Y, ::Val{nitems}) where {nitems}
        (; x1, x2, x3, x4) = X
        (; y1, y2, y3, y4) = Y
        i =
            CUDA.threadIdx().x +
            (CUDA.blockIdx().x - Int32(1)) * CUDA.blockDim().x
        @inbounds begin
            if i ≤ nitems
                y1[i] = x1[i] + x2[i] + x3[i] + x4[i]
                y2[i] = x1[i] * x2[i] * x3[i] * x4[i]
                y3[i] = x1[i] + x2[i] - x3[i] + x4[i]
                y4[i] = x1[i] * x2[i] + x3[i] * x4[i]
            end
        end
        return nothing
    end
end

# ===========================================

function perf_kernel_unfused!(X, Y)
    (; x1, x2, x3, x4) = X
    (; y1, y2, y3, y4) = Y
    # 4 writes; 4 unique reads
    # 4 writes; 16 reads including redundant ones
    @. y1 = x1 + x2 + x3 + x4
    @. y2 = x1 * x2 * x3 * x4
    @. y3 = x1 + x2 - x3 + x4
    @. y4 = x1 * x2 + x3 * x4
    return nothing
end

function perf_kernel_fused!(X, Y)
    (; x1, x2, x3, x4) = X
    (; y1, y2, y3, y4) = Y
    MBF.@fused_direct begin
        @. y1 = x1 + x2 + x3 + x4
        @. y2 = x1 * x2 * x3 * x4
        @. y3 = x1 + x2 - x3 + x4
        @. y4 = x1 * x2 + x3 * x4
    end
end

@static get(ENV, "USE_CUDA", nothing) == "true" && using CUDA
use_cuda = @isdefined(CUDA) && CUDA.has_cuda() # will be true if you first run `using CUDA`
AType = use_cuda ? CUDA.CuArray : Array
device_name = use_cuda ? CUDA.name(CUDA.device()) : "CPU"
bm = Benchmark(; device_name, float_type = Float32)

problem_size = (50, 5, 5, 6, 5400)

array_size = problem_size # array
X = get_arrays(:x, AType, bm.float_type, array_size)
Y = get_arrays(:y, AType, bm.float_type, array_size)
test_kernel!(
    use_cuda;
    fused! = perf_kernel_fused!,
    unfused! = perf_kernel_unfused!,
    X,
    Y,
)
use_cuda && test_kernel!(
    use_cuda;
    fused! = perf_kernel_hard_coded!,
    unfused! = perf_kernel_unfused!,
    X,
    Y,
)
push_benchmark!(
    bm,
    use_cuda,
    perf_kernel_unfused!,
    X,
    Y;
    n_reads_writes = 4 + 4,
    problem_size = array_size,
)
push_benchmark!(
    bm,
    use_cuda,
    perf_kernel_fused!,
    X,
    Y;
    n_reads_writes = 4 + 4,
    problem_size = array_size,
)
use_cuda && push_benchmark!(
    bm,
    use_cuda,
    perf_kernel_hard_coded!,
    X,
    Y;
    n_reads_writes = 4 + 4,
    problem_size = array_size,
)

array_size = (prod(problem_size),) # vector
X = get_arrays(:x, AType, bm.float_type, array_size)
Y = get_arrays(:y, AType, bm.float_type, array_size)
test_kernel!(
    use_cuda;
    fused! = perf_kernel_fused!,
    unfused! = perf_kernel_unfused!,
    X,
    Y,
)
use_cuda && test_kernel!(
    use_cuda;
    fused! = perf_kernel_hard_coded!,
    unfused! = perf_kernel_unfused!,
    X,
    Y,
)
push_benchmark!(
    bm,
    use_cuda,
    perf_kernel_unfused!,
    X,
    Y;
    n_reads_writes = 4 + 4,
    problem_size = array_size,
)
push_benchmark!(
    bm,
    use_cuda,
    perf_kernel_fused!,
    X,
    Y;
    n_reads_writes = 4 + 4,
    problem_size = array_size,
)
use_cuda && push_benchmark!(
    bm,
    use_cuda,
    perf_kernel_hard_coded!,
    X,
    Y;
    n_reads_writes = 4 + 4,
    problem_size = array_size,
)


tabulate_benchmark(bm)

nothing
