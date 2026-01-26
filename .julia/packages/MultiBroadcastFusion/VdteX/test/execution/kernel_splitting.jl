#=
using TestEnv
TestEnv.activate()
using CUDA # (optional)
using Revise; include(joinpath("test", "execution", "kernel_splitting.jl"))
=#

include("utils_test.jl")
include("utils_setup.jl")
include("utils_benchmark.jl")

import MultiBroadcastFusion as MBF

#! format: off
function perf_kernel_shared_reads_fused!(X, Y)
    (; x1, x2, x3, x4) = X
    (; y1, y2, y3, y4) = Y
    # TODO: can we write this more compactly with `@fused_assemble`?

    # Let's make sure that every broadcasted object is different,
    # so that we use up a lot of parameter memory:
    MBF.@fused_direct begin
        @. y1 = x1
        @. y2 = x1 + x2
        @. y3 = x1 + x2 + x3
        @. y4 = x1 * x2 + x3 + x4
        @. y1 = x1 * x2 + x3 + x4 + x1
        @. y2 = x1 * x2 + x3 + x4 + x1 + x2
        @. y3 = x1 * x2 + x3 + x4 + x1 + x2 + x3
        @. y4 = x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4
        @. y1 = x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1
        @. y2 = x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 + x2
        @. y3 = x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 + x2 + x3
        @. y4 = x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4
        @. y1 = x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1
        @. y2 = x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 + x2
        @. y3 = x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 + x2 + x3
        @. y4 = x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4
        @. y1 = x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1
        @. y2 = x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 + x2
        @. y3 = x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 + x2 + x3
        @. y4 = x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4
        @. y1 = x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1
        @. y2 = x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 + x2
        @. y3 = x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 + x2 + x3
        @. y4 = x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4
        @. y1 = x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1
        @. y2 = x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 + x2
        @. y3 = x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 + x2 + x3
        @. y4 = x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4
        @. y1 = x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1
        @. y2 = x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 * x2 + x3 + x4 + x1 + x2 # breaks on A100 due to too much parameter memory
    end
end
#! format: on

@static get(ENV, "USE_CUDA", nothing) == "true" && using CUDA
use_cuda = @isdefined(CUDA) && CUDA.has_cuda() # will be true if you first run `using CUDA`
AType = use_cuda ? CUDA.CuArray : Array
device_name = use_cuda ? CUDA.name(CUDA.device()) : "CPU"
bm = Benchmark(; device_name, float_type = Float32)
problem_size = (50, 5, 5, 6, 5400)

array_size = problem_size # array
X = get_arrays(:x, AType, bm.float_type, array_size)
Y = get_arrays(:y, AType, bm.float_type, array_size)
@testset "Test kernel splitting with too much parameter memory" begin
    use_cuda && perf_kernel_shared_reads_fused!(X, Y)
end

nothing
