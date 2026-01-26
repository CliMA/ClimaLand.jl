#=
using TestEnv
TestEnv.activate()
using CUDA # (optional)
using Revise; include(joinpath("test", "execution", "measure_parameter_memory.jl"))
=#

include("utils_test.jl")
include("utils_setup.jl")
include("utils_benchmark.jl")

@static get(ENV, "USE_CUDA", nothing) == "true" && using CUDA
use_cuda = @isdefined(CUDA) && CUDA.has_cuda() # will be true if you first run `using CUDA`

import MultiBroadcastFusion as MBF
@static if use_cuda
    const MBFExt = Base.get_extension(MBF, :MultiBroadcastFusionCUDAExt)
end

@static if use_cuda
    AType = use_cuda ? CUDA.CuArray : Array
    device_name = use_cuda ? CUDA.name(CUDA.device()) : "CPU"
    bm = Benchmark(; device_name, float_type = Float32)
    problem_size = (50, 5, 5, 6, 5400)

    array_size = problem_size # array
    X = get_arrays(:x, AType, bm.float_type, array_size)
    Y = get_arrays(:y, AType, bm.float_type, array_size)

    function perf_kernel_shared_reads_fused!(X, Y)
        (; x1, x2, x3, x4) = X
        (; y1, y2, y3, y4) = Y
        fmb = MBF.@get_fused_direct begin
            @. y1 = x1
            @. y2 = x1 + x2
            @. y3 = x1 + x2 + x3
        end
        @test MBFExt.param_usage_args(MBFExt.fused_multibroadcast_args(fmb)) ==
              672
        MBFExt.@rprint_parameter_memory(fmb)
    end
    @testset "Test measuring parameter memory" begin
        perf_kernel_shared_reads_fused!(X, Y)
    end
end

nothing
