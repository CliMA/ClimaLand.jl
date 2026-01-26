#=
using Revise; include(joinpath("test", "execution", "runtests.jl"))
=#

#! format: off
@safetestset "fused_shared_reads" begin; @time include("bm_fused_shared_reads.jl"); end
@safetestset "fused_shared_reads_writes" begin; @time include("bm_fused_shared_reads_writes.jl"); end
@safetestset "bm_fused_reads_vs_hard_coded" begin; @time include("bm_fused_reads_vs_hard_coded.jl"); end
@safetestset "measure_parameter_memory" begin; @time include("measure_parameter_memory.jl"); end
@safetestset "kernel_splitting" begin; @time include("kernel_splitting.jl"); end
#! format: on
