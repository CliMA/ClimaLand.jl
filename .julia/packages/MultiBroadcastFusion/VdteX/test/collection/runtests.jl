#=
using Revise; include(joinpath("test", "collection", "runtests.jl"))
=#
using Test
using SafeTestsets

#! format: off
@safetestset "expr_code_lowered_single_expression" begin; @time include("expr_code_lowered_single_expression.jl"); end
@safetestset "expr_materialize_args" begin; @time include("expr_materialize_args.jl"); end
@safetestset "expr_fused_direct" begin; @time include("expr_fused_direct.jl"); end
@safetestset "expr_fused_assemble" begin; @time include("expr_fused_assemble.jl"); end
@safetestset "expr_errors_and_edge_cases" begin; @time include("expr_errors_and_edge_cases.jl"); end
#! format: on
