#=
julia --project
using Revise; include(joinpath("test", "runtests.jl"))
=#
using Test
using SafeTestsets

#! format: off
@safetestset "expr_code_lowered_single_expression" begin; @time include("expr_code_lowered_single_expression.jl"); end
@safetestset "lazy_broadcast" begin; @time include("lazy_broadcast.jl"); end
#! format: on
