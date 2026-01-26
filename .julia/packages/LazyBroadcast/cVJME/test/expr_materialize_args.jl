#=
using Revise; include(joinpath("test", "expr_materialize_args.jl"))
=#
using Test
import LazyBroadcast as LB

@testset "materialize_args" begin
    expr_in = :(Base.materialize!(y1, Base.broadcasted(+, x1, x2, x3, x4)))
    tuple_out = (:(y1), :(Base.broadcasted(+, x1, x2, x3, x4)))
    @test LB.materialize_args(expr_in) == tuple_out

    expr_in = :(Base.materialize(Base.broadcasted(+, x1, x2, x3, x4)))
    entry_out = :(Base.broadcasted(+, x1, x2, x3, x4))
    @test LB.materialize_args(expr_in) == (entry_out, entry_out)

    expr_in = :(foo(Base.broadcasted(+, x1, x2, x3, x4)))
    @test expr_in == :(foo(Base.broadcasted(+, x1, x2, x3, x4)))
end
