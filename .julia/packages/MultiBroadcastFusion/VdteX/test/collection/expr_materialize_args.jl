#=
using Revise; include(joinpath("test", "collection", "expr_materialize_args.jl"))
=#
using Test
import MultiBroadcastFusion as MBF

@testset "materialize_args" begin
    expr_in = :(Base.materialize!(y1, Base.broadcasted(+, x1, x2, x3, x4)))
    tuple_out = (:(y1), :(Base.broadcasted(+, x1, x2, x3, x4)))
    @test MBF.materialize_args(expr_in) == tuple_out
end
