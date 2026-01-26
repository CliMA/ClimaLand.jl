#=
using Revise; include(joinpath("test", "collection", "expr_fused_direct.jl"))
=#
using Test
import MultiBroadcastFusion as MBF

@testset "fused_direct" begin
    expr_in = quote
        @. y1 = x1 + x2 + x3 + x4
        @. y2 = x2 + x3 + x4 + x5
    end

    expr_out = :(tuple(
        Pair(y1, Base.broadcasted(+, x1, x2, x3, x4)),
        Pair(y2, Base.broadcasted(+, x2, x3, x4, x5)),
    ))
    @test MBF.fused_direct(expr_in) == expr_out
end
