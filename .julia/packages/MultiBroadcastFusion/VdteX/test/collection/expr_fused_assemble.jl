#=
using Revise; include(joinpath("test", "collection", "expr_fused_assemble.jl"))
=#
using Test
import MultiBroadcastFusion as MBF

@testset "fused_assemble - simple sequential" begin
    expr_in = quote
        @. y1 = x1 + x2 + x3 + x4
        @. y2 = x2 + x3 + x4 + x5
    end

    expr_out = quote
        tup = ()
        tup = (tup..., Pair(y1, Base.broadcasted(+, x1, x2, x3, x4)))
        tup = (tup..., Pair(y2, Base.broadcasted(+, x2, x3, x4, x5)))
        tup
    end

    @test MBF.linefilter!(MBF.fused_assemble(expr_in, :tup)) ==
          MBF.linefilter!(expr_out)
    @test MBF.fused_assemble(expr_in, :tup) == expr_out
end


@testset "fused_assemble - loop" begin
    expr_in = quote
        for i in 1:10
            @. y1 = x1 + x2 + x3 + x4
            @. y2 = x2 + x3 + x4 + x5
        end
    end

    expr_out = quote
        tup = ()
        for i in 1:10
            tup = (tup..., Pair(y1, Base.broadcasted(+, x1, x2, x3, x4)))
            tup = (tup..., Pair(y2, Base.broadcasted(+, x2, x3, x4, x5)))
        end
        tup
    end

    @test MBF.linefilter!(MBF.fused_assemble(expr_in, :tup)) ==
          MBF.linefilter!(expr_out)
    @test MBF.fused_assemble(expr_in, :tup) == expr_out
end

@testset "fused_assemble - loop with @inbounds" begin
    expr_in = quote
        @inbounds for i in 1:10
            @. y1 = x1 + x2 + x3 + x4
            @. y2 = x2 + x3 + x4 + x5
        end
    end

    expr_out = quote
        tup = ()
        @inbounds for i in 1:10
            tup = (tup..., Pair(y1, Base.broadcasted(+, x1, x2, x3, x4)))
            tup = (tup..., Pair(y2, Base.broadcasted(+, x2, x3, x4, x5)))
        end
        tup
    end
    @test MBF.linefilter!(MBF.fused_assemble(expr_in, :tup)) ==
          MBF.linefilter!(expr_out)
    @test MBF.fused_assemble(expr_in, :tup) == expr_out
end

@testset "fused_assemble - if" begin
    expr_in = quote
        if a && B || something(x, y, z)
            @. y1 = x1 + x2 + x3 + x4
            @. y2 = x2 + x3 + x4 + x5
        end
    end

    expr_out = quote
        tup = ()
        if a && B || something(x, y, z)
            tup = (tup..., Pair(y1, Base.broadcasted(+, x1, x2, x3, x4)))
            tup = (tup..., Pair(y2, Base.broadcasted(+, x2, x3, x4, x5)))
        end
        tup
    end
    @test MBF.linefilter!(MBF.fused_assemble(expr_in, :tup)) ==
          MBF.linefilter!(expr_out)
    @test MBF.fused_assemble(expr_in, :tup) == expr_out
end
