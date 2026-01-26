#=
julia --project
using Revise; include(joinpath("test", "lazy_broadcast.jl"))
=#
using Test
import LazyBroadcast as LB

a = rand(3, 3)
b = rand(3, 3)

@testset "lazy_broadcast" begin
    bco = LB.@lazy_broadcast @. a + b
    @test bco == Base.Broadcast.instantiate(Base.broadcasted(+, a, b))

    bco = LB.@lazy_broadcast a .+ b
    @test bco == Base.Broadcast.instantiate(Base.broadcasted(+, a, b))

    bco = LB.lazy_broadcast.(a .+ b)
    @test bco == Base.Broadcast.instantiate(Base.broadcasted(+, a, b))
end
