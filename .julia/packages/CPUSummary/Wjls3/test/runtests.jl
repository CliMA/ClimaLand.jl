using CPUSummary
using Test

@testset "CPUSummary.jl" begin
  @test @inferred(CPUSummary.sys_threads()) == Sys.CPU_THREADS::Int
end
