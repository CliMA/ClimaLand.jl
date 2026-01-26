using FunctionWrappersWrappers
using Test

@testset "FunctionWrappersWrappers.jl" begin
  
  fwplus = FunctionWrappersWrapper(+, (Tuple{Float64,Float64}, Tuple{Int,Int}), (Float64,Int));
  @test fwplus(4.0, 8.0) === 12.0
  @test fwplus(4, 8) === 12

  fwexp2 = FunctionWrappersWrapper(exp2, (Tuple{Float64}, Tuple{Float32}, Tuple{Int}), (Float64,Float32,Float64));
  @test fwexp2(4.0) === 16.0
  @test fwexp2(4f0) === 16f0
  @test fwexp2(4) === 16.0

end
