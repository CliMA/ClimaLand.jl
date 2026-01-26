using SIMDTypes
using Test

@testset "SIMDTypes.jl" begin
  @test SIMDTypes.Bit(true).data
  @test !SIMDTypes.Bit(false).data
  @test SIMDTypes.Bit(true) isa SIMDTypes.NativeTypes
  @test !(SIMDTypes.Bit(true) isa SIMDTypes.NativeTypesExceptBit)
  @test SIMDTypes.Bit(true) isa SIMDTypes.NativeTypesExceptFloat16
  @test !(SIMDTypes.Bit(true) isa SIMDTypes.NativeTypesExceptBitandFloat16)
  @test Float16(1) isa SIMDTypes.NativeTypes
  @test Float16(1) isa SIMDTypes.NativeTypesExceptBit
  @test !(Float16(1) isa SIMDTypes.NativeTypesExceptFloat16)
  @test !(Float16(1) isa SIMDTypes.NativeTypesExceptBitandFloat16)
  @test 1f0 isa SIMDTypes.NativeTypesExceptBitandFloat16
  
  @test 1.3 isa SIMDTypes.FloatingTypes
  # @test !(1.3 isa SIMDTypes.IntegerTypes)
  @test !(1.3 isa SIMDTypes.IntegerTypesHW)
  # @test 1 isa SIMDTypes.IntegerTypes
  @test 1 isa SIMDTypes.SignedHW
  @test !(1 isa SIMDTypes.UnsignedHW)
  @test 1 isa SIMDTypes.IntegerTypesHW
  # @test one(UInt) isa SIMDTypes.IntegerTypes
  @test !(one(UInt) isa SIMDTypes.SignedHW)
  @test one(UInt) isa SIMDTypes.UnsignedHW
  @test one(UInt) isa SIMDTypes.IntegerTypesHW
  

  @test ntuple(Core.VecElement, 2) isa SIMDTypes._Vec
end
