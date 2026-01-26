using ManualMemory: MemoryBuffer, load, store!, LazyPreserve, preserve, PseudoPtr, Reference
using Test

@testset "ManualMemory.jl" begin

  @test_throws AssertionError MemoryBuffer{4,String}(undef)
  m = MemoryBuffer{4,Float64}(undef);
  store!(pointer(m), 1.23)
  @test load(pointer(m)) == 1.23
  str = "Hello world!"
  GC.@preserve str begin
    store!(Base.unsafe_convert(Ptr{String}, m), str)
    @test load(Base.unsafe_convert(Ptr{String}, m)) === str
  end

  x = [0 0; 0 0];
  preserve(store!, LazyPreserve(x), 1)
  @test x[1] === 1
  p = PseudoPtr(x, 1)
  @test load(p) === 1
  p += 1
  store!(p, 2)
  @test load(p) === 2
  p = 1 + p
  store!(p, 3)
  @test load(p) === 3

  x = Reference{Int}()
  y = Reference(1)
  GC.@preserve x y begin
    store!(pointer(x), 1)
    @test load(pointer(x)) === 1 === load(pointer(y))
  end

  # Test construction with existing data
  s = (1,2,3,4,5)
  @test MemoryBuffer(s).data === s

  # Test equality and inequality
  @test MemoryBuffer((1,2,3,4,5)) == MemoryBuffer((1,2,3,4,5))
  @test MemoryBuffer((1,2,3,4,6)) != MemoryBuffer((1,2,3,4,5))
  @test MemoryBuffer((0x01,0x02)) == MemoryBuffer((1,2))
  @test MemoryBuffer((0x01,0x02)) != MemoryBuffer((0x01,0x02,0x03))

end

using ThreadingUtilities
include(joinpath(pkgdir(ThreadingUtilities), "test", "runtests.jl"))
