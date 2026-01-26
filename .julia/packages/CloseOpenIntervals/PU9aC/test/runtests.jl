using CloseOpenIntervals
using Test

function mysum2(X = CartesianIndices((SafeCloseOpen(10), SafeCloseOpen(10))))
  s = 0
  for I in X
    s += sum(Tuple(I))
  end
  s
end
mysum2_allocated() = @allocated(mysum2())
const ArrayInterface = CloseOpenIntervals.StaticArrayInterface
@testset "CloseOpenIntervals.jl" begin
  function mysum(x, N)
    s = zero(eltype(x))
    @inbounds @fastmath for i ∈ CloseOpen(N)
      s += x[i+1]
    end
    s
  end
  function mysafesum(x, N)
    s = zero(eltype(x))
    @inbounds @fastmath for i ∈ SafeCloseOpen(N)
      s += x[i+1]
    end
    s
  end
  x = rand(128)
  for n ∈ 1:128
    @test mysum(x, n) ≈ mysafesum(x, n) ≈ sum(view(x, 1:n))
    @test length(CloseOpen(n)) == n
  end
  @test @inferred(mysum(x, 0)) == first(x)
  @test @inferred(mysafesum(x, 0)) == 0.0
  @test @inferred(mysum(x, CloseOpenIntervals.StaticInt(128))) ≈ sum(x)
  @test @inferred(
    ArrayInterface.static_length(CloseOpen(CloseOpenIntervals.StaticInt(128)))
  ) === CloseOpenIntervals.StaticInt(128)
  @test @inferred(eltype(CloseOpen(7))) === Int
  @test ArrayInterface.known_length(CloseOpen(CloseOpenIntervals.StaticInt(128))) == 128
  @test @inferred(mysum2()) == sum(0:9) * 2 * length(0:9)
  @test mysum2_allocated() == 0
end
