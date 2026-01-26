using LinearAlgebra
using SparseArrays
using Test
using Random

using AMD
import AMD: SS_Int, AMD_STATUS, AMD_OK, COLAMD_STATUS, COLAMD_OK

Random.seed!(666)

for n in [10, 20, 30]
  for density in [0.25, 0.75, 1.0]
    for T in [Cint, SS_Int]
      A = convert(SparseMatrixCSC{Float64, T}, sprand(n, n, density))
      @test amd_valid(A)

      meta = Amd()
      p = amd(A, meta)
      @test meta.info[AMD_STATUS] == AMD_OK
      @test minimum(p) == 1
      @test maximum(p) == n
      @test isperm(p)

      q = amd(A)
      @test all(p .== q)

      # test that amd(A) with A unsymmetric is amd(A + A')
      A = A + A'
      pA = amd(A)
      @test minimum(pA) == 1
      @test maximum(pA) == n
      @test isperm(pA)
    end
  end
end

for n in [10, 20, 30]
  for m in [10, 20, 30]
    for density in [0.25, 0.75, 1.0]
      for T in [Cint, SS_Int]
        A = convert(SparseMatrixCSC{Float64, T}, sprand(n, m, density))

        meta = Colamd{T}()
        p = colamd(A, meta)
        @test meta.stats[COLAMD_STATUS] == COLAMD_OK
        @test minimum(p) == 1
        @test maximum(p) == m
        @test isperm(p)
        q = colamd(A)
        @test all(p .== q)

        if n == m
          A = A * A'
          meta = Colamd{T}()
          p = symamd(A, meta)
          @test meta.stats[COLAMD_STATUS] == COLAMD_OK
          @test minimum(p) == 1
          @test maximum(p) == m
          @test isperm(p)

          q = symamd(A)
          @test all(p .== q)
        end
      end
    end
  end
end

# Test Symmetric and Hermitian wrappers
for ordering in (:amd, :symamd, :colamd)
  for wrapper in (:Symmetric, :Hermitian)
    for T in [Cint, SS_Int]
      @eval begin
        A = convert(SparseMatrixCSC{Float64, $T}, rand(10, 10))
        A = A * A'
        p = $ordering(A)

        L = $wrapper(A, :L)
        pL = $ordering(L)
        @test all(p .== pL)

        U = $wrapper(A, :U)
        pU = $ordering(U)
        @test all(p .== pU)
      end
    end
  end
end

# For coverage.
meta = Amd()
A = convert(SparseMatrixCSC{Float64, Cint}, sprand(10, 10, 0.5))
p = amd(A, meta)
show(meta)
print(meta)

meta = Colamd{Cint}()
p = colamd(A, meta)
show(meta)
print(meta)
