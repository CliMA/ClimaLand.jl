# StrideArraysCore

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaSIMD.github.io/StrideArraysCore.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaSIMD.github.io/StrideArraysCore.jl/dev)
[![Build Status](https://github.com/JuliaSIMD/StrideArraysCore.jl/workflows/CI/badge.svg)](https://github.com/JuliaSIMD/StrideArraysCore.jl/actions)
[![Coverage](https://codecov.io/gh/JuliaSIMD/StrideArraysCore.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaSIMD/StrideArraysCore.jl)

Defines the core `PtrArray` type so that some libraries can make use of it internally without the need for circular dependencies.
[StrideArrays](https://github.com/chriselrod/StrideArrays.jl) extends this type with many methods and functionality. It is
recommended you depend on and use `StrideArrays` instead.


Example initialization:
```julia
julia> using StrideArraysCore

julia> A = StrideArray(undef, 2, 5) # `Float64` is default eltype
2×5 StrideArray{Tuple{Int64, Int64}, (true, true), Float64, 2, 1, 0, (1, 2), Tuple{StaticInt{8}, Int64}, Tuple{StaticInt{1}, StaticInt{1}}, Vector{Float64}}:
 1.5e-323  -2.0092e-237  4.94e-321   1.0e-322  1.6e-322
 1.6e-322   4.94e-322    1.235e-321  1.5e-323  2.0e-323

julia> B = StrideArray{Float32}(zero, 2, 5) # can specify initialization function; function must have 1-arg method accepting eltype as argument
2×5 StrideArray{Tuple{Int64, Int64}, (true, true), Float32, 2, 1, 0, (1, 2), Tuple{StaticInt{4}, Int64}, Tuple{StaticInt{1}, StaticInt{1}}, Matrix{Float32}}:
 0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0

julia> C = StrideArray{Float32}(one, static(2), 5) # sizes may optionally be static
2×5 StrideArray{Tuple{StaticInt{2}, Int64}, (true, true), Float32, 2, 1, 0, (1, 2), Tuple{StaticInt{4}, StaticInt{8}}, Tuple{StaticInt{1}, StaticInt{1}}, Vector{Float32}} with indices 1:1:2×Base.OneTo(5):
 1.0  1.0  1.0  1.0  1.0
 1.0  1.0  1.0  1.0  1.0

julia> D = StrideArray{Float32}(one, static(2), static(5)) # all sizes being static will allow the compiler to elide the allocation if the array does not escape.
2×5 StrideArraysCore.StaticStrideArray{Tuple{StaticInt{2}, StaticInt{5}}, (true, true), Float32, 2, 1, 0, (1, 2), Tuple{StaticInt{4}, StaticInt{8}}, Tuple{StaticInt{1}, StaticInt{1}}, 10} with indices 1:1:2×1:1:5:
 1.0  1.0  1.0  1.0  1.0
 1.0  1.0  1.0  1.0  1.0

julia> using StrideArraysCore, BenchmarkTools

julia> @inline function alloctest()
           D = StrideArray(one, static(2), static(5))
           s = 0.0
           for i in eachindex(D)
               s += D[i]
           end
           s
       end
alloctest (generic function with 1 method)

julia> @btime alloctest()
  1.214 ns (0 allocations: 0 bytes)
10.0

julia> @code_llvm debuginfo=:none alloctest() # compiler compiled-away function
define double @julia_alloctest_1199() #0 {
L67.9:
  ret double 1.000000e+01
}
```



