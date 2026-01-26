# ManualMemory

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaSIMD.github.io/ManualMemory.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaSIMD.github.io/ManualMemory.jl/dev)
[![Build Status](https://github.com/JuliaSIMD/ManualMemory.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/JuliaSIMD/ManualMemory.jl/actions/workflows/CI.yml)
[![Coverage](https://codecov.io/gh/JuliaSIMD/ManualMemory.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaSIMD/ManualMemory.jl)

Manually managed memory buffers backed by NTuples

### Examples

```julia
julia> using ManualMemory: MemoryBuffer, load, store!, LazyPreserve, preserve, PseudoPtr, Reference

julia> m = MemoryBuffer{4,Float64}(undef)
MemoryBuffer{4, Float64}((2.283825594e-314, 2.2157350003e-314, 2.216358792e-314, 2.08e-322))

julia> store!(pointer(m), 1.23)

julia> load(pointer(m))
1.23
```
Specifying an existing `NTuple` of data:
```julia
julia> s = (1,2,3,4,5);

julia> m = MemoryBuffer(s)
MemoryBuffer{5, Int64}((1, 2, 3, 4, 5))

julia> load(p)
1

julia> load(p+sizeof(Int64))
2

julia> load(p+sizeof(Int64)*2)
3
```
