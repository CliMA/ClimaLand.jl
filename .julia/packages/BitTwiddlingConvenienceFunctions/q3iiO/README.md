# BitTwiddlingConvenienceFunctions

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaSIMD.github.io/BitTwiddlingConvenienceFunctions.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaSIMD.github.io/BitTwiddlingConvenienceFunctions.jl/dev)
[![Build Status](https://github.com/JuliaSIMD/BitTwiddlingConvenienceFunctions.jl/workflows/CI/badge.svg)](https://github.com/JuliaSIMD/BitTwiddlingConvenienceFunctions.jl/actions)
[![Coverage](https://codecov.io/gh/JuliaSIMD/BitTwiddlingConvenienceFunctions.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaSIMD/BitTwiddlingConvenienceFunctions.jl)



Useful for going to next/previous mask size, or calculating corresponding shifts:
```julia
julia> using BitTwiddlingConvenienceFunctions: prevpow2, nextpow2, intlog2

julia> prevpow2.(7:9)'
1×3 adjoint(::Vector{Int64}) with eltype Int64:
 4  8  8

julia> nextpow2.(7:9)'
1×3 adjoint(::Vector{Int64}) with eltype Int64:
 8  8  16

julia> intlog2.(7:9)' # truncated
1×3 adjoint(::Vector{Int64}) with eltype Int64:
 2  3  3
 ```

