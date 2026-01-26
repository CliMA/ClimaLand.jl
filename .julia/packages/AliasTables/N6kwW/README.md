# AliasTables

[![docs](https://img.shields.io/badge/docs-v1-blue.svg)](https://aliastables.lilithhafner.com/dev)
[![Build Status](https://github.com/LilithHafner/AliasTables.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/LilithHafner/AliasTables.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/LilithHafner/AliasTables.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/LilithHafner/AliasTables.jl)
[![PkgEval](https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/A/AliasTables.svg)](https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/A/AliasTables.html)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

AliasTables provides the `AliasTable` type, which is an object that defines a probability
distribution over `1:n` for some `n`. They are efficient to construct and very efficient to
sample from.

An alias table can be combined with a dense vector of values to create a discrete
distribution over anything.

Internally, AliasTables define a mapping from an unsigned integer type to the sampling
domain. To get a random sample according to the AliasTable's distribution, one must provide
a random unsigned integer uniformly at random. One can also provide a `Random.AbstractRNG`
object instead and a random unsigned integer will be generated using that rng. When using
the random API, this latter approach is taken.

```julia
julia> using AliasTables

julia> at = AliasTable([5,10,1])
AliasTable([0x5000000000000000, 0xa000000000000000, 0x1000000000000000])

julia> rand(at, 10)
10-element Vector{Int64}:
 2
 1
 2
 2
 2
 2
 1
 1
 3
 2

julia> using Chairmarks

julia> @b at rand
2.898 ns

julia> @b rand(UInt)
2.738 ns

julia> @b rand(1000) AliasTable
9.167 μs (2 allocs: 16.031 KiB)

julia> @b AliasTable(rand(1000)) rand(_, 1000)
1.506 μs (3 allocs: 7.875 KiB)

julia> @b AliasTable(rand(1000)), rand(1000) AliasTables.set_weights!(_...)
8.427 μs

julia> using StatsBase

julia> at = AliasTable{UInt16}([5,10,1])
AliasTable{UInt16}([0x5000, 0xa000, 0x1000])

julia> countmap(AliasTables.sample(x, at) for x in typemin(UInt16):typemax(UInt16))
Dict{Any, Int64} with 3 entries:
  2 => 40960
  3 => 4096
  1 => 20480
```
