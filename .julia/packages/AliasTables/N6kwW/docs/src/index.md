```@meta
CurrentModule = AliasTables
```

# AliasTables

[AliasTables](https://github.com/LilithHafner/AliasTables.jl) provides the
[`AliasTable`](@ref) type, which is an object that defines a probability distribution over
`1:n` for some `n`. They are efficient to construct and very efficient to sample from.

## Constructing an AliasTable

Construct an AliasTable by calling [`AliasTable(probabilities)`](@ref) for some collection
of probabilities. For example, to create a table with a 30% chance of returning 1, and a
70% chance of returning 2, you would call `AliasTable([0.3, 0.7])`.

`probabilities` must be an abstract vector of real numbers. The sum need not be 1 as the
input will be automatically normalized.

See [`AliasTables.set_weights!`](@ref) for a way to change the weights of an existing alias
table without GC-managed allocations.

## Sampling from an AliasTable

Sample from an `AliasTable` the same way you would sample from any sampleable object using
the Random API. For example, to draw a single sample, call `rand(at::AliasTable)`, to draw
`n` samples, call `rand(at::AliasTable, n)`, to sample using a specific random number
generator, call `rand(rng::Random.AbstractRNG, at::AliasTable)`, and to populate an
existing array, call `rand!(x, at::AliasTable)`.

## Example

```jldoctest; filter=[r" [1-3]"]
julia> using AliasTables

julia> at = AliasTable([5,10,1])
AliasTable([0x5000000000000000, 0xa000000000000000, 0x1000000000000000])

julia> rand(at, 8)
8-element Vector{Int64}:
 2
 2
 1
 2
 3
 1
 2
 2

julia> length(at)
3

julia> AliasTables.probabilities(float, at)
3-element Vector{Float64}:
 0.3125
 0.625
 0.0625
```

## Implementation details

Alias tables are composed of a list of (acceptance probability, alias) pairs. To sample from
an alias table, first pick an element `(p, alias)` from that list uniformly at random. Then,
with probability `p`, return the index of that element and with probability `1-p`, return
`alias`. For more information, see the
[wikipedia article](https://en.wikipedia.org/wiki/Alias_method), or a publication by the
original author [Walker, A. J. "An Efficient
Method for Generating Discrete Random Variables with General Distributions." _ACM
Transactions on Mathematical Software_ 3 (3): 253, 1977.](https://lilithhafner.com/An-Efficient-Method-for-Generating-Discrete-Random-Variables-with-General-Distributions.pdf)

---

While this package does follow the general structure of the algorithms described in the
above articles, it makes some departures for performance reasons.

Conventional alias tables map an integer in the range `1:n` and a real number to a number in
the range `1:n`. This package's alias tables, however, map an integer in the range `0:2^k-1`
to a number in the range `1:n`. Where `k = 64` by default. This results in increased
precision and increased performance.

To apply this mapping to `x` using an alias table with `2^b` elements, we use the most
significant `b` bits of `x` to index into the table, retrieving a pair `(p, alias)`. We then
compare the least significant `k` bits of `x` to `p`, and if `x`'s low bits are less than
`p`, we return the aliased index, otherwise we return the index itself. Any finite
probability distribution can be thought of as a distribution over `2^b` elements by simply
appending zeros to the end of the distribution.

This whole process uses integer arithmetic which allows both very fast sampling and exact
construction.

We can count exactly how many inputs map to a given output as follows.

For a given output `m ∈ 1:n`, drawn from an `AliasTable{UIntK}` with a `k`-bit domain and a
range of `1:2^b`, the inputs that produce `m` come from two disjoint sets
- The integers between `m * 2^(k-b) + p`, inclusive, and `(m+1) * 2^(k-b)`, exclusive where
  `p` is the alias probability of the `m`th table entry. This set has `2^(k-b) - p`
  elements.
- The integers between `j * 2^(k-b)`, inclusive and `j + 2^(k-b) + p_j`, exclusive for all `j`
  whose table entry aliases to `m` where `p_j` is alias probability of the `j`th table
  entry. This set has size `sum(p_j for j in 1:n if alias_j == m)`.

The default constructors in this package utilize those formulae to produce alias tables that
can exactly represent any distribution where all probabilities are of the form `p/2^k` for
some integer `p`.

## Alternate sampling API

You can bypass the Random API and sample directly from an alias table `at::AliasTable{T}`
using the public [`AliasTables.sample`](@ref) function which is branchless, deterministic,
and not pseudorandom. If given an input drawn uniformly at random from the domain of `T`,
this method will produce a sample drawn from the distribution represented by `at`.

```@docs; canonical=false
AliasTables.sample
```

## Performance characteristics

Constructing an `AliasTable{T}` is O(n) in the number of elements in the input collection,
with a low constant factor. Sampling itself is O(1) with a very low constant factor. It is
branchless, involves one random array read, and takes about 20 instructions more than
`rand(T)`.

```julia
julia> using Chairmarks

julia> @b rand(1000) AliasTable
13.250 μs (5 allocs: 23.906 KiB)

julia> @b AliasTable(rand(1000)) rand
3.059 ns

julia> @b rand(UInt64)
2.891 ns

julia> @b AliasTable(rand(1000)) rand(_, 1000)
1.588 μs (3 allocs: 7.875 KiB)

julia> @b rand(UInt64, 1000)
606.870 ns (3 allocs: 7.875 KiB)
```

Bulk generation of UInt64 has hand written llvm instructions to support SIMD, while alias
tables don't SIMD as nicely and have not been as aggressively optimized; hence the
difference in bulk generation time while scalar generation time is similar.

## Docstrings

The docstring for the [`AliasTable`](@ref) constructor defines the API for constructing
`AliasTable`s, the [`AliasTables.probabilities`](@ref) function allows recovery of the exact
sampling probabilities provided by an `AliasTable`, and the [`AliasTables.sample`](@ref)
function provides an alternative API for sampling from `AliasTable`s.

However, the primary sampling API is the Random API. `AliasTable`s may be used as a sampling
domain according to the specifications layed out by the Random stdlib.

```@docs
AliasTable
AliasTables.sample
AliasTables.probabilities
AliasTables.set_weights!
length(::AliasTable)
```
