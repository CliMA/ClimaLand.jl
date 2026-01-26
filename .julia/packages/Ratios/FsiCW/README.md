# Ratios

[![CI](https://github.com/timholy/Ratios.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/timholy/Ratios.jl/actions/workflows/ci.yml) [![Coverage](https://codecov.io/gh/timholy/Ratios.jl/branch/master/graph/badge.svg?token=ZVcLnVyTBB)](https://codecov.io/gh/timholy/Ratios.jl)

This package provides types similar to Julia's `Rational` type, which make some sacrifices but have better computational performance at the risk of greater risk of overflow.

Currently the only type provided is `SimpleRatio(num, den)` for two integers `num` and `den`.

Demo:

```julia
julia> x, y, z = SimpleRatio(1, 8), SimpleRatio(1, 4), SimpleRatio(2, 8)
(SimpleRatio{Int}(1, 8), SimpleRatio{Int}(1, 4), SimpleRatio{Int}(2, 8))

julia> x+y
SimpleRatio{Int}(12, 32)

julia> x+z
SimpleRatio{Int}(3, 8)
```

`y` and `z` both represent the rational number `1//4`, but when performing arithmetic with `x`
`z` is preferred because it has the same denominator and is less likely to overflow.

To detect overflow, [SaferIntegers.jl](https://github.com/JeffreySarnoff/SaferIntegers.jl) is recommended:

```julia
julia> using Ratios, SaferIntegers

julia> x, y = SimpleRatio{SafeInt8}(1, 20), SimpleRatio{SafeInt8}(1, 21)
(SimpleRatio{SafeInt8}(1, 20), SimpleRatio{SafeInt8}(1, 21))

julia> x + y
ERROR: OverflowError: 20 * 21 overflowed for type Int8
Stacktrace:
[...]
```

[FastRationals](https://github.com/JeffreySarnoff/FastRationals.jl) is another package with safety and performance characteristics that lies somewhere between `SimpleRatio` and `Rational`:

```julia
julia> @btime x + y setup=((x, y) = (SimpleRatio(rand(-20:20), rand(2:20)), SimpleRatio(rand(-20:20), rand(2:20))));
  1.969 ns (0 allocations: 0 bytes)

julia> @btime x + y setup=((x, y) = (FastRational(rand(-20:20), rand(2:20)), FastRational(rand(-20:20), rand(2:20))));
  3.192 ns (0 allocations: 0 bytes)

julia> @btime x + y setup=((x, y) = (Rational(rand(-20:20), rand(2:20)), Rational(rand(-20:20), rand(2:20))));
  23.065 ns (0 allocations: 0 bytes)
```
