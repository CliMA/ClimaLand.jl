# Showoff

[![Build Status](https://travis-ci.org/JuliaGraphics/Showoff.jl.svg?branch=master)](https://travis-ci.org/JuliaGraphics/Showoff.jl)
[![Coverage Status](https://coveralls.io/repos/github/JuliaGraphics/Showoff.jl/badge.svg?branch=master)](https://coveralls.io/github/JuliaGraphics/Showoff.jl?branch=master)

Showoff provides an interface for consistently formatting an array of n things,
e.g. numbers, dates, unitful values. It's used in Gadfly, Plots and Makie to
label axes and keys.

It defines a function called `showoff` that takes an `AbstractArray` of some
type, and returns an array of strings of the same length.

If you want your type to look nice when plotted, just define a `showoff`
function. Here's an example.

```julia
using Showoff

struct Percent
    value::Float64
end

function Showoff.showoff(xs::AbstractArray{Percent})
    return [string(x, "%") for x in showoff([x.value for x in xs])]
end
```

Now we (and more importantly, Gadfly) can print percentages like:

```julia
map(println, showoff([Percent(100 * rand()) for _ in 1:20]))
```
```julia
60.505943%
73.255897%
97.477079%
43.330976%
69.023165%
52.580184%
13.011683%
22.718034%
93.843776%
29.875979%
64.110999%
91.203653%
91.534161%
80.684188%
81.674362%
11.530227%
30.498260%
38.876922%
35.444115%
8.857208%
```

Notice, that compared to `show`, these all have the same number of digits
trailing the `.`, and look nice when right-aligned.

When no specialized `showoff` is defined, it falls back on the `show` function.

This package was originally written by [Daniel C. Jones](https://github.com/dcjones).
