# FreeType.jl

[![Build Status](https://travis-ci.org/JuliaGraphics/FreeType.jl.svg?branch=master)](https://travis-ci.org/JuliaGraphics/FreeType.jl)
[![codecov](https://codecov.io/gh/JuliaGraphics/FreeType.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaGraphics/FreeType.jl)

[FreeType](http://www.freetype.org/) bindings for [Julia](http://julialang.org/).

## Example

```julia
using FreeType

library = Vector{FT_Library}(undef, 1)
error = FT_Init_FreeType(library)
@assert error == 0
```
