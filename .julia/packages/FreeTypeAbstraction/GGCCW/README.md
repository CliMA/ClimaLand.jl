[![codecov](https://codecov.io/gh/JuliaGraphics/FreeTypeAbstraction.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaGraphics/FreeTypeAbstraction.jl)

[![Build Status](https://travis-ci.org/JuliaGraphics/FreeTypeAbstraction.jl.svg?branch=master)](https://travis-ci.org/JuliaGraphics/FreeTypeAbstraction.jl)

# FreeTypeAbstraction

Draw text into a Matrix.

```Julia

using FreeTypeAbstraction

# search for a font on your system
# Choose the name of a font you have installed
face = findfont("arial") # if this returns `nothing`, that means it could not find the font

# render a character
img, metric = renderface(face, 'C', 100)

# render a string into an existing matrix
myarray = zeros(UInt8, 100, 100)
pixelsize = 10
x0, y0 = 90, 10
renderstring!(myarray, "hello", face, pixelsize, x0, y0, halign=:hright)
```

credits to @aaalexandrov from whom most of the early code comes.
