# OpenEXR.jl

[![Build Status](https://github.com/twadleigh/OpenEXR.jl/workflows/CI/badge.svg)](https://github.com/twadleigh/OpenEXR.jl/actions?query=workflow%3A%22CI%22+branch%3Amaster)
[![codecov.io](http://codecov.io/github/twadleigh/OpenEXR.jl/coverage.svg?branch=master)](http://codecov.io/github/twadleigh/OpenEXR.jl?branch=master)

Saving and loading of OpenEXR files.

# Load an EXR file

```jl
using OpenEXR

myimage = OpenEXR.load("myimage.exr")
```

`myimage` has type `Array{T,2}` where `T` is one of:

- RGBA{Float16}
- RGB{Float16}
- GrayA{Float16}
- Gray{Float16}

depending on which channels are present in the file.

# Save an image as an EXR file

```jl
OpenEXR.save("myimage2.exr", myimage)
```

`myimage` can be any subtype of `AbstractArray{C,2}` where `C` is any color type defined in
[ColorTypes.jl](https://github.com/JuliaGraphics/ColorTypes.jl).
