[![Build Status](https://travis-ci.org/JuliaGraphics/SignedDistanceFields.jl.svg?branch=master)](https://travis-ci.org/JuliaGraphics/SignedDistanceFields.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/ans0khpf1j3hktrs/branch/master?svg=true)](https://ci.appveyor.com/project/SimonDanisch/signeddistancefields-jl/branch/master)
[![codecov.io](http://codecov.io/github/JuliaGraphics/SignedDistanceFields.jl/coverage.svg?branch=master)](http://codecov.io/github/JuliaGraphics/SignedDistanceFields.jl?branch=master)
[![Coverage Status](https://coveralls.io/repos/JuliaGraphics/SignedDistanceFields.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/JuliaGraphics/SignedDistanceFields.jl?branch=master)

# SignedDistanceFields

This package implements an algorithm for quickly computing signed distance fields in 2D.

`result = sdf(img)`, where `img` is a 2-dimensional boolean array with `true` indicating the foreground and `false` indicating the background.

If you want a downsampled distance field, you can call `sdf(img, xsize[, ysize=xsize])`, and pass it your desired output dimensions. Right now we only support scaling by integer amounts (e.g. turning a 1024x1024 image into a 64x64 SDF or a 500x500 image into a 125x125 SDF).

There is also an `edf` function with identical usage that calculates the Euclidean distance transform—the distance from every background pixel to the closest pixel in the foreground.

Signed distance fields are a useful representation for rendering glyphs and other shapes with crisp edges using OpenGL. By allowing the graphics hardware to interpolate between distance samples and testing for the shape boundary in the vertex shader, you can get surprisingly high-resolution output from smaller SDFs.

Here's a [paper from Valve](http://www.valvesoftware.com/publications/2007/SIGGRAPH2007_AlphaTestedMagnification.pdf) introducing the idea and here's a [good blog post](https://www.mapbox.com/blog/text-signed-distance-fields/) about using signed distance fields for rendering text.

[![Build Status](https://travis-ci.org/yurivish/SignedDistanceFields.jl.svg?branch=master)](https://travis-ci.org/yurivish/SignedDistanceFields.jl)
