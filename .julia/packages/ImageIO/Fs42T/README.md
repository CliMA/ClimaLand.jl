# ImageIO.jl

FileIO.jl integration for image files

![Julia version](https://img.shields.io/badge/julia-%3E%3D%201.6-blue)
[![CI](https://github.com/JuliaIO/ImageIO.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/JuliaIO/ImageIO.jl/actions/workflows/CI.yml)
[![Codecov](https://codecov.io/gh/JuliaIO/ImageIO.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaIO/ImageIO.jl)

| Format | Extensions | Provider | Implementation | Comment |
| ------- | ---------- | -------- | ---- | ----------- |
| JPEG | `.jpg`, `.jpeg` | [JpegTurbo.jl](https://github.com/johnnychen94/JpegTurbo.jl) | Julia wrapper of [libjpeg-turbo](https://github.com/libjpeg-turbo/libjpeg-turbo) | [Benchmark results against other backends](https://github.com/johnnychen94/JpegTurbo.jl/issues/15) |
| [OpenEXR](https://www.openexr.com/) | `.exr` | [OpenEXR.jl](https://github.com/twadleigh/OpenEXR.jl) | Julia wrapper of [OpenEXR](https://github.com/AcademySoftwareFoundation/openexr) | |
| Portable Bitmap formats | `.pbm`, `.pgm`, `.ppm` | [Netpbm.jl](https://github.com/JuliaIO/Netpbm.jl) | pure Julia | |
| PNG (Portable Network Graphics) | `.png` | [PNGFiles.jl](https://github.com/JuliaIO/PNGFiles.jl) | Julia wrapper of [libpng](https://github.com/glennrp/libpng) | [Benchmark vs. ImageMagick & QuartzImageIO](https://github.com/JuliaIO/PNGFiles.jl/issues/1#issuecomment-586749654) |
| [QOI (Quite Okay Image)](https://qoiformat.org/) format | `.qoi` | [QOI.jl](https://github.com/KristofferC/QOI.jl) | pure Julia | |
| DEC SIXEL (six-pixels) graphics | `.six`, `.sixel` | [Sixel.jl](https://github.com/johnnychen94/Sixel.jl) | Julia wrapper of [libsixel](https://github.com/libsixel/libsixel) | |
| TIFF (Tag Image File Format) | `.tiff`, `.tif` | [TiffImages.jl](https://github.com/tlnagy/TiffImages.jl) | pure Julia | check [OMETIFF.jl](https://github.com/tlnagy/OMETIFF.jl) for OMETIFF support |
| WebP | `.webp` | [WebP.jl](https://github.com/stemann/WebP.jl) | Julia wrapper of [libwebp](https://developers.google.com/speed/webp/docs/api) | |

## Installation

Install with Pkg:

```jl
pkg> add ImageIO  # Press ']' to enter te Pkg REPL mode
```

## Usage

```jl
using FileIO
save("test.png", rand(Gray, 100, 100))
load("test.png")
save("test.ppm", rand(RGB, 100, 100))
load("test.ppm")
save("test.tiff", rand(RGB, 100, 100))
load("test.tiff")
```

## Canonicalization

Some image loaders may return custom AbstractArray types. By default, this package "canonicalizes" the returned type to be either `Array` or [`IndirectArray`](https://github.com/JuliaArrays/IndirectArrays.jl).
An exception is for calls like `load(filename; mmap=true)` where the image data will be "lazily" loaded using [memory-mapped IO](https://en.wikipedia.org/wiki/Memory-mapped_I/O), in which case the default is to allow the lower-level I/O package to return whatever AbstractArray type it chooses.

You can manually control canonicalization with `load(filename; canonicalize=tf)` where `tf` is `true` or `false`.

## Compatibility

If you're using old Julia versions (`VERSION < v"1.3"`), a dummy ImageIO version v0.0.1 with no real function will be installed.
In this case, you still need to install [ImageMagick.jl] to make `FileIO.save`/`FileIO.load` work.

[ImageMagick.jl]: https://github.com/JuliaIO/ImageMagick.jl
