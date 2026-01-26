# JpegTurbo

[JpegTurbo.jl](https://github.com/johnnychen94/JpegTurbo.jl) is a Julia wrapper of the C library
[libjpeg-turbo](https://github.com/libjpeg-turbo/libjpeg-turbo) that provides IO support for the
JPEG image format. This package also backs the JPEG IO part of
[ImageIO](https://github.com/JuliaIO/ImageIO.jl) and [FileIO](https://github.com/JuliaIO/FileIO.jl).

For benchmark results against other image IO backends, please check
[here](https://github.com/johnnychen94/JpegTurbo.jl/issues/15).

## Usage

There are two different usages for this package:

- (convenient) via the [FileIO](https://github.com/JuliaIO/FileIO.jl): `save`/`load`
- (powerful) via the JpegTurbo.jl interfaces: `jpeg_encode`/`jpeg_decode`

### FileIO interface: `save`/`load`

[FileIO](https://github.com/JuliaIO/FileIO.jl) is an IO frontend with various IO backends;
[ImageIO](https://github.com/JuliaIO/ImageIO.jl) is the default IO backend provided by the
JuliaImages ecosystem. When JpegTurbo (and/or ImageIO) are available in `DEPOT_PATH`, FileIO will
uses JpegTurbo to load and save the JPEG images:

```julia
using FileIO
img = rand(64, 64)
save("test.jpg", img)
load("test.jpg")
```

Note that you do not necessarily need to install them in your project environments. For instance,
you can do `(@v1.8) pkg> add JpegTurbo` or `(@v1.8) pkg> add ImageIO` and it should work for your
local setup.

### JpegTurbo interface: `jpeg_encode`/`jpeg_decode`

`jpeg_encode` is used to compress 2D colorant matrix as JPEG image.

```julia
jpeg_encode(filename::AbstractString, img; kwargs...) -> Int
jpeg_encode(io::IO, img; kwargs...) -> Int
jpeg_encode(img; kwargs...) -> Vector{UInt8}
```

`jpeg_decode` is used to decompress JPEG image as 2D colorant matrix.

```julia
jpeg_decode([T,] filename::AbstractString; kwargs...) -> Matrix{T}
jpeg_decode([T,] io::IO; kwargs...) -> Matrix{T}
jpeg_decode([T,] data::Vector{UInt8}; kwargs...) -> Matrix{T}
```

### Advanced: in-memory encode/decode

For some applications, it can be faster to do encoding/decoding without the need
to read/write disk:

```julia
using JpegTurbo
img = rand(64, 64)
bytes = jpeg_encode(img) # Vector{UInt8}
img_saveload = jpeg_decode(bytes) # size: 64x64
```

### Advanced: preview optimization

One can request a single-component output or a downsampled output so that fewer calculation is
needed during the decompression. This can be particularly useful to accelerate image preview.

```julia
using BenchmarkTools, TestImages, JpegTurbo
filename = testimage("earth", download_only=true)
# full decompression
@btime jpeg_decode(filename); # 224.760 ms (7 allocations: 51.54 MiB)
# only decompress luminance component
@btime jpeg_decode(Gray, filename); # 91.157 ms (6 allocations: 17.18 MiB)
# process only a few pixels
@btime jpeg_decode(filename; scale_ratio=0.25); # 77.254 ms (8 allocations: 3.23 MiB)
# process only a few pixels for luminance component
@btime jpeg_decode(Gray, filename; scale_ratio=0.25); # 63.119 ms (6 allocations: 1.08 MiB)
```

An exclusive alternative to `scale_ratio` is `preferred_size`:

```julia
# minimal `scale_ratio` that output size is greater than or equal to (512, 512)
jpeg_decode(filename; preferred_size=(512, 512)) # size: (751, 750)
# maximal `scale_ratio` that output size is less than or equal to (512, 512)
jpeg_decode(filename; preferred_size=(<=, (512, 512))) # size: (376, 375)
```
