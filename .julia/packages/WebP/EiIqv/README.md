# WebP

[![Build Status](https://github.com/stemann/WebP.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/stemann/WebP.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Coverage](https://codecov.io/gh/stemann/WebP.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/stemann/WebP.jl)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

WebP.jl is a Julia library for handling [WebP](https://developers.google.com/speed/webp) images. WebP provides both lossy and lossless compression of images, and may offer smaller file sizes compared to JPEG and PNG.

The core functionality of this package is supported by the [libwebp](https://developers.google.com/speed/webp/docs/api) C library.

## Usage

This package provides functions for reading and writing WebP image files,

* `WebP.read_webp`
* `WebP.write_webp`

as well as functions for decoding and encoding WebP image data,

* `WebP.decode`
* `WebP.encode`

### Reading and writing

An image may be written,
```julia
using TestImages
using WebP

image = testimage("lighthouse")
WebP.write_webp("lighthouse.webp", image)
```

and subsequently read,
```julia
image = WebP.read_webp("lighthouse.webp")
```

### Decoding and encoding

An image may be encoded,
```julia
using TestImages
using WebP

image = testimage("lighthouse")
data = WebP.encode(image) # data is a Vector{UInt8}
```

and subsequently decoded,
```julia
image = WebP.decode(data)
```
