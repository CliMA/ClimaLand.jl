# TiffImages.jl

*Pure-Julia TIFF reader and writer with a focus on correctness* 🧐

TIFF (Tagged Image File Format) is a notoriously flexible file format that is
very difficult to support properly so why not just link `libtiff` and call it
a day? Because [Julia developers are
greedy](https://julialang.org/blog/2012/02/why-we-created-julia/). I wanted to
design a clean, minimal, and standards-compliant TIFF reader and writer that can
have the speed and compliance of `libtiff` while adding modern features like
streaming, out-of-memory support, and fancy color support. I wanted to design it
to be extensible such that packages like
[`OMETIFF.jl`](https://github.com/tlnagy/OMETIFF.jl) can hook right in with
minimal overhead. I wanted to leverage the wonderful Julia Arrays ecosystem to
do as much lazily and flexibly as possible.

## Features

TiffImages.jl supports:

- The TIFF 6.0 baseline spec
- Thorough testing
- HDR images stored as 32bit or 64bit floats
- BigTIFF standard (TIFFs larger than 4 GB)
- Memory-mapped and lazy loading/writing
- Arbitrary bit depths (e.g. 12 or 14 bit cameras)
- Common compression algorithms like LZW and Packbits

## Usage

Check out the examples to see how to use `TiffImages.jl`

```@contents
Pages = ["examples/reading.md", "examples/writing.md", "examples/mmap_lazyio.md"]
Depth = 1
```
