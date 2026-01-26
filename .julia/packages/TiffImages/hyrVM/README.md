<img src="https://raw.githubusercontent.com/tlnagy/TiffImages.jl/b4f946b75caae33992eb818551230ded7a9aa9de/docs/src/assets/fulllogo.svg" width="500">

_"Don't get into a tiff with your images"_

| **Stable release** | **Documentation**                 | **Build Status**                                              |
|:------------------------------------------------------|:-------------------------------------------------------------------------|:--------------------------------------------------------------|
| ![](https://juliahub.com/docs/TiffImages/version.svg) | [![][docs-stable-img]][docs-stable-url][![][docs-dev-img]][docs-dev-url] | [![][status-img]][status-url] [![][ci-img]][ci-url] [![][codecov-img]][codecov-url] [![][aqua-img]][aqua-url] |

This package aims to be a fast, minimal, and correct TIFF reader and writer
written in Julia.

## Features

- Fast reading and writing of many common TIFFs
- Extensible core for other TIFF packages to build on
- Native integration with `Colors.jl` and the Julia Array ecosystem
- Memory-mapping for loading images too large to fit in memory
- BigTIFF standard (TIFFs larger than 4 GB)
- Arbitrary bit depths (e.g. 12 or 14 bit cameras)
- Common compression algorithms like LZW and Packbits

## Installation

`TiffImages.jl` is available through Julia's general repository. You can install
it by running the following commands in the Julia REPL:

```julia
using Pkg
Pkg.install("TiffImages")
```

Please see the documentation above for usage details and examples

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://tamasnagy.com/TiffImages.jl/stable

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://tamasnagy.com/TiffImages.jl/dev

[ci-img]: https://github.com/tlnagy/TiffImages.jl/workflows/CI/badge.svg
[ci-url]: https://github.com/tlnagy/TiffImages.jl/actions

[codecov-img]: https://codecov.io/gh/tlnagy/TiffImages.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/tlnagy/TiffImages.jl

[status-img]: https://www.repostatus.org/badges/latest/active.svg
[status-url]: https://www.repostatus.org/#active

[aqua-img]: https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg
[aqua-url]: https://github.com/JuliaTesting/Aqua.jl
