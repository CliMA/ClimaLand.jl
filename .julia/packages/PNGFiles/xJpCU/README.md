# PNGFiles.jl

FileIO.jl integration for PNG files

[![Codecov](https://codecov.io/gh/JuliaIO/PNGFiles.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaIO/PNGFiles.jl)

## Installation

Installation is recommended via the `ImageIO` thin IO wrapper for `FileIO`:

```jl
pkg> add ImageIO  # Press ']' to enter te Pkg REPL mode.
```

## Usage

Once `ImageIO` is installed, usage is as simple as:

```jl
using FileIO
save("img.png", rand(Gray, 100, 100))
img = load("img.png")
```

Or direct usage, if `PNGFiles` has been directly installed:
```jl
using PNGFiles
PNGFiles.save("img.png", rand(Gray, 100, 100))
img = PNGFiles.load("img.png")
```
