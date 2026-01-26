# Isoband

[![Build Status](https://travis-ci.com/jkrumbiegel/Isoband.jl.svg?branch=master)](https://travis-ci.com/jkrumbiegel/Isoband.jl)
[![Coverage](https://codecov.io/gh/jkrumbiegel/Isoband.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jkrumbiegel/Isoband.jl)


`Isoband.jl` wraps `isoband_jll`, which gives access to [wilkelab's isoband package](https://github.com/wilkelab/isoband), which powers contour plots in ggplot.


# Installation

```
] add Isoband
```

# Use

```julia
xs = 0:10 # correspond to columns of zs
ys = 0:20 # correspond to rows of zs
zs = ys * xs'

# there are single band / line methods...

one_band = isobands(xs, ys, zs, 0, 10)
# access polygon data
one_band.x
one_band.y
one_band.id

one_line = isolines(xs, ys, zs, 5)


# ...and multi band / line methods

multiple_bands = isobands(xs, ys, zs, [0, 1, 10], [1, 10, 100])

multiple_lines = isolines(xs, ys, zs, [1, 10, 20])
# access first line data
multiple_lines[1].x
multiple_lines[1].y
multiple_lines[1].id

```