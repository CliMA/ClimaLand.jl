# ColorBrewer.jl

[![](https://travis-ci.org/timothyrenner/ColorBrewer.jl.svg?branch=master)](https://travis-ci.org/timothyrenner/ColorBrewer.jl)

A small library for generating color brewer color palettes.
For more information on color brewer, see http://colorbrewer2.org/.

## About

The color brewer palette collection was originally designed for maps, but is also useful for quickly generating nice colors for any kind of plot.

## Setup

This package is now registered. Add it with

```julia
Pkg.add("ColorBrewer");
```

You can also add it with

```julia
Pkg.clone("git://github.com/timothyrenner/ColorBrewer.jl.git");
```

## Usage

The module contains one function, `palette`, which takes a scheme name and an integer for the number of colors in the palette.
It returns an array of `RGB` values, which are defined in the [ColorTypes.jl](https://github.com/JuliaGraphics/ColorTypes.jl) package.

For example, to obtain a five color sequential green scheme, call

```julia
    using ColorBrewer

    greens = palette("Greens", 5);
```
Pretty simple stuff.

It can drop right into [Gadfly](https://github.com/dcjones/Gadfly.jl) plots.

```julia
    using Gadfly

    p = plot([x -> sin(n*x) for n in 1:3], 0, 25, 
             Scale.color_discrete_manual(greens[3:5]...))
```

![Example](images/color_example.png)

If the `palette` function is passed an invalid scheme or number of colors, it throws an `ArgumentError`.

The schemes are broken up into three categories based on the type of data you want to show: sequential, diverging, and qualitative.

All sequential color schemes can have between 3 and 9 colors.
The available sequential color schemes are:

| Name    | Example                            |
|---------|------------------------------------|
| Blues   | ![Blues](images/Blues.png "Blues") |
| Oranges | ![Oranges](images/Oranges.png "Oranges") |
| Greens  | ![Greens](images/Greens.png "Greens") |
| Reds    | ![Reds](images/Reds.png "Reds") |
| Purples | ![Purples](images/Purples.png "Purples") |
| Greys   | ![Greys](images/Greys.png "Greys") |
| OrRd    | ![OrRd](images/OrRd.png "OrRd") |
| GnBu    | ![GnBu](images/GnBu.png "GnBu") |
| PuBu    | ![PuBu](images/PuBu.png "PuBu") |
| PuRd    | ![PuRd](images/PuRd.png "PuRd") |
| BuPu    | ![BuPu](images/BuPu.png "BuPu") |
| BuGn    | ![BuGn](images/BuGn.png "BuGn") |
| YlGn    | ![YlGn](images/YlGn.png "YlGn") |
| RdPu    | ![RdPu](images/RdPu.png "RdPu") |
| YlOrBr  | ![YlOrBr](images/YlOrBr.png "YlOrBr") |
| YlGnBu  | ![YlGnBu](images/YlGnBu.png "YlGnBu") |
| YlOrRd  | ![YlOrRd](images/YlOrRd.png "YlOrRd") |
| PuBuGn  | ![PuBuGn](images/PuBuGn.png "PuBuGn") |

All diverging color schemes can have between 3 and 11 colors.
The available diverging color schemes are:

| Name     | Example                                     |
| -------- | ------------------------------------------- |
| Spectral | ![Spectral](images/Spectral.png "Spectral") |
| RdYlGn   | ![RdYlGn](images/RdYlGn.png "RdYlGn")       |
| RdBu     | ![RdBu](images/RdBu.png "RdBu")             |
| PiYG     | ![PiYG](images/PiYG.png "PiYG")             |
| PRGn     | ![PRGn](images/PRGn.png "PRGn")             |
| RdYlBu   | ![RdYlBu](images/RdYlBu.png "RdYlBu")       |
| BrBG     | ![BrBG](images/BrBG.png "BrBG")             |
| RdGy     | ![RdGy](images/RdGy.png "RdGy")             |
| PuOr     | ![PuOr](images/PuOr.png "PuOr")             |

The number of colors a qualitative color scheme can have depends on the scheme.
The available qualitative color schemes are:

| Name | Example                         |
| ---- | ------------------------------- |
| Set1 | ![Set1](images/Set1.png "Set1") |
| Set2 | ![Set2](images/Set2.png "Set2") |
| Set3 | ![Set3](images/Set3.png "Set3") |
| Accent | ![Accent](images/Accent.png "Accent") |
| Dark2 | ![Dark2](images/Dark2.png "Dark2") |
| Paired | ![Paired](images/Paired.png "Paired") |
| Pastel1 | ![Pastel1](images/Pastel1.png "Pastel1") |
| Pastel2 | ![Pastel2](images/Pastel2.png "Pastel2") |

The `Paired` and `Set3` schemes can have between 3 and 12 colors.
The `Pastel1` and `Set1` schemes can have between 3 and 9 colors.
The `Set2`, `Accent`, `Dark2`, and `Pastel2` schemes can have between 3 and 8 colors.


## Thanks

Huge thanks to ikirill for providing the color scheme previews!
