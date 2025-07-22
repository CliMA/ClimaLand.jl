```@eval
using Pkg
io = IOBuffer()
v = filter(x-> x.second.name == "ClimaLand", Pkg.dependencies()) |> x -> first(x)[2].version
print(io, """
    # ClimaLand.jl Documentation (v$(v))

    """)
import Markdown
Markdown.parse(String(take!(io)))
```

## Introduction

ClimaLand is the Land model of the Climate Modeling Alliance (CliMA) Earth System Model, which
also contains other components ([atmosphere](https://github.com/CliMA/ClimaAtmos.jl), [ocean](https://github.com/CliMA/ClimaOcean.jl), [sea-ice](https://github.com/CliMA/ClimaSeaIce.jl)).
Details about the CliMA project can be found on the [project website](https://clima.caltech.edu/).

ClimaLand can be run coupled ("online") with these other components via [ClimaCoupler](https://github.com/CliMA/ClimaCoupler.jl),
or it can be run as a standalone model, via prescribed meteorological data ("offline").

The ClimaLand library, described in this documentation, is written in the Julia
programming language. It aims to be fast and have a clear syntax. ClimaLand
can run on CPU or GPU, it has a modular design, and is flexible in many ways. This documentation will expand on each of these elements.

## Documentation for Users and Developers

ClimaLand has documentation for both users and developers. The documentation for users is aimed at scientists who wants
to run simulations using ClimaLand, whereas the documentation for developers is aimed at contributors of the ClimaLand
code library. As such, users can skip reading the docs for developers, and vice-versa.

## This documentation includes information about:
- How to run your first ClimaLand simulation
- A series of tutorials explaining ClimaLand models and physics
- The structure of ClimaLand models and code
- How to analyze model output, calibrate models using CliMA's pipeline, and restart simulations
- and more!

## Important Links

- [CliMA GitHub Organisation](https://github.com/CliMA)
- [Julia Homepage](https://julialang.org)
- [Julia Manual](https://docs.julialang.org/en/v1/manual/getting-started/)
