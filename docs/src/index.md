```@eval
using Pkg
io = IOBuffer()
v = Pkg.installed()["ClimaLand"]
print(io, """
    # ClimaLand.jl Documentation (v$(v))

    """)
import Markdown
Markdown.parse(String(take!(io)))
```

## Introduction

ClimaLand is the Land model of the Climate Modeling Alliance (CliMA) Earth System Model, which
also contains other components ([atmosphere](https://github.com/CliMA/ClimaAtmos.jl), [ocean](https://github.com/CliMA/ClimaOcean.jl), [sea-ice](https://github.com/CliMA/ClimaSeaIce.jl)).

ClimaLand can be run coupled (or "online") with these other components via [ClimaCoupler](https://github.com/CliMA/ClimaCoupler.jl),
or it can be run as a standalone, via prescribed meteorological data ("offline").

ClimaLand library, described in this documentation, is written in the Julia
programming language. It aims to be fast and have a clear syntax. ClimaLand
can run on CPU or GPU, it has a modular design, and is flexible in many ways. This documentation will expand on each of these elements.

## Important Links

- [CliMA Homepage](https://clima.caltech.edu/)
- [CliMA GitHub Organisation](https://github.com/CliMA)
- [ClimaCoupler](https://github.com/CliMA/ClimaCoupler.jl)
- [ClimaAnalysis](https://github.com/CliMA/ClimaAnalysis.jl)
- [Julia Homepage](https://julialang.org)

## Documentation for Users and Developers

ClimaLand has documentation for both users and developers. The documentation for users is aimed at scientists who wants
to run simulations using ClimaLand, whereas the documentation for developers is aimed at contributors of the ClimaLand
code library. As such, users can skip reading the docs for developers, and vice-versa.

## Physical units

Note that CliMA, in all its repositories, uses Standard Units, reminded below

| Quantity              | Unit Name | SI Symbol | SI Unit Equivalent  |
|:----------------------|:----------|:----------|:--------------------|
| Length                | Meter     | m         | 1 m                 |
| Mass                  | Kilogram  | kg        | 1 kg                |
| Time                  | Second    | s         | 1 s                 |
| Temperature           | Kelvin    | K         | 1 K                 |
| Amount of Substance   | Mole      | mol       | 1 mol               |
| Energy                | Joule     | J         | 1 J = 1 N·m         |
| Power                 | Watt      | W         | 1 W = 1 J/s         |
| Pressure              | Pascal    | Pa        | 1 Pa = 1 N/m²       |
| Frequency             | Hertz     | Hz        | 1 Hz = 1 s⁻¹        |

