# RoundingEmulator.jl

Emulate directed rounding using only the default rounding mode. 

[![Build Status](https://travis-ci.com/matsueushi/RoundingEmulator.jl.svg?branch=master)](https://travis-ci.com/matsueushi/RoundingEmulator.jl) [![Dev Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://matsueushi.github.io/RoundingEmulator.jl/dev/)

This package is meant to produce the exact same results of `Rounding.setrounding` ([deprecated](https://github.com/JuliaLang/julia/pull/27166)) without switching rounding modes.

## Requirements 
 - Julia 1.3 or higher
 - `Base.Rounding.get_zero_subnormals() == true`. (See [Base.Rounding.get_zero_subnormals](https://docs.julialang.org/en/v1/base/numbers/#Base.Rounding.get_zero_subnormals))
