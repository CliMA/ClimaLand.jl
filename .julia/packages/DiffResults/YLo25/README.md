# DiffResults

![CI](https://github.com/JuliaDiff/DiffResults.jl/workflows/CI/badge.svg)
[![Coverage Status](https://coveralls.io/repos/github/JuliaDiff/DiffResults.jl/badge.svg?branch=master)](https://coveralls.io/github/JuliaDiff/DiffResults.jl?branch=master)

[![](https://img.shields.io/badge/docs-stable-blue.svg)](http://www.juliadiff.org/DiffResults.jl/stable)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](http://www.juliadiff.org/DiffResults.jl/latest)

Many differentiation techniques can calculate primal values and multiple orders of
derivatives simultaneously. In other words, there are techniques for computing `f(x)`,
`∇f(x)` and `H(f(x))` in one fell swoop!

For this purpose, DiffResults provides the `DiffResult` type, which can be passed
to in-place differentiation methods instead of an output buffer. The method
then loads all computed results into the given `DiffResult`, which the user
can then query afterwards.
