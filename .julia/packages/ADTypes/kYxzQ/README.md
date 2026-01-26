# ADTypes.jl: Multi-Valued Logic for Automatic Differentiation Choices

[![Docs stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://docs.sciml.ai/ADTypes/stable/)
[![Docs dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://docs.sciml.ai/ADTypes/dev/)
[![Build Status](https://github.com/SciML/ADTypes.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/SciML/ADTypes.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/SciML/ADTypes.jl/branch/main/graph/badge.svg)](https://app.codecov.io/gh/SciML/ADTypes.jl)

ADTypes.jl is a multi-valued logic system to choose an automatic differentiation (AD) package and specify its parameters.

## Which AD libraries are supported?

See the API reference in the documentation.
If a given package is missing, feel free to open an issue or pull request.

## Why should AD users adopt this standard?

A natural approach is to use a keyword argument with e.g. `Bool` or `Symbol` values.
Let's see a few examples to understand why this is not enough:

  - `autodiff = true`: ambiguous, we don't know which AD package should be used
  - `autodiff = :forward`: ambiguous, there are several AD packages implementing both forward and reverse mode (and there are other modes beyond that)
  - `autodiff = :Enzyme`: ambiguous, some AD packages can work both in forward and reverse mode
  - `autodiff = (:Enzyme, :forward)`: not too bad, but many AD packages require additional configuration (number of chunks, tape compilation, etc.)

A more involved struct is thus required, with package-specific parameters.
If every AD user develops their own version of said struct, it will ruin interoperability.
This is why ADTypes.jl provides a single set of shared types for this task, as an extremely lightweight dependency.
They are types and not enums because we need AD choice information statically to use it for dispatch.
