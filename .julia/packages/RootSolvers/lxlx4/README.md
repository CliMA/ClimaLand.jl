<div align="center">
  <img src="docs/src/assets/logo.svg" alt="RootSolvers.jl Logo" width="128" height="128">
</div>

# RootSolvers.jl

A high-performance root solver package with GPU support and broadcasting across abstract types

RootSolvers.jl provides robust, efficient numerical methods for finding roots of nonlinear equations. It supports broadcasting across abstract types including GPU arrays and custom field types, making it ideal for high-performance computing applications in climate modeling, machine learning, and scientific computing.

|                           |                                                                          |
|--------------------------:|:-------------------------------------------------------------------------|
| **Stable Release**        | [![stable][stable-img]][stable-url] [![docs-stable][docs-stable-img]][docs-stable-url] |
| **Latest Documentation**  | [![dev][docs-latest-img]][docs-latest-url]                                |
| **Unit Tests**            | [![unit tests][unit-tests-img]][unit-tests-url] [![codecov][codecov-img]][codecov-url] |
| **Downloads**             | [![Downloads][dlt-img]][dlt-url]                                          |

[docs-latest-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-latest-url]: https://CliMA.github.io/RootSolvers.jl/dev/

[stable-img]: https://img.shields.io/github/v/release/CliMA/RootSolvers.jl?label=stable
[stable-url]: https://github.com/CliMA/RootSolvers.jl/releases/latest

[docs-stable-img]: https://img.shields.io/badge/docs-stable-green.svg
[docs-stable-url]: https://CliMA.github.io/RootSolvers.jl/stable/

[unit-tests-img]: https://github.com/CliMA/RootSolvers.jl/actions/workflows/OS-UnitTests.yml/badge.svg
[unit-tests-url]: https://github.com/CliMA/RootSolvers.jl/actions/workflows/OS-UnitTests.yml

[codecov-img]: https://codecov.io/gh/CliMA/RootSolvers.jl/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/CliMA/RootSolvers.jl

[dlt-img]: https://img.shields.io/badge/dynamic/json?url=http%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Ftotal_downloads%2FRootSolvers&query=total_requests&label=Downloads
[dlt-url]: https://juliapkgstats.com/pkg/RootSolvers

## Features

- **Multiple Root-Finding Methods**: Secant, Regula Falsi, Brent's method, Newton's method with automatic differentiation
- **GPU Support**: Full GPU acceleration with CUDA.jl and other GPU array types
- **Broadcasting**: Allows broadcasting over distributed arrays and custom field types
- **Dual Number Support**: Compatible with automatic differentiation frameworks, allowing integration into differentiable models
- **Flexible Convergence Criteria**: Multiple tolerance types for different applications
- **High-Performance**: Optimized for large-scale parallel processing

## Quick Example

```julia
using RootSolvers

# Simple scalar root finding
sol = find_zero(x -> x^2 - 4, SecantMethod(0.0, 3.0))

# Broadcasting over arrays
x0 = rand(100, 100)
x1 = rand(100, 100)
f(x) = x.^2 .- 2.0
sol = find_zero.(f, SecantMethod(x0, x1), CompactSolution())
```
