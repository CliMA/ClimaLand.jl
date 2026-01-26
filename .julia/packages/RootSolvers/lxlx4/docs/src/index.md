# RootSolvers.jl

RootSolvers.jl is a Julia package for finding roots of nonlinear equations using robust, efficient, and GPU-capable numerical methods. It provides a simple, unified interface for a variety of classic root-finding algorithms, with flexible convergence criteria and solution reporting. The package supports dual numbers for automatic differentiation, making it suitable for integration into differentiable models and optimization problems.

- [Getting Started](GettingStarted.md): Installation, quick start, and how-to guide
- [API Reference](API.md): Full documentation of all methods and types

## Quick Example
See the [Getting Started](GettingStarted.md) page for more details and examples.

Install stable release:
```julia
using Pkg
Pkg.add("RootSolvers")
```

Find a root of a quadratic equation:
```@example howto
using RootSolvers

# Find the root of x^2 - 100^2 using the secant method
sol = find_zero(x -> x^2 - 100^2, SecantMethod(0.0, 1000.0))
```

Or use Brent's method for robust bracketing
```@example howto
sol = find_zero(x -> x^2 - 100^2, BrentsMethod(-200.0, 0.0))
```

## Documentation
- [Getting Started](GettingStarted.md)
- [API Reference](API.md)
- [Developer Documentation](DeveloperDocs.md)
