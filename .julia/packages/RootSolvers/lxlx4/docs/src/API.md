# API Reference

```@meta
CurrentModule = RootSolvers
```

This page provides the complete API reference for RootSolvers.jl. For a more narrative introduction, see the Getting Started guide.

## Module Overview

```@docs
RootSolvers
```

## Main Function

The primary entry point for the package is the `find_zero` function.

```@docs
find_zero
```

---

## Numerical Methods

The following structs are used to select the root-finding algorithm.

| Method              | Requirements                        | Best For                        |
|:--------------------|:------------------------------------|:--------------------------------|
| [`SecantMethod`](@ref)      | 2 initial guesses                   | No derivatives, **fast** convergence|
| [`RegulaFalsiMethod`](@ref) | Bracketing interval (sign change)   | **Guaranteed** convergence      |
| [`BisectionMethod`](@ref)   | Bracketing interval (sign change)   | **Guaranteed** convergence, simple |
| [`BrentsMethod`](@ref)      | Bracketing interval (sign change)   | **Superlinear** convergence, robust |
| [`NewtonsMethodAD`](@ref)   | 1 initial guess, differentiable f   | **Fastest**, uses autodiff, robust step control |
| [`NewtonsMethod`](@ref)     | 1 initial guess, f and f' provided  | **Analytical** derivatives, robust step control |

```@docs
SecantMethod
RegulaFalsiMethod
BisectionMethod
BrentsMethod
NewtonsMethodAD
NewtonsMethod
```

---

## Solution Types

These types control the level of detail in the output returned by `find_zero`.

| Solution Type      | Features                              | Best For                        |
|:-------------------|:--------------------------------------|:--------------------------------|
| [`CompactSolution`](@ref)  | Minimal output, GPU-friendly          | **High-performance**, GPU, memory efficiency |
| [`VerboseSolution`](@ref)  | Full diagnostics, iteration history   | **Debugging**, analysis, CPU    |

```@docs
CompactSolution
VerboseSolution
```

---

## Tolerance Types

Tolerance types define the convergence criteria for the solver.

| Tolerance Type                        | Criterion                        | Best For                                 |
|:--------------------------------------|:---------------------------------|:-----------------------------------------|
| [`SolutionTolerance`](@ref)           )                   | `abs(x₂ - x₁)`                   | When you want iterates to **stabilize** |
| [`ResidualTolerance`](@ref)                   | `abs(f(x))`                      | When you want the function value near **zero** |
| [`RelativeSolutionTolerance`](@ref)           | `abs((x₂ - x₁)/x₁)`              | When root magnitude **varies widely** |
| [`RelativeOrAbsoluteSolutionTolerance`](@ref) | Relative or Absolute             | **Robust** for both small and large roots    |

```@docs
AbstractTolerance
SolutionTolerance
ResidualTolerance
RelativeSolutionTolerance
RelativeOrAbsoluteSolutionTolerance
```

---

## Broadcasting and High-Performance Computing

RootSolvers.jl is designed for high-performance computing applications and supports broadcasting to solve many problems in parallel. This is especially useful for GPU arrays or custom field types used in scientific modeling.

The custom broadcasting rule unpacks initial guesses from the `method` struct while treating all other arguments as scalars. This enables a clean API for batch-solving.

For more information about broadcasting, see the examples in the `find_zero` documentation.

---

## Developer Documentation

For information about internal methods, extending RootSolvers.jl, and developer-focused functionality, see the [Developer Documentation](DeveloperDocs.md).
