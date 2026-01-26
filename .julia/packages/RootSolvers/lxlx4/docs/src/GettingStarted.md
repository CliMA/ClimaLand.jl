# Getting Started

RootSolvers.jl is a Julia package for finding roots of nonlinear equations using robust, efficient, and GPU-capable numerical methods. It provides a simple, unified interface for a variety of classic root-finding algorithms, with flexible convergence criteria and solution reporting.

---

## Installation

The package is registered in the Julia General registry.

**Stable Release:**
```julia
using Pkg
Pkg.add("RootSolvers")
```

---

## Quick Start Example
```@example howto
using RootSolvers
# Find the root of x^2 - 100^2 using the secant method
sol = find_zero(x -> x^2 - 100^2, SecantMethod(0.0, 1000.0))
sol
```
The numerical value of the root is contained in `sol.root`:
```@example howto
sol.root
```

---

## How-to Guide

This guide shows the basic steps for solving a root-finding problem.

### General Workflow

#### 1. Define Your Function
Write your function as a Julia callable.
```@example howto
f(x) = x^3 - 2x - 5
nothing # hide
```

#### 2. Choose a Root-Finding Method
Pick a method and provide initial guesses. The type parameter (e.g., `Float64`) is often inferred automatically.
```@example howto
# For SecantMethod, provide two initial guesses
method = SecantMethod(1.0, 3.0)
nothing # hide
```

#### 3. (Optional) Set Tolerance and Solution Type
Customize the convergence criteria and the level of detail in the output.
```@example howto
# Stop when iterates are closer than 1e-6
tol = SolutionTolerance(1e-6)

# Request detailed output for debugging
soltype = VerboseSolution()
nothing # hide
```

#### 4. Call `find_zero`
All arguments after `method` are optional.
```@example howto
sol = find_zero(f, method, soltype, tol)
```

#### 5. Interpret Results
- `sol.converged`: `true` if a root was found.
- `sol.root`: The root value.
- `sol.err`, `sol.iter_performed`, `sol.root_history` (available with [`VerboseSolution`](@ref)).


### Specific Example: Newton's Method with a Provided Derivative

When using [`NewtonsMethod`](@ref), you must provide a function that returns both the value `f(x)` and its derivative `f'(x)` as a tuple. This avoids the overhead of automatic differentiation and is highly efficient if you can provide an analytical derivative.

#### 1. Define Function and Derivative
```@example howto
# This function finds the root of f(x) = x^2 - 4.
# It returns the tuple (f(x), f'(x)).
f_with_deriv(x) = (x^2 - 4, 2x)
nothing # hide
```

#### 2. Choose the Method and Call `find_zero`
```@example howto
# Provide a single initial guess for Newton's method
method = NewtonsMethod(1.0)

# The function f_with_deriv is passed to find_zero
sol = find_zero(f_with_deriv, method)
```

### Specific Example: Brent's Method for Robust Root Finding

Brent's method combines the bisection method, secant method, and inverse quadratic interpolation. It provides superlinear convergence while maintaining the robustness of bracketing methods.

#### 1. Define Your Function
```@example howto
# This function finds the root of f(x) = x^3 - 2.
f(x) = x^3 - 2
nothing # hide
```

#### 2. Choose the Method and Call `find_zero`
```@example howto
# Provide a bracketing interval where f(x0) and f(x1) have opposite signs
method = BrentsMethod(-1.0, 2.0)  # f(-1) = -3, f(2) = 6

# Solve the root-finding problem
sol = find_zero(f, method)
```

---

## Automatic Differentiation and Dual Number Support 🔄

RootSolvers.jl is fully compatible with automatic differentiation frameworks, making it suitable for integration into differentiable models and optimization problems. The package supports dual numbers (from ForwardDiff.jl and other AD packages) as input arguments, allowing gradients to flow through root-finding computations. Dual number support works on GPU arrays when using compatible AD frameworks.

### Using RootSolvers in Differentiable Models

When your function `f(x)` accepts dual numbers, RootSolvers.jl can be used within larger differentiable computations:

```@example
using RootSolvers, ForwardDiff

# Create a function that uses root finding
function solve_and_evaluate(θ)
    # θ is a parameter that affects the root-finding problem
    f(x) = x^3 - θ * x - 5
    sol = find_zero(f, SecantMethod(1.0, 3.0))
    return sol.root
end

# Compute the derivative with respect to θ
θ = 2.0
deriv = ForwardDiff.derivative(solve_and_evaluate, θ)
println("Derivative: ", deriv)
```
This enables integration, for example, with derivative-based optimization algorithms, when an objective function may include a root finding problem. 

## High-Performance and GPU Computing 🚀

RootSolvers.jl is designed for high-performance computing, supporting broadcasting over custom data structures and GPU acceleration. This makes it ideal for solving many problems in parallel.

### Broadcasting with Abstract Types
The package works seamlessly with any abstract type that supports broadcasting, making it well-suited for scientific domains like climate modeling.

**Example: Solving over a custom field type**
```@example howto
using RootSolvers

# Example using regular arrays to represent a field grid
x0 = rand(10, 10)  # A 10x10 field of initial guesses
x1 = x0 .+ 1       # A second field of initial guesses

# Define a function that operates element-wise on the field
f(x) = x^2 - 2

# Solve the root-finding problem across the entire field
method = SecantMethod(x0, x1)
sol = find_zero.(f, method, CompactSolution()) # sol is an Array of structs
```

Use `getproperty.()` to extract the fields from each struct in the array:
```@example howto
converged_field = getproperty.(sol, :converged)
root_field = getproperty.(sol, :root)

println("All converged: ", all(converged_field))
println("Root field shape: ", size(root_field))
```

### GPU Acceleration for Batch Processing
You can achieve significant speedups by running large batches of problems on a GPU.
!!! note "GPU Backends"
    The following examples use ['CUDA.jl`](https://cuda.juliagpu.org/stable/), but similar results\
    can be achieved for different GPU backends with [`KernelAbstractions.jl`](https://juliagpu.github.io/KernelAbstractions.jl/stable/).



**GPU Usage Tips:**
- **Use[`CompactSolution`](@ref):** Only [`CompactSolution`](@ref) is GPU-friendly. [`VerboseSolution`](@ref) is for CPU debugging only.
- **GPU-Compatible Function:** Ensure your function `f(x)` uses only GPU-supported operations.
- **Minimize Data Transfer:** Keep initial guesses and results on the GPU.

**Broadcasting Example: 1 Million problems on the GPU**
```julia
using RootSolvers, CUDA

# Create GPU arrays for batch processing
x0 = CUDA.fill(1.0f0, 1000, 1000)  # 1M initial guesses on GPU
x1 = CUDA.fill(2.0f0, 1000, 1000)  # Second initial guesses

# Define GPU-compatible function
f(x) = x^3 - x - 2

# Solve all problems in parallel using broadcasting
method = SecantMethod(x0, x1) # method = SecantMethod.(x0, x1) is also supported
sol = find_zero.(f, method, CompactSolution()) # broadcast launches kernel

# Results are on the GPU as an array of CompactSolutions
converged_field = map(sol_i -> sol_i.converged, sol)
root_field = map(sol_i -> sol_i.root, sol)

println("All converged: ", all(converged_field)) # Ouput: "All converged: true"
println("Root field shape: ", size(root_field)) # Output "Root field shape: (1000, 1000)"
```

**Map Example: 1 Million problems on the GPU**
```julia
using RootSolvers, CUDA

# Create GPU arrays for batch processing
x0 = CUDA.fill(1.0f0, 1000, 1000)  # 1M initial guesses on GPU
x1 = CUDA.fill(2.0f0, 1000, 1000)  # Second initial guesses

# Define GPU-compatible function
f(x) = x^3 - x - 2

# Solve all problems in parallel using map
const METHOD = SecantMethod
sol = map(x0, x1) do x0, x1 # map launches kernel
    find_zero(f, METHOD(x0, x1), CompactSolution())
end

# Results are on the GPU as an array of CompactSolutions
converged_field = map(sol_i -> sol_i.converged, sol)
root_field = map(sol_i -> sol_i.root, sol)

println("All converged: ", all(converged_field)) # Ouput: "All converged: true"
println("Root field shape: ", size(root_field)) # Output "Root field shape: (1000, 1000)"
```

---

## Reference Tables

### Available Root-Finding Methods

| Method | Requirements | Best For |
| :--- | :--- | :--- |
| [`SecantMethod`](@ref) | 2 initial guesses | No derivatives, **fast** convergence|
| [`RegulaFalsiMethod`](@ref) | Bracketing interval | **Guaranteed** convergence |
| [`BisectionMethod`](@ref) | Bracketing interval | **Guaranteed** convergence, simple |
| [`BrentsMethod`](@ref) | Bracketing interval | **Superlinear** convergence, robust |
| [`NewtonsMethodAD`](@ref) | 1 initial guess, differentiable `f` | **Fastest**, uses autodiff, robust step control |
| [`NewtonsMethod`](@ref) | 1 initial guess, `f` and `f'` provided | **Analytical** derivatives, robust step control |

### Available Tolerance Types

| Tolerance Type | Criterion | Best For |
| :--- | :--- | :--- |
| [`SolutionTolerance`](@ref) | `abs(x₂ - x₁)` | When you want iterates to **stabilize** |
| [`ResidualTolerance`](@ref) | `abs(f(x))` | When you want the function value near **zero** |
| [`RelativeSolutionTolerance`](@ref) | `abs((x₂ - x₁)/x₁)` | When root magnitude **varies widely** |
| [`RelativeOrAbsolute...`](@ref RelativeOrAbsoluteSolutionTolerance)| Relative or Absolute | **Robust** for both small and large roots |

### Available Solution Types

| Solution Type | Features | Best For |
| :--- | :--- | :--- |
| [`CompactSolution`](@ref) | Minimal output, GPU-friendly | **High-performance**, GPU, memory efficiency |
| [`VerboseSolution`](@ref) | Full diagnostics, iteration history | **Debugging**, analysis, CPU |

### Advanced Features

| Feature | Description | Use Cases |
| :--- | :--- | :--- |
| **Dual Number Support** | Compatible with automatic differentiation | **Differentiable models**, optimization, gradient-based learning |
| **GPU Acceleration** | Full CUDA.jl support with broadcasting | **Large-scale parallel processing**, batch computations |
| **Custom Field Types** | Works with any broadcastable type | **Scientific computing**, climate modeling, custom data structures |

---

## Troubleshooting
- If not converging, try different initial guesses or a bracketing method such as [`BrentsMethod`](@ref).
- Use [`VerboseSolution()`](@ref) to inspect the iteration history and diagnose issues.
- Adjust the tolerance for stricter or looser convergence criteria.

## Extending RootSolvers.jl

If you want to add custom root-finding methods, tolerance types, or solution formats, see the [Developer Documentation](DeveloperDocs.md) for detailed guidance on extending the package.
