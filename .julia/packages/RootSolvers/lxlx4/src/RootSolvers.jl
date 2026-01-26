"""
RootSolvers.jl

A Julia package for solving roots of non-linear equations using various numerical methods.
Contains functions for finding zeros of scalar functions using robust iterative algorithms.

The main entry point is [`find_zero`](@ref), which supports multiple root-finding methods
and tolerance criteria.

## Supported Methods
- **Secant Method**: Requires two initial guesses, uses linear interpolation
- **Bisection Method**: Requires bracketing interval with sign change, converges linearly
- **Regula Falsi Method**: Requires bracketing interval with sign change
- **Brent's Method**: Requires bracketing interval, combines bisection, secant, and inverse quadratic interpolation
- **Newton's Method with AD**: Requires one initial guess, uses automatic differentiation
- **Newton's Method**: Requires one initial guess and user-provided derivative

## Method Selection Guide
- **Bracketing methods** (Bisection, Regula Falsi, Brent's): Use when you know an interval containing the root
- **Brent's method**: Recommended for robust, guaranteed convergence with superlinear rate
- **Secant method**: Good when you have two guesses but no bracketing interval
- **Newton's methods**: Fastest convergence when you have a good initial guess

## Example

```julia
using RootSolvers

# Find the square root of a quadratic equation using the secant method
sol = find_zero(x -> x^2 - 100^2,
               SecantMethod{Float64}(0.0, 1000.0),
               CompactSolution());

println(sol)
# CompactSolutionResults{Float64}:
# ├── Status: converged
# └── Root: 99.99999999994358

# Access the root value
root_value = sol.root  # 99.99999999994358

# Use Brent's method for robust convergence with bracketing interval
sol_brent = find_zero(x -> x^3 - 8,
                     BrentsMethod{Float64}(-1.0, 3.0),
                     CompactSolution());

# Use Newton's method with automatic differentiation for faster convergence
sol_newton = find_zero(x -> x^3 - 27,
                      NewtonsMethodAD{Float64}(2.0),
                      VerboseSolution());
```

"""
module RootSolvers

export find_zero,
    RootSolvingMethod,
    SecantMethod,
    RegulaFalsiMethod,
    BrentsMethod,
    BisectionMethod,
    NewtonsMethodAD,
    NewtonsMethod
export SolutionType, CompactSolution, VerboseSolution
export ResidualTolerance,
    SolutionTolerance,
    RelativeSolutionTolerance,
    RelativeOrAbsoluteSolutionTolerance
export method_args, value_deriv, default_tol

import ForwardDiff
import Printf: @printf

base_type(::Type{FT}) where {FT} = FT
base_type(::Type{FT}) where {T, FT <: ForwardDiff.Dual{<:Any, T}} = base_type(T)
base_type(::Type{FT}) where {T, FT <: AbstractArray{T}} = base_type(T)

# Input types
const FTypes = Union{Real, AbstractArray}

"""
    RootSolvingMethod{FT} <: AbstractType

Abstract type for root-finding methods in RootSolvers.jl.

This is the base type for all numerical methods used to find roots of scalar functions.
Each concrete method type should implement the `method_args` function to extract
initial guesses and parameters for the solver.

## Type Parameters
- `FT`: The floating-point type for computations (e.g., `Float64`, `Float32`)

## Concrete Implementations
- [`SecantMethod{FT}`](@ref): Linear interpolation between two points
- [`RegulaFalsiMethod{FT}`](@ref): Bracketing method with linear interpolation
- [`BisectionMethod{FT}`](@ref): Simple bracketing method with guaranteed convergence
- [`BrentsMethod{FT}`](@ref): Robust bracketing method combining multiple techniques
- [`NewtonsMethodAD{FT}`](@ref): Newton's method with automatic differentiation
- [`NewtonsMethod{FT}`](@ref): Newton's method with user-provided derivatives

## Interface Requirements
All concrete subtypes must implement:
- `method_args(method)`: Return initial guesses as a tuple
- `find_zero(f, method, args..., soltype, tol, maxiters)`: Main solver implementation

## Example
```julia
# Define a custom method
struct MyCustomMethod{FT} <: RootSolvingMethod{FT}
    x0::FT
end

# Implement required interface
method_args(method::MyCustomMethod) = (method.x0,)

# Use with find_zero
sol = find_zero(x -> x^2 - 4, MyCustomMethod{Float64}(1.0))
```
"""
abstract type RootSolvingMethod{FT <: FTypes} end

"""
    SecantMethod{FT} <: RootSolvingMethod{FT}

The secant method for root finding, which uses linear interpolation between two points
to approximate the derivative. This method requires two initial guesses but does not
require the function to be differentiable or the guesses to bracket a root.

The method uses the recurrence relation:
```math
x_{n+1} = x_n - f(x_n) \\frac{x_n - x_{n-1}}{f(x_n) - f(x_{n-1})}
```

## Convergence
- **Order**: Approximately 1.618 (superlinear)
- **Requirements**: Two initial guesses, continuous function
- **Advantages**: No derivative required, fast convergence
- **Disadvantages**: May not converge if initial guesses are poor

## Fields
- `x0::FT`: First initial guess
- `x1::FT`: Second initial guess

## Example
```julia
method = SecantMethod{Float64}(0.0, 2.0)
sol = find_zero(x -> x^3 - 8, method)
```
"""
struct SecantMethod{FT} <: RootSolvingMethod{FT}
    x0::FT
    x1::FT
end

"""
    BisectionMethod{FT} <: RootSolvingMethod{FT}

The bisection method for root finding. This is a simple bracketing method that
divides the interval in half at each iteration. The method requires that `f(x0)`
and `f(x1)` have opposite signs, ensuring that a root exists in the interval `[x0, x1]`.

The method uses the recurrence relation:
```math
x_{n+1} = \\frac{x_0 + x_1}{2}
```
where the interval `[x0, x1]` is updated based on the sign of `f(x_{n+1})`.

## Convergence
- **Order**: Linear (slower than Newton's method)
- **Requirements**: Bracketing interval with `f(x0) * f(x1) < 0`
- **Advantages**: Guaranteed convergence, simple implementation
- **Disadvantages**: Slower convergence than other methods

## Fields
- `x0::FT`: Lower bound of bracketing interval
- `x1::FT`: Upper bound of bracketing interval

## Example
```julia
# Find root of x^3 - 2 in interval [-1, 2]
method = BisectionMethod{Float64}(-1.0, 2.0)
sol = find_zero(x -> x^3 - 2, method)
```
"""
struct BisectionMethod{FT} <: RootSolvingMethod{FT}
    x0::FT
    x1::FT
end

"""
    RegulaFalsiMethod{FT} <: RootSolvingMethod{FT}

The Regula Falsi (false position) method for root finding. This is a bracketing method
that maintains the sign change property and uses linear interpolation to find the root.

The method requires that `f(x0)` and `f(x1)` have opposite signs, ensuring that
a root exists in the interval `[x0, x1]`.

## Convergence
- **Order**: Linear (slower than Newton's method)
- **Requirements**: Bracketing interval with `f(x0) * f(x1) < 0`
- **Advantages**: Guaranteed convergence, robust
- **Disadvantages**: Slower convergence than other methods

## Fields
- `x0::FT`: Lower bound of bracketing interval
- `x1::FT`: Upper bound of bracketing interval

## Example
```julia
# Find root of x^3 - 2 in interval [-1, 2]
method = RegulaFalsiMethod{Float64}(-1.0, 2.0)
sol = find_zero(x -> x^3 - 2, method)
```
"""
struct RegulaFalsiMethod{FT} <: RootSolvingMethod{FT}
    x0::FT
    x1::FT
end

"""
    BrentsMethod{FT} <: RootSolvingMethod{FT}

Brent's method for root finding, which combines the bisection method, secant method, and inverse quadratic interpolation.
This is a bracketing method that maintains the sign change property and provides superlinear convergence.

The method requires that `f(x0)` and `f(x1)` have opposite signs, ensuring that
a root exists in the interval `[x0, x1]`.

## Convergence
- **Order**: Superlinear (faster than Regula Falsi)
- **Requirements**: Bracketing interval with `f(x0) * f(x1) < 0`
- **Advantages**: Guaranteed convergence, fast convergence, robust
- **Disadvantages**: More complex than simpler bracketing methods

## Fields
- `x0::FT`: Lower bound of bracketing interval
- `x1::FT`: Upper bound of bracketing interval

## Example
```julia
# Find root of x^3 - 2 in interval [-1, 2]
method = BrentsMethod{Float64}(-1.0, 2.0)
sol = find_zero(x -> x^3 - 2, method)
```
"""
struct BrentsMethod{FT} <: RootSolvingMethod{FT}
    x0::FT
    x1::FT
end

"""
    NewtonsMethodAD{FT} <: RootSolvingMethod{FT}

Newton's method for root finding using automatic differentiation to compute derivatives.
This method provides quadratic convergence when close to the root and the derivative
is non-zero. The implementation includes step size limiting and backtracking line search
for robustness.

The method uses the iteration
```math
x_{n+1} = x_n - \\frac{f(x_n)}{f'(x_n)}
```
where the derivative `f'(x_n)` is computed using ForwardDiff.jl.

## Convergence
- **Order**: Quadratic (very fast near the root)
- **Requirements**: Differentiable function, good initial guess
- **Advantages**: Fast convergence, automatic derivative computation, robust step size control
- **Disadvantages**: May not converge if initial guess is poor or derivative is zero

## Fields
- `x0::FT`: Initial guess for the root

## Example
```julia
# Find cube root of 27
method = NewtonsMethodAD{Float64}(2.0)
sol = find_zero(x -> x^3 - 27, method)
```
"""
struct NewtonsMethodAD{FT} <: RootSolvingMethod{FT}
    x0::FT
end

"""
    NewtonsMethod{FT} <: RootSolvingMethod{FT}

Newton's method for root finding where the user provides both the function and its derivative.
This method provides quadratic convergence when close to the root. The implementation includes
step size limiting and backtracking line search for robustness.

The method uses the iteration
```math
x_{n+1} = x_n - \\frac{f(x_n)}{f'(x_n)}
```

## Convergence
- **Order**: Quadratic (very fast near the root)
- **Requirements**: Function and derivative, good initial guess
- **Advantages**: Fast convergence, no automatic differentiation overhead, robust step size control
- **Disadvantages**: Requires manual derivative computation

## Fields
- `x0::FT`: Initial guess for the root

## Note
When using this method, your function `f` should return a tuple `(f(x), f'(x))` containing
both the function value and its derivative at `x`.

## Example
```julia
# Find root of x^2 - 4, providing both function and derivative
f_and_df(x) = (x^2 - 4, 2x)
method = NewtonsMethod{Float64}(1.0)
sol = find_zero(f_and_df, method)
```
"""
struct NewtonsMethod{FT} <: RootSolvingMethod{FT}
    x0::FT
end

"""
    SolutionType <: AbstractType

Abstract type for solution formats in RootSolvers.jl.

This is the base type for all solution types that control the level of detail
returned by `find_zero`. Each concrete solution type should implement the
`SolutionResults` constructor and history management functions.

## Concrete Implementations
- [`CompactSolution`](@ref): Memory-efficient, GPU-compatible solution with minimal output
- [`VerboseSolution`](@ref): CPU-only solution with detailed diagnostics and iteration history

## Interface Requirements
All concrete subtypes must implement:
- `SolutionResults(soltype, args...)`: Constructor for the corresponding results type
- `init_history(soltype, x)`: Initialize history storage for the solution type
- `push_history!(history, x, soltype)`: Add a value to the history

## Example
```julia
# Define a custom solution type
struct MySolution <: SolutionType end

# Implement required interface
SolutionResults(soltype::MySolution, args...) = MySolutionResults(args...)
init_history(::MySolution, x) = [x]
push_history!(history, x, ::MySolution) = push!(history, x)

# Use with find_zero
sol = find_zero(x -> x^2 - 4, SecantMethod(1.0, 3.0), MySolution())
```
"""
abstract type SolutionType end
Base.broadcastable(soltype::SolutionType) = Ref(soltype)

"""
    VerboseSolution <: SolutionType

A solution type that returns detailed information about the root-finding process,
including iteration history and convergence diagnostics.

When used with `find_zero`, returns a `VerboseSolutionResults` object
containing the root, convergence status, error information, and complete iteration history.

## Accessing Results
The returned `VerboseSolutionResults` object contains the following fields:
- `sol.root`: The found root value
- `sol.converged`: Boolean indicating if the method converged
- `sol.err`: Final error value (function value at the root)
- `sol.iter_performed`: Number of iterations performed
- `sol.root_history`: Vector of all root values during iteration
- `sol.err_history`: Vector of all error values during iteration

## Note
This solution type stores iteration history and is primarily intended for CPU computations.
For GPU computations or when memory usage is a concern, use `CompactSolution` instead.

## Example
```julia
sol = find_zero(x -> x^2 - 4, SecantMethod(1.0, 3.0), VerboseSolution())

# Access the root
root_value = sol.root

# Check convergence and diagnostics
if sol.converged
    println("Root found: ", root_value)
    println("Converged in ", sol.iter_performed, " iterations")
    println("Final error: ", sol.err)
else
    println("Method failed to converge")
end

# Access iteration history
println("First iteration root: ", sol.root_history[1])
println("Last iteration root: ", sol.root_history[end])
```
"""
struct VerboseSolution <: SolutionType end

abstract type AbstractSolutionResults{Real} end

struct VerboseSolutionResults{FT} <: AbstractSolutionResults{FT}
    "solution ``x^*`` of the root of the equation ``f(x^*) = 0``"
    root::FT
    "indicates convergence"
    converged::Bool
    "error of the root of the equation ``f(x^*) = 0``"
    err::FT
    "number of iterations performed"
    iter_performed::Int
    "solution per iteration"
    root_history::Vector{FT}
    "error of the root of the equation ``f(x^*) = 0`` per iteration"
    err_history::Vector{FT}
end
SolutionResults(soltype::VerboseSolution, args...) =
    VerboseSolutionResults(args...)

function Base.show(io::IO, sol::VerboseSolutionResults{FT}) where {FT}
    color = sol.converged ? :green : :red
    status = sol.converged ? "converged" : "failed to converge"
    println(io, "VerboseSolutionResults{$FT}:")
    print(io, "├── Status: ")
    printstyled(io, status; color)
    println(io)
    println(io, "├── Root: ", sol.root)
    println(io, "├── Error: ", sol.err)
    println(io, "├── Iterations: ", sol.iter_performed)
    println(io, "└── History:")
    n_iters = length(sol.root_history)
    for i in 1:n_iters
        if n_iters > 20 && 9 < i < n_iters - 9
            i == 11 && println(io, "    ⋮            ⋮                ⋮")
            i == 12 && println(io, "    ⋮            ⋮                ⋮")
            continue
        end
        prefix = i == n_iters ? "    └──" : "    ├──"
        @printf(
            io,
            "%s iter %2d: x = %8.5g, err = %.4g\n",
            prefix,
            i,
            sol.root_history[i],
            sol.err_history[i]
        )
    end
end

"""
    CompactSolution <: SolutionType

A memory-efficient solution type that returns only essential information about the root.

When used with `find_zero`, returns a `CompactSolutionResults` object
containing only the root value and convergence status. This solution type is GPU-compatible
and suitable for high-performance applications where memory usage is critical.

## Accessing Results
The returned `CompactSolutionResults` object contains the following fields:
- `sol.root`: The found root value
- `sol.converged`: Boolean indicating if the method converged

## Example
```julia
sol = find_zero(x -> x^2 - 4, 
               SecantMethod{Float64}(0.0, 3.0), 
               CompactSolution())

# Access the root
println("Root: \$(sol.root)")

# Check convergence
if sol.converged
    println("Root found successfully!")
else
    println("Method failed to converge")
end
```
"""
struct CompactSolution <: SolutionType end

struct CompactSolutionResults{FT} <: AbstractSolutionResults{FT}
    "solution ``x^*`` of the root of the equation ``f(x^*) = 0``"
    root::FT
    "indicates convergence"
    converged::Bool
end
SolutionResults(soltype::CompactSolution, root, converged, args...) =
    CompactSolutionResults(root, converged)

function Base.show(io::IO, sol::CompactSolutionResults{FT}) where {FT}
    color = sol.converged ? :green : :red
    if get(io, :compact, false)
        # compact printing within a vector or similar
        printstyled(io, sol.root; color)
    else
        # standalone printing
        status = sol.converged ? "converged" : "failed to converge"
        println(io, "CompactSolutionResults{$FT}:")
        print(io, "├── Status: ")
        printstyled(io, status; color)
        println(io)
        print(io, "└── Root: ", sol.root)
    end
end

init_history(::VerboseSolution, x::FT) where {FT <: Real} = FT[x]
init_history(::CompactSolution, x) = nothing
init_history(::VerboseSolution, ::Type{FT}) where {FT <: Real} = FT[]
init_history(::CompactSolution, ::Type{FT}) where {FT <: Real} = nothing

function push_history!(
    history::Vector{FT},
    x::FT,
    ::VerboseSolution,
) where {FT <: Real}
    push!(history, x)
end
function push_history!(
    history::Nothing,
    x::FT,
    ::CompactSolution,
) where {FT <: Real}
    nothing
end

"""
    AbstractTolerance{FT} <: AbstractType

Abstract type for tolerance criteria in RootSolvers.jl.

This is the base type for all tolerance types that define convergence criteria
for root-finding algorithms. Each concrete tolerance type should implement
the callable interface `(tol)(x1, x2, y)` for convergence checking.

## Type Parameters
- `FT`: The floating-point type for tolerance values (e.g., `Float64`, `Float32`)

## Concrete Implementations
- [`ResidualTolerance{FT}`](@ref): Based on `|f(x)|`
- [`SolutionTolerance{FT}`](@ref): Based on `|x_{n+1} - x_n|`
- [`RelativeSolutionTolerance{FT}`](@ref): Based on `|(x_{n+1} - x_n)/x_n|`
- [`RelativeOrAbsoluteSolutionTolerance{FT}`](@ref): Combined relative and absolute tolerance

## Interface Requirements
All concrete subtypes must implement:
- `(tol)(x1, x2, y)`: Check convergence using three arguments (previous iterate, current iterate, function value)

## Example
```julia
# Define a custom tolerance type
struct MyTolerance{FT} <: AbstractTolerance{FT}
    threshold::FT
end

# Implement the required interface
(tol::MyTolerance)(x1, x2, y) = abs(x2 - x1) < tol.threshold || abs(y) < eps(typeof(y))

# Use with find_zero
sol = find_zero(x -> x^2 - 4, SecantMethod(1.0, 3.0), CompactSolution(), MyTolerance(1e-6))
```
"""
abstract type AbstractTolerance{FT <: FTypes} end
Base.broadcastable(tol::AbstractTolerance) = Ref(tol)

"""
    ResidualTolerance{FT} 

A convergence criterion based on the absolute value of the residual (function value).
The iteration stops when `|f(x)| < tol`, where `tol` is the specified tolerance (limited to 
the machine epsilon of the function value type).

This tolerance is appropriate when you want to ensure that the function value is
sufficiently close to zero, regardless of how close consecutive iterates are.

## Fields
- `tol::FT`: Tolerance threshold for `|f(x)|`

## Example
```julia
tol = ResidualTolerance(1e-10)
sol = find_zero(x -> x^2 - 4, 
               NewtonsMethodAD{Float64}(1.0),
               CompactSolution(),
               tol)
```
"""
struct ResidualTolerance{FT} <: AbstractTolerance{FT}
    tol::FT
end

"""
    (tol::ResidualTolerance)(x1, x2, y)

Evaluates residual tolerance, based on ``|f(x)|``, limiting the residual to the machine 
epsilon of the function value type.
"""
(tol::ResidualTolerance)(x1, x2, y) =
    abs(y) < tol.tol || abs(y) < eps(typeof(y))


"""
    SolutionTolerance{FT} 

A convergence criterion based on the absolute difference between consecutive iterates.
The iteration stops when `|x_{n+1} - x_n| < tol`, where `tol` is the specified tolerance.
Convergence is also triggered if |f(x)| is smaller than the machine epsilon for the value type.

This tolerance is appropriate when you want to ensure that consecutive iterates are
sufficiently close, indicating that the solution has stabilized.

## Fields
- `tol::FT`: Tolerance threshold for `|x_{n+1} - x_n|`

## Example
```julia
tol = SolutionTolerance(1e-8)
sol = find_zero(x -> x^3 - 8, 
               SecantMethod{Float64}(1.0, 3.0),
               CompactSolution(),
               tol)
```
"""
struct SolutionTolerance{FT} <: AbstractTolerance{FT}
    tol::FT
end

"""
    (tol::SolutionTolerance)(x1, x2, y)

Evaluates solution tolerance, based on ``|x2-x1|``
"""
(tol::SolutionTolerance)(x1, x2, y) =
    abs(x2 - x1) < tol.tol || abs(y) < eps(typeof(y))

"""
    RelativeSolutionTolerance{FT} 

A convergence criterion based on the relative difference between consecutive iterates.
The iteration stops when `|(x_{n+1} - x_n)/x_n| < tol`, where `tol` is the specified tolerance. 
Convergence is also triggered if |f(x)| is smaller than the machine epsilon for the value type.

This tolerance is appropriate when you want to convergence relative to the magnitude of the solution,
which is useful when the root value might be very large or very small.

## Fields
- `tol::FT`: Relative tolerance threshold

## Warning
This tolerance criterion can fail if `x_n ≈ 0` during iteration, as it involves division by `x_n`.
Consider using `RelativeOrAbsoluteSolutionTolerance` for more robust behavior.

## Example
```julia
tol = RelativeSolutionTolerance(1e-6)
sol = find_zero(x -> x^2 - 1e6, 
               NewtonsMethodAD{Float64}(500.0),
               CompactSolution(),
               tol)
```
"""
struct RelativeSolutionTolerance{FT} <: AbstractTolerance{FT}
    tol::FT
end

"""
    (tol::RelativeSolutionTolerance)(x1, x2, y)

Evaluates solution tolerance, based on ``|(x2-x1)/x1|``
"""
(tol::RelativeSolutionTolerance)(x1, x2, y) =
    abs((x2 - x1) / x1) < tol.tol || abs(y) < eps(typeof(y))

"""
    RelativeOrAbsoluteSolutionTolerance{FT} 

A robust convergence criterion combining both relative and absolute tolerances.
The iteration stops when either `|(x_{n+1} - x_n)/x_n| < rtol` OR `|x_{n+1} - x_n| < atol`.

This tolerance provides robust behavior across different scales of root values:
- The relative tolerance `rtol` ensures accuracy for large roots
- The absolute tolerance `atol` ensures convergence when the root is near zero

## Fields
- `rtol::FT`: Relative tolerance threshold
- `atol::FT`: Absolute tolerance threshold

## Example
```julia
# Use relative tolerance of 1e-6 and absolute tolerance of 1e-10
tol = RelativeOrAbsoluteSolutionTolerance(1e-6, 1e-10)
sol = find_zero(x -> x^2 - 1e-8, 
               NewtonsMethodAD{Float64}(1e-3),
               CompactSolution(),
               tol)
```
"""
struct RelativeOrAbsoluteSolutionTolerance{FT} <: AbstractTolerance{FT}
    rtol::FT
    atol::FT
end

"""
    (tol::RelativeOrAbsoluteSolutionTolerance)(x1, x2, y)

Evaluates combined relative and absolute tolerance, based
on ``|(x2-x1)/x1| || |x2-x1|``
"""
(tol::RelativeOrAbsoluteSolutionTolerance)(x1, x2, y) =
    abs((x2 - x1) / x1) < tol.rtol ||
    abs(x2 - x1) < tol.atol ||
    abs(y) < eps(typeof(y))

"""
    find_zero(f, method, soltype=CompactSolution(), tol=nothing, maxiters=1_000)

Find a root of the scalar function `f` using the specified numerical method.

This is the main entry point for root finding in RootSolvers.jl. Given a function `f`,
it finds a value `x` such that `f(x) ≈ 0` using iterative numerical methods. The function
supports various root-finding algorithms, tolerance criteria, and solution formats.

## Arguments
- `f::Function`: The function for which to find a root. Should take a scalar input and return a scalar output.
- `method::RootSolvingMethod`: The numerical method to use. Available methods:
    - [`BisectionMethod`](@ref): Bracketing method maintaining sign change (linear convergence, guaranteed)
    - [`SecantMethod`](@ref): Uses linear interpolation between two points (superlinear convergence)
    - [`RegulaFalsiMethod`](@ref): Bracketing method maintaining sign change (linear convergence, guaranteed)
    - [`NewtonsMethodAD`](@ref): Newton's method with automatic differentiation (quadratic convergence)
    - [`NewtonsMethod`](@ref): Newton's method with user-provided derivative (quadratic convergence)
- `soltype::`[`SolutionType`](@ref): Format of the returned solution (default: [`CompactSolution()`](@ref)):
    - [`CompactSolution`](@ref): Returns only root and convergence status (GPU-compatible)
    - [`VerboseSolution`](@ref): Returns detailed diagnostics and iteration history (CPU-only)
- `tol::Union{Nothing, AbstractTolerance}`: Convergence criterion (default: [`SolutionTolerance(1e-3)`](@ref)):
    - [`ResidualTolerance`](@ref): Based on `|f(x)|`
    - [`SolutionTolerance`](@ref): Based on `|x_{n+1} - x_n|`
    - [`RelativeSolutionTolerance`](@ref): Based on `|(x_{n+1} - x_n)/x_n|`
    - [`RelativeOrAbsoluteSolutionTolerance`](@ref): Combined relative and absolute tolerance
- `maxiters::Int`: Maximum number of iterations allowed (default: 10,000)

## Returns
- `AbstractSolutionResults`: Solution object containing the root and convergence information.
  The exact type depends on the `soltype` parameter:
  - `CompactSolutionResults`: Contains `root` and `converged` fields
  - `VerboseSolutionResults`: Additionally contains `err`, `iter_performed`, and iteration history

## Examples

```julia
using RootSolvers

# Find square root of 2 using secant method
sol = find_zero(x -> x^2 - 2, SecantMethod{Float64}(1.0, 2.0))
println("√2 ≈ \$(sol.root)")  # √2 ≈ 1.4142135623730951

# Use Newton's method with automatic differentiation for faster convergence
sol = find_zero(x -> x^3 - 27, NewtonsMethodAD{Float64}(2.0))
println("∛27 = \$(sol.root)")  # ∛27 = 3.0

# Get detailed iteration history
sol = find_zero(x -> exp(x) - 2, 
               NewtonsMethodAD{Float64}(0.5), 
               VerboseSolution())
println("ln(2) ≈ \$(sol.root) found in \$(sol.iter_performed) iterations")

# Use custom tolerance
tol = RelativeOrAbsoluteSolutionTolerance(1e-12, 1e-15)
sol = find_zero(x -> cos(x), 
               NewtonsMethodAD{Float64}(1.0), 
               CompactSolution(), 
               tol)
println("π/2 ≈ \$(sol.root)")

# Robust bracketing method for difficult functions
sol = find_zero(x -> x^3 - 2x - 5, RegulaFalsiMethod{Float64}(2.0, 3.0))
```

## Batch and GPU Root-Finding (Broadcasting)

You can broadcast `find_zero` over arrays of methods or initial guesses to solve many root-finding problems in parallel, including on the GPU:

```julia
using CUDA, RootSolvers
x0 = CUDA.fill(1.0, 1000)  # 1000 initial guesses on the GPU
method = SecantMethod.(x0, x0 .+ 1)
# f should be broadcastable over arrays
sol = find_zero.(x -> x.^2 .- 2, method, CompactSolution())
```

This is especially useful for large-scale or batched root-finding on GPUs. Only [`CompactSolution`](@ref) is GPU-compatible.

## Method Selection Guide
- **BisectionMethod**: Simple general-purpose bracketing method, slow but guaranteed convergence
- **SecantMethod**: Good general-purpose method, no derivatives needed
- **RegulaFalsiMethod**: Use when you need guaranteed convergence with a bracketing interval
- **NewtonsMethodAD**: Fastest convergence when derivatives are available via autodiff
- **NewtonsMethod**: Use when you can provide analytical derivatives efficiently

## See Also
- [`BisectionMethod`](@ref), [`SecantMethod`](@ref), [`RegulaFalsiMethod`](@ref), [`NewtonsMethodAD`](@ref), [`NewtonsMethod`](@ref)
- [`CompactSolution`](@ref), [`VerboseSolution`](@ref)
- [`ResidualTolerance`](@ref), [`SolutionTolerance`](@ref)
"""
function find_zero end

# Helper to get the default tolerance for a given type
"""
    default_tol(FT)

Returns the default tolerance for a given type `FT`.
This is a helper function to provide a consistent default tolerance
for different numerical types.

## Arguments
- `FT`: The type of the numerical value (e.g., `Float64`, `ComplexF64`).

## Returns
- `AbstractTolerance`: A default tolerance object.

## Example
```julia
using RootSolvers

# Find the default tolerance for Float64
tol = default_tol(Float64)
println("Default tolerance for Float64: ", tol)
# Default tolerance for Float64: SolutionTolerance{Float64}(1e-4)
```
"""
function default_tol(::Type{Float64})
    return SolutionTolerance{Float64}(1e-4)
end

function default_tol(::Type{FT}) where {FT}
    return SolutionTolerance{base_type(FT)}(1e-3)
end

# Main entry point: Dispatch to specific method
function find_zero(
    f::F,
    method::M,
    soltype::SolutionType = CompactSolution(),
    tol::Union{Nothing, AbstractTolerance} = nothing,
    maxiters::Int = 1_000,
) where {FT <: FTypes, F <: Function, M <: RootSolvingMethod{FT}}
    if tol === nothing
        tol = default_tol(FT)
    end
    return find_zero(f, M, method_args(method)..., soltype, tol, maxiters)
end

function Broadcast.broadcasted(
    ::typeof(find_zero),
    f::F,
    method::M,
    soltype::SolutionType = CompactSolution(),
    tol::Union{Nothing, AbstractTolerance} = nothing,
    maxiters::Int = 1_000,
) where {FT <: FTypes, F, M <: RootSolvingMethod{FT}}
    if tol === nothing
        tol = default_tol(FT)
    end
    return Broadcast.broadcasted(
        find_zero,
        f,
        M,
        method_args(method)...,
        soltype,
        tol,
        maxiters,
    )
end

####
#### Numerical methods
####

@inline method_args(method::SecantMethod) = (method.x0, method.x1)
function find_zero(
    f::F,
    ::Type{M},
    x0::FT,
    x1::FT,
    soltype::SolutionType,
    tol::AbstractTolerance,
    maxiters::Int,
) where {F <: Function, FT <: FTypes, M <: SecantMethod}
    return _find_zero_secant(f, x0, x1, soltype, tol, maxiters)
end

@inline method_args(method::BisectionMethod) = (method.x0, method.x1)
function find_zero(
    f::F,
    ::Type{M},
    x0::FT,
    x1::FT,
    soltype::SolutionType,
    tol::AbstractTolerance,
    maxiters::Int,
) where {F <: Function, FT, M <: BisectionMethod}
    return _find_zero_bracketed(
        f,
        _bisection_rule,
        _bisection_y_update,
        x0,
        x1,
        soltype,
        tol,
        maxiters,
    )
end

@inline method_args(method::RegulaFalsiMethod) = (method.x0, method.x1)
function find_zero(
    f::F,
    ::Type{M},
    x0::FT,
    x1::FT,
    soltype::SolutionType,
    tol::AbstractTolerance,
    maxiters::Int,
) where {F <: Function, FT, M <: RegulaFalsiMethod}
    return _find_zero_bracketed(
        f,
        _regula_falsi_rule,
        _regula_falsi_y_update,
        x0,
        x1,
        soltype,
        tol,
        maxiters,
    )
end

@inline method_args(method::BrentsMethod) = (method.x0, method.x1)
function find_zero(
    f::F,
    ::Type{M},
    x0::FT,
    x1::FT,
    soltype::SolutionType,
    tol::AbstractTolerance,
    maxiters::Int,
) where {F <: Function, FT, M <: BrentsMethod}
    return _find_zero_brent(f, x0, x1, soltype, tol, maxiters)
end

@inline method_args(method::NewtonsMethodAD) = (method.x0,)
function find_zero(
    f::F,
    ::Type{M},
    x0::FT,
    soltype::SolutionType,
    tol::AbstractTolerance,
    maxiters::Int,
) where {F <: Function, FT, M <: NewtonsMethodAD}
    return _find_zero_newton(
        Base.Fix1(value_deriv, f),
        f,
        x0,
        soltype,
        tol,
        maxiters,
    )
end

@inline method_args(method::NewtonsMethod) = (method.x0,)
function find_zero(
    f::F,
    ::Type{M},
    x0::FT,
    soltype::SolutionType,
    tol::AbstractTolerance,
    maxiters::Int,
) where {F <: Function, FT, M <: NewtonsMethod}
    return _find_zero_newton(f, x -> f(x)[1], x0, soltype, tol, maxiters)
end

####
#### Method-specific helper functions
####

# --- Rules for calculating the next `x` value ---

# Update rule for Bisection method: `x` is the midpoint.
@inline _bisection_rule(x0, y0, x1, y1) = x0 + (x1 - x0) / 2

# Update rule for Regula Falsi method: `x` is the secant intercept.
@inline _regula_falsi_rule(x0, y0, x1, y1) = ifelse(
    abs(y1 - y0) < 100 * eps(y1),
    x0 + (x1 - x0) / 2, # Fallback to bisection if function is flat
    (x0 * y1 - x1 * y0) / (y1 - y0),
)

# --- Rules for updating the endpoint `y` values ---

# For Bisection, `y` values are simply replaced.
@inline _bisection_y_update(y0, y1, y, is_neg, lastside) =
    (ifelse(is_neg, y0, y), ifelse(is_neg, y, y1))

# For Regula Falsi, apply the Illinois method stagnation fix.
@inline _regula_falsi_y_update(y0, y1, y, is_neg, lastside) = (
    ifelse(is_neg, ifelse(lastside == +1, y0 / 2, y0), y),
    ifelse(is_neg, y, ifelse(lastside == -1, y1 / 2, y1)),
)


# Internal function implementing the core bracketing algorithm used by bisection, 
# Regula Falsi, and other bracketing methods.
@inline function _find_zero_bracketed(
    f,
    update_x_rule,
    update_y_rule,
    x0,
    x1,
    soltype,
    tol,
    maxiters,
)
    FT = typeof(x0)
    if !isfinite(x0) || !isfinite(x1)
        y = FT(Inf)
        x_history = init_history(soltype, FT)
        y_history = init_history(soltype, FT)
        return SolutionResults(soltype, x0, false, y, 0, x_history, y_history)
    end

    y0 = f(x0)
    y1 = f(x1)
    if y0 * y1 >= 0
        x_history = init_history(soltype, x0)
        y_history = init_history(soltype, y0)
        return SolutionResults(soltype, x0, false, y0, 0, x_history, y_history)
    end

    x_history = init_history(soltype, x0)
    y_history = init_history(soltype, y0)
    lastside = 0

    local x, y
    for i in 1:maxiters
        x = update_x_rule(x0, y0, x1, y1)
        y = f(x)

        push_history!(x_history, x, soltype)
        push_history!(y_history, y, soltype)

        # Check for convergence first (including exact root)
        if tol(x0, x1, y)
            return SolutionResults(soltype, x, true, y, i, x_history, y_history)
        end

        # Update x and y values using the provided policy rules
        is_neg = y * y0 < 0
        x0_new = ifelse(is_neg, x0, x)
        x1_new = ifelse(is_neg, x, x1)
        y0_new, y1_new = update_y_rule(y0, y1, y, is_neg, lastside)
        x0, x1, y0, y1 = x0_new, x1_new, y0_new, y1_new

        lastside = ifelse(is_neg, +1, -1)
    end

    return SolutionResults(soltype, x, false, y, maxiters, x_history, y_history)
end


# Internal function implementing Brent's method, which combines bisection, secant, and 
# inverse quadratic interpolation.
@inline function _find_zero_brent(f, x0, x1, soltype, tol, maxiters)
    FT = typeof(x0)
    if !isfinite(x0) || !isfinite(x1)
        y = FT(Inf)
        x_history = init_history(soltype, FT)
        y_history = init_history(soltype, FT)
        return SolutionResults(soltype, x0, false, y, 0, x_history, y_history)
    end

    a, b = x0, x1
    fa, fb = f(a), f(b)

    if fa * fb >= 0
        # Return failed solution instead of error for GPU compatibility
        x_history = init_history(soltype, a)
        y_history = init_history(soltype, fa)
        return SolutionResults(soltype, a, false, fa, 0, x_history, y_history)
    end

    if abs(fa) < abs(fb)
        a, b = b, a
        fa, fb = fb, fa
    end

    x_history = init_history(soltype, a)
    y_history = init_history(soltype, fa)
    push_history!(x_history, b, soltype)
    push_history!(y_history, fb, soltype)

    c = a
    fc = fa
    d = b - a
    e = d

    for i in 1:maxiters
        # Convergence check
        if (tol(a, b, fb)) || (fb == 0)
            return SolutionResults(
                soltype,
                b,
                true,
                fb,
                i,
                x_history,
                y_history,
            )
        end

        # On GPUs, `ifelse` is often more performant than `if-else` blocks
        # by avoiding branch divergence.
        s = ifelse(
            (fa != fc) & (fb != fc),
            # Inverse Quadratic Interpolation
            a * fb * fc / ((fa - fb) * (fa - fc)) +
            b * fa * fc / ((fb - fa) * (fb - fc)) +
            c * fa * fb / ((fc - fa) * (fc - fb)),
            # Secant Method
            b - fb * (b - a) / (fb - fa),
        )

        # Check if interpolation result is acceptable
        m = (3 * a + b) / 4
        s_is_between = ifelse(a < b, (s > m) & (s < b), (s < m) & (s > b))

        # If step is not acceptable, fall back to bisection
        use_bisection =
            (!s_is_between) || (abs(s - b) >= abs(e / 2)) || (abs(e) < eps(FT))

        s_new = ifelse(use_bisection, (a + b) / 2, s)
        d_new = s_new - b
        e_new = ifelse(use_bisection, d_new, d)
        s, d, e = s_new, d_new, e_new

        fs = f(s)
        push_history!(x_history, s, soltype)
        push_history!(y_history, fs, soltype)

        c, fc = b, fb

        # Update bracket
        cond_bracket = fa * fs < 0
        a_new = ifelse(cond_bracket, a, s)
        fa_new = ifelse(cond_bracket, fa, fs)
        b_new = ifelse(cond_bracket, s, b)
        fb_new = ifelse(cond_bracket, fs, fb)
        a, fa, b, fb = a_new, fa_new, b_new, fb_new

        # Ensure abs(fb) <= abs(fa)
        cond_swap = abs(fa) < abs(fb)
        a_new_swap = ifelse(cond_swap, b, a)
        b_new_swap = ifelse(cond_swap, a, b)
        fa_new_swap = ifelse(cond_swap, fb, fa)
        fb_new_swap = ifelse(cond_swap, fa, fb)
        a, fa, b, fb = a_new_swap, fa_new_swap, b_new_swap, fb_new_swap
    end

    return SolutionResults(
        soltype,
        b,
        false,
        fb,
        maxiters,
        x_history,
        y_history,
    )
end

# Internal function implementing the Secant method, an open, two-point method.
@inline function _find_zero_secant(f, x0, x1, soltype, tol, maxiters)
    FT = typeof(x0)
    if !isfinite(x0) || !isfinite(x1)
        y = FT(Inf)
        x_history = init_history(soltype, FT)
        y_history = init_history(soltype, FT)
        return SolutionResults(soltype, x0, false, y, 0, x_history, y_history)
    end

    y0 = f(x0)
    y1 = f(x1)

    x_history = init_history(soltype, x0)
    y_history = init_history(soltype, y0)

    for i in 1:maxiters
        Δx = x1 - x0
        Δy = y1 - y0

        if abs(Δy) <= 100 * eps(y1)
            # Exiting because the function is flat. This is a stall.
            # The method has only "converged" if the residual is already small.
            converged = tol(x0, x1, y1) # Check for convergence
            return SolutionResults(
                soltype,
                x1,
                converged,
                y1,
                i,
                x_history,
                y_history,
            )
        end

        # Update x0, y0 to the "previous" state for the next iteration
        # and for the tolerance check below.
        x0, y0 = x1, y1

        # Update x1 in-place using the now-stored previous state in y0.
        x1 -= y0 * Δx / Δy
        y1 = f(x1)

        push_history!(x_history, x1, soltype)
        push_history!(y_history, y1, soltype)

        # Check for convergence
        if tol(x0, x1, y1)
            return SolutionResults(
                soltype,
                x1,
                true,
                y1,
                i,
                x_history,
                y_history,
            )
        end
    end

    return SolutionResults(
        soltype,
        x1,
        false,
        y1,
        maxiters,
        x_history,
        y_history,
    )
end


# Internal function implementing the main iteration loop for Newton methods (with AD 
# or analytical derivatives)
@inline function _find_zero_newton(
    f_value_and_deriv,
    f_value_only,
    x0,
    soltype,
    tol,
    maxiters,
)
    FT = typeof(x0)
    x_history = init_history(soltype, FT)
    y_history = init_history(soltype, FT)

    if !isfinite(x0)
        y = FT(Inf)
        return SolutionResults(soltype, x0, false, y, 0, x_history, y_history)
    end

    # Initial state
    x = x0
    y, y′ = f_value_and_deriv(x)

    # Log the initial state
    push_history!(x_history, x, soltype)
    push_history!(y_history, y, soltype)

    c = FT(1e-4) # Conservative Armijo constant
    for i in 1:maxiters

        # Fallback to secant method when derivative is too small
        if abs(y′) <= FT(1e-5)
            x_pert = x + (iszero(x) ? sqrt(eps(FT)) : x * sqrt(eps(FT)))
            # Note: History is lost on fallback for VerboseSolution 
            return _find_zero_secant(
                f_value_only,
                x,
                x_pert,
                soltype,
                tol,
                maxiters - i,
            )
        end

        Δx = y / y′ # Full Newton step

        # --- Step Limiting for Robustness ---
        max_step = FT(20) * (abs(x) + sqrt(eps(FT)))
        Δx = sign(Δx) * min(abs(Δx), max_step)

        # --- Backtracking Line Search ---
        α = one(FT)
        x_new = x - α * Δx
        y_new = f_value_only(x_new)
        max_backtrack = 5
        j = 0
        while j < max_backtrack && (
            isfinite(y_new) &
            (abs(y_new) > abs(y) * (FT(1) - c * α)) &
            (α > eps(FT))
        )
            α /= 2
            x_new = x - α * Δx
            y_new = f_value_only(x_new)
            j += 1
        end

        # If the line search failed to find a finite point,
        # the original point `x` is the best we can do. Abort.
        if !isfinite(y_new)
            return SolutionResults(
                soltype,
                x,
                false,
                y,
                i,
                x_history,
                y_history,
            )
        end

        # Log the accepted new point
        push_history!(x_history, x_new, soltype)
        push_history!(y_history, y_new, soltype)

        # Check for convergence 
        if tol(x, x_new, y_new)
            return SolutionResults(
                soltype,
                x_new,
                true,
                y_new,
                i,
                x_history,
                y_history,
            )
        end

        # Update for next iteration
        x = x_new
        y, y′ = f_value_and_deriv(x)
    end

    y_final = f_value_only(x)
    return SolutionResults(
        soltype,
        x,
        false,
        y_final,
        maxiters,
        x_history,
        y_history,
    )
end

"""
    method_args(method::RootSolvingMethod)

Extract the initial guess(es) for a root-solving method for internal dispatch.

This function is used internally to unpack method parameters for passing to the
appropriate `find_zero` implementation.

## Arguments
- `method::RootSolvingMethod`: The root-solving method instance

## Returns
- `Tuple`: Initial guess(es) specific to the method type

## Example
```julia
method = SecantMethod{Float64}(0.0, 1.0)
args = method_args(method)  # Returns (0.0, 1.0)
```
"""
function method_args end

"""
    value_deriv(f, x)

Compute both the function value and its derivative at point `x` using automatic differentiation.

This function uses ForwardDiff.jl to simultaneously compute `f(x)` and `f'(x)`, which is
more efficient than computing them separately when both are needed (as in Newton's method). It
is used internally by `NewtonsMethodAD`.

## Arguments
- `f`: Function to evaluate
- `x::FT`: Point at which to evaluate the function and derivative

## Returns
- `Tuple{FT, FT}`: `(f(x), f'(x))` where the second element is the derivative

## Example
```julia
f(x) = x^3 - 2x + 1
val, deriv = value_deriv(f, 1.5)
# val ≈ 1.375, deriv ≈ 4.75
```
"""
function value_deriv(f, x::FT) where {FT}
    tag = typeof(ForwardDiff.Tag(f, FT))
    y = f(ForwardDiff.Dual{tag}(x, one(x)))
    ForwardDiff.value(tag, y), ForwardDiff.partials(tag, y, 1)
end

end
