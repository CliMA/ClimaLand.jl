# Simple Demonstration of a Symbolic System Structure

In this tutorial we will show how to implement a system structure type for defining the
symbolic indexing of a domain-specific language. This tutorial will show how the
`SymbolCache` type is defined to take in arrays of symbols for its independent, dependent,
and parameter variable names and uses that to define the symbolic indexing interface.

## Defining the ODE

For this example, we will use the Robertson equations:

```math
\begin{aligned}
\frac{dy_1}{dt} &= -0.04y₁ + 10^4 y_2 y_3 \\
\frac{dy_2}{dt} &= 0.04 y_1 - 10^4 y_2 y_3 - 3*10^7 y_{2}^2 \\
\frac{dy_3}{dt} &= 3*10^7 y_{2}^2 \\
\end{aligned}
```

The in-place function for this ODE system can be defined as:

```@example symbolcache
function rober!(du, u, p, t)
    y₁, y₂, y₃ = u
    k₁, k₂, k₃ = p
    du[1] = -k₁ * y₁ + k₃ * y₂ * y₃
    du[2] = k₁ * y₁ - k₂ * y₂^2 - k₃ * y₂ * y₃
    du[3] = k₂ * y₂^2
    nothing
end
```

To add symbolic names for the states in this example, a [`SymbolCache`](@ref) can be
created and passed as the `sys` keyword argument to the `ODEFunction` constructor,
as shown below:

```@example symbolcache
using OrdinaryDiffEq, SymbolicIndexingInterface

sys = SymbolCache([:y₁, :y₂, :y₃])
odefun = ODEFunction(rober!; sys = sys)
nothing # hide
```

This is then used to create and solve the `ODEProblem`

```@example symbolcache
prob = ODEProblem(odefun, [1.0, 0.0, 0.0], (0.0, 1e5), [0.04, 3e7, 1e4])
sol = solve(prob, Rosenbrock23())
```

The solution can now be indexed symbolically:

```@example symbolcache
sol[:y₁]
```

```@example symbolcache
sol(1e3, idxs=:y₁)
```

However, we did not give names to the parameters or the independent variables. They can
be specified using `SymbolCache` as well:

```@example symbolcache
sys = SymbolCache([:y₁, :y₂, :y₃], [:k₁, :k₂, :k₃], :t)
odefun = ODEFunction(rober!; sys = sys)
prob = ODEProblem(odefun, [1.0, 0.0, 0.0], (0.0, 1e5), [0.04, 3e7, 1e4])
sol = solve(prob, Rosenbrock23())
getk1 = getp(sys, :k₁)

getk1(prob)
```

```@example symbolcache
getk1(sol)
```

```@example symbolcache
sol[:t]
```
