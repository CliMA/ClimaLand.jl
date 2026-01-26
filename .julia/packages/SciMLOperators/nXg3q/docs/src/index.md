# SciMLOperators.jl: Unified operator interface for Julia and SciML

`SciMLOperators` is a package for managing linear, nonlinear,
time-dependent, and parameter dependent operators acting on vectors,
(or column-vectors of matrices). We provide wrappers for matrix-free
operators, fast tensor-product evaluations, pre-cached mutating
evaluations, as well as `Zygote`-compatible non-mutating evaluations.

The lazily implemented operator algebra allows the user to update the
operator state by passing in an update function that accepts arbitrary
parameter objects. Further, our operators behave like `AbstractMatrix` types
thanks to overloads defined for methods in `Base`, and `LinearAlgebra`.

Therefore, an `AbstractSciMLOperator` can be passed to `LinearSolve.jl`,
or `NonlinearSolve.jl` as a linear or nonlinear operator, or to
`OrdinaryDiffEq.jl` as an `ODEFunction`. Examples of usage within the
`SciML` ecosystem are provided in the documentation.

## Installation

To install SciMLOperators.jl, use the Julia package manager:

```julia
using Pkg
Pkg.add("SciMLOperators")
```

## Why `SciMLOperators`?

Many functions, from linear solvers to differential equations, require
the use of matrix-free operators to achieve maximum performance in
many scenarios. `SciMLOperators.jl` defines the abstract interface for how
operators in the SciML ecosystem are supposed to be defined. It gives the
common set of functions and traits that solvers can rely on for properly
performing their tasks. Along with that, `SciMLOperators.jl` provides
definitions for the basic standard operators that are used as building
blocks for most tasks, simplifying the use of operators while also
demonstrating to users how such operators can be built and used in practice.

`SciMLOperators.jl` has the design that is required to be used in
all scenarios of equation solvers. For example, Magnus integrators for
differential equations require defining an operator ``u' = A(t) u``, while
Munthe-Kaas methods require defining operators of the form ``u' = A(u) u``.
Thus, the operators need some form of time and state dependence, which the
solvers can update and query when they are non-constant
(`update_coefficients!`). Additionally, the operators need the ability to
act like “normal” functions for equation solvers. For example, if `A(v,u,p,t)`
has the same operation as `update_coefficients(A, u, p, t); A * v`, then `A`
can be used in any place where a differential equation definition
`(u,p,t) -> A(u, u, p, t)` is used without requiring the user or solver to do any extra
work.

Another example is state-dependent mass matrices. `M(u,p,t)*u' = f(u,p,t)`.
When solving such an equation, the solver must understand how to "update M"
during operations, and thus the ability to update the state of `M` is a required
function in the interface. This is also required for the definition of Jacobians
`J(u,p,t)` in order to be properly used with Krylov methods inside of ODE solves
without reconstructing the matrix-free operator at each step.

Thus while previous good efforts for matrix-free operators have existed
in the Julia ecosystem, such as
[LinearMaps.jl](https://github.com/JuliaLinearAlgebra/LinearMaps.jl), those
operator interfaces lack these aspects to actually be fully seamless
with downstream equation solvers. This necessitates the definition and use of
an extended operator interface with all of these properties, hence the
`AbstractSciMLOperator` interface.

!!! warn
    
    This means that LinearMaps.jl is fundamentally lacking and is incompatible
    with many of the tools in the SciML ecosystem, except for the specific cases
    where the matrix-free operator is a constant!

## Features

  - Matrix-free operators with `FunctionOperator`
  - Fast tensor product evaluation with `TensorProductOperator`
  - Lazy algebra: addition, subtraction, multiplication, inverse, adjoint, and transpose
  - Couple fast methods for operator evaluation with inversion via `InvertibleOperator`
  - One-line API to update operator state depending on arbitrary parameters.
  - Mutating and nonmutating update behavior (Zygote compatible)
  - One-line API for pre-caching operators for in-place operator evaluations

## Contributing

  - Please refer to the
    [SciML ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://github.com/SciML/ColPrac/blob/master/README.md)
    for guidance on PRs, issues, and other matters relating to contributing to SciML.

  - See the [SciML Style Guide](https://github.com/SciML/SciMLStyle) for common coding practices and other style decisions.
  - There are a few community forums:
    
      + The #diffeq-bridged and #sciml-bridged channels in the
        [Julia Slack](https://julialang.org/slack/)
      + The #diffeq-bridged and #sciml-bridged channels in the
        [Julia Zulip](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged)
      + On the [Julia Discourse forums](https://discourse.julialang.org)
      + See also [SciML Community page](https://sciml.ai/community/)

## Reproducibility

```@raw html
<details><summary>The documentation of this SciML package was built using these direct dependencies,</summary>
```

```@example
using Pkg # hide
Pkg.status() # hide
```

```@raw html
</details>
```

```@raw html
<details><summary>and using this machine and Julia version.</summary>
```

```@example
using InteractiveUtils # hide
versioninfo() # hide
```

```@raw html
</details>
```

```@raw html
<details><summary>A more complete overview of all dependencies and their versions is also provided.</summary>
```

```@example
using Pkg # hide
Pkg.status(; mode = PKGMODE_MANIFEST) # hide
```

```@raw html
</details>
```

```@eval
using TOML
using Markdown
version = TOML.parse(read("../../Project.toml", String))["version"]
name = TOML.parse(read("../../Project.toml", String))["name"]
link_manifest = "https://github.com/SciML/" * name * ".jl/tree/gh-pages/v" * version *
                "/assets/Manifest.toml"
link_project = "https://github.com/SciML/" * name * ".jl/tree/gh-pages/v" * version *
               "/assets/Project.toml"
Markdown.parse("""You can also download the
[manifest]($link_manifest)
file and the
[project]($link_project)
file.
""")
```
