# CommonSolve

[![Join the chat at https://julialang.zulipchat.com #sciml-bridged](https://img.shields.io/static/v1?label=Zulip&message=chat&color=9558b2&labelColor=389826)](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged)
[![Global Docs](https://img.shields.io/badge/docs-SciML-blue.svg)](https://docs.sciml.ai/CommonSolve/stable)

[![codecov](https://codecov.io/gh/SciML/CommonSolve.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/SciML/CommonSolve.jl)
[![Build Status](https://github.com/SciML/CommonSolve.jl/workflows/CI/badge.svg)](https://github.com/SciML/CommonSolve.jl/actions?query=workflow%3ACI)

[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)

This holds the common `solve`, `init`, `step!`, and `solve!` commands. By using the same definition,
solver libraries from other completely different ecosystems can extend the functions and thus
not clash with SciML if both ecosystems export the `solve` command. The rules are that 
you must dispatch on one of your own types. That's it. No pirates.

## General recommendation

`solve` function has the default definition

```julia
solve(args...; kwargs...) = solve!(init(args...; kwargs...))
```

So, we recommend defining

```julia
init(::ProblemType, args...; kwargs...) :: SolverType
solve!(::SolverType) :: SolutionType
```

where `ProblemType`, `SolverType`, and `SolutionType` are the types defined in
your package. 

In many cases, the `SolverType` is an object that is iteratively progressed to achieve the solution. 
In such cases, the `step!` function can be used:

```julia
step!(::SolverType, args...; kwargs...)
```

To avoid method ambiguity, the first argument of `solve`, `solve!`, `step!`, and `init`
_must_ be dispatched on the type defined in your package.  For example, do
_not_ define a method such as

```julia
init(::AbstractVector, ::AlgorithmType)
```
