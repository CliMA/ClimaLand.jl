# CommonSolve.jl: The Common Solve Definition and Interface

This holds the common `solve`, `init`, `solve!`, and `step!` commands. By using the same definition,
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

In many cases, the `SolverType` is an object that is iteratively progressed to achieve the solution. In such cases, the `step!` function can be used:


```julia
step!(::SolverType, args...; kwargs...)
```

To avoid method ambiguity, the first argument of `solve`, `solve!`, `step!`, and `init`
_must_ be dispatched on the type defined in your package.  For example, do
_not_ define a method such as

```julia
init(::AbstractVector, ::AlgorithmType)
```

## API

```@docs
CommonSolve.init
CommonSolve.solve
CommonSolve.solve!
CommonSolve.step!
```

## Contributing

- Please refer to the
  [SciML ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://github.com/SciML/ColPrac/blob/master/README.md)
  for guidance on PRs, issues, and other matters relating to contributing to SciML.
- There are a few community forums:
    - The #diffeq-bridged and #sciml-bridged channels in the
      [Julia Slack](https://julialang.org/slack/)
    - [JuliaDiffEq](https://gitter.im/JuliaDiffEq/Lobby) on Gitter
    - On the Julia Discourse forums (look for the [modelingtoolkit tag](https://discourse.julialang.org/tag/modelingtoolkit)
    - See also [SciML Community page](https://sciml.ai/community/)

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
Pkg.status(;mode = PKGMODE_MANIFEST) # hide
```
```@raw html
</details>
```
```@raw html
You can also download the
<a href="
```
```@eval
using TOML
version = TOML.parse(read("../../Project.toml",String))["version"]
name = TOML.parse(read("../../Project.toml",String))["name"]
link = "https://github.com/SciML/"*name*".jl/tree/gh-pages/v"*version*"/assets/Manifest.toml"
```
```@raw html
">manifest</a> file and the
<a href="
```
```@eval
using TOML
version = TOML.parse(read("../../Project.toml",String))["version"]
name = TOML.parse(read("../../Project.toml",String))["name"]
link = "https://github.com/SciML/"*name*".jl/tree/gh-pages/v"*version*"/assets/Project.toml"
```
```@raw html
">project</a> file.
```
