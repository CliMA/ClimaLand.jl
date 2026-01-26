# SymbolicIndexingInterface.jl: Standardized Symbolic Indexing of Julia

SymbolicIndexingInterface.jl is a set of interface functions for handling containers
of symbolic variables.

## Installation

To install SymbolicIndexingInterface.jl, use the Julia package manager:

```julia
using Pkg
Pkg.add("SymbolicIndexingInterface")
```

## Introduction

The symbolic indexing interface has 2 levels:

1. The user level. At the user level, a modeler or engineer simply uses terms from a
   domain-specific language (DSL) inside of SciML functionality and will receive the requested
   values. For example, if a DSL defines a symbol `x`, then `sol[x]` returns the solution
   value(s) for `x`.
2. The DSL system structure level. This is the structure which defines the symbolic indexing
   for a given problem/solution. DSLs can tag a constructed problem/solution with this
   object in order to endow the SciML tools with the ability to index symbolically according
   to the definitions the DSL writer wants.

## Contributing

- Please refer to the
  [SciML ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://github.com/SciML/ColPrac/blob/master/README.md)
  for guidance on PRs, issues, and other matters relating to contributing to SciML.
- There are a few community forums:
    - the #diffeq-bridged channel in the [Julia Slack](https://julialang.org/slack/)
    - [JuliaDiffEq](https://gitter.im/JuliaDiffEq/Lobby) on Gitter
    - on the [Julia Discourse forums](https://discourse.julialang.org)
    - see also [SciML Community page](https://sciml.ai/community/)

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
