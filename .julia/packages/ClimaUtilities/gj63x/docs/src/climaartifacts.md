# `ClimaArtifacts`

[Julia artifacts](https://pkgdocs.julialang.org/v1/artifacts/) are pieces
of data that can be distributed alongside a package. Julia artifacts were
developed to distribute application binaries (e.g., compiled libraries). In
[`CliMA`](https://github.com/CliMA), we use them to distribute data required to
perform our simulations (e.g., input data).

The `ClimaArtifacts` module extends the Julia artifact system to solve two
issues:
1. Ensuring that artifacts can be gracefully acquired by parallel runs;
2. _Tagging_ artifacts that are accessed during a simulation.

We will examine these two problems below. In the meantime, this entire
documentation page can be summarized in a one short directive for package
developers:

>  Instead of accessing artifacts with
>  [ArtifactWrappers.jl](https://github.com/CliMA/ArtifactWrappers.jl) or using
>  Julia artifacts directly, use the `@clima_artifact` macro instead.

Also, keep in mind that

!!! note

    Julia artifacts are always entire folders, never single files!

### Julia artifacts and MPI

Package developers can specify one of two modes for any given artifact: greedy
(default) or lazy download. Artifact that are not marked as `lazy` are
automatically downloaded by Julia when the package is instantiated. On the other
hand, `lazy` artifacts are downloaded the first time they are accessed.

CliMA packages can distribute tens of artifacts that are relevant for very
different types of simulations, and it is a good idea to mark artifacts as
`lazy` unless they are strictly required for the operation of the package (e.g.,
the orbital parameters in
[Insolation.jl](https://github.com/CliMA/Insolation.jl)).

Lazy artifacts are incompatible with MPI. In parallel runs, each process tries
to download the same file, resulting in a race condition and corrupted files
(not to mention tens of processing pinging the same server at the same time).
`ClimaArtifacts` implements a new macro, `@clima_artifact`, to solve this
problem.

`@clima_artifact` is a near drop-in replacement for the
[`@artifact_str`](https://docs.julialang.org/en/v1/stdlib/Artifacts/#Artifacts.@artifact_str)
Julia macro.

For greedy artifacts and non-MPI runs, it is possible to simply call
`@clima_artifact(artifact_name)`. This will return the path of the artifact
folder. This macro will fail for lazy artifacts. In that case, one has to also
pass the [ClimaComms.jl](https://github.com/CliMA/ClimaComms.jl) `context`. The
context is required because `@clima_artifact` needs to synchronize different MPI
processes.

#### Example

> This extension is loaded when loading `ClimaComms`

Assume `socrates` is a lazy artifact, we can access the `socrates` artifact folder as in
```julia
using ClimaUtilities.ClimaArtifacts
# If the artifact is lazy, we also need LazyArtifacts
using LazyArtifacts

import ClimaComms
# When loading ClimaComms, a Julia extension for ClimaUtilities will be loaded

my_mpi_context = ClimaComms.context()

socrates_path_folder = @clima_artifact("socrates", context)
```

The `@clima_artifact` macro is executed at parse time when the argument is a
literal string (e.g., `@clima_artifact("socrates")`), and at runtime when it is
a variable `@clima_artifact(artifact_name)`.

!!! note

    `context` is a positional argument, not a keyword one. Calling
    `@clima_artifact("socrates"; context)` will fail. (This is due to how Julia
    macros handle keyword arguments)

### Tagging artifacts

A full climate simulation requires lots of external input data. Most of this
data comes from scientific experiments that have to be properly acknowledged.
`ClimaArtifacts` allows users to know what artifacts were used in a given
simulation. As long as artifacts are being accessed with `@clima_artifacts`, the
`ClimaArtifacts` keeps track of what is being used. The set of artifacts
accessed can be obtained with `ClimaArtifacts.accessed_artifacts`.

#### Example
```julia
using ClimaUtilities.ClimaArtifacts
art1 = @clima_artifact("socrates")
art2 = @clima_artifact("zeno")

ClimaArtifacts.accessed_artifacts()
# Set(["socrates", "zeno"])
```

## API

```@docs
ClimaUtilities.ClimaArtifacts.@clima_artifact
ClimaUtilities.ClimaArtifacts.accessed_artifacts
```
