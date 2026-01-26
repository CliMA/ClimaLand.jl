<h1 align="center">
  <img src="logo.svg" width="180px"> <br>
ClimaUtilities.jl

[![CI](https://github.com/CliMA/ClimaUtilities.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/CliMA/ClimaUtilities.jl/actions/workflows/ci.yml)
[![Downstream](https://github.com/CliMA/ClimaUtilities.jl/actions/workflows/downstream.yml/badge.svg)](https://github.com/CliMA/ClimaUtilities.jl/actions/workflows/downstream.yml)
[![codecov](https://codecov.io/gh/CliMA/ClimaUtilities.jl/graph/badge.svg?token=sbYkr5ydzp)](https://codecov.io/gh/CliMA/ClimaUtilities.jl)
[![Docs](https://img.shields.io/badge/docs_are_here-click_me!-blue.svg)](https://clima.github.io/ClimaUtilities.jl/dev/)
[![Downloads](https://img.shields.io/badge/dynamic/json?url=http%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Ftotal_downloads%2FClimaUtilities&query=total_requests&suffix=%2Ftotal&label=Downloads)](http://juliapkgstats.com/pkg/ClimaUtilities)

</h1>

`ClimaUtilities` is collection of general-purpose tools and types for CliMA packages.

`ClimaUtilities.jl` contains:
- [`ClimaArtifacts`](https://clima.github.io/ClimaUtilities.jl/dev/climaartifacts/),
  a module that provides an MPI-safe (with `ClimaComms`) way to lazily download
  artifacts and to tag artifacts that are being accessed in a given simulation.
- [`TimeManager`](https://clima.github.io/ClimaUtilities.jl/dev/timemanager/) to
  handle time and dates. `TimeManager` implements `ITime`, the time type used in
  CliMA simulations.
- [`SpaceVaryingInputs` and
  `TimeVaryingInputs`](https://clima.github.io/ClimaUtilities.jl/dev/inputs/) to
  work with external input data.
- [`DataStructures`](https://clima.github.io/ClimaUtilities.jl/dev/datastrctures/)
  with useful general purpose data structures.
- [`FileReaders`](https://clima.github.io/ClimaUtilities.jl/dev/filereaders/),
  [`DataHandling`](https://clima.github.io/ClimaUtilities.jl/dev/datahandling/),
  and [`Regridders`](https://clima.github.io/ClimaUtilities.jl/dev/regridders/)
  to process input data and remap it onto the simulation grid.
- [`OnlineLogging`](https://clima.github.io/ClimaUtilities.jl/dev/onlinelogging/)
  to add output information while a simulation is running.
- [`OutputPathGenerator`](https://clima.github.io/ClimaUtilities.jl/dev/outputpathgenerator/)
  to prepare the output directory structure of a simulation.

## ClimaUtilities.jl Developer Guidelines

These guidelines aim to ensure consistent code quality, maintainability, and a
smooth collaborative workflow for `ClimaUtilities.jl`. Please, read these
guidelines even if you are familiar with other CliMA packages as there may be
some differences.

### A key design principle

`ClimaUtilities.jl` is meant to be a foundational package, so it should have
near-zero-cost for features that are not explicitly requested. This means that
it should not depend on any package outside of the standard library and
`ClimaComms` and should have negligible import costs. To accomplish this,
everything is implemented in Julia extensions.

Extensions are organized in the following way. The extensions defined in the
`Project.toml` are defined in terms of the packages they require to be loaded.
This avoids circular dependencies among extensions. Each of these extensions
`include`s modules that match what is defined in `src`.

For example, `ClimaUtilitiesClimaCoreNCDatasetsExt` is loaded when `ClimaCore`
and `NCDatasets` are loaded. Internally, `ClimaUtilitiesClimaCoreNCDatasetsExt`
loads `DataHandlingExt.jl`, `SpaceVaryingInputsExt.jl`, and
`TimeVaryingInputsExt.jl`. Each of them implements methods defined in
`src/DataHandling.jl`, `src/SpaceVaryingInputs.jl`, and
`src/TimeVaryingInputs.jl` respectively

### Tests and environments

We prioritize well-tested code to guarantee `ClimaUtilities.jl` functions
reliably. Here are some principles we follow:

#### Tests are collected in the `test` folder and are exclusively there

This means that:
- All the tests can be run with `Pkg.test()` or `julia --project=test tests/runtests.jl`;
- There are no special tests in the Buildkite pipeline;
- In fact, we use Buildkite only to run the tests with MPI/GPUs.

#### There are no checked-in `Manifest.toml` files

While checking in `Manifest.toml` files ensures reproducibility, it also
introduces some nuisance, including:
- lot of git/repository noise just for "up deps";
- multiple environments that have to be managed;
- busywork to keep the manifests updated.

In this repository, we have four environments:
- project,
- documentation,
- test,
- buildkite.

The buildkite environment, in `.buildkite`, contains everything that is required
to test `ClimaUtilities.jl` in its most general configuration (ie, including
`CUDA.jl` and `MPI.jl`). As the name suggests, this is the (only) environment
used by our Buildkite pipeline.

> :note: Please, open an issue if you find workflow problems/friction with this
> system.

#### Running tests

There are two equivalent ways to run tests.

First, Start a Julia session in the `ClimaUtilities` directory:
``` sh
julia --project
```
Enter `Pkg` mode by typing `]`. This will change the prompt. Run `test`.

Equivalently, you can run
``` sh
julia --project=test tests/runtests.jl
```

> :note: Please, open an issue if you find workflow problems/friction with this
> system.

#### Code Formatting with `JuliaFormatter.jl`

One of the tests consists in checking that the code is uniformly formatted. We
use [JuliaFormatter.jl](https://github.com/domluna/JuliaFormatter.jl) to achieve
consistent formatting. Here's how to use it:

You can either install in your base environment with
``` sh
julia -e 'using Pkg; Pkg.add("JuliaFormatter")'
```
or use it from within the `TestEnv` or the `.buildkite` environments (see previous section).

Then, you can format the package running:
``` julia
using JuliaFormatter; format(".")
```
or just with `format(".")` if the package is already imported.

The rules for formatting are defined in the `.JuliaFormatter.toml`.

If you are used to formatting from the command line instead of the REPL, you can
install `JuliaFormatter` in your base environment and call
``` sh
julia -e 'using JuliaFormatter; format(".")'
```
You could also define a shell alias
``` sh
alias julia_format_here="julia -e 'using JuliaFormatter; format(\".\")'"
```

> :note: Please, open an issue if you find workflow problems/friction with this
> system.

### Documentation

Documentation is generated with
[Documenter.jl](https://documenter.juliadocs.org/stable/). We strive to have
complete and up-to-date information.

To generate documentation, run
``` sh
julia --project=docs docs/make.jl
```

Please, update the documentation if you add new features or change the behavior
of existing ones.

### Pull Request (PR) and commits

Here's how to structure your contributions effectively:

- Descriptive Title: Briefly summarize the changes your PR introduces. Commit
  titles should preferably be under 50 characters, start with a capital latter,
  and use imperative verbs (e.g., "Remove superfluous function call").
- Detailed Description: Explain the purpose of your changes. Focus on the
  intent.
- Breaking Changes: If your PR introduces breaking changes, highlight them
  clearly in the description.

Your pull request can contain one or multiple commits. In either cases, it is
important that each commit is atomic (meaning that each commit represents a
single logical change).

Please, squash commits that represent a single logical change (e.g., do not have
two commits when the second just fixes the first).

### Credits

The `Space` and `TimeVaryingInputs` modules were initially developed in the
context of [`ClimaLand`](https://github.com/CliMA/ClimaLand.jl), the
`TempestRegridder` and `TimeManager` ones were initially developed in
[`ClimaCoupler`](https://github.com/CliMA/ClimaCoupler.jl).
