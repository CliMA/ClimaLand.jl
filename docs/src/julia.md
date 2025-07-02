# Julia Background

New to Julia? We're glad you've found us, and hope you'll like it as much as we do!

The [Julia manual](https://docs.julialang.org/en/v1/manual/getting-started/) is a very
helpful resource, and we recommend reading through it while getting started with the language.

If you don't have it already, Julia can be downloaded from [https://julialang.org/downloads/](https://julialang.org/downloads/).

## Julia REPL

The Julia REPL provides a convenient interface to interactively run Julia code,
similar in spirit to a Jupyter notebook. To start a Julia REPL, simply run the command
`julia` from your terminal.

Many of our tutorials are designed to be run within the REPL, to see live output
as you set up and run the simulation (though they can also be run from the command line).

## Julia environments

Julia environments are specified by a Project.toml and a Manifest.toml.
The Project.toml lists all packages required to setup the environment,
as well as a list of package version requirements (under `compat`),
and additional optional information.
The Manifest.toml is automatically generated from the Project.toml and
contains more detailed information about the specific package versions
to be used when setting up the environment.

Additional details about Julia environments can be found on the
[Code Loading](https://docs.julialang.org/en/v1/manual/code-loading/#Environments-1) manual page.

The ClimaLand.jl package contains three Julia environments: the top-level ClimaLand.jl environment,
a `docs/` environment, and a `.buildkite/` environment.
When running simulations (e.g. tests, experiments, or other simulations), the `.buildkite`
environment should be used.
You can enter this environment using the following command from the terminal:
```
julia --project=.buildkite
```

## Handy Julia tools

### Revise.jl

One downside to Julia is that it requires (often time-consuming) precompilation.
During local development, any source code that has been changed will need to be precompiled again
to incorporate the new changes. Given that our package depends on many others, this can take some time.

[Revise.jl](https://github.com/timholy/Revise.jl) has been developed to minimize this latency, and is very useful in doing so.
Their [documentation](https://timholy.github.io/Revise.jl/stable/) provides more information
about how the package works.

To install Revise and load it whenever you start Julia, first add it to your base Julia environment.
Enter your base Julia environment by simply running `julia` from the terminal, then run the following:
```
julia> using Pkg
julia> Pkg.add(Revise)
```

Next, create the file `~/.julia/config/startup.jl` and add the following line to it: `using Revise`.
This will load Revise anytime you start a Julia session in any local Julia environment.

## Juliaup

Juliaup is a useful tool to manage installed Julia versions. It allows you to easily download new versions,
remove old ones, set a default version, and switch between multiple installed versions.

Please see the [documentation](https://github.com/JuliaLang/juliaup?tab=readme-ov-file#juliaup---julia-version-manager)
for more information and details about installing juliaup.
