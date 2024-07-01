## ClimaLand Experiments

These experiment scripts are primarily intended for testing models via
simulation, often comparing to other model output or to data. These run
as part of CI on the Caltech cluster, and must run without error in order
to pass. Many create plots and these plots should be visually assessed before
merging a PR. These experiments do not take the place of unit tests or
tutorials.

In our CI pipeline we use the environment defined by ClimaLand.jl/.buildkite/Manifest.toml and associate Project.toml.  From ClimaLand.jl, run

```
julia --project=.buildkite path_to_experiment_script
```