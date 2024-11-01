# ClimaLand Experiments

These experiment scripts are primarily intended for testing models via
simulation, often comparing to other model output or to data. These run
as part of CI on the Caltech cluster, and must run without error in order
to pass. Many create plots and these plots should be visually assessed before
merging a PR. These experiments do not take the place of unit tests or
tutorials.

In our CI pipeline we use the environment defined by ClimaLand.jl/.buildkite/Manifest.toml and associate Project.toml.  From ClimaLand.jl, run

``` bash
julia --project=.buildkite path_to_experiment_script
```

## Experiment Output Paths

Many of the outputs from the experiments are saved using a non-destructive approach called
`ActiveLinkStyle`. Each experiment has some unique save directory. For example, an experiment might have an
output directory of `example_outdir`. When run, if the directory does
not exist relative to the working directory, then `example_outdir` is created. The results of the latest experiment
run will be saved into `example_outdir/output_xxx`, where xxxx is an increasing counter. The lower indices store the outputs of previous runs. The latest results are also linked to
in `example_outdir/output_active`, which can be assumed to always contain the most recent output. More details on this style of output directory handling
can be found [here](https://clima.github.io/ClimaUtilities.jl/dev/outputpathgenerator/#ActiveLinkStyle-(Non-Destructive))

### benchmarks

The experiments in experiments/benchmarks generate and save files relative to the working directory.
When a benchmark is run, its outputs are saved using `ActiveLinkStyle` into (benchmark_name)\_benchmark\_(device_suffix).

### long_runs

The experiments in experiments/long_runs save its image outputs to (experiment_name)\_longrun\_(device_suffix).
All other outputs are saved using `ActiveLinkStyle` to (experiment_name)\_longrun\_(device_suffix)/global_diagnostics.

### integrated and standalone

The experiments in experiments/long_runs and experiments/standalone have save paths that use `ActiveLinkStyle` and assume the experiments are being run from
the ClimaLand.jl directory. If run from the ClimaLand.jl directory, then the outputs of an experiment will be saved into
the same directory the experiment script is found in. Otherwise the outputs will be saved to
the path of the script relative to the ClimaLand.jl directory.
