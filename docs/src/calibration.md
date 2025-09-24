# Calibration framework of the CliMA land model

## Introduction

In general, models attempt to reproduce real world observations. Calibration is
the process of finding the parameter set that will best reproduce real world
observation.

The key ingredients in calibration are:
- the model we want to run;
- the parameters we want to tune along with their prior distributions;
- the observational data we want to reproduce and how that data is represented
  in the model;
- the noise associated to such observational data.

The process of calibrating consists of optimizing how different parameters match
the given observations within the given noise. In ClimaLand, we use
`EnsembleKalmanProcesses.jl` (EKP) to perform automatic calibration. Before
discussing the details of EKP, we introduce the following terminology to mirror
the key ingredients introduced above:

- a `forward_model` that runs the model for a given parameter set,
- an `observation` vector that contains the data or statistics of data we want
  the model to reproduce,
- an `observation_map` which maps the equivalent of the observation vector, but
  from the model output,
- `priors` which contains the parameters we want to calibrate (find the value
  that makes the model match the observation best), priors gives which
  parameters, but also their distribution
- a covariance matrix that defines the observational error and correlations

[EnsembleKalmanProcesses.jl](https://github.com/CliMA/EnsembleKalmanProcesses.jl)
(EKP) is at the center of CliMA's calibration efforts. EKP implements a suite of
Ensemble Kalman methods to find a (locally) optimal parameter set `U` for a
model `G` to fit noisy `Γ` observational data `Y`. These methods are optimized
for problems where the model `G` is computationally expensive and no analtyic
derivatives are available, as in the case of weather forcasting, where Ensemble
Kalman techniques have a long history of success.

Large calibration campaigns often require supercomputers and while direct use of
EKP.jl is possible, CliMA's preferred approach is using
[ClimaCalibrate.jl](https://github.com/CliMA/ClimaCalibrate.jl), a package
optimized for running on compute clusters. `ClimaCalibrate` handles efficient
job orchestration and abstracts the details of the underlying system, providing
a simpler user experience. Consult the
[ClimaCalibrate documentation](https://clima.github.io/ClimaCalibrate.jl/stable/)
for further information.

## Calibrate a land model

In this tutorial, we will perform a calibration using `ClimaCalibrate`.
`ClimaCalibrate` provides an interface to `EnsembleKalmanProcesses.jl` that is
more optimized for use on supercomputers. The
[tutorial to calibrate a single site latent heat flux](https://clima.github.io/ClimaLand.jl/stable/generated/calibration/minimal_working_example_obs/)
shows how to perform a calibration using `EKP` directly.

The `calibrate` function is at the heart of performing a calibration with
`ClimaCalibrate`:

```julia
import ClimaCalibrate

ClimaCalibrate.calibrate(
    ClimaCalibrate.WorkerBackend,
    utki,
    n_iterations,
    prior,
    output_dir,
)
```

where the `utki` object defines your EKP configurations, for example, the
default is:

```julia
EKP.EnsembleKalmanProcess(
    obs_series,
    EKP.TransformUnscented(prior, impose_prior = true),
    verbose = true,
    rng = rng,
    scheduler = EKP.DataMisfitController(terminate_at = 100),
)
```

where `obs_series` is the "truth" you want to calibrate your model on. It can
take many forms. For example, you may want to calibrate monthly global average
of latent heat flux, monthly average at 100 random locations on land, or the
annual amplitude and phase. You will create `obs_series` from some data (for
example ERA5), as a vector.

Note that `obs_series` object contains the covariance matrix of the noise, which
informs the uncertainties in space and time of your targeted "truth". It can be
set, for example, to the inter-annual variance of a variable, or a percentage
(e.g., 5%) of the average of the variable, or to a flat noise (e.g., 5 W m-2 for
latent heat). This will inform the EKP algorithm that if the model is within the
target plus or minus the noise at specific space and time, the goal is reached.

- `TransformUnscented.Unscented` is a method in EKP, that requires `2p + 1`
  ensemble members for each iteration, where `p` is the number of parameters.
  For more information, read the
  [EKP documentation for that method](https://clima.github.io/EnsembleKalmanProcesses.jl/stable/unscented_kalman_inversion/).

- `verbose = true` is a setting that writes information about your calibration
  run to a log file.

- `rng` is a set random seed.

- `Scheduler` is a EKP setting for timestepping, please read
  [EKP schedulers documentation](https://clima.github.io/EnsembleKalmanProcesses.jl/stable/learning_rate_scheduler/).

- `ClimaCalibrate.WorkerBackend` defines how to interact with the underlying
  compute system. For other possible backends (for example, `JuliaBackend`,
  `ClimaGPUBackend`, or `DerechoBackend`), see the
  [backend documentation in ClimaCalibrate](https://clima.github.io/ClimaCalibrate.jl/stable/backends/).

Each backend is optimized for specific use cases and computing resources. The
backend system is implemented through Julia's multiple dispatch, so that code
written for one environment can seamlessly be ported to a new/different
environments.

- `prior` is the distribution of the parameters you want to calibrate. For
  example, if you want to calibrate two parameters called `sc` and `pc`, you
  would define your priors like this, for example:
```julia
prior_sc = EKP.constrained_gaussian("sc", 5e-6, 5e-4, 0, Inf);
prior_pc = EKP.constrained_gaussian("pc", -2e6, 1e6, -Inf, Inf);
prior = EKP.combine_distributions([prior_sc, prior_pc]);
```
For more documentation about prior distribution, see [this EKP documentation page](https://clima.github.io/EnsembleKalmanProcesses.jl/stable/parameter_distributions/).

- `n_iterations` is the number of times your priors distribution will be
  updated, at each iteration your model is run for  the number of
  `ensemble_size`. So in total, your model will be run `ensemble_size` *
  `n_iterations`.

- `output_dir` is the path to your calibration output directory. Inside this
  folder, the parameter set of each ensemble member for each iteration is
  stored, as well as the output of your model simulations. For example, if you
  ran a calibration with 1 iteration and 2 members, output_dir would be
  structured like this:

```
.
├── iteration_000
│   ├── eki_file.jld2
│   ├── G_ensemble.jld2
│   ├── member_001
│   │   ├── global_diagnostics
│   │   │   ├── output_0000
│   │   │   │   ├── lhf_1M_average.nc
│   │   │   │   ├── lwu_1M_average.nc
│   │   │   │   ├── shf_1M_average.nc
│   │   │   │   └── swu_1M_average.nc
│   │   │   └── output_active -> output_0000
│   │   └── parameters.toml
│   ├── member_002
│   │   ├── global_diagnostics
│   │   │   ├── output_0000
│   │   │   │   ├── lhf_1M_average.nc
│   │   │   │   ├── lwu_1M_average.nc
│   │   │   │   ├── shf_1M_average.nc
│   │   │   │   └── swu_1M_average.nc
│   │   │   └── output_active -> output_0000
│   │   └── parameters.toml
├── iteration_001
│   ├── eki_file.jld2
│   ├── G_ensemble.jld2
│   ├── member_001
│   │   ├── global_diagnostics
│   │   │   ├── output_0000
│   │   │   │   ├── lhf_1M_average.nc
│   │   │   │   ├── lwu_1M_average.nc
│   │   │   │   ├── shf_1M_average.nc
│   │   │   │   └── swu_1M_average.nc
│   │   │   └── output_active -> output_0000
│   │   └── parameters.toml
│   ├── member_002
│   │   ├── global_diagnostics
│   │   │   ├── output_0000
│   │   │   │   ├── lhf_1M_average.nc
│   │   │   │   ├── lwu_1M_average.nc
│   │   │   │   ├── shf_1M_average.nc
│   │   │   │   └── swu_1M_average.nc
│   │   │   └── output_active -> output_0000
│   │   └── parameters.toml
```
Each iteration contains directories for each member, inside which you can find
the parameters value inside `parameters.toml`, and model outputs inside
`global_diagnostics`.

Two additional functions need to be defined in order to run
`ClimaCalibrate.calibrate`. `ClimaCalibrate.forward_model(iteration, member)`
and `ClimaCalibrate.observation_map(iteration)`. The
`ClimaCalibrate.forward_model(iteration, member)` needs to generate your model
output for a specific iteration and member. The
`ClimaCalibrate.observation_map(iteration)` needs to return your loss, a vector
of the same format as `observations` but created with your model output (for
example, monthly average of latent heat flux), for all members. To make this
easier, it can be useful to implement a `process_member_data(root_path)`
function that generates one member output from your model output path.

Once you have defined `ClimaCalibrate.forward_model`,
`ClimaCalibrate.observation_map`, `output_dir`, `noise`, `observations`,
`n_iterations`, `ensemble_size`, and your backend, you can call
`ClimaCalibrate.calibrate`!

## Job script

A calibration job will likely take hours to complete, so you will probably have
to submit a job with a job scheduler. Below is an example job .pbs script (for
PBS, e.g., Derecho):

```bash
#!/bin/bash
#PBS -N derecho_calibration
#PBS -o output.txt
#PBS -e error.txt
#PBS -l walltime=10:00:00
#PBS -l select=1:ncpus=64:ngpus=4

## Account number for CliMA
#PBS -A UCIT0011
#PBS -q main

module load climacommon
# Run your command
export CLIMACOMMS_DEVICE="CUDA"
export CLIMACOMMS_CONTEXT="SINGLETON"
julia --project=.buildkite -e 'using Pkg; Pkg.instantiate(;verbose=true)'

julia --project=.buildkite/ experiments/calibration/run_calibration.jl
```

and below is a example job .sh script (for Slurm, e.g., central or clima):

```bash
#!/bin/bash
#SBATCH --job-name=slurm_calibration
#SBATCH --output=output.txt
#SBATCH --error=error.txt
#SBATCH --time=12:00:00
#SBATCH --ntasks=5
#SBATCH --cpus-per-task=4
#SBATCH --gpus-per-task=1

# Set environment variables for CliMA
export CLIMACOMMS_DEVICE="CUDA"
export CLIMACOMMS_CONTEXT="SINGLETON"

# Build and run the Julia code
module load climacommon
julia --project=.buildkite -e 'using Pkg; Pkg.instantiate(;verbose=true)'
julia --project=.buildkite/ experiments/calibration/run_calibration.jl
```

where `run_calibration.jl` is a script that set up the calibration and call
`ClimaCalibrate.calibrate`. You would start the job with a command such as `qsub
name_of_job_script` for PBS or `sbatch name_of_job_script` for Slurm, and a few
hours later, you would get a calibrated parameter set. You can check the status
of your job with `qstat -u username` of PBS or `squeue -u username` on Slurm.

Note that with the default EKP configuration, UTKI, the number of ensemble is
set by the number of parameters, as explained in the documentation above. The
number of workers (if you use the worker backend) is automatically set to that
numbers, so that all members are run in parallel for each iteration.

## Configure your land calibration

For configuring the land calibration, you can configure the optimization method
used by EnsembleKalmanProcesses.jl, the observations and how they are
preprocessed, and the simulation itself. For most settings, you can modify
the `CalibrateConfig` struct. See the example below.

```julia
CalibrateConfig(;
    short_names = ["swu"],
    minibatch_size = 2,
    n_iterations = 10,
    sample_date_ranges = [("2007-12-1", "2008-9-1"),
                          ("2008-12-1", "2009-9-1"),
                          ("2009-12-1", "2010-9-1"),
                          ("2010-12-1", "2011-9-1")],
    extend = Dates.Month(3),
    spinup = Dates.Month(3),
    nelements = (101, 15),
    output_dir = "experiments/calibration/land_model",
    rng_seed = 42,
)
```

With the configuration above, a calibration is being done using the `swu`
observation. The calibration will run for 10 iterations, and each iteration will
use a minibatch size of 2. The start and end dates of the simulation is
automatically determined by `sample_date_ranges`, `spinup`, and `extend`. The
amount to spin up the simulation for is three months and the amount to run the
simulation for past the dates in `sample_date_ranges` is also three months.
Since the minibatch size is 2, the first iteration of the calibration will run
from 1 September 2007 to 1 December 2009 and the second iteration of the
calibration is 1 September 2009 to 1 December 2011. Afterward, the iterations
will repeat, so the third iteration will be the same as the first iteration, the
fourth iteration will be the same as the second iteration, and so on.

The period chosen for `extend` is ensure that all the data is gathered. If you
are calibrating against seasonal averages, then `extend` should be 3 months and
if you are calibrating against monthly averages, then `extend` should be 1
month. For more information, see the sections [Data pipelines](@ref) and
[Simulation settings](@ref).

The number of horizontal and vertical elements to use for the simulation is
determined by `nelements`. The calibration is saved at `output_dir` and the
random number generator is seeded by `rng_seed`.

### EKP settings

In `CalibrateConfig`, you can modify the number of iterations via `n_iterations`
and the size of the minibatch by `minibatch_size`. For reproducibility, you can
pass in an integer for `rng_seed` which is used internally by
EnsembleKalmanProcesses.jl to seed the random number generator.

For the land calibration, the default optimization method is
`EnsembleKalmanProcesses.TransformUnscented` and the default scheduler is
`EKP.DataMisfitController`. To change this optimization method or change the
scheduler, you need to go to the `experiments/calibration/run_calibration.jl`
file and modify it there. For more information about the different optimization
methods, see
[EnsembleKalmanProcesses.jl](https://clima.github.io/EnsembleKalmanProcesses.jl/stable/)
for more information.

### Observation settings

In `CalibrateConfig`, you can modify which observations that are used via
`short_names`. As of now, the currently supported observations for calibration
are seasonal averages of `lhf`, `shf`, `lwu`, and `swu`. See the
[Data pipelines](@ref) on how to add more variables.

Observations are generated with `generate_observations.jl`. To generate
observations, you can run

```julia
julia --project=.buildkite experiments/calibration/generate_observations.jl`
```

!!! note "When do I need to generate observations?"
    Whenever `nelements`, `sample_date_ranges`, `short_names`, or the
    preprocessing of any of the variables change, then you must regenerate the
    observations!

Each observation contains time series data of the variables from ERA5 and a
covariance matrix. You can adjust `sample_date_ranges` if you want to calibrate
over a particular year or over multiple years for example.

!!! note "What can be passed for `sample_date_ranges`?"
    The sample date ranges should be times from the time series data of the
    observations. For example, when using seasonal averages, the times passed
    must be the first day of December, March, June, or September. The seasons
    are December to February (DJF), March to May (MAM), June to August (JJA),
    and September to November (SON). The convention for time is use the first
    time of the time reduction.

To change the covariance matrix, you need to go to `generate_observations.jl`
and modify the covariance matrix that is passed in. See
[ClimaCalibrate.jl documentation](https://clima.github.io/ClimaCalibrate.jl/stable/api/#Observation-Recipe-Interface)
for a list of different covariance matrices that can be used.

#### Data pipelines

!!! note "Data preprocessing"
    The data preprocessing uses ClimaAnalysis.jl. See the
    [ClimaAnalysis.jl documentation](https://clima.github.io/ClimaAnalysis.jl/stable/api/)
    for functions you can use for preprocessing.

To add additional variables or change how a variable is preprocessed, you must
add the variable to the data sources and specify how the variable should be
preprocessed. The name of the variable should match the name in the diagnostics.
As of now, each variable from the ERA5 data is loaded by the function
`get_era5_obs_var_dict` in `ext/land_sim_vis/leaderboard/data_sources.jl`. See
the documentation of `get_obs_var_dict` for how to add a new variable. In
addition, you should also add the same variable to `get_sim_var_dict` in the
same file as before.

Further preprocessing is done in the function `preprocess_single_era5_var` in
`generate_observations.jl`. After adding the variable, you should specify any
additional preprocessing. For example, if you want seasonal averages, you should
use `ClimaAnalysis.average_season_across_time`. Furthermore, you should check
the units of the variable, ensure that the grid and mask are the same between
the simulation and observational data, and ensure the conventions for times are
the same. Finally, the same preprocessing should be applied between both
observational and simulation data. Preprocessing the observational data is done
in `preprocess_single_era5_var` in
`experiments/calibration/generate_observations.jl` and preprocessing the
simulation data is done in `process_member_data` in
`experiments/calibration/observation_map.jl`.

!!! note "Time conventions"
    Different data sources have different conventions for time. For example,
    data sources use January 1, January 15, and Feburary 1 for the monthly
    average of January. For the diagnostics saved from a CliMA simulation, the
    current convention is to save the monthly average on the Feburary 1. You
    must ensure that the time conventions are the same. For calibration, we
    choose the start of the time reduction, so for this example, the time
    associated with the monthly average of January is January 1. For seasonal
    averages, the dates would be December 1, March 1, June 1, and September 1.

### Simulation settings

In `CalibrateConfig`, you can modify the resolution of the model via
`nelements`. To change the directory of where the calibration starts, you can
change `output_dir`, whose default is `experiments/calibration/land_model`.

For adjusting the start and end dates of the simulation, you can change `spinup`
and `extend`.

!!! note "How are the start and end dates of the simulation chosen?"
    The start and end dates of a simulation is automatically inferred from the
    `sample_date_ranges`. For example, suppose that `sample_date_ranges =
    [("2007-12-1", "2008-9-1"), ("2008-12-1", "2008-9-1)]` and the minibatch
    size is 2. The start and end dates are `2007-12-1` to `2008-9-1`. This
    should cover the time ranges of the observational data. However, observe
    that the simulation output may not be realistic at the beginning of the
    simulation while the land model is equilibrating and the end date is not
    correct, because of the time convention we chosen earlier. To solve the
    first problem, you can specify how the long the model should spin up for
    with `spinup` in the `CalibrateConfig`. For the second problem, you can make
    the simulation extend past `2008-9-1` with `extend`. For calibrating with
    monthly averages, `extend` should be `Dates.Month(1)` and for calibrating
    with seasonal averages, `extend` should be `Dates.Month(3)`.

## Which parameters can I calibrate?
ClimaLand.jl provides a full list of spatially-constant parameters that may be
calibrated in `toml/default_parameters.toml`. For each parameter, this file
includes the parameter name, type, default value, units, and model or
parameterization it is used for.
Please see that file to see which parameters could be calibrated
for your model of interest.
