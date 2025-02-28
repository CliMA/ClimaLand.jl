# Calibration framework of the CliMA land model

## Introduction

In general, models attempt to reproduce real world observations.
Calibration is the process of finding the parameter set that will best reproduce real world observation.

The key ingredients in calibration are:
- the model we want to run;
- the parameters we want to tune along with their prior distributions;
- the observational data we want to reproduce and how that data is represented in the model;
- the noise associated to such observational data.

The process of calibrating consists of optimizing how different parameters match the given observations within the given noise. In ClimaLand, we use `EnsembleKalmanProcesses.jl` (EKP) to perform automatic calibration.
Before discussing the details of EKP, we introduce the following terminology to mirror the key ingredients introduced above:

- a `forward_model` that runs the model for a given parameter set,
- an `observation` vector that contains the data or statistics of data we want the model to reproduce,
- an `observation_map` which maps the equivalent of the observation vector, but from the model output,
- `priors` which contains the parameters we want to calibrate (find the value that makes the model match the observation best), priors gives which parameters, but also their distribution
- a covariance matrix, that defines the observational error and correlations

[EnsembleKalmanProcesses.jl](https://github.com/CliMA/EnsembleKalmanProcesses.jl) (EKP) is at the center of CliMA's calibration efforts. EKP implements a suite of Ensemble Kalman methods to find a (locally) optimal parameter set `U` for a model `G` to fit noisy `Γ` observational data `Y`. These methods are optimized for problems where the model `G` is computationally expensive and no analtyic derivatives are available, as in the case of weather forcasting, where Ensemble Kalman techniques have a long history of success.

Large calibration campaigns often require supercomputers and while direct use of EKP.jl is possible, CliMA's preferred approach is using [ClimaCalibrate.jl](https://github.com/CliMA/ClimaCalibrate.jl), a package optimized for running on compute clusters. `ClimaCalibrate` handles efficient job orchestration and abstracts the details of the underlying system, providing a simpler user experience. Consult the [ClimaCalibrate documentation](https://clima.github.io/ClimaCalibrate.jl/dev/) for further information.

## Calibrate a land model

In this tutorial, we will perform a calibration using `ClimaCalibrate`. `ClimaCalibrate` provides an interface to `EnsembleKalmanProcesses.jl` that is more optimized for use on supercomputers. The [tutorial to calibrate a single site latent heat flux](https://clima.github.io/ClimaLand.jl/stable/generated/calibration/minimal_working_example_obs/) shows how to perform a calibration using `EKP` directly.

The `calibrate` function is at the heart of performing a calibration with `ClimaCalibrate`:

```julia
import ClimaCalibrate as CAL

CAL.calibrate(
    CAL.WorkerBackend,
    utki,
    n_iterations,
    prior,
    caldir,
)
```

where the utki object defines your EKP configurations, for example, the default is:

```julia
    observationseries,
    EKP.TransformUnscented(prior, impose_prior = true);
    verbose = true,
    rng,
    scheduler = EKP.DataMisfitController(terminate_at = 100),
```

where `observationseries` is the "truth" you want to calibrate your model on. It can take many forms.
For example, you may want to calibrate your land model latent heat flux (lhf), the
observations could be monthly global average of lhf, or monthly average at 100 random locations on land, or the annual amplitude and phase...
You will create `observationsseries` from some data (for example ERA5), as a vector.

note that `observationseries` object contains the covariance matrix of the noise, which informs the uncertainties in space and time of your targeted "truth". It can be set, for example, to the inter-annual variance of a variable, or to the average of the variable times a % (e.g., 5%), or to a flat noise (for example, 5 W m-2 for latent heat). This will inform the EKP algorithm that if the model is within the target +- noise at specific space and time, the goal is reached.

`TransformUnscented.Unscented` is a method in EKP, that requires `2 x number of parameters + 1` ensemble members (`ensembe_size`, the number of parameter set drawn for your prior distribution tested at each iteration). For more information, read the [EKP documentation for that method](https://clima.github.io/EnsembleKalmanProcesses.jl/dev/unscented_kalman_inversion/).

`verbose = true` is a setting that writes information about your calibration run to a log file.

`rng` is a set random seed.

`Scheduler` is a EKP setting for timestepping, please read [EKP schedulers documentations](https://clima.github.io/EnsembleKalmanProcesses.jl/dev/learning_rate_scheduler/).

`CAL.WorkerBackend` defines how to interact with the underlying compute system. For other possible backends (for example, `JuliaBackend`, `ClimaGPUBackend`, or `DerechoBackend`),
see [this doc page](https://clima.github.io/ClimaCalibrate.jl/dev/backends/).

Each backend is optimized for specific use cases and computing resources. The backend system is implemented through Julia's multiple dispatch,
so that code written for one computer can seamlessly be ported to a new/different environments.

`prior` is the distribution of the parameters you want to calibrate. For example, if you want to calibrate two parameters called `sc` and `pc`,
you would define your priors like this, for example:
```julia
prior_sc = EKP.constrained_gaussian("sc", 5e-6, 5e-4, 0, Inf);
prior_pc = EKP.constrained_gaussian("pc", -2e6, 1e6, -Inf, Inf);
prior = EKP.combine_distributions([prior_sc, prior_pc]);
```
For more documentation about prior distribution, see [this EKP documentation page](https://clima.github.io/EnsembleKalmanProcesses.jl/dev/parameter_distributions/).

`n_iterations` is the number of times your priors distribution will be updated, at each iteration your model is run for the number of `ensemble_size`.
So in total, your model will be run `ensemble_size` * `n_iterations`.

`caldir` is the path to your calibration output directory. For example `calibration_output`. Inside this folder, the parameter set of each iteration * member will
be stored, as well as the output of your model simulations. For example, if you ran a calibration with 1 iteration and 2 members, caldir would be structured
like this:

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
each iteration contains folders for each member, inside which you can find the parameters value inside `parameters.toml`, and model outputs inside `global_diagnostics`.

Two additional functions need to be defined in order to run `CAL.calibrate`. `CAL.forward_model(iteration, member)` and `observation_map(iteration)`.
The `CAL.forward_model(iteration, member)` needs to generate your model output for a specific iteration and member. The `observation_map(iteration)`
needs to return your loss, a vector of the same format as `observations` but created with your model output (for example, monthly average of latent heat flux),
for all members. To make this easier, it can be useful to implement a `process_member_data(root_path)` function that generates one member output from your
model output path.

Once you have defined `CAL.forward_model`, `CAL.observation_map`, `caldir`, `noise`, `observations`, `n_iterations`, `ensemble_size`, and your backend, you can
call `CAL.calibrate`!

## job script

A calibration job will likely take hours to complete, so you will probably have to submit a job with a job scheduler.
Below is an example job .pbs script, slurm would be slightly different.

```
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

julia --project=.buildkite/ experiments/calibration/global_land/calibrate_land.jl
```

where `calibrate_land.jl` is a script that generates all the arguments needed and eventually calls `CAL.calibrate`.
You would start the job with a command such as `qsub name_of_job_script`, and a few hours later, you would get a calibrated parameter set.

Note that with the default EKP configuration, UTKI, the number of ensemble is set by the number of parameters, as explained in the documentation above. The number of workers (if you use the worker backend) is automatically set to that numbers, so that all members are run in parallel for each iteration.

## Configure your calibration job

To run a calibration, you can modify the following objects in `calibrate_land.jl`:
- include `forward_model.jl` or `forward_model_bucket.jl`, depending on which model you want to calibrate
- `variable_list`: choose the variables you want to calibrate. For example, `["swu"]` or `["lhf", "shf"]`
- `n_iterations`: how many iterations you want to run
- `spinup_period`: the time length of your spinup
- `start_date`: the start date of the calibration
- `nelements`: to adjust the model resolution (n_horizontal elements, n_vertical elements)
- `caldir`: you can change the name and path of your calibration output directory
- `training_locations`: the default is all locations on land, but you can change this to a subset of land coordinates

You should also modify the files:
- `priors.jl`, to give your parameters and their distribution
- `forward_mode.jl` or `forward_model_bucket.jl`, to ensure it reads the parameters from `priors.jl` and uses them

other things to consider:
- `observationseries_era5.jl`: currently the target data is era5. This could be changed to another dataset, or a perfect model target
- `observationseries_era5.jl`: the noise is currently set to 5^2 (flat noise, in W m^-2), this should be changed to desired noise
- you could change the EnsembleKalmanProcesses settings, set in `EKP.EnsembleKalmanProcesses()` in `calibrate_land.jl`
- you can adjust the ClimaCalibrate backend, see the [documentation](https://clima.github.io/ClimaCalibrate.jl/dev/backends/)
- Because the variable names are currently different in era5 file, you may have to add variables to map in observationseries_era5.jl (for example, "lhf" => "mslhf")

Also, note that currently:
- the temporal resolution of observation_maps are seasonal (3 months average)
- the length of calibration is 1 year (data used after spinup)
- the entire globe is used
These could also be changed in the code, but would currently requires significant changes.

Finally, note that the HPC job script and command will slightly differ between slurm (for example central or clima) and pbs (for example derecho).
