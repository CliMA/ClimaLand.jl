# Calibration framework of the CliMA land model

## Introduction

In general, models attempt to reproduce real world observations.
Calibration is the process of finding the parameter set that will best reproduce real worl observation.

To put it simply, a calibration framework needs observation data (for example, monthly latent heat flux from ERA5, for n locations, and m years) that it will attempt to reproduce via the model with the best parameters. In our framework, we need to define:
- a `forward_model` that runs the model for a given parameter set,
- an `observation` vector that contains the data or statistics of data we want the model to reproduce,
- an `observation_map` which maps the equivalent of the observation vector, but from the model output,
- `priors` which contains the parameters we want to calibrate (find the value that makes the model match the observation best), priors gives which parameters, but also their distribution.

[EnsembleKalmanProcesses.jl](https://github.com/CliMA/EnsembleKalmanProcesses.jl) (EKP) is at the center of CliMA's calibration efforts. EKP implements a suite of Ensemble Kalman methods to find a (locally) optimal parameter set `U` for a model `G` to fit noisy `Γ` observational data `Y`. These methods are optimized for problems where the model `G` is computationally expensive and no analtyic derivatives are available, as in the case of weather forcasting, where Ensemble Kalman techniques have a long history of success.

Large calibration campaigns often require supercomputers and while direct use of EKP.jl is possible, CliMA's preferred approach is using [ClimaCalibrate.jl](https://github.com/CliMA/ClimaCalibrate.jl), a package optimized for running on compute clusters. `ClimaCalibrate` handles efficient job orchestration and abstracts the details of the underlying system, providing a simpler user experience. Consult the [ClimaCalibrate documentation](https://clima.github.io/ClimaCalibrate.jl/dev/) for further information.

## Calibrate a land model

In this tutorial, we will perform a calibration using `ClimaCalibrate`. `ClimaCalibrate` provides an interface to `EnsembleKalmanProcesses.jl` that is more optimized for use on supercomputers. The [tutorial to calibrate a single site latent heat flux](https://clima.github.io/ClimaLand.jl/stable/generated/calibration/minimal_working_example_obs/) shows how to perform a calibration using `EKP` directly.

The `calibrate` function is at the heart of performing a calibration with `ClimaCalibrate`:

```julia
import ClimaCalibrate as CAL

CAL.calibrate(
    CAL.WorkerBackend,
    ensemble_size,
    n_iterations,
    observations,
    noise,
    prior,
    caldir,
)
```
where `CAL.WorkerBackend` defines how to interact with the underlying compute system. For other possible backends (for example, `JuliaBackend`, `ClimaGPUBackend, or `DerechoBackend`),
see [this doc page](https://clima.github.io/ClimaCalibrate.jl/dev/backends/).

Each backend is optimized for specific use cases and computing resources. The backend system is implemented through Julia's multiple dispatch,
so that code written for one computer can seamlessly be ported to a new/different environments.

`prior` is the distribution of the parameters you want to calibrate. For example, if you want to calibrate two parameters called sc and pc,
you would define your priors like this, for example:
```julia
prior_sc = EKP.constrained_gaussian("sc", 5e-6, 5e-4, 0, Inf);
prior_pc = EKP.constrained_gaussian("pc", -2e6, 1e6, -Inf, Inf);
prior = EKP.combine_distributions([prior_sc, prior_pc]);
```
For more documentation about prior distribution, see [this EKP documentation page](https://clima.github.io/EnsembleKalmanProcesses.jl/dev/parameter_distributions/).

`ensemble_size` is the number of parameter set drawn for your prior distribution tested at each iteration. For example, if you set ensemble_size to 5,
your model will be run 5 times (5 member) at each iteration, using 5 different parameter set drawn from your priors distribution. The distribution is
updated at each iteration by the EKP algorithm.

`n_iterations` is the number of times your priors distribution will be updated, at each iteration your model is run for the number of ensemble\_size.
So in total, your model will be run `ensemble_size` * `n_iterations`.

`observations` is the "truth" you want to calibrate your model on. It can take many forms.
For example, you may want to calibrate your land model latent heat flux (lhf), the
observations could be monthly global average of lhf, or monthly average at 100 random locations on land, or the annual amplitude and phase...
You will create `observations` from some data (for example ERA5), as a vector.

`noise` is the variance of each element of your `observations` vector. For example, if `observations` is monthly average of latent heat flux in a specific year,
your noise has to be the standard deviation of each month (over multiple years).

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

## Results exploration

After running a calibration, you can explore the results via CliCal dashboard, which allows users to navigate between:
- calibration runs,
- variables,
- iterations,
- ensemble member,
- and seasons.

For each combination of these options, the dashboard will provide:
- a table of the parameter values and how they changed relative to initial values,
- EKP error for each iterations and how this error changed relative to initial error,
- RMSE overall (all variables) and RMSE for the selection (variable, ensemble), for each iteration.
- A figure (global map) or Y (e.g., era5 data) and G (ClimaLand output),
- A figure of the seasonal average of Y and G,
- A figure of anomalies G - Y.

To generate such dashboard, the user can run the script `make_dashboard.jl` directly from the HPC where the calibrations ran, if they connected via port forwarding, e.g., ssh -L 8888:localhost:8888 your_user@hpc-address, then the user can open the following URL on their local computer: http://localhost:8888/browser-display

Such dashboard can also be served on the internet, to be accessed anywhere anytime via a URL, for example [clickme](https://clima.westus3.cloudapp.azure.com/jsserve/calibration_dash):

```@raw html
<iframe src="https://clima.westus3.cloudapp.azure.com/jsserve/calibration_dash"
   style="height:1200px;width:90%;">
</iframe>
```
