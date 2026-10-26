# Calibration on ERA5 columns

This experiment calibrates the land model of the global calibration
(`experiments/calibration/models/snowy_land.jl`: soil, canopy with the P-model
and Weibull moisture stress, prescribed MODIS LAI, snow with the neural depth
and albedo models and the topographic snow cover fraction, slab lake, soil CO2)
forced by ERA5, with the calibrated variables, noise, prior, and number of
iterations of the global configuration `configs/energy_fluxes.jl`. Instead of
the global grid, it runs a `ColumnEnsemble` of 200 columns chosen to represent
global land, which makes each forward run cheap.

- **Observations:** seasonal means from March 2008 to February 2010 of the
  variables of `energy_fluxes.jl`, from the datasets of the global calibration:
  ERA5 latent and sensible heat fluxes and upwelling shortwave and longwave
  radiation, FLUXCOM GPP (from the `ilamb_data` artifact), and the MODIS snow
  cover fraction of ILAMB.
- **Simulations:** March 2007 to March 2010, with the first year as spinup. The
  dates of `energy_fluxes.jl` only apply to the global calibration.
- **Parameters:** the prior of `energy_fluxes.jl`.
- **Algorithm:** TransformUnscented EKI (`2p + 1` members for `p` parameters),
  with a diagonal noise covariance per variable scaled by the area weight of
  each column, so that the loss approximates an area-weighted global misfit.
  The standard deviation of each variable (`NOISE_STD`) is the square root of
  its variance in `NOISE_SCALARS` of `energy_fluxes.jl`.

## Data that are not artifacts

- `topographic_data.nc`, at the root of the repository: the statistics of the
  topography in 0.5° cells, made from ETOPO 2022 by `scf_topography_artifact.jl`.
  `TopoSnowCoverFractionModel` reads it from the working directory, so the
  calibration runs from the root of the repository.
- The monthly 0.5° MODIS snow cover fraction of ILAMB (2001-2021), read from
  `ILAMB_SCF_FILE`, or by default from `scf_0.5x0.5_ILAMB_filled.nc` at the
  root of the repository. The Buildkite step reads it from a group directory
  on the Caltech cluster.

## Running on Buildkite

Add the `calibrate ERA5` label (together with `Launch Buildkite`) to a pull
request, or start a manual build with `CALIBRATE_ERA5=""`. The step drives
the calibration on a CPU and submits every forward run as a Slurm job with one
H100 GPU (the neural snow depth model does not run on P100s), so the members of
an iteration run at once when enough GPUs are free. It uploads `era5_calibration/results` as Buildkite artifacts and
summarizes it in a build annotation:

- `era5_parameters.toml`: the calibrated parameters.
- `parameters.png`: the parameter ensemble and the loss over iterations.
- `rmse.png`: the area-weighted RMSE of each variable with the default and the
  calibrated parameters.
- `columns.png`: the columns, colored by how the calibration changed their
  error.
- `columns.csv`: the columns and their area weights.
- `summary.md`: the tables of the annotation.

To use the calibrated parameters, copy `era5_parameters.toml` to
`toml/era5_parameters.toml` and pass it to
`LP.create_toml_dict(FT; override_files = ["toml/era5_parameters.toml"])`.

## Running elsewhere

```bash
julia --project=.buildkite experiments/calibration/era5_columns/run_calibration.jl [OUTPUT_DIR]
```

The hourly 1° ERA5 forcing (`forty_yrs_era5_land_forcing_data`) is only
available on some clusters. With `ERA5_LOWRES_FORCING` set, the columns are
forced by the 8° ERA5 of 2008, repeated every year, which only tests the
pipeline. With `CALIBRATION_BACKEND=caltech`, every forward run is a Slurm job
on the GPU partition of the Caltech cluster (see `CALTECH_DIRECTIVES` in
`run_calibration.jl`); otherwise, the forward runs follow one another in the
driver process, on the device selected by `CLIMACOMMS_DEVICE`. Rerunning with
the same `OUTPUT_DIR` resumes an interrupted calibration (Caltech backend
only).

## Files

- `config.jl`: dates and number of columns, and the calibrated variables,
  noise, prior, and number of iterations of `energy_fluxes.jl`.
- `observations.jl`: reads the observations on the points of the ERA5 1° grid,
  computes seasonal means, and builds the `EKP.Observation`. The columns are
  ERA5 points, so they are forced by and compared to ERA5 at the same place.
  The other datasets are averaged over the four cells around each point: the
  four 0.5° cells of FLUXCOM and MODIS tile the 1° cell centered on the point.
- `columns.jl`: selects the columns. The candidates are the points inside land
  where all the calibrated variables are observed, with GOSIF GPP (from the
  `inversion_nee` artifact) in place of FLUXCOM GPP, which is missing in most
  deserts and ends in 2013. GOSIF does not cover the ice sheets, so they are
  excluded. Area-weighted k-means groups the candidates by their observed
  2002-2020 seasonal climatology, and each cluster contributes its most
  central point, weighted by the cluster's area. The driver writes them to
  `results/columns.csv`; to preview them, run
  `include("experiments/calibration/era5_columns/model_interface.jl"); write_columns("columns.csv", select_columns())`.
- `model_interface.jl`: the forward model and observation map.
- `analysis.jl`: the calibrated parameter file, metrics, figures, and summary.
- `run_calibration.jl`: the driver.
