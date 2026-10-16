# Calibration on ERA5 columns

This experiment calibrates the land model of the default long run (soil,
canopy with the P-model, prescribed MODIS LAI, and spatially varying canopy
height, snow, slab lake, soil CO2) forced by ERA5. Instead of the global grid,
it runs a `ColumnEnsemble` of 200 columns chosen to represent global land,
which makes each forward run cheap.

- **Observations:** seasonal means from March 2008 to February 2010 of ERA5
  latent and sensible heat fluxes and upwelling shortwave and longwave
  radiation, and GOSIF GPP (from the `inversion_nee` artifact).
- **Simulations:** March 2007 to March 2010, with the first year as spinup.
- **Parameters:** 16 parameters of the P-model, moisture stress, stem and root
  area indices, canopy turbulence and heat capacity, turbulence below the
  canopy, and radiation (see `get_calibration_prior` in `config.jl`). A single
  factor scales the PAR and NIR reflectances of stems.
- **Algorithm:** 5 iterations of TransformUnscented EKI (33 members), with a
  diagonal noise covariance per variable (`NOISE_STD`) scaled by the area
  weight of each column, so that the loss approximates an area-weighted global
  misfit.

## Running on Buildkite

Add the `calibrate ERA5` label (together with `Launch Buildkite`) to a pull
request, or start a manual build with `CALIBRATE_ERA5=""`. The step drives
the calibration on a CPU and submits every forward run as a Slurm job with one
P100 GPU, so all the members of an iteration run at once when enough GPUs are
free. It uploads `era5_calibration/results` as Buildkite artifacts and
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

- `config.jl`: dates, calibrated variables, noise, and prior.
- `observations.jl`: reads ERA5 and GOSIF on the points of the ERA5 1° grid,
  computes seasonal means, and builds the `EKP.Observation`. The columns are
  ERA5 points, so they are forced by and compared to ERA5 at the same place;
  GOSIF, whose 1° cells are centered half a degree away, is averaged over the
  four cells around each point.
- `columns.jl`: selects the columns. The candidates are the points inside land
  where all five variables are observed (which excludes the ice sheets).
  Area-weighted k-means groups them by their observed 2002-2020 seasonal
  climatology, and each cluster contributes its most central point, weighted by
  the cluster's area. The driver writes them to `results/columns.csv`; to
  preview them, run
  `include("experiments/calibration/era5_columns/model_interface.jl"); write_columns("columns.csv", select_columns())`.
- `model_interface.jl`: the forward model and observation map.
- `analysis.jl`: the calibrated parameter file, metrics, figures, and summary.
- `run_calibration.jl`: the driver.
