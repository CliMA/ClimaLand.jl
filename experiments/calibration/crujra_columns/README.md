# Calibration on CRUJRA columns

This experiment calibrates the land model of the default long run (soil,
canopy with the P-model and prescribed MODIS LAI, snow, slab lake, soil CO2)
forced by CRUJRA. Instead of the global grid, it runs a `ColumnEnsemble` of 200
columns chosen to represent global land, which makes each forward run cheap.

- **Observations:** seasonal means from March 2008 to February 2010 of FLUXCOM
  latent heat flux, sensible heat flux, and GPP, and CERES EBAF ed4.2 upwelling
  shortwave and longwave radiation, in the 1° cell of each column.
- **Simulations:** March 2007 to March 2010, with the first year as spinup.
- **Parameters:** 10 P-model, moisture stress, canopy turbulence, and longwave
  parameters (see `get_calibration_prior` in `config.jl`).
- **Algorithm:** 5 iterations of TransformUnscented EKI (21 members), with a
  diagonal noise covariance per variable (`NOISE_STD`) scaled by the area
  weight of each column, so that the loss approximates an area-weighted global
  misfit.

## Running on Buildkite

Add the `calibrate CRUJRA` label (together with `Launch Buildkite`) to a pull
request, or start a manual build with `CALIBRATE_CRUJRA=""`. The step uploads
`crujra_calibration/results` as Buildkite artifacts and summarizes it in a
build annotation:

- `crujra_parameters.toml`: the calibrated parameters. Copy it to
  `toml/crujra_parameters.toml` and use it with
  `LP.create_toml_dict(FT; override_files = ["toml/crujra_parameters.toml"])`.
- `parameters.png`: the parameter ensemble and the loss over iterations.
- `rmse.png`: the area-weighted RMSE of each variable with the default and the
  calibrated parameters.
- `columns.png`: the columns, colored by how the calibration changed their
  error.
- `summary.md`: the tables of the annotation.

## Running elsewhere

```bash
CLIMACOMMS_DEVICE=CUDA CALIBRATION_N_WORKERS=7 \
    julia --project=.buildkite experiments/calibration/crujra_columns/run_calibration.jl [OUTPUT_DIR]
```

CRUJRA is only available on the CliMA clusters. With
`CALIBRATION_N_WORKERS = 0` (the default), the ensemble members run one after
another in a single process. Otherwise, they run concurrently on that many
local worker processes on the device selected by `CLIMACOMMS_DEVICE`. When
`CUDA_VISIBLE_DEVICES` lists several GPUs, the workers are spread over them, one
GPU per worker in turn.

## Files

- `config.jl`: dates, calibrated variables, noise, and prior.
- `select_columns.jl`: selects the columns and writes `columns.csv`. The
  candidates are 1° land cells, in two strata that get columns in proportion
  to their area: cells where all five variables are observed, and interior
  land cells equatorward of 60° without FLUXCOM coverage (mostly deserts),
  which only constrain SWU and LWU. Within each stratum, area-weighted k-means
  groups the cells by their observed 2001-2013 seasonal climatology, and each
  cluster contributes its most central cell, weighted by the cluster's area.
- `observations.jl`: reads the observations at the columns and builds the
  `EKP.Observation`.
- `model_interface.jl`: the forward model and observation map.
- `analysis.jl`: the calibrated parameter file, metrics, figures, and summary.
- `run_calibration.jl`: the driver.
- `tmp_artifacts/ilamb_fluxcom_energy`: a 2008-2010 subset of the FLUXCOM
  latent and sensible heat fluxes, used until the `ilamb_fluxcom_energy`
  artifact (see `ClimaArtifacts/ilamb_fluxcom_energy`) is available.
