# Calibration on FLUXNET sites

This experiment calibrates the land model of the default long run (soil,
canopy with the P-model and prescribed MODIS LAI, snow, slab lake, soil CO2)
on the FLUXNET2015 sites. The sites run
together as a `ColumnEnsemble`, one column per site, and each column is forced
by the meteorology measured at its tower.

- **Observations:** monthly means of the observed latent and sensible heat
  fluxes (`LE_CORR`, `H_CORR`), upwelling shortwave and longwave radiation
  (`SW_OUT`, `LW_OUT`), and GPP (`GPP_DT_VUT_REF`). A month is compared if at
  least `MIN_GOOD_FRACTION` of its half-hours are measured or good quality
  gap-filled (for the gap-filled fluxes) and `MIN_COVERAGE` have a value.
- **Simulations:** four years, the first of which is spinup; the next three
  (`CALIBRATION_YEARS`) are compared to observations.
- **Parameters:** 18 parameters of the P-model, moisture stress, stem and root
  area indices, canopy turbulence and heat capacity, turbulence below the
  canopy, bare soil, and radiation (see `get_calibration_prior` in
  `config.jl`). A single factor scales the PAR and NIR reflectances of stems.
- **Algorithm:** 5 iterations of TransformUnscented EKI (37 members), with a
  diagonal noise covariance per variable (`NOISE_STD`). Every site-month has the
  same weight.

## One calendar for all sites

A `ColumnEnsemble` has one clock, but the FLUXNET2015 records cover different
years, between 1991 and 2014. Each site contributes up to three consecutive
calendar years in which every forcing variable but CO2 is observed or
gap-filled for at least `MIN_FORCING_COVERAGE` of the half-hours, choosing the
years with the most observed monthly means. These years are simulated as the
last of `CALIBRATION_YEARS` (by month, day, and time of day, with February 29
mapped to February 28 of a common year), and the earlier years, starting with
the spinup year, cycle through them. Only the site's own years are compared to
observations. A site with one complete year thus repeats it four times, and
only the last repetition is compared.

The forcing is linearly interpolated in time to half-hours in UTC (FLUXNET
timestamps are in local standard time), and precipitation is split into rain
and snow as in `FluxnetSimulations.prescribed_forcing_fluxnet`. A site that
measures CO2 for less than half of its years is forced by the annual mean at
Mauna Loa instead. Since the years of a site differ from the simulated years
and FLUXNET2015 mostly predates MODIS, the LAI is the MODIS climatology.

The canopy height of each column is that of its site: the height in the
site's configuration in `FluxnetSimulations` for US-MOz, US-NR1, US-Ha1 and
US-Var, or else in the FLUXNET2015 metadata. Where a site has none, it is the
CLM canopy height of the default model, which, at the scale of a grid cell,
can be far from the vegetation around a tower (0.5 m at the US-MOz forest).
The forcing height of each column is the height of its tower above the ground,
or `MIN_CANOPY_CLEARANCE` above the canopy where the tower is not above it.

## Running on Buildkite

Add the `calibrate fluxnet` label (together with `Launch Buildkite`) to a pull
request, or start a manual build with `CALIBRATE_FLUXNET=""`. The step reads
the sites of the `fluxnet2015` artifact, drives the calibration on a CPU, and
submits every forward run as a Slurm job with one P100 GPU, so all the members
of an iteration run at once when enough GPUs are free. It uploads
`fluxnet_calibration/results` as Buildkite artifacts and summarizes it in a
build annotation:

- `fluxnet_parameters.toml`: the calibrated parameters.
- `parameters.png`: the parameter ensemble and the loss over iterations.
- `rmse.png`: the RMSE of each variable with the default and the calibrated
  parameters.
- `seasonal_cycle.png`: the mean seasonal cycle over the sites of each
  hemisphere, observed and simulated.
- `sites.png`: the sites, colored by how the calibration changed their error.
- `sites.csv`: the sites, their tower and canopy heights (`NaN` where the
  CLM height is used), and their simulated years.
- `summary.md`: the tables of the annotation.

To use the calibrated parameters, copy `fluxnet_parameters.toml` to
`toml/fluxnet_parameters.toml` and pass it to
`LP.create_toml_dict(FT; override_files = ["toml/fluxnet_parameters.toml"])`.

## Running elsewhere

```bash
julia --project=.buildkite experiments/calibration/fluxnet_sites/run_calibration.jl [OUTPUT_DIR]
```

The `fluxnet2015` artifact is too large to be downloadable (see the
`fluxnet2015` entry of
[ClimaArtifacts](https://github.com/CliMA/ClimaArtifacts/tree/main/fluxnet2015)).
To calibrate on a few sites instead, list them in `FLUXNET_SITES`, for example
the four sites of the downloadable `fluxnet_sites` artifact:

```bash
FLUXNET_SITES=US-MOz,US-NR1,US-Ha1,US-Var julia --project=.buildkite experiments/calibration/fluxnet_sites/run_calibration.jl
```

With `CALIBRATION_BACKEND=caltech`, every forward run is a Slurm job on the GPU
partition of the Caltech cluster (see `CALTECH_DIRECTIVES` in
`run_calibration.jl`); otherwise, the forward runs follow one another in the
driver process, on the device selected by `CLIMACOMMS_DEVICE`. Rerunning with
the same `OUTPUT_DIR` resumes an interrupted calibration (Caltech backend
only), and reuses the site data saved in `OUTPUT_DIR/site_data.jld2`.

## Files

- `config.jl`: dates, calibrated variables, data selection, noise, and prior.
- `sites.jl`: reads the FLUXNET files, selects the years of each site, and
  saves the forcing on the simulated calendar and the observed monthly means
  (`site_data.jld2`), and the table of the sites (`results/sites.csv`).
- `model_interface.jl`: the forward model and observation map.
- `analysis.jl`: the calibrated parameter file, metrics, figures, and summary.
- `run_calibration.jl`: the driver.
