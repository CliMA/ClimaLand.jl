# Default model at any FLUXNET2015 site

`run_generic_site.jl` runs the default ClimaLand `LandModel` at any FLUXNET2015
site: parameters are the spatially varying defaults at the site coordinates,
LAI is the MODIS climatology, the forcing is the site's observed meteorology,
and the initial conditions come from the site observations. Unlike
`experiments/integrated/fluxnet/run_fluxnet.jl`, it needs no site-specific
configuration.

## Data

ClimaLand can read FLUXNET data from two artifacts:

- `fluxnet_sites` is downloadable and contains one year of data at four sites
  (US-MOz, US-NR1, US-Ha1, US-Var). Their coordinates and tower heights are
  set in `ext/fluxnet_simulations/<site>.jl`.
- `fluxnet2015` is the full FLUXNET2015 dataset (200+ sites) and its metadata
  table. It is too large to be downloadable: to run any other site, download
  the data from fluxnet.org (free account required) and set it up as described
  in the `fluxnet2015` entry of
  [ClimaArtifacts](https://github.com/CliMA/ClimaArtifacts/tree/main/fluxnet2015).
  `FluxnetSimulations.get_location` and `get_fluxtower_height` then read the
  site coordinates and tower height from the metadata table.

## Files

| File | Purpose |
| --- | --- |
| `run_generic_site.jl` | Runs one site and plots model vs observations to `out/<SITE_ID>/`. |
| `list_fluxnet_sites.jl` | Prints every site ID in the `fluxnet2015` artifact, one per line. |
| `fluxnet_ilamb_rmse_sites.jl` | Runs sites in series and writes their monthly LE, H, SWup and LWup to `out/ilamb_rmse/sites/`. |
| `fluxnet_ilamb_rmse_plot.jl` | Computes the ILAMB RMSE from those files and plots it against the ILAMB land-hist models. |
| `ilamb_land_hist_site_rmse.jl` | Writes `ilamb_land_hist_site_rmse.csv`, the per-site RMSE of the land-hist models. |

## Run a site

```bash
julia --project=.buildkite experiments/integrated/generic_site/run_generic_site.jl US-MOz

# A site that needs the `fluxnet2015` artifact:
julia --project=.buildkite experiments/integrated/generic_site/run_generic_site.jl DE-Tha
```

The run lasts 7 days (`duration` in the script) from the first time step at
which every forcing variable is observed. The timestep and soil layers match the
global runs, and diagnostics are averaged at the resolution of the site data
(half-hourly or hourly) so they line up with the observations.

The script builds each piece (domain, forcing, LAI, `LandModel`, initial
conditions, diagnostics) explicitly, so any of them can be swapped, for example
by passing `canopy = ...` to `LandModel` as in the
[ERA5 single-column tutorial](../../../docs/src/tutorials/integrated/snowy_land_era5_tutorial.jl).

## RMSE against the ILAMB land-hist models

`boxplot_rmse_fluxnet2015.png` compares ClimaLand's FLUXNET2015 RMSE of LE, H,
SWup and LWup with the ILAMB land-hist models (CLM, ISBA-CTRIP and JSBACH,
each forced by CRUJRA, GSWP3 and Princeton;
[dashboard](https://www.ilamb.org/land-hist/)), computed the way ILAMB does
(`AnalysisMeanStateSites` in
[ILAMB](https://github.com/rubisco-sfa/ILAMB)): at each site, the RMSE of the
monthly means against the FLUXNET2015 monthly files (`LE_F_MDS`, `H_F_MDS`,
`SW_OUT`, `LW_OUT`, the columns of ILAMB's benchmark), then the mean over
sites. All models are averaged over the same sites: those where ClimaLand and
every land-hist model have a value.

- `fluxnet_ilamb_rmse_sites.jl` runs every site of
  `ilamb_land_hist_site_rmse.csv` in the `fluxnet2015` artifact, except the
  four sites of `fluxnet_sites`, which resolve to that single-year artifact.
  Each site spins up for one year from the first time step with all forcing
  variables, up to the next month boundary in local standard time, and is
  scored on the following 12 calendar months; sites with a shorter record are
  skipped. It also accepts a list of site IDs.
- `fluxnet_ilamb_rmse_plot.jl` writes the figure, `rmse_summary.csv` (mean over
  sites per model) and `rmse_per_site.csv` (ClimaLand) to `out/ilamb_rmse/`.

ClimaLand is forced by the tower meteorology and scored on one year, while the
land-hist models are global runs forced by gridded reanalyses, sampled at the
sites and scored over each site's whole record.

```bash
julia --project=.buildkite experiments/integrated/generic_site/fluxnet_ilamb_rmse_sites.jl DE-Tha BE-Vie
julia --project=.buildkite experiments/integrated/generic_site/fluxnet_ilamb_rmse_plot.jl
```

## CI

The buildkite step `generic_site BE-Vie run` runs `run_generic_site.jl` at
BE-Vie, which is not in `fluxnet_sites`, so it exercises the `fluxnet2015`
metadata and data paths.

The group `FLUXNET2015 RMSE vs ILAMB land-hist` runs `fluxnet_ilamb_rmse_sites.jl`
split over 24 parallel CPU jobs, then `fluxnet_ilamb_rmse_plot.jl`, and keeps
the figure and CSVs as artifacts.
