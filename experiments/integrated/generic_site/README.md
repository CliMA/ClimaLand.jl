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
  table. It is too large to be downloadable: to run any other site, obtain the
  data yourself as described in the `fluxnet2015` entry of
  [ClimaArtifacts](https://github.com/CliMA/ClimaArtifacts/tree/main/fluxnet2015).
  `FluxnetSimulations.get_location` and `get_fluxtower_height` then read the
  site coordinates and tower height from the metadata table.

## Files

| File | Purpose |
| --- | --- |
| `run_generic_site.jl` | Runs one site and plots model vs observations to `out/<SITE_ID>/`. |
| `list_fluxnet_sites.jl` | Prints every site ID in the `fluxnet2015` artifact, one per line. |

## Run a site

```bash
julia --project=.buildkite experiments/integrated/generic_site/run_generic_site.jl US-MOz

# 30 days at a site that needs the `fluxnet2015` artifact:
SIM_DURATION_DAYS=30 julia --project=.buildkite \
    experiments/integrated/generic_site/run_generic_site.jl DE-Tha
```

The script builds each piece (domain, forcing, LAI, `LandModel`, initial
conditions, diagnostics) explicitly, so any of them can be swapped, for example
by passing `canopy = ...` to `LandModel` as in the
[ERA5 single-column tutorial](../../../docs/src/tutorials/integrated/snowy_land_era5_tutorial.jl).

## CI

The buildkite step `generic_site BE-Vie run` runs `run_generic_site.jl` at
BE-Vie, which is not in `fluxnet_sites`, so it exercises the `fluxnet2015`
metadata and data paths.
