# Generic Single-Site Fluxnet Calibration

Per-site calibration of ClimaLand parameters against daily means of fluxnet
observations. Each site is an independent `TransformUnscented` calibration
with a small ensemble; sites are embarrassingly parallel and run as one
Slurm array task each on HPC.

The driver builds an `EKP.ObservationSeries` with up to `max_n_samples`
year-long sample windows per site. EKP draws `minibatch_size` windows per
iteration; each window is preceded by a `spinup_days` spinup so deep soil
state has time to settle before observations are scored. Sites with only one
year of data degenerate cleanly (single window, repeated each iteration).

## Files

| File | Purpose |
| --- | --- |
| `calibrate_site.jl` | One-site driver. Loads a config, builds an ObservationSeries, runs `ClimaCalibrate.calibrate` on `JuliaBackend`. Appends a row to `calibration_out/site_results.csv` when done. |
| `site_drivers.jl` | Aggregates fluxnet forcing + spatial soil/LAI/IGBP data at the site for the per-site CSV. |
| `calibrate_all_fluxnet.jl` | Loops over every fluxnet site (or a subset) sequentially in one Julia process. Used by CI. |
| `calibrate_all_fluxnet.sbatch` | Slurm array wrapper — one array task per site. **Production path on HPC.** |
| `list_fluxnet_sites.jl` | Prints every FLUXNET2015 site ID in the `fluxnet2015` artifact, one per line. Used to size the Slurm array. |
| `configs/default.jl` | Production GPP + LHF config: 365-day window, 5 EKP iterations, no spinup, no minibatching. |
| `configs/smoke_test.jl` | GPP + LHF, 90-day window, 2 iterations — buildkite CI smoke. |
| `configs/er.jl` | Production Reco/DAMM config: 4 DAMM kinetic parameters, 365-day windows, 90-day spinup, up to 5 yearly samples per site, `minibatch_size=1`. |
| `configs/er_smoke.jl` | Reco config trimmed for buildkite CI (1 sample, 30-day window, 15-day spinup, 1 iteration). |

## Configs

Configs live in `configs/` and define three names that `calibrate_site.jl`
loads via `include`:

- `CALIBRATE_CONFIG` — NamedTuple with `short_names`, `n_iterations`,
  `cal_window_days`, `spinup_days`, `minibatch_size`, `max_n_samples`, and
  `rng_seed`. To preserve today's "single window, no spinup, no minibatching"
  behavior, set `spinup_days = 0`, `max_n_samples = 1`, `minibatch_size = 1`.
- `NOISE_VARIANCES` — `Dict` mapping each entry of `short_names` to a per-bin
  scalar variance used to build the diagonal noise covariance. Variances are
  in the obs vector's units (see calibrate_site.jl `_obs_unit_factor` for the
  per-variable conversion — e.g. `er` is in `g C m⁻² day⁻¹`).
- `get_calibration_prior()` — returns an `EKP.combine_distributions(...)` prior.

Select a config at run time with `CALIBRATION_CONFIG=<filename> julia ...`.
Defaults to `default.jl`. To add a new experiment, copy an existing config,
tweak, and point `CALIBRATION_CONFIG` at the new file — no other code change.

`N_ITERS` env var overrides the config's value when set; useful for one-off
probes without editing files.

## Reco / DAMM calibration

`configs/er.jl` calibrates four DAMM kinetic parameters
(`soilCO2_pre_exponential_factor`, `soilCO2_activation_energy`,
`michaelis_constant`, `O2_michaelis_constant`) against the fluxnet
`RECO_NT_VUT_REF` column (g C m⁻² day⁻¹). The driver compares **daily-mean**
fluxnet RECO to **daily-mean** model `er` diagnostics on each sample window.

Per-site outputs land in
`calibration_out/<SITE_ID>_<cal_window_days>d/`. In addition, every site
appends a row to a single shared CSV at
`calibration_out/site_results.csv` containing the calibrated DAMM means and
SDs alongside per-site climate / soil drivers (mean TA / VPD / SW / TS,
annual cumulative precip, MODIS max LAI, soil porosity, IGBP class). Use it
to study what controls parameter variation across sites.

## Run all fluxnet sites for Reco (Slurm, central HPC)

The `fluxnet2015` artifact is HPC-only, so this only works on a machine where
it has been pre-staged.

```bash
tmux new -s reco_cal
module load climacommon/2026_02_18

# 1. Generate the site list once on the login node:
julia --project=.buildkite \
    experiments/integrated/generic_site/list_fluxnet_sites.jl \
    > experiments/integrated/generic_site/calibration_out/site_ids.txt

# 2. Submit the array job with the Reco config:
N=$(wc -l < experiments/integrated/generic_site/calibration_out/site_ids.txt)
sbatch --array=0-$((N-1))%50 \
    --export=ALL,CALIBRATION_CONFIG=er.jl \
    experiments/integrated/generic_site/calibrate_all_fluxnet.sbatch
```

Each site runs in its own array task (serial Julia process per task — no MPI
or GPU). The shared `site_results.csv` is written under a `mkpidlock` so
concurrent appends are safe.

## Run all fluxnet sites (Slurm, default GPP+LHF config)

```bash
module load climacommon/2026_02_18
julia --project=.buildkite \
    experiments/integrated/generic_site/list_fluxnet_sites.jl \
    > experiments/integrated/generic_site/calibration_out/site_ids.txt
N=$(wc -l < experiments/integrated/generic_site/calibration_out/site_ids.txt)
sbatch --array=0-$((N-1))%50 \
    experiments/integrated/generic_site/calibrate_all_fluxnet.sbatch
```

To pick a different config or override the iterations, use `--export`:

```bash
sbatch --array=0-$((N-1))%50 \
    --export=ALL,CALIBRATION_CONFIG=er.jl,N_ITERS=10 \
    experiments/integrated/generic_site/calibrate_all_fluxnet.sbatch
```

Per-site outputs land in
`experiments/integrated/generic_site/calibration_out/<SITE_ID>_<window>d/`,
slurm logs in `calibration_out/slurm_logs/`, the shared per-site row in
`calibration_out/site_results.csv`.

## Run all fluxnet sites (single Julia process)

For smoke tests or a small subset, run sequentially on a login node /
workstation. Slower but simpler — no scheduler involvement.

```bash
# All sites, default config:
julia --project=.buildkite \
    experiments/integrated/generic_site/calibrate_all_fluxnet.jl

# A 2-site subset with the smoke-test config:
SITES="AU-Cum,BE-Vie" CALIBRATION_CONFIG=smoke_test.jl \
    julia --project=.buildkite \
    experiments/integrated/generic_site/calibrate_all_fluxnet.jl
```

## Run a single site

```bash
SITE_ID=US-MOz \
    julia --color=yes --project=.buildkite \
    experiments/integrated/generic_site/calibrate_site.jl
```

The four bundled fluxnet sites (`US-MOz`, `US-Var`, `US-NR1`, `US-Ha1`) work
with default `LOCAL=true`. For any other FLUXNET2015 site, set `LOCAL=false` so
coordinates are auto-resolved from the FLUXNET2015 metadata table.

## CI smoke tests (buildkite)

Five steps in `.buildkite/pipeline.yml` exercise this pipeline on every
buildkite run:

- `generic_site_calibrate_er AU-Cum` — Reco / DAMM smoke at AU-Cum.
- `generic_site_calibrate_er BE-Vie` — Reco / DAMM smoke at BE-Vie. Both
  publish `site_results.csv` as a buildkite artifact for inspection.
- `generic_site_calibrate_gpp_lhf AU-Cum` — original GPP+LHF smoke retained.
- `generic_site_calibrate_gpp_lhf BE-Vie` — same, for BE-Vie.
- `generic_site_calibrate_er_all_fluxnet (smoke, 2 sites)` — exercises
  `calibrate_all_fluxnet.jl` with `SITES="AU-Cum,BE-Vie"` and the Reco config.
  This is the same driver `calibrate_all_fluxnet.sbatch` invokes per array
  task, so a green step here means the all-sites HPC submission is ready to
  fan out. The shared `site_results.csv` is published as an artifact and
  contains one row per site, which is the form most useful for cross-site
  analysis.

All steps are `soft_fail: true` while the pipeline beds in.
