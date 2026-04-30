# Generic Single-Site Fluxnet Calibration

Per-site calibration of ClimaLand parameters against monthly means of fluxnet
observations (default: GPP and LHF). Each site is an independent
`TransformUnscented` calibration with a small ensemble; sites are embarrassingly
parallel and run as one Slurm array task each on HPC.

This pipeline is the single-site cousin of `experiments/calibration/`: same
priors style, but no minibatching, no global model, no MPI/GPU.

## Files

| File | Purpose |
| --- | --- |
| `calibrate_site.jl` | One-site driver. Loads a config, builds observations, runs `ClimaCalibrate.calibrate` on `JuliaBackend`. |
| `calibrate_all_fluxnet.jl` | Loops over every fluxnet site (or a subset) sequentially in one Julia process. Used by CI. |
| `calibrate_all_fluxnet.sbatch` | Slurm array wrapper — one array task per site. **Production path on HPC.** |
| `list_fluxnet_sites.jl` | Prints every FLUXNET2015 site ID in the `fluxnet2015` artifact, one per line. Used to size the Slurm array. |
| `configs/default.jl` | Production config: GPP + LHF, 365-day window, 5 EKP iterations. |
| `configs/smoke_test.jl` | 90-day window, 2 iterations — used by the buildkite CI smoke tests. |

## Configs

Configs live in `configs/` and define three names that `calibrate_site.jl`
loads via `include`:

- `CALIBRATE_CONFIG` — NamedTuple with `short_names`, `n_iterations`,
  `cal_duration_days`, and `rng_seed`.
- `NOISE_VARIANCES` — `Dict` mapping each entry of `short_names` to a per-bin
  scalar variance used to build the diagonal noise covariance.
- `get_calibration_prior()` — returns an `EKP.combine_distributions(...)` prior.

Select a config at run time with `CALIBRATION_CONFIG=<filename> julia ...`.
Defaults to `default.jl`. To add a new experiment, copy `default.jl`, tweak,
and point `CALIBRATION_CONFIG` at the new file — no other code change.

`CAL_DURATION_DAYS` and `N_ITERS` env vars override the config's values when
set; useful for one-off probes without editing files.

## Run all fluxnet sites (Slurm, central HPC)

The `fluxnet2015` artifact is HPC-only (no download URL), so this only works
on a machine where it has been pre-staged.

```bash
module load climacommon/2026_02_18

# 1. Generate the site list once on the login node:
julia --project=.buildkite \
    experiments/integrated/generic_site/list_fluxnet_sites.jl \
    > experiments/integrated/generic_site/calibration_out/site_ids.txt

# 2. Submit the array job. `%50` caps concurrent tasks at 50 — tune to fairshare.
N=$(wc -l < experiments/integrated/generic_site/calibration_out/site_ids.txt)
sbatch --array=0-$((N-1))%50 \
    experiments/integrated/generic_site/calibrate_all_fluxnet.sbatch
```

To pick a different config or override the window, use `--export`:

```bash
sbatch --array=0-$((N-1))%50 \
    --export=ALL,CALIBRATION_CONFIG=my_experiment.jl,CAL_DURATION_DAYS=730 \
    experiments/integrated/generic_site/calibrate_all_fluxnet.sbatch
```

Per-site outputs land in
`experiments/integrated/generic_site/calibration_out/<SITE_ID>_<DAYS>d/`,
slurm logs in `calibration_out/slurm_logs/`.

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

Three steps in `.buildkite/pipeline.yml` exercise this pipeline on every
buildkite run:

- `generic_site_calibrate AU-Cum` — single-site path, `smoke_test.jl` config.
- `generic_site_calibrate BE-Vie` — single-site path, southern + northern coverage.
- `generic_site_calibrate_all_fluxnet (smoke, 2 sites)` — exercises
  `calibrate_all_fluxnet.jl` with `SITES="AU-Cum,BE-Vie"`. This is the same
  driver `calibrate_all_fluxnet.sbatch` invokes per array task, so a green
  step here means the all-sites HPC submission is ready to fan out.

All three are `soft_fail: true` while the pipeline beds in.
