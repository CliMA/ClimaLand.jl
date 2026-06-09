# NEON calibration pipeline

A config-driven pipeline for the NEON soil-CO₂ DAMM calibration. One TOML config
drives one driver, which loops over runs and writes a single joint results CSV.

It replaces the old hand-edited workflow (`set_Station.jl` comment-toggling,
priors hard-coded in `run_calibration*.jl`, the `start_calibration_pipeline_mult`
scripts) with: **config → driver → CSV**.

## Quick start

```bash
cd /home/evametz/Github/ClimaLand/ClimaLand.jl

# run everything in a config (single run or a batch):
julia --project=.buildkite experiments/calibrate_neon_pipeline/run_pipeline.jl \
    experiments/calibrate_neon_pipeline/config/cper_allyears.toml

# run just one (site, year) out of a batch config:
julia --project=.buildkite experiments/calibrate_neon_pipeline/run_pipeline.jl \
    experiments/calibrate_neon_pipeline/config/cper_allyears.toml --site NEON-cper --year 2021
```

Driver options: `--site S`, `--year Y`, `--run-identifier ID`, `--stop-on-error`.

## What a run does

Per `[[runs]]` entry the driver runs the enabled `[steps]`:

1. **generate_observations** — builds `observations.jld2` for the site/period.
2. **calibrate** — UKI calibration; writes `prior_values.txt`,
   `final_parameter_means.txt`, the `iteration_*/eki_file.jld2`, and a
   `priors.toml` provenance copy.
3. **plot_eki_diagnostics** — timeseries / RMSE / 1-to-1 figures.
4. **run_prior_mean** — a forward run at the *optimized* parameters, with its own
   diagnostic figures (SWC, soil T, O₂, SOC profile, CO₂ budget…).

Then it appends one row (priors + posteriors) to the master CSV.

## The config (TOML)

See [`config/default.toml`](config/default.toml) (template) and
[`config/cper_allyears.toml`](config/cper_allyears.toml) (full-featured example).

- `[settings]` — spinup, iterations, `cal_depth` (0.02 or 0.06 m), `settingsdesc`,
  `dt`, `output_root`, `results_csv`, optional `run_identifier`.
- `[steps]` — toggle each step on/off.
- `[priors.<name>]` — `mean, std, lower, upper` for each calibrated parameter.
  **The parameters present = the parameters calibrated.** Bounds may be the
  strings `"Inf"`/`"-Inf"`.
- `[[runs]]` — one calibration each. Use `start`/`stop`, or `years = [...]` sugar
  (expands to full-year ranges). Per-run `[runs.priors.<name>]` overrides only
  the fields you set. Optional `restart_from` and `settingsdesc` per run.

### Calibrating the labile parameter (or not)

`labile_depth_scale` (k, 1/m) scales SOC by `exp(k·z)`. It is the optional 5th
parameter:

- **Present** in priors → calibrated (5 params).
- **Absent** → not calibrated; the model runs with `k = 0` (`exp(0·z) = 1`), i.e.
  the plain non-labile behaviour. The pipeline always uses the wLabile model
  interface and pins `labile_depth_scale = 0.0` when it isn't calibrated, so a
  single code path covers both modes.

### Restarting from a previous posterior

A run with `restart_from = "<run_identifier>"` seeds its prior **means** from that
previous run's posterior (`final_parameter_means.txt`), keeping std/bounds from
the current config. Parameters absent from the previous posterior keep their
config prior — so you can *add* a parameter (e.g. switch on labile) on restart.
The previous run is located via the CSV (fast path) or a scan of `output_root`
for `output_<id>`.

## Output layout

```
<output_root>/<site>/<site>_<start>_<stop>/SpinUP-<spinup>d/CalDepth-<depthM>/<n_iter>-It/<settingsdesc>/
    observations.jld2                 # shared by the run(s) at this base path
    output_<run_identifier>/          # the run folder (replaces output_N)
        prior_values.txt
        priors.toml                   # resolved priors (provenance)
        final_parameter_means.txt
        iteration_*/eki_file.jld2
        figures_eki_diagnostics/
        prior_mean_optimized/         # forward-run figures
        model_scripts/                # snapshot of the model code used
```

The master CSV (`<output_root>/<results_csv>`) has one row per run with metadata,
`run_identifier`, `restart_from`, all prior columns (4 fields × 5 params), and the
posterior columns. Uncalibrated parameters have `missing` priors; labile's
posterior cell is `0.0` when it is off.

## Design

Every step is a **pure function** that takes the `RunConfig` (and prior step
results) and RETURNS its results. There is **no ENV interface between steps** and
no top-level `const`s, so a whole batch runs safely in **one Julia session** — no
subprocess-per-run needed.

The science is rewritten as self-contained function copies in `pipeline/` (the
originals in `../calibrate_neon` are left untouched as a reference):

| `pipeline/` file | function | from original |
|---|---|---|
| `Observations.jl` | `generate_observations(run; base_dir)` | `generate_observations.jl` |
| `Calibration.jl` | `run_calibration(run; obs_filepath, output_dir)` | `run_calibration_wLabile.jl` |
| `Diagnostics.jl` | `plot_eki_diagnostics(run; output_dir, eki_path, obs_filepath)` | `plot_eki_diagnostics.jl` |
| `ForwardRun.jl` | `forward_run(run; output_dir, params)` | `run_prior_mean_wLabile.jl` |

`Config.jl` / `ResultsTable.jl` / `Pipeline.jl` are the config / CSV /
orchestration layer. `run_pipeline.jl` loads `Pipeline.jl`, builds + filters the
config, and calls `run_pipeline(cfg)`.

The driver threads results between steps (no ENV):

```julia
obs   = generate_observations(run; base_dir)
calib = run_calibration(run; obs_filepath = obs, output_dir)     # -> final_params, eki_path
diag  = plot_eki_diagnostics(run; output_dir, eki_path = calib.eki_path, obs_filepath = obs)
forward_run(run; output_dir = .../prior_mean_optimized, params = calib.final_params)
```

### Calibration workers

The calibration step still uses `Distributed` workers (inherent to
ClimaCalibrate — one ensemble member per worker). `Calibration.jl` broadcasts the
config to workers and `include`s the **original** `model_interface_wLabile.jl`
plus `labile_pin_shim.jl` (which makes one interface serve both labile modes by
injecting `labile_depth_scale = 0.0` when it isn't calibrated). The model
interface reads a few module globals on each worker — that is the model
interface's own contract, internal to the calibration step, not a pipeline-level
ENV interface.

The master CSV is appended under a file lock, so it is also safe if you ever run
several pipelines concurrently.

See `PIPELINE_PLAN.md` in `../calibrate_neon` for the full design rationale.
