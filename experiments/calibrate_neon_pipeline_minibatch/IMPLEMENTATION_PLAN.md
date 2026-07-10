# `calibrate_neon_pipeline_minibatch` — Implementation Plan

Status: **planning / not yet implemented.** This document captures the design so
we can resume work later. Written 2026-06-25; revised 2026-07-09 to base on the
`AllDepth` pipeline, keep the `Unscented(α_reg=1.0)` estimator, and default to a
20-day (adaptable) spinup.

## Goal

Set up a new folder `experiments/calibrate_neon_pipeline_minibatch/` that
calibrates the NEON soil-CO₂ (DAMM) parameters for **one site using all of its
years jointly**, instead of the current pipeline's one-calibration-per-year
approach. The mechanism mirrors ClimaLand PR
[#1732](https://github.com/CliMA/ClimaLand.jl/pull/1732)
(`experiments/integrated/generic_site/`): an EKP `ObservationSeries` of yearly
windows + a minibatcher.

## Reference points

- **Base folder (copy from this):**
  `experiments/calibrate_neon_pipeline_AllDepth/`
  - Entry: `run_pipeline.jl`
  - `pipeline/Config.jl`, `Observations.jl`, `Observation_flag.jl` (the cap-mask
    variant actually wired in by `Pipeline.jl`), `Calibration.jl`, `ForwardRun.jl`,
    `Diagnostics.jl`, `ResultsTable.jl`, `Pipeline.jl`, `labile_pin_shim.jl`
  - **Self-contained `src/`** (this is a key AllDepth difference from the original
    single-depth pipeline — the model code is copied here, NOT `include`d from
    `../calibrate_neon`): `src/model_interface_wLabile.jl`,
    `src/site_metadata.jl`, `src/neon_depth_lookup.jl`.
- **PR #1732 minibatch reference:** `experiments/integrated/generic_site/calibrate_site.jl`
  and `configs/er.jl` (the DAMM-kinetics config, closest analog to NEON soil-CO₂).
  **Note (2026-07-09):** this folder is NOT present in the current local checkout.
  Treat the summary below as the working reference; pull the actual PR files
  before implementing the forward-model minibatch loop if any detail is unclear.

---

## How PR #1732 does it (summary)

The core is `ObservationSeries` + minibatching, all in `calibrate_site.jl`:

1. `_sample_windows()` chops the site's fluxnet record into up to
   `max_n_samples` non-overlapping `cal_window_days` (365) windows, each
   preceded by a `spinup_days` spinup.
2. `_build_window_observation()` → one `EKP.Observation` per window
   (daily-mean target vector + diagonal noise covariance).
3. `generate_observation_series()` wraps the per-year `Observation`s with
   `ClimaCalibrate.minibatcher_over_samples(length(obs_vec), MINIBATCH_SIZE)`
   into one `EKP.ObservationSeries`. Saves `samples.jld2` (the
   `(window_start, window_stop)` list) so the forward model maps minibatch
   index → dates.
4. Each EKP **iteration** draws `minibatch_size` windows
   (`EKP.get_current_minibatch(ekp)`). With `minibatch_size = 1`, each iteration
   runs one year, cycling through years across iterations.
5. `forward_model` reads the current minibatch; for each window index runs a
   fresh `spinup + 365 d` simulation, drops the spinup prefix, extracts daily
   diagnostics, and **concatenates** them into the member's `g` vector (layout
   matches the concatenated observation for that minibatch).
6. Estimator (in PR): `EKP.TransformUnscented(prior; impose_prior = true)`.
   **This plan intentionally diverges — see the estimator note below.**
7. Scale-out: **sites in parallel** as a Slurm job array (one array task per
   site); ensemble members within a site run on `JuliaBackend`.

`er.jl` config: 4 DAMM params, `cal_window_days=365`, `spinup_days=90`,
`max_n_samples=5`, `minibatch_size=1`, `n_iterations=5`, `rng_seed=42`.

---

## Requested design choices for THIS pipeline (2026-07-09)

These are decided, not open:

1. **Base on `calibrate_neon_pipeline_AllDepth`** (multi-depth, self-contained
   `src/`, cap-mask obs), not the original single-depth `calibrate_neon_pipeline`.
2. **Estimator: keep `EKP.Unscented(prior; α_reg = 1.0, update_freq = 1)`** —
   exactly as `Calibration.jl:145-149` does today. Do **not** switch to
   `TransformUnscented`. Minibatching is a property of the obs + minibatcher, not
   the estimator, so `Unscented` works unchanged with an `ObservationSeries`.
3. **Spinup: default 20 days, kept adaptable.** `spinup_days` stays a per-run /
   per-settings TOML knob (as it already is in `Config.jl`), just with a default
   of 20. `run_from_env` already defaults to 20 (`Config.jl:357`); mirror that as
   the `[settings]` default so no config change is needed to get 20-day spinup,
   and any TOML can override it.

---

## Differences from the AllDepth pipeline (what this plan changes)

| | AllDepth `calibrate_neon_pipeline_AllDepth` | New `..._minibatch` (this plan) |
|---|---|---|
| Unit of calibration | one **year** per `RunConfig` (`years=[…]` fans out) | one **site**, all years jointly (one `RunConfig`) |
| Multi-year handling | independent runs, separate posteriors | one `ObservationSeries`, one posterior |
| Years enter EKP via | — (separate runs) | **minibatcher** (draws `minibatch_size`/iter) |
| Observation object | single stacked-over-depth `y_obs` + `noise_cov` | vector of per-window `EKP.Observation`s (each itself stacked over depth) |
| Estimator | `EKP.Unscented(α_reg=1.0, update_freq=1)` | **same** — `EKP.Unscented(α_reg=1.0)` (kept, per request) |
| Forward model per member | one sim over fixed `NEON_START/STOP_DATE` ENV | loop over current minibatch windows, concat `g` |
| Parallelism | ensemble members across workers (`WorkerBackend`) | keep `WorkerBackend` (orthogonal to minibatch) |
| Spinup | `spinup_days` before window, then trimmed off obs | same semantics; **default 20 d**, still adaptable |
| Depth handling | per-depth concat over `cal_depth_codes` (501/502/503) | unchanged — nests INSIDE each window's observation |

The multi-depth stacking is the AllDepth-specific complication the original plan
did not account for: each window's observation is itself a per-depth concatenated
vector. So the member `g` for one minibatch is
`[window₁: depth501; depth502; depth503]; [window₂: …]; …`.

---

## Version check (re-verified 2026-07-09)

`.buildkite` environment of `ClimaLand.jl`:
- **EnsembleKalmanProcesses 2.7.1** (installed at
  `~/.julia/packages/EnsembleKalmanProcesses/1ArOR`, `Project.toml` version 2.7.1).
  `ObservationSeries`, `FixedMinibatcher`, `get_current_minibatch`,
  `update_minibatch!`, `get_obs`/`get_minibatch` all present
  (`src/Observations.jl`, `src/EnsembleKalmanProcess.jl`). `Unscented` builds
  fine on an EKP constructed from an `ObservationSeries`.
- **ClimaCalibrate** (installed at `~/.julia/packages/ClimaCalibrate/V6O6K`) has
  `minibatcher_over_samples(n_samples, batch_size)` and
  `observation_series_from_samples(samples, batch_size, names)`
  (`src/ekp_utils.jl`). `observation_series_from_samples` is the convenient
  one-call builder; it warns and drops the remainder when `n_samples` is not
  divisible by `batch_size`.

So keeping `Unscented` requires **no estimator-side change** — only the obs
object handed to `EnsembleKalmanProcess` changes from `(y_obs, noise_cov)` to an
`ObservationSeries`.

---

## What gets reused vs. changed

### Reused essentially unchanged
- `src/model_interface_wLabile.jl`, `src/site_metadata.jl`,
  `src/neon_depth_lookup.jl` — heavy model + labile + per-depth layer lookup.
  **Copied into this folder's `src/`** (as AllDepth already does). The per-depth
  layer-selection logic in `forward_model` (`argmin` against the live grid) and
  the per-depth `observation_map` alignment stay as-is; minibatching wraps around
  them.
- `labile_pin_shim.jl` — labile pinning wrapper; `NeonPipelineInterface` dispatch
  tag. The forward-model minibatch loop is added here (see change 3).
- `ResultsTable.jl`, `Diagnostics.jl`, `ForwardRun.jl` (optimized/prior-mean run),
  provenance snapshotting — minor path/metadata tweaks only.
- The cap mask in `Observation_flag.jl` (spring-peak `is_capped_day`, union over
  501/502/503, all-years 25%-quantile cold-soil baseline) — reused per window
  (see the open question about per-year vs. all-years baseline).
- TOML-driven `Config` design, variable-length priors, depth handling.

### Changed

**1. `Config.jl` — `RunConfig` carries multiple windows**
- Today `_date_ranges` *fans out* `years=[…]` into separate `RunConfig`s
  (`Config.jl:169-176`, `327-333`). Now **one `RunConfig` holds the full
  `Vector{(start,stop)}`** for the site (add a `windows::Vector{Tuple{Date,Date}}`
  field; `start_date`/`stop_date` become the span min/max for naming).
- New fields: `minibatch_size`, `max_n_samples`, `rng_seed`. Keep interpreting
  `spinup_days` as "spinup *before each window*" (already the AllDepth semantics —
  `Observation_flag.jl:210-211` computes `spinup_date = start + Day(spinup_days)`).
- New `[settings]`/`[runs]` TOML keys: `minibatch_size` (default 1),
  `max_n_samples`, `rng_seed`. **`spinup_days` default = 20** (mirror
  `run_from_env`'s existing 20-day default at `Config.jl:357`); still overridable
  per-run/per-settings.
- `cal_depth_codes` handling is unchanged (still the ordered 501/502/503 list).

**2. `Observation_flag.jl` — build an `ObservationSeries`**
- Refactor the stacked-over-depth builder into a **per-window** function
  `build_window_observation(run, window_start, window_stop) -> EKP.Observation`
  that runs the existing per-depth loop (`daily_sco2_for_depth` +
  `capped`), scoped to one window, and returns one stacked `EKP.Observation`
  (this is basically the body of today's `generate_observations`,
  `Observation_flag.jl:236-272`, parameterized by window).
- New `generate_observation_series(run)` loops the run's windows, builds one
  `Observation` each (dropping windows that fail the per-depth data-coverage
  check instead of `error`ing — the current code `error`s on an empty depth at
  `Observation_flag.jl:249`; under minibatching a thin window should be skipped,
  not fatal), then calls
  `ClimaCalibrate.observation_series_from_samples(obs_vec, run.minibatch_size,
  window_names)` (or `EKP.ObservationSeries(obs_vec, minibatcher_over_samples(…))`).
- Save `observation_series.jld2` **and** `samples.jld2` — the ordered
  `[(window_start, window_stop), …]` list AND the per-window, per-depth
  `obs_dates_per_code` / `depth_codes` / `z_obs_m` metadata the forward model +
  observation_map need to reconstruct the stacking order for any minibatch.
- The cap mask stays; open question (below) is whether the cold-soil baseline is
  computed once over all years (current) or per window.

**3. `Calibration.jl` + `labile_pin_shim.jl` + `model_interface_wLabile.jl` —
series + minibatch-aware forward model / observation_map**
- Build EKP from the series instead of `(y_obs, noise_cov)`, **keeping
  `Unscented`**:
  ```julia
  obs_series = JLD2.load(obs_filepath)["observation_series"]
  ekp = EKP.EnsembleKalmanProcess(
      obs_series,
      EKP.Unscented(prior; α_reg = 1.0, update_freq = 1);
      rng,
  )
  ```
  (Only the first argument changes from `y_obs, noise_cov`; the estimator line is
  identical to `Calibration.jl:147`.)
- **Forward model** (`labile_pin_shim.jl`'s `forward_model`, delegating into
  `model_interface_wLabile.jl`): read the EKP's current minibatch, map each window
  index → `(start,stop)` via `samples.jld2`, and for each window run a fresh
  `spinup + window` simulation, drop the spinup prefix, extract the per-depth daily
  diagnostics, and **concatenate** across depths then across windows into one
  member `g`. Today `forward_model` reads a single fixed date range from
  `NEON_START_DATE`/`NEON_STOP_DATE` ENV (`model_interface_wLabile.jl:82-83`) — that
  per-member ENV date pair is replaced by a per-window loop. The per-depth
  layer-selection (`argmin` against grid, `:120-133`) runs unchanged inside the
  loop.
- **`observation_map`** (`model_interface_wLabile.jl:523+`) currently rebuilds the
  stacked G from a single fixed `OBS_FILEPATH`'s `obs_dates_per_code`. It must
  instead read the **current minibatch** from the EKP it already loads
  (`ekp = JLD2.load_object(ekp_path(...))`, `:524`), pull that minibatch's window
  list + per-window/per-depth obs dates from `samples.jld2`, and align + concat in
  the SAME window×depth order the forward model wrote. This is the second
  genuinely new model-side piece (the first being the forward-model loop).
- How does the forward model know the current minibatch? Two options, pick during
  implementation: (a) `forward_model` loads the EKP at
  `ekp_path(OUTPUT_DIR, iteration)` (same file `observation_map` loads) and calls
  `EKP.get_current_minibatch(ekp)`; (b) the driver writes the per-iteration
  minibatch window list to a small file the workers read. (a) mirrors the PR and
  avoids extra I/O contracts — **lean (a)**.
- Keep `WorkerBackend` (Slurm-workers-per-ensemble) — minibatching is orthogonal
  to the backend, so no change there (`Calibration.jl:224-231`).

**4. Output layout / `output_dir_for`**
- Current per-run dir embeds `<site>_<start>_<stop>` (`Config.jl:220-235`).
  Replace the date span with a site+span id like
  `<site>_<firstyear>-<lastyear>_minibatch` so one site = one output dir = one
  posterior. Keep the `CalDepth-501-502-503` depth tag and `SpinUP-<n>d` folders.

**5. Driver (`run_pipeline.jl` / `Pipeline.jl`)**
- `--year` filter loses meaning (a run is multi-year); keep `--site`.
- Step sequence unchanged (generate obs → calibrate → diagnostics →
  optimized-param forward run); each step now operates on a series.
- `Pipeline.jl` currently `include`s `Observation_flag.jl` (not `Observations.jl`)
  at `:25` — keep wiring the cap-mask variant, now exposing
  `generate_observation_series`.

---

## Open decisions (resolve before coding)

1. **Window definition.** PR uses fixed 365-day windows from
   `data_start + spinup`. For NEON, choose:
   - (a) **Calendar years** (Jan 1–Dec 31, matching the AllDepth `years=[…]`
     fan-out) — smaller change, continuity with existing per-year results.
     **Leaning this.**
   - (b) PR's rolling 365-day-from-record-start windows.
2. **Cap-mask baseline scope under windowing.** `Observation_flag.jl` computes the
   cold-soil baseline **once over all years** (`capped_dates` on the full
   `obs_df`, `:229-233`). Under per-window observations, decide: keep the
   all-years baseline (compute `capped` once, apply per window — simplest, keeps
   the cross-year definition) vs. per-window baseline. **Leaning all-years**
   (compute the capped-date set once, then drop within each window). This is the
   open `is_capped_day` per-year-baseline decision from memory, now scoped to the
   minibatch design.
3. **Thin-window policy.** A window with too few valid days at some depth: skip
   the whole window, or skip just that depth for that window (ragged stacking)?
   **Leaning skip-the-window** (keeps every window's stacking order identical, so
   `observation_map` alignment stays simple). Log every dropped window/depth so a
   silently short series is never mistaken for full coverage.

Estimator is **decided** (keep `Unscented(α_reg=1.0)`), so it is no longer an open
decision.

---

## Suggested first cut

Thin, validate-then-generalize:
1. **Calendar-year windows + `Unscented(α_reg=1.0)` + `minibatch_size = 1` +
   `spinup_days = 20` + all-years cap baseline.**
2. Validate on **one site** (e.g. CPER — already has a `cper_allyears.toml`
   config in the AllDepth pipeline to port from).
3. Then generalize (other sites, tune `minibatch_size`/`max_n_samples`, wire up
   the results CSV and diagnostics for the multi-window posterior).

## Next actions when we resume
- [ ] Confirm the three open decisions above.
- [ ] Create folder `experiments/calibrate_neon_pipeline_minibatch/` by copying
      `calibrate_neon_pipeline_AllDepth/` wholesale (including `src/`).
- [ ] Pull PR #1732 `generic_site/calibrate_site.jl` for the forward-model
      minibatch-loop reference (not in the local checkout).
- [ ] Implement the 5 changes (Config → Observation series → Calibration +
      forward_model/observation_map minibatch → output layout → driver),
      in that order.
- [ ] Make/port a `config/cper_minibatch.toml` (span the site's years,
      `minibatch_size=1`, `spinup_days=20`, `cal_depth_codes=["501","502","503"]`).
- [ ] Smoke-test on CPER with small `n_iterations`.
