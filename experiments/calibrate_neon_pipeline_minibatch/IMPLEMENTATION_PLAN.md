# `calibrate_neon_pipeline_minibatch` — Implementation Plan

Status: **planning / not yet implemented.** This document captures the design so
we can resume work later. Written 2026-06-25.

## Goal

Set up a new folder `experiments/calibrate_neon_pipeline_minibatch/` that
calibrates the NEON soil-CO₂ (DAMM) parameters for **one site using all of its
years jointly**, instead of the current pipeline's one-calibration-per-year
approach. The mechanism mirrors ClimaLand PR
[#1732](https://github.com/CliMA/ClimaLand.jl/pull/1732)
(`experiments/integrated/generic_site/`): an EKP `ObservationSeries` of yearly
windows + a minibatcher.

## Reference points

- **Current pipeline:** `experiments/calibrate_neon_pipeline/`
  - Entry: `run_pipeline.jl`
  - `pipeline/Config.jl`, `Observations.jl`, `Calibration.jl`, `ForwardRun.jl`,
    `Diagnostics.jl`, `ResultsTable.jl`, `Pipeline.jl`
  - Heavy model (reused, not duplicated): `../calibrate_neon/model_interface_wLabile.jl`,
    `site_metadata.jl`, and `pipeline/labile_pin_shim.jl`
- **PR #1732 minibatch reference:** `experiments/integrated/generic_site/calibrate_site.jl`
  and `configs/er.jl` (the DAMM-kinetics config, closest analog to NEON soil-CO₂).

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
6. Estimator: `EKP.TransformUnscented(prior; impose_prior = true)` via
   `ClimaCalibrate.calibrate(JuliaBackend(), ...)`.
7. Scale-out: **sites in parallel** as a Slurm job array (one array task per
   site); ensemble members within a site run on `JuliaBackend`.

`er.jl` config: 4 DAMM params, `cal_window_days=365`, `spinup_days=90`,
`max_n_samples=5`, `minibatch_size=1`, `n_iterations=5`, `rng_seed=42`.

---

## Differences from the current NEON pipeline

| | Current `calibrate_neon_pipeline` | New `..._minibatch` (this plan) |
|---|---|---|
| Unit of calibration | one **year** (separate calibration each) | one **site**, all years jointly |
| Multi-year handling | independent runs, separate posteriors | one `ObservationSeries`, one posterior |
| Years enter EKP via | — (separate calibrations) | **minibatcher** (draws `minibatch_size`/iter) |
| Observation object | single `y_obs` + `noise_cov` | vector of per-year `EKP.Observation`s |
| Estimator | `EKP.Unscented(α_reg=1.0)` | `EKP.TransformUnscented(impose_prior=true)` |
| Forward model per member | one sim of the year | loop over current minibatch windows, concat `g` |
| Parallelism | ensemble members across workers (`WorkerBackend`) | keep `WorkerBackend` (orthogonal to minibatch) |
| Spinup | trims `spinup_days` off front of obs | run `spinup_days` *before* each window, then drop |

---

## Version check (done 2026-06-25)

`.buildkite` environment of `ClimaLand.jl`:
- **EnsembleKalmanProcesses 2.7.1** — `TransformUnscented(prior; impose_prior=true)`
  exists (`UnscentedKalmanInversion.jl:257`; `impose_prior` kw at line 81,
  default `false`). `Unscented(prior; ...)` is the same package → swapping
  estimators is one line.
- **ClimaCalibrate** has `minibatcher_over_samples(n_samples, batch_size)`
  (`ekp_utils.jl`).

So both estimator options are viable; plan defaults to `TransformUnscented`.

### Why `TransformUnscented` (not strictly required)

Minibatching is a property of the EKP (obs + minibatcher), **not** the
estimator — `Unscented` would also work. We prefer `TransformUnscented`
because:
- `impose_prior = true` regularizes the posterior toward the prior, which
  matters more under minibatching (each iter sees only a subset → noisier
  updates). `er.jl` explicitly widened a prior to stop EKI leaning against a
  bound "on single-site data" — the instability `impose_prior` addresses.
- The `Transform*` family is the current EKP-recommended/most-tested path.
- Same derivative-free deterministic sigma-point method as `Unscented`, same
  ensemble size → low-cost switch.

**Fallback:** keep `Unscented` for a minimal-diff first cut if we want to
validate mechanics before adopting the new estimator.

---

## What gets reused vs. changed

### Reused essentially unchanged
- `model_interface_wLabile.jl`, `site_metadata.jl`, `labile_pin_shim.jl` — heavy
  model + labile pinning. **Not duplicated**; `include`d from `../calibrate_neon`
  exactly as the current pipeline does.
- `ResultsTable.jl`, `Diagnostics.jl`, `ForwardRun.jl` (prior-mean run),
  provenance snapshotting — minor path/metadata tweaks only.
- TOML-driven `Config` design, variable-length priors, depth handling
  (501 → 0.02 m, 502 → 0.06 m).

### Changed

**1. `Config.jl` — `RunConfig` carries multiple windows**
- Today `_date_ranges` *fans out* `years=[…]` into separate `RunConfig`s.
  Now **one `RunConfig` holds the full `Vector{(start,stop)}`** for the site.
- New fields: `minibatch_size`, `max_n_samples`, `rng_seed`; reinterpret
  `spinup_days` as "spinup *before each window*" (PR semantics), not
  "trim off the front".
- New `[settings]`/`[runs]` TOML keys: `minibatch_size`, `max_n_samples`,
  `rng_seed`. Default `minibatch_size = 1` (mirrors `er.jl`).

**2. `Observations.jl` — build an `ObservationSeries`**
- Refactor the single-vector builder into a **per-window** function
  `build_window_observation(run, window_start, window_stop) -> EKP.Observation`
  (existing daily-mean soil-CO₂ + inter-sensor-variance noise, scoped to one
  window).
- New `generate_observation_series(run)` loops the run's windows, builds one
  `Observation` each, wraps with
  `minibatcher_over_samples(n_windows, run.minibatch_size)` into an
  `EKP.ObservationSeries`. Saves `observation_series.jld2` **and** `samples.jld2`
  (the `[(window_start, window_stop), …]` list) for the forward model to map a
  minibatch index → dates.
- Drop windows failing the NEON data-coverage check (too few valid days),
  paralleling the PR's per-window validity guard.

**3. `Calibration.jl` — series + minibatch-aware forward model**
- Build EKP from the series instead of `(y_obs, noise_cov)`:
  ```julia
  ekp = EKP.EnsembleKalmanProcess(
      obs_series,
      EKP.TransformUnscented(prior; impose_prior = true);
      rng, verbose = true,
  )
  ```
- **Forward model** (in `labile_pin_shim.jl`'s `NeonPipelineInterface`, or a
  minibatch variant): read `EKP.get_current_minibatch(ekp)`; for each window
  index set the dates to that window, run the sim, **concatenate** each window's
  daily `g` into one member vector (layout matches the concatenated obs for the
  minibatch). This is the one genuinely new model-side piece — today the forward
  model runs a single fixed date range from `NEON_START_DATE`/`NEON_STOP_DATE`.
- Per-window date now comes from the minibatch + `samples.jld2`, not fixed ENV.
- Keep `WorkerBackend` (Slurm-workers-per-ensemble) — minibatching is orthogonal
  to the backend, so no change there.

**4. Output layout / `output_dir_for`**
- Current per-run dir embeds `<start>_<stop>`. Replace with a span/site id like
  `<site>_<firstyear>-<lastyear>_minibatch` so one site = one output dir = one
  posterior.

**5. Driver (`run_pipeline.jl` / `Pipeline.jl`)**
- `--year` filter loses meaning (a run is multi-year); keep `--site`.
- Step sequence unchanged (generate obs → calibrate → diagnostics →
  prior-mean run); each step now operates on a series.

---

## Open decisions (resolve before coding)

1. **Window definition.** PR uses fixed 365-day windows from
   `data_start + spinup`. For NEON, choose:
   - (a) **Calendar years** (Jan 1–Dec 31, matching current `years=[…]`) —
     smaller change, continuity with existing per-year results. **Leaning this.**
   - (b) PR's rolling 365-day-from-record-start windows.
2. **Estimator.** `TransformUnscented(impose_prior=true)` (PR-matching,
   recommended) vs. `Unscented` (minimal diff first). **Leaning
   TransformUnscented.**

Also note (NEON spring-peak cap mask, from memory): there's an open question
about an `is_capped_day` mask for NEON soil-CO₂ in `Observations.jl` and a
per-year-baseline decision — revisit how that interacts with per-window obs.

---

## Suggested first cut

Thin, validate-then-generalize:
1. **Calendar-year windows + `TransformUnscented(impose_prior=true)` +
   `minibatch_size = 1`.**
2. Validate on **one site** (e.g. CPER — already has a `cper_allyears.toml`
   config in the current pipeline).
3. Then generalize (other sites, tune `minibatch_size`/`max_n_samples`,
   wire up the results CSV and diagnostics).

## Next actions when we resume
- [ ] Confirm the two open decisions above.
- [ ] Create folder `experiments/calibrate_neon_pipeline_minibatch/` by copying
      the current pipeline structure.
- [ ] Implement the 5 changes (Config → Observations → Calibration → output
      layout → driver), in that order.
- [ ] Make/port a `config/cper_minibatch.toml`.
- [ ] Smoke-test on CPER with small `n_iterations`.
