# CalLMIP Phase 1b — Work Log

Append-only journal. Newest entries at the bottom. Every session: date, what was done,
what was checked, what broke, decisions made. Companion to [PLAN.md](PLAN.md) and
[DATA_MANIFEST.md](DATA_MANIFEST.md).

---

## 2026-08-18 — Phase 0 (session 1)

### Context gathered (read-only)

- Read CalLMIP protocol v2 (checklist, §3.1 scenarios, §4 forcing, §5 flux obs,
  §6 outputs, §8 spin-up, Table 1/2 site lists); key points captured in PLAN.md.
- Reviewed PR #1784 (`rb/DK-Sor`, Phase 1a pipeline — our template), PR #1826
  (ColumnEnsemble), ClimaUtilities #248 (matrix TVI, branch `kp/point-cloud-intp`),
  ClimaCore #2525 (MERGED, v0.15.1), Kevin's `kp/diff-multi-cols` demo branch,
  issue #1723 (GPU multi-site design).
- Kevin's Slack caveats logged in PLAN.md ("Key design constraints").
- Machine: 2× Tesla V100-PCIE-16GB, Julia 1.12.6, CUDA driver 13.0.

### Data verification (see DATA_MANIFEST.md)

- Both callmip artifact tree hashes recomputed and match published values.
- **Found stale flux artifact**: AU-How present / US-Var missing (upstream swap
  2026-07-17 post-dates artifact creation 2026-06-30). Built refreshed artifact
  `9b4bd0d8…` at `/net/sampo/data1/renatob/callmip_artifacts/callmip_phase1_20260818`
  from callmip-org/Phase1 @ `d6d36ad`; bumped `Artifacts.toml` on this branch.
  Follow-up: upstream ClimaArtifacts refresh.
- Overrides.toml: added 3 entries (backup at `Overrides.toml.bak-20260818`).
- Confirmed met-window intersection is 2009–2010 only → padded-union design (D3).
- Confirmed met files carry no UTC-offset attribute (checked US-SRG) → site table needed.

### Branch setup

- Committed Phase 1a WIP snapshot on `callmip-phase1a-dksor-ces` (`026039141`) to
  clear the working tree. Untracked local analysis files left in place (131 files,
  no collisions with target branch — verified).
- Created `rb/callmip-phase1b` from `origin/kp/diff-multi-cols` (`f829d7cfd`).
- Brought in from `origin/rb/DK-Sor` (PR #1784), applied as 3-way diffs vs its
  merge-base with main (`25d1af070`):
  - `experiments/callmip_dksor/` (whole pipeline, 10 files) + `ext/fluxnet_simulations/DK-Sor.jl` (plain checkout, new files);
  - `src/diagnostics/{default_diagnostics,land_compute_methods}.jl` (+7 lines: soilrn/soillhf/soilshf for LandModel) — applied cleanly;
  - `ext/FluxnetSimulationsExt.jl`, `ext/fluxnet_simulations/forcing.jl`, `src/ext/FluxnetSimulations.jl` — applied cleanly;
  - `src/standalone/Vegetation/spatially_varying_parameters.jl` — **conflicted; dropped**:
    the CLM-map hole-patching (`zeros_to_val`) is already on main (positional-arg form),
    identical values; kept HEAD.
  - Deliberately NOT taken: `src/Artifacts.jl` (kp branch version is newer — rosetta
    rename; callmip accessors functionally identical), `Artifacts.toml` rosetta entries,
    `.github/`, `docs/` manifests, `.buildkite/Manifest` (kp's pins ClimaUtilities
    `kp/point-cloud-intp` — required).
- Net diff vs `kp/diff-multi-cols`: +2,672 lines, 16 files (verified `git diff --stat`).

### Site metadata extraction + UTC-offset verification (DATA_MANIFEST)

- Extracted lat/lon/z_ref/elev/time-units from all 21 met files (netCDF4/python).
- Met lengths validate exactly against year spans (e.g. DK-Sor 18 yr = 6574 d × 48).
- **UTC offsets verified empirically** from SWdown diurnal cycle:
  - centroid method flagged US-NR1 (−7.79) and US-MMS (−5.56) — turned out to be
    afternoon-convection bias, NOT wrong offsets;
  - sunrise/sunset-midpoint method (unbiased) confirms FLUXNET political offsets at
    **all 21 sites within ±0.15 h**. Naive `round(lon/15)` wrong at FR-Pue, NL-Loo
    (+1), RU-Fyo (+3), US-MMS (−5) — do NOT derive offsets from longitude.
- Verified all 21 refreshed flux files: `{NEE,Qle,Qh}_daily` (+`_uc_daily`), NEE in
  gC/m2/d. **Coverage finding:** NEE valid days very sparse at 5 sites (US-Ton 0.3%,
  CH-Dav/DE-Hai/IT-Lav/US-Var ≤ 1.8%); Qle/Qh much better (26–97%). Totals:
  NEE 16.6k / Qle 60.5k / Qh 67.4k valid daily obs. → per-month masks mandatory;
  misfit energy-dominated at NEE-poor sites (noted for Phase 4).

### Environment

- `.buildkite` instantiate + precompile started in background
  (Julia 1.12.6; Manifest pins ClimaCore 0.15.1 + ClimaUtilities `kp/point-cloud-intp`).

### Housekeeping (same day)

- Renato: **never add `Co-Authored-By: Claude` trailers** — rewrote all 3 local commits
  (filter-branch / amend); hashes cited here are the rewritten ones.
- Formatter: CI pins JuliaFormatter **2.10.1** (`.github/workflows/JuliaFormatter.yml`).
  Running a newer version dirtied 60 unrelated files (reverted). With 2.10.1 the repo
  base is clean; committed real fixes for our 11 changeset files (`770b26986`).
  Kevin's 4 demo files intentionally left unformatted (minimize divergence vs his
  branch for future rebases).

### Phase 0 gate checklist

- [x] Data verified (artifacts, hashes, windows, offsets, flux structure/coverage).
- [x] Branch `rb/callmip-phase1b` assembled and committed (`1ac8e4e81`).
- [x] Docs scaffolding (this file, PLAN.md, DATA_MANIFEST.md).
- [x] `.buildkite` project instantiated + precompiled clean (728 pkgs, exit 0).
- [x] Gate 1/4: `callmip_forcing_demo.jl` CPU — PASS (2026-08-18). Shared UTC axis
  = exact DK-Sor∩US-NR1 overlap (280,496 half-hours); per-column drivers physically
  sane at t=0 and t=+12 h (night/day, winter T, LAI 0.21 vs 1.38); per-column
  `atmos_h` Field = [57, 26] m in correct site order.
- [x] Gate 2/4: `callmip_run_demo.jl` CPU — PASS (2026-08-18). 5-day, 2-column
  LandModel, Δt=450 s: **~780 SYPD** on CPU after compile (1.5 s integration;
  97% of 70 s wall = one-time compilation). Final states physically consistent
  (DK-Sor 01:00 local: night, canopy < air; US-NR1 17:00 local: SW=379 W/m²,
  canopy > air). Observed: US-NR1 top-soil T=272.3 K in July = default cold IC
  still visible after only 5 days (expected; equivalence scrutiny in Phase 1,
  spin-up in real runs). Note: run scripts with ABSOLUTE paths (shell cwd resets
  broke a first attempt with a relative path).
- [x] Gate 3/4: `callmip_forcing_demo.jl` CUDA — PASS (2026-08-18, after CUDA-runtime
  fix below). Driver values bit-identical to CPU; per-column heights [57, 26] correct.
- [x] Gate 4/4: `callmip_run_demo.jl` CUDA — PASS (2026-08-18, after atmos_h device
  fix below). Final states match CPU to ~11 significant digits (US-NR1 T_canopy equal
  to last printed digit). GPU ≈ 125 SYPD at 2 columns vs CPU ≈ 780 SYPD
  (launch-latency-bound at tiny N — Phase 3 must find the crossover). CPU rerun after
  the fix: **bit-identical** to pre-fix CPU run (no regression).

**PHASE 0 CLOSED 2026-08-18.** All four gates pass; environment, data, branch, and
docs verified. Upstream findings recorded in UPSTREAM_NOTES.md (atmos_h GPU bug for
Kevin; stale ClimaArtifacts flux artifact).

### GPU environment incident + fix (2026-08-18) — READ IF GPU BREAKS

First CUDA attempt failed at `ColumnEnsemble(...)` space construction:
`ERROR: Compute capability 7.0.0 is not supported by CUDA 13.2.0`.
**Cause:** our Tesla V100s are Volta (sm_70); NVIDIA dropped Volta in CUDA 13.x, and
CUDA.jl (5.11.3) auto-picked the 13.2 runtime to match the v580 driver.
**Failed first fix (instructive):** `CUDA.set_runtime_version!(v"12.9")` in `.buildkite`
writes `LocalPreferences.toml` **and adds `CUDA_Runtime_jll` to `[extras]` in
`Project.toml`** — the extras entry is required for the preference to be honored.
Reverting the Project.toml churn (`git checkout`) silently disabled the preference →
same crash on the next run (Preferences.jl only honors prefs for packages listed in a
load-path project).

**Working fix (machine-global, zero repo churn):** put the preference in the default
`v1.12` environment, which is on every project's load path:
- `~/.julia/environments/v1.12/Project.toml`: add `[extras] CUDA_Runtime_jll = "76a88914-d11a-5bdc-97e0-2f5a05c973a2"`
- `~/.julia/environments/v1.12/LocalPreferences.toml`: `[CUDA_Runtime_jll]` / `version = "12.9"`

(backup: `Project.toml.bak-20260818`). Verified: kernel launch OK on sm_70,
`CUDA.runtime_version() = 12.9.0`, both V100s visible (~15.5 GiB free each); repo
tracked files stay clean; the redundant `.buildkite/LocalPreferences.toml` was removed.
If Julia is upgraded to 1.13+, replicate the two entries in the new default env.

### GPU bug found by gate 4 + fix (2026-08-18) — per-column atmos_h was CPU-resident

With CUDA working, `callmip_run_demo.jl` crashed in the snow model's surface-temp
broadcast: `ThisHost and ThisThreadPool do not overlap` — a **host-resident field mixed
into a GPU broadcast**. Culprit: `prescribed_forcing_callmip` built the per-column
tower-height field as `array2field(heights::Vector, surface_space)`, wrapping a CPU
`Array` even when the space lives on the GPU. Exactly the scalar→Field class of issue
Kevin warned about; invisible on CPU, fatal on GPU.
Fix in `experiments/integrated/fluxnet/callmip_forcing.jl`: convert with
`ClimaComms.array_type(ClimaComms.device(surface_space))(heights)` before
`array2field`. **TODO: report upstream to Kevin (kp/diff-multi-cols / PR #1826)** —
any other `array2field`/`parent(...) .=` per-column constructions must use the space's
device array type. Grep pattern for future scalar→Field audits:
`ThisHost.*ThisThreadPool` in errors = host/device mix.


---

## 2026-08-18 — Phase 1 (session 1, continued)

### Correctness harness

- era5 `column_ensemble_comparison.jl` (Kevin's duplicated-column test) CPU:
  **7/7 pass, 15.7 s** (duplicates exactly identical; single-vs-ensemble ≤ 1e-9;
  diagnostics agree at ~1e-13 rel). GPU repeat: running.
- NEW `experiments/callmip_phase1b/verify_single_vs_ensemble.jl` — covers what the
  era5 test cannot: **different forcing per column** (matrix TVI), cross-column
  contamination, row-order robustness, forcing fidelity. CPU results:
  - Single DK-Sor `Column` vs DK-Sor as col 1 of [DK-Sor, US-NR1] ensemble:
    **9/9 pass**; worst err 4.1e-12 (`Y.soil.θ_i`, scale 2e-5), most fields ≤ 8e-14
    (tolerance 1e-9, era5 precedent). → no cross-column contamination.
  - **Column-order flip [US-NR1, DK-Sor]: EXACTLY 0 error** for both sites' state and
    cache — solution is bitwise independent of column position.
  - Forcing fidelity: **63/63** — 7 drivers × 9 raw met sample times reproduce the
    file values per column (rtol 1e-12) → row order and interpolation verified
    end-to-end.
- GPU repeats: first attempts were SIGTERM-killed by a session restart (leftover
  processes killed by hand — check `ps aux | grep julia` if GPUs look busy with no
  progress). Relaunched:
  - `verify_single_vs_ensemble.jl` GPU (V100, CUDA 12.9): **9/9 + 63/63 pass**;
    order-flip exactly 0 on GPU as well; worst single-vs-ensemble err 3.2e-14.
  - era5 `column_ensemble_comparison.jl` GPU: **7/7 pass, 16.3 s** (worst err 1.7e-12).

### Phase 1 gate — CLOSED 2026-08-18

- [x] era5 duplicated-column test: CPU 7/7, GPU 7/7.
- [x] DK-Sor single-Column vs 2-site ensemble (different forcing per column):
  CPU 9/9 (worst 4.1e-12), GPU 9/9 (worst 3.2e-14) — tolerance 1e-9.
- [x] Column-order flip: EXACTLY 0 difference, CPU and GPU.
- [x] Forcing row-order/fidelity: 63/63 on both devices (drivers reproduce raw met
  values per column at sample times, rtol 1e-12).
- [x] scalar→Field gaps: atmos_h fixed in Phase 0; none new surfaced at N=2.
  (Larger N with the default LandModel may surface more — watched in Phase 2.)

**PHASE 1 CLOSED.** Multi-forcing ColumnEnsemble machinery is verified correct at
N=2 on CPU and GPU. Next: Phase 2 (21-site forcing, union-axis padding).


---

## 2026-08-18 — Phase 2 (session 1, continued)

### D7: LAI source

All 21 met files carry BOTH `LAI` and `LAI_alternative` at 100% validity
(PLUMBER2 climatology-fills outside product coverage) → availability does not
constrain the choice. Decision D7 (PLAN.md): use `LAI` (the PLUMBER2-preferred
product per protocol §4) for Phase 1b, `lai_var` override available. Differs from
Phase 1a DK-Sor (`LAI_alternative`/MODIS) — documented.

### Implementation

- `sites.jl` — verified 21-site table (coords, UTC offsets, z_ref, met years) +
  validation-site ID list.
- `forcing_phase1b.jl` — `load_site` (table-checked reads), `union_axis`,
  `pad_to_axis` (interior-year recycling, Feb-29→28 fallback, `native` mask),
  `load_calibration_sites(; window)`.
- `read_callmip_met` gained a `lai_var` kwarg (default unchanged =
  `LAI_alternative`, so demos/Phase-1 scripts are untouched).
- `verify_21site_forcing.jl` — Phase 2 gate script: padding fidelity (native
  untouched, recycled = mapped interior-year sample, all finite), 21-column
  driver order/fidelity, tower-height order, 5-day 21-column forward run with
  per-column finiteness + physical canopy-T range. CPU run in progress.


### Gate iterations (findings while making verify_21site_forcing.jl pass)

1. `US-SRG tower height mismatch` — table said 3.2 (rounded print); file says
   **3.25**. Fixed table + manifest. Also: heights are Float32 in files → all
   equality checks use atol=1e-3.
   Process lesson: a `sed` fix silently no-opped (pattern predated formatter
   reflow; sed exits 0 on zero matches) → **only assert-guarded python edits**.
2. `KeyError: 2010-12-31T22:30` in padding — **US-MMS met is HOURLY** (only such
   site; 5,844 d × 24). Added `resample_to_30min` (linear midpoints ≡ TVI linear
   interp; physically a no-op) wired into `load_site`. Kevin's `align_sites`
   would also have rejected the mixed grid (its equispaced/shared-dt assert).
3. D7 revised per Renato mid-session: **MODIS everywhere**. Implemented by the
   `source` attribute via `modis_lai_var` — variable names flip per site
   (13/21 keep MODIS under `LAI`), so name-based selection would silently pick
   Copernicus at 13 sites.


### Phase 2 gate — CPU PASS (2026-08-18)

- Padding correctness (full union axis 1995-12-31T22:00 → 2015-01-01T07:30,
  333,140 half-hours): **840/840** — native samples untouched, recycled samples
  equal the independently recomputed interior-year mapping, everything finite.
- 21-column driver order/fidelity: **22/22** (7 sample times × 3 drivers vs raw
  values per column; tower heights in column order at atol 1e-3).
- 21-column × 5-day forward run (default LandModel): **85/85** finiteness checks;
  July canopy T spans 278–304 K across the ensemble (boreal → desert). Wall 54 s
  including compile → the 21.9 "SYPD" is a floor, not a throughput number
  (Phase 3 measures properly).
- Committed `344ed7932`. GPU repeat running.
- Data-location decision (Renato): artifacts resolve from **sampo**; the kiwi copy
  at /home/renatob/data/callmip_data is a spare (its flux subset is the STALE set).
- Validation download list sent to Renato: 10 met files from ME-org "Phase 1b
  Validation" (site IDs + Table-2 record lengths in the chat + DATA_MANIFEST).


### Validation-site met verified + 31-site artifact (2026-08-18)

- Renato downloaded the 10 ME-org validation met files (zip →
  /net/sampo/data1/renatob/callmip_validation_met/). Verified: full driver set,
  UTC offsets empirical-confirmed (AU = +9.5 half-hour zone!), 3 HOURLY sites
  (AU-Tum, US-Ha1, US-UMB), provenance mix (OzFlux/LaThuile/FLUXNET2015).
- Built 31-site forcing artifact `a597745c…`; Artifacts.toml bumped, override added.
- sites.jl: verified VALIDATION_SITES table; site_info searches both tables.
- Bug found by the check: my start-year cross-check used `Dates.Hour(offset)` →
  `InexactError` on +9.5; fixed with `FSExt.hour_offset_to_period` (fraction-aware).
- End-to-end load of all 10 validation sites: CLEAN (30-min grids after resample,
  MODIS by attribute, zero non-finite).
- Per Renato: data stays on sampo (kiwi copy is a spare only).


### Phase 2 gate — GPU PASS → PHASE 2 CLOSED (2026-08-18)

- GPU (V100, CUDA 12.9): padding 840/840, driver order 22/22, forward run 85/85;
  final canopy-T extrema IDENTICAL to CPU to all printed digits
  (278.1357825374…, 304.0405517522…). Wall 99.6 s incl. compile.
- Phase 2 deliverables: sites.jl (31 verified sites), forcing_phase1b.jl
  (padding + resampling + MODIS-by-attribute), 31-site artifact `a597745c…`.
  Commits `344ed7932`, `08addda41`.

## 2026-08-18 — Phase 3 (session 1, continued)

Benchmark: SYPD vs column count, CPU vs V100, default LandModel, Δt=450 s,
10-day timed window (1-day warmup discarded), duplicated sites cycled to reach
N ∈ {1, 2, 5, 10, 21, 30, 100, 300}. Purpose: pick the EKI member schedule
(GPU vs CPU crossover) + post numbers to PR #1826 (kmdeck's request).


### Phase 3 gate — CLOSED (2026-08-19)

Benchmark ran detached overnight, both sweeps exit 0
(benchmark_logs/bench_{gpu,cpu}.log; N=1 rows are compile-contaminated, ignore).

| N cols | V100 SYPD | V100 col-SYPD | CPU-1core SYPD | CPU col-SYPD |
|---|---|---|---|---|
| 2 | 145 | 291 | 1099 | 2199 |
| 5 | 141 | 703 | 784 | 3921 |
| 10 | 150 | 1497 | 505 | 5049 |
| 21 | **151** | 3178 | **289** | 6070 |
| 30 | 147 | 4421 | 212 | 6368 |
| 100 | 151 | 15087 | 69 | 6862 |
| 300 | 146 | 43659 | 24 | 7183 |

- V100 is launch-latency-bound: SYPD flat ~150 from N=2 to N=300; column-throughput
  linear and UNSATURATED at 43.7k col-SYPD (N=300). CPU saturates ~7k col-SYPD.
  Crossover ≈ 50 columns/process.
- **Decision: calibration engine = CPU member-parallel** (one 21-col ColumnEnsemble
  per member per core; machine has 96 logical CPUs / 754 GB). Wall-time model:
  3-sim-yr member run ≈ 15 min → EKI iteration (30 members parallel) ≈ 15 min →
  10 iterations ≈ **2.5 h**. GPUs reserved for the sites×members mega-batch
  (630 cols ≈ 29 min/iteration on ONE V100 — viable fallback if per-column params
  ever get implemented) and for PR #1826 evidence.
- Numbers added to UPSTREAM_NOTES.md for kmdeck's #1826 question.


## 2026-08-19 — Phase 4 (session 2)

### Implementation (mirrors Phase 1a conventions exactly)

- `observations_phase1b.jl` — per-site Phase 1a recipe (monthly means,
  MIN_VALID_DAYS=5, SIGMA2_MISS=1e12, noise = FLUXNET σ²/n + inter-annual floor,
  valid-dates threading), site-major 756-entry windows. Built 2003–2014:
  **12 windows, 5,708/9,072 unmasked (63%)**. NEE contributes at 16/21 sites
  (0 months at CH-Dav/DE-Hai/US-Ton/US-Var, 1 at IT-Lav); LHF/SHF everywhere.
- `forward_model_phase1b.jl` — `run_member(θ, years)`: TOML override → default
  LandModel (D1) on 21-col ColumnEnsemble → leading spin-up year (D5) → daily
  nee/lhf/shf via DictWriter → monthly all-day means per column (GMATCH=allday),
  NEE ×12×86400 → gC/m²/d; year-major site-major layout; NaN on failure.
  Padded sites cached once per process (`sites_full()`), windows cropped
  in-memory (`crop_sites` added to forcing_phase1b.jl).
- `run_calibration_phase1b.jl` — UTKI (TransformUnscented, impose_prior,
  2p+1=23 members), biennial consecutive minibatches (G dim 1512), pmap on 23
  CPU workers, DataMisfitController(100), JLD2 checkpoint/resume, --test mode
  (serial 2 members × 1 iteration).
- Priors: reusing Phase 1a `priors.jl` verbatim (D2) — note `pmodel_β_c4` is now
  genuinely constrained (3 C4-grass sites).


### Phase 4 gate — CLOSED (2026-08-19)

- G smoke: 5/5 (756 finite entries, plausible ranges, 12.5 min for 2 sim-yr).
- Gate iteration attempt 1: all 23 members NaN — missing `using Statistics` in
  forward_model_phase1b.jl (smoke test masked it by importing before include).
  Fixed (`776f2042f`), stale vacuous checkpoint deleted, rerun.
- Gate iteration (real, 23 workers, window [2005,2006]): **NaN members 0**,
  G_ens (1512×23), wall **19.6 min**, real UTKI update (T=0.061; β_c3 159→42,
  cstar 0.36→0.51; β_c4 ~unchanged as expected from one window).
- **PHASE 4 CLOSED.** Phase 5 = resume checkpoint to iteration 10
  (CALLMIP1B_N_ITER=10 via output_calibration/run_iter.sh, detached).


### Upstream feedback (2026-08-19)

- Renato reported the atmos_h GPU bug + ordering verification + timings to Kevin
  (Slack) and the scaling numbers toward kmdeck's #1826 question.
- Kevin confirmed: (a) the shared-time-axis requirement is a known limitation of
  the ClimaUtilities matrix TimeVaryingInput (by design, simplifies the input);
  (b) his own benchmarking found the same qualitative scaling (CPU linear-ish,
  GPU not saturating until large N) — independent corroboration of Phase 3.


## 2026-08-19 — Phase 6 prep (while Phase 5 runs)

- Phase 5 progress: checkpoint iteration 5/10, ~20 min/iter, healthy. (Master
  stdout is buffered until exit — read progress from the checkpoint's
  `iteration` field, not the log.)
- `write_callmip_netcdf_phase1b.jl` — parametrized port of the Phase 1a writer
  (identical variable mapping, units, hfg sign convention, template time axis
  "days since <y0-1>-12-31"; per-site lat/lon/window; Phase1b/Scen1 attrs;
  Cal/Val token).
- `deliverable_runs_phase1b.jl` — per-site prior/posterior runner using THE
  SAME model as the calibration engine (default LandModel, MODIS LAI, 1-site
  ColumnEnsemble): min(5, window−1)-yr spin-up cycling the site's own met,
  FULL-state handoff (recursive parent-copy), production run over the native
  window (day 1 = fill since the UTC start must sit inside all offsets' records;
  tail padded so Dec 31 closes), daily CalLMIP diagnostics → Phase 1a-format
  JLD2 → NetCDF. CLI: `<site_id> <Prior|Posterior>` (Posterior reads the
  UTKI ϕ_mean).
- Smoke: US-SRG Prior (11 sim-yr) running.


### ⚠️ MAJOR CATCH: default LandModel omits soilco2 — Phase 5 restarted (2026-08-19)

The deliverable smoke (US-SRG) failed the diagnostics whitelist assertion for
`hr`/`soc` — which exposed that the convenience constructor
`LandModel{FT}(forcing, LAI, toml_dict, domain, Δt)` defaults to
`prognostic_land_components = (:canopy, :snow, :soil)` — **no SoilCO2 model**
(type parameter MM = Nothing). Consequences of the model we had been calibrating:
- 4 of the 11 calibrated parameters (all DAMM heterotrophic-respiration params)
  had ZERO effect;
- modeled NEE omitted heterotrophic respiration entirely (er = autotrophic only);
- protocol variable cSoil impossible.
The calibration G-map only requested nee/lhf/shf, so nothing errored — silent
scientific wrongness caught by the wider deliverable diagnostics. Kevin's demo
and our Phases 0–3 verification are UNAFFECTED (they test forcing/space
machinery, not carbon completeness); the Phase 3 SYPD numbers were measured on
the no-soilco2 model and will shift somewhat.

Actions:
- **Phase 5 stopped at iteration ~6; checkpoint quarantined**
  (`ekp_checkpoint_NOSOILCO2_INVALID.jld2`, log `calibration_iter_NOSOILCO2.log`).
- `prognostic_land_components = (:canopy, :snow, :soil, :soilco2)` added to the
  LandModel calls in forward_model_phase1b.jl and deliverable_runs_phase1b.jl.
- Re-verification ladder restarted: G smoke → deliverable smoke → Phase 4 gate
  iteration → fresh Phase 5.
Lesson recorded: when using a convenience constructor, ASSERT the components you
scientifically require (a `@assert !isnothing(land.soilco2)` now belongs in the
forward model).


### Re-verification after soilco2 fix + Phase 5 relaunch (2026-08-19)

- G smoke with soilco2: **5/5**; NEE source range doubles (+0.8 → +1.8 gC/m²/d),
  LHF/SHF bit-identical (soilco2 has no energy feedback — a good consistency
  check); member cost 14.5 min (+16%).
- Deliverable smoke US-SRG Prior: **full chain green** — first Phase 1b NetCDF
  written and validated against protocol (dims/attrs/units/time axis exact;
  99.9% valid = the 2 by-design boundary fill days; physics sane for semi-arid
  grassland). Science note: prior reco floor ~1 gC/m²/d looks high for desert
  winter — prior bias, calibration's job.
- Fix committed `b3c0f2200`. **Fresh Phase 5 launched** (10 iterations from
  scratch, soilco2 active, ~23 min/iter ≈ 4 h). Iteration 1 doubles as the
  Phase 4 gate re-run — verify at first checkpoint.


### All-site prior runs + fig5 (2026-08-25)

- All 21 Cal_Prior deliverable runs completed (every log exit=0) — a third of the
  final deliverables done, and the prior-skill baseline for the Phase 5 gate.
- figures/fig5{a,b,c}_prior_vs_obs_all_sites.png: prior vs obs climatology, all sites.
- Prior-skill read: (1) forest NEE seasonality is in phase, amplitude mixed
  (over-uptake at DE-Gri/RU-Fyo/US-Ton, under at IT-Lav); (2) the semi-arid
  cluster (US-SRG/SRM/Whs/Wkg) NEVER shows model uptake while obs have monsoon
  uptake — the calibration's clearest target (moisture_stress_c already ↑ in the
  posterior trajectory); (3) systematic energy bias at forests: LE over, H under
  (Bowen ratio too low) — carbon params touch this only via moisture stress; if
  energy misfit dominates post-calibration, D2 (no energy params) needs revisiting.


### Phase 5/6 skill verdict (2026-08-25) — figs 6+7, skill_{prior,posterior}.csv

Calibration extended to 15 iterations; noise-weighted misfit 10.5 → 0.35–0.66
plateau (~25×). Window oscillation on β_c3/msc → posterior = epoch average of
iteration means 10–15 (epoch_posterior.jld2). All 21 posterior deliverable runs
completed clean; per-site monthly-RMSE gate:
- **LE: 20/21 improved (95%)** — the LE-over bias largely fixed.
- **NEE: 10/16 improved (62%)** — drylands much better (US-SRG R² +0.30) but
  **DE-Tha 0.76→1.41 and DK-Sor 1.30→1.85 gC/m²/d degrade >30%** (the global
  β_c3=21.5 serves drylands at the temperate forests' expense).
- **H: 4/21 improved (19%)** — small ubiquitous tax (~1–4 W/m²).
Pre-registered ≥⅔ gate: LE passes, NEE just misses, H fails → NOT a clean pass.
Decision needed (Renato): accept Scenario-1 posterior as-is (defensible;
objective improved 25×; document trade-offs) vs PFT-split β_c3 (protocol allows
PFT-specific params; directly targets the forest/dryland tension; ~5 h rerun).
