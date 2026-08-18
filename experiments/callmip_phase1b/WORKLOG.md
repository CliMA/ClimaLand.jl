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
