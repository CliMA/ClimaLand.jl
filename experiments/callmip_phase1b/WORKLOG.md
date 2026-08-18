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

### Phase 0 gate checklist

- [x] Data verified (artifacts, hashes, windows, offsets, flux structure/coverage).
- [x] Branch `rb/callmip-phase1b` assembled and committed (`1ac8e4e81`).
- [x] Docs scaffolding (this file, PLAN.md, DATA_MANIFEST.md).
- [ ] `.buildkite` project instantiated + precompiled clean.
- [ ] Gate: `callmip_forcing_demo.jl` + `callmip_run_demo.jl` on CPU.
- [ ] Gate: same demos on CUDA (`CLIMACOMMS_DEVICE=CUDA`), values sane per column.
