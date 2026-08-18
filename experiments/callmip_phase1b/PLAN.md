# CalLMIP Phase 1b — Multisite Calibration Plan (ClimaLand)

**Branch:** `rb/callmip-phase1b` (based on `kp/diff-multi-cols` = PR #1826 ColumnEnsemble + demos)
**Started:** 2026-08-18 · **Owner:** Renato Braghiere · **Protocol:** CalLMIP_Phase1_Protocol_v2.pdf
**Companion docs:** [WORKLOG.md](WORKLOG.md) (append-only journal), [DATA_MANIFEST.md](DATA_MANIFEST.md) (data verification)

## Goal

CalLMIP Phase 1b, Calibration Scenario 1 (required): calibrate **one operational
(globally shared) ClimaLand parameter set** against daily **NEE, Qle (LE), Qh (H)** at
**21 FLUXNET calibration sites**, then deliver prior + posterior daily NetCDFs for the
21 calibration sites (`Cal`) **and** 10 validation sites (`Val`) to modelevaluation.org
(ME-org): **62 files** total (21×2 Cal + 10×2 Val), named
`ClimaLand.<version>_Phase1b_Scen1_<Site>_<Cal|Val>_<Prior|Posterior>.nc`.

## Architecture (decided 2026-08-18, approved by Renato)

1. **Calibration engine (GPU):** the EKI forward map `G(θ)` = **one `ColumnEnsemble`
   run with 21 columns (one per site) per ensemble member**. Calibrated parameters are
   global scalars injected via ClimaParams TOML override (as in Phase 1a
   `run_calibration.jl`) — all columns share θ within a run, so no scalar→Field surgery
   for calibrated params. 30 members split across the 2 local V100 GPUs.
2. **Deliverable runs (per site, single column):** prior/posterior continuous runs over
   each site's full native met window with spin-up, reusing the Phase 1a
   `run_callmip_continuous.jl` machinery generalized per site. Keeps per-site time axes
   protocol-clean and avoids padding artifacts in outputs.
3. **Stretch goal (only if profiling demands):** sites×members mega-batch
   (21×30 = 630 columns, one run per EKI iteration) — requires per-column parameter
   Fields (issue #1723 step 3). Not attempted initially.

### Key design constraints (from Kevin Phan, 2026-08-17, and code inspection)

- Multi-column `TimeVaryingInput` requires **one shared time axis** for all columns
  (matrix of size n_sites × n_times). Site met windows differ (intersection is only
  2009–2010) → **pad each site's row with its own recycled met climatology** on a union
  axis; compute obs misfit only inside each site's real obs window. Padding-by-recycling
  matches the protocol Section 8 spin-up/transient convention (cycled site met).
- Very little error handling in ClimaUtilities for the matrix input → **verify forcing
  row order == column order** explicitly (assertions + spot checks vs ncdump).
- Diagnostics: **only `DictWriter`** works on `ColumnEnsemble` (NetCDFWriter broken);
  postprocessing must be updated for the extra columns.
- Some scalars must become per-column `ClimaCore.Fields.Field`s where sites differ
  (`atmos_h` already done on the kp branch; expect more — e.g. site parameters).
- Compare single-`Column` vs `ColumnEnsemble` runs with
  `experiments/integrated/era5/comparison_utils.jl`; expect small FP differences
  (Cartesian vs Spherical geometry).
- UTC offsets are **not** in the PLUMBER2 met files — derive from longitude
  (`round(lon/15)`) and curate in a site table (verify vs known sites, e.g. DK-Sor=+1,
  US-NR1=−7).

### Decisions taken (with rationale)

| # | Decision | Rationale |
|---|----------|-----------|
| D1 | Prior model = **default `LandModel`** with global-map (CLM/SoilGrids) parameters at each tower's coordinates | Honest "out-of-the-box prior" per protocol §6.1; scales to 31 sites; per-site Table-1 params would require scalar→Field work ×21. Documented deviation: Phase 1a DK-Sor used curated site params, so Phase 1b DK-Sor prior ≠ Phase 1a prior. |
| D2 | Calibrated parameters = the **11 carbon params from Phase 1a `priors.jl`** | Comparability with Phase 1a; they demonstrably fix NEE (RMSE 10.9→0.97 at DK-Sor). Revisit after prior-skill assessment if grass/savanna sites are pathological. |
| D3 | Calibration window = **union axis with per-site recycled-met padding**, obs masked to each site's real window | Intersection (2009–2010) wastes ~80% of obs. Padding = protocol-consistent cycling. |
| D4 | Calibration objective = **monthly means** of NEE/Qle/Qh per site (site-major blocks), minibatched over year-windows | Matches Phase 1a (worked); daily misfit is noisy and 100× larger G vector. |
| D5 | Spin-up: calibration loop uses short water/T spin-up + prescribed carbon ICs (protocol §8 Note 1 allows prescribed ICs); deliverable runs use the 5-yr water/T spin-up + full transient | Full 1850 carbon equilibrium per member×iteration is computationally prohibitive; documented deviation, same treatment for prior and posterior. |
| D7 (revised by Renato 2026-08-18) | Phase 1b LAI = **MODIS at every site**, resolved per file via the `source` attribute (`modis_lai_var`) — NOT by variable name, since 13/21 sites keep MODIS under `LAI` and 8 under `LAI_alternative` | Renato's call ("I rather keep MODIS everywhere"): one consistent sensor across sites and consistency with Phase 1a DK-Sor. Both products are 100% valid everywhere, so no availability cost. |
| D6 | Flux obs = **refreshed artifact** `9b4bd0d8…` (callmip-org/Phase1 @ d6d36ad) | Published artifact `c8014d3…` is stale: contains AU-How (validation site), missing US-Var (calibration site). See DATA_MANIFEST.md. |

## Dependency state (checked 2026-08-18)

| Dependency | State | Action |
|---|---|---|
| ClimaCore #2525 (multi-column space) | MERGED, v0.15.1+ (latest 0.15.3) | pinned via kp branch Project.toml (`= 0.15.1`) |
| ClimaUtilities #248 (matrix TVI) | **OPEN** | `.buildkite` Manifest pins branch `kp/point-cloud-intp` |
| ClimaLand #1826 (ColumnEnsemble) | **OPEN** | our branch is based on it (`kp/diff-multi-cols`) |
| ClimaArtifacts #172 (callmip artifacts) | MERGED but flux artifact stale | refreshed copy + Overrides (D6); upstream refresh = follow-up |

Risk: both OPEN PRs may be rebased/changed on merge → periodically diff upstream;
keep our changes additive and isolated to `experiments/` where possible.

## Phases and gates

Do **not** start phase N+1 until phase N's gate is checked and logged in WORKLOG.md.

- [ ] **Phase 0 — Environment + data + branch.** Verify data (done, see DATA_MANIFEST),
  branch `rb/callmip-phase1b` with Phase 1a pipeline + diagnostics merged in,
  artifact overrides, `.buildkite` env instantiated (CPU+CUDA), **gate:** Kevin's
  2-site demos (`callmip_forcing_demo.jl`, `callmip_run_demo.jl`) run clean on CPU
  and on GPU with sane per-column values.
- [ ] **Phase 1 — Correctness harness.** era5 `column_ensemble_comparison.jl` passes
  locally; DK-Sor: single `Column` vs 2-col ensemble compared via `comparison_utils.jl`;
  forcing row-order spot checks (driver values per column vs ncdump at known datetimes).
  **Gate:** documented tolerance report in WORKLOG; any scalar→Field gaps found are fixed
  or explicitly deferred with justification.
- [ ] **Phase 2 — 21-site forcing.** Site table (coords, UTC offset, met path, windows);
  union-axis padding with recycled climatology; per-site LAI (`LAI_alternative` → fallback
  `LAI`); TRENDY CO2 splice for padded segments. **Gate:** 21-column ensemble builds and
  integrates ≥5 days on GPU; per-column drivers spot-checked at 3 dates × 5 sites.
- [ ] **Phase 3 — Benchmark.** SYPD vs #columns (1/2/5/10/21; duplicated 30/100/300)
  on V100 vs CPU. **Gate:** wall-time model for a 30-member × 10-iteration calibration
  ≤ a few days on 2 GPUs, or fallback plan (CPU pmap / fewer members) chosen. Numbers
  posted to PR #1826 (kmdeck asked).
- [ ] **Phase 4 — Multisite calibration wiring.** `generate_observations.jl` → 21-site
  site-major monthly blocks with masks; G(θ) extraction from DictWriter per column;
  `run_calibration.jl` → member scheduling on 2 GPUs. **Gate:** 2-member × 1-iteration
  × 1-year smoke run: finite G, plausible misfit, EKP update step executes.
- [ ] **Phase 5 — Production calibration.** UTKI, 30 members, ≤10 iterations,
  checkpointed each iteration. **Gate:** converged (DataMisfitController or plateau),
  posterior physically plausible, misfit ↓ vs prior at ≥ ~2/3 of sites.
- [ ] **Phase 6 — Calibration deliverables.** Per-site continuous prior/posterior runs
  (full native window + spin-up), Phase 1b NetCDF writer + format tests.
  **Gate:** 42 files pass `test_callmip_output.jl` (extended for Phase 1b) vs template.
- [ ] **Phase 7 — Validation sites.** ME-org met download for 10 sites (**needs Renato's
  ME-org login**), site table extension, prior+posterior runs. **Gate:** 20 Val files pass
  format tests.
- [ ] **Phase 8 — Upload + methods.** ME-org uploads (Calibration and Validation are
  separate experiments!); METHODS document (all deviations: prescribed carbon ICs,
  padding scheme, LAI source per site, error model, D1–D6).

## Site lists

**Calibration (21):** CA-Qfo, CH-Dav, DE-Gri, DE-Hai, DE-Tha, DK-Sor, FI-Hyy, FR-Pue,
IT-Lav, IT-MBo, IT-Noe, NL-Loo, RU-Fyo, US-MMS, US-NR1, US-SRG, US-SRM, US-Ton, US-Var,
US-Whs, US-Wkg.
**Validation (10):** AU-ASM, AU-How, AU-Stp, AU-Tum, CH-Cha, US-FPe, US-GLE, US-Ha1,
US-Me2, US-UMB. (Met only via ME-org "Phase 1b Validation" experiment; no flux files.)

## Protocol quick-reference (v2)

- Output: 14 daily variables (Table 3, CMIP names; `nep` = −NEE, positive into land,
  kgC m-2 s-1), dims `(time, lat, lon)`, lat/lon scalar variables, time
  `days since <met-start-year>-01-01 00:00:00`, steps 1,2,3,… or 0.5,1.5,….
- Global attrs: `Model=ClimaLand.<version>`, `CalLMIP_Phase=Phase1b`,
  `Calibration_Scenario=Scen1`, `Calibration_stage=Prior|Posterior`,
  `Cal_Val=Calibration|Validation`.
- New in Phase 1b: top-10 cm soil moisture (`mrsos`).
- Optional: ensemble outputs / uncertainties for optimized variables (we do EKI ensemble).
- Spin-up (§8): cycle site met to 1850 carbon equilibrium (CO2 287.42 ppm), transient
  1851→first met year w/ TRENDY CO2; **prescribing ICs is explicitly allowed** (Note 1).
- CO2: TRENDY file `CO2_1700_2024_TRENDYv2025.txt` (in `callmip_phase1` artifact).
- Flux obs: daily, non-gapfilled, with `_uc` uncertainties; days >25% gapfilled are NaN;
  combine measurement + model error ourselves.
