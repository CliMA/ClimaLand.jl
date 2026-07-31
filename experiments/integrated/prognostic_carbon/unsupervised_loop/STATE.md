# Prognostic carbon — pipeline state

Single source of truth for where the work is. The loop reads this first every
iteration and commits it every iteration. See `loop_prompt.md` for what each
stage means and how it is gated, and `../MODEL.md` for the physics.

Status values: `not_started` | `in_progress` | `done` | `failed` | `blocked`

## Fixed configuration

- branch `ar/prognostic_carbon`, PR #1834, base `ar/climate_responsive_lai_inputs`
- PBS: account UCIT0011, queue `develop`, job names prefixed `pc_`
- scratch: `/glade/derecho/scratch/arenchon/claude/prognostic_carbon`
- LAI: prescribed MODIS for stages 1–3, prognostic Zhou from stage 4
- spinup: recycle ERA5 1979–1989 with CO₂ held at ~340 ppm
- `gh` requires `module load gh`

## `~/.julia/artifacts/Overrides.toml` — last observed (2026-07-31)

Shared with the other loop; this loop must never write it. Re-read and compare
every iteration.

```
2ce5e5e05e4f17c86b07eadff1f9dd551e779524 = "/glade/derecho/scratch/arenchon/artifacts_local/optimal_lai_inputs"
f269a0b057b9f438b4caafdef17da73746310787 = "/glade/campaign/univ/ucit0011/ClimaArtifacts2/artifacts/forty_yrs_era5_land_forcing_data/forty_yrs_era5_land_forcing_data_artifact"
```

---

## Stage 0 — test harness (single-pixel ERA5 battery)

- **status:** in_progress
- **ported from:** `origin/ar/derecho_loop` (`single_site.jl`, `run_battery.pbs`,
  `test_sites.csv`)
- **lives at:** `experiments/integrated/prognostic_carbon/harness/`
  (`site_driver.jl`, `run_battery.pbs`, `test_sites.csv`) — deliberately NOT
  `experiments/integrated/era5/single_site.jl`, which is this branch's
  hard-coded desert experiment and is left untouched
- **PBS job ids:** 6967513 (4-site smoke test, submitted 2026-07-31)
- **baseline recorded (GPP / LAI / Ra / Rh per site):** no — pending 6967513,
  then a full 20-site run
- **notes:** driver is ENV-parametrized (`SITE_LON/SITE_LAT/SITE_NAME/START/
  STOP/DT/SITE_OUTDIR/LAI_MODE/SPINUP_YEARS`). `LAI_MODE=prescribed` (MODIS,
  default, stages 1–3) or `prognostic` (Zhou, stage 4+). Writes
  `carbon_metrics.txt` per site; `run_battery.pbs` collects
  `baseline_summary.tsv`. GPP/Ra/Rh reported in g C m⁻² day⁻¹.

## Stage 1 — carbon pools in `biomass.jl`

- **status:** not_started
- **pools:** C_sugar, C_leaf, C_stem, C_root
- **conservation test written:** no
- **GPP/LAI unchanged vs stage-0 baseline:** —
- **PBS job ids:** —
- **notes:** —

## Stage 2 — pool-based autotrophic respiration

- **status:** not_started
- **annual NPP/GPP by site:** —
- **PBS job ids:** —
- **notes:** —

## Stage 3 — prognostic SOC

- **status:** not_started
- **`dY.soilco2.SOC` wired to litter − Sm:** no
- **soil C conservation test passes:** no
- **Rh vs stage-0 baseline:** —
- **PBS job ids:** —
- **notes:** —

## Stage 4 — offline spinup to steady state

- **status:** not_started
- **driver dump path (scratch):** —
- **offline integrator validated against coupled run:** no
- **IC file path:** —
- **artifact override verified to resolve:** no
- **notes:** —

## Stage 5 — global run, comparison, tuning

- **status:** not_started
- **reference biomass dataset:** — (TBD, see open questions)
- **PBS job ids:** —
- **cVeg RMSE / bias by zone:** —
- **calibrated parameters:** —
- **committed to `toml/default_parameters.toml`:** no
- **notes:** —

---

## Job table

| job_id | stage | purpose | submitted | status | result | output path |
|---|---|---|---|---|---|---|
| 6967513 | 0 | 4-site smoke test of the ported harness (amazon_central, sahara, ozark_us, us_great_plains), 1 yr, prescribed LAI | 2026-07-31 | submitted | pending | `/glade/derecho/scratch/arenchon/claude/prognostic_carbon/battery_6967513.desched1/` |

Jobs seen in `qstat -u arenchon` that are NOT `pc_*` belong to the other loop —
never reconcile or `qdel` them. Observed 2026-07-31: `clima_cal*` (6967486).

## Decisions & findings

- **2026-07-31, harness location.** The ported harness lives in
  `experiments/integrated/prognostic_carbon/harness/`, not in
  `experiments/integrated/era5/single_site.jl`. That file on this branch is the
  hard-coded desert experiment, not the parametrized driver; the parametrized
  one only exists on `ar/derecho_loop`. Neither CI nor the docs reference it, so
  a separate directory keeps the carbon harness self-contained and leaves the
  existing experiment working.
- **2026-07-31, dropped ENV switches.** `BETA_IN_A0`, `ONLINE_F0`, `F0_SCALE`,
  `ONLINE_VPD_GS`, `ONLINE_GSL` were removed from the port: the corresponding
  `optimal_lai_*` parameters are absent from `toml/default_parameters.toml` on
  this branch (that behaviour is now unconditional). Keeping them would write
  overrides for parameters nothing reads — inert, and silently so. Surviving
  overrides: `optimal_lai_z`, `z_c4`, `sigma`, `sigma_c4`, `alpha`, `f0`,
  `z_a0`, `online_c3c4`.
- **2026-07-31, `module load gh` must be run bare.** `module` is a shell
  function; piping it into anything runs it in a subshell and the `PATH` update
  is lost, which presents exactly as "gh: command not found". Run bare, `gh` is
  present (2.74.2) and authenticated as AlexisRenchon.
- **2026-07-31, diagnostic name matching.** `output_short_name` is
  `"<short>_<schedule>_<suffix>"` (e.g. `gpp_1d_average`). The driver anchors
  its lookup with `startswith(name, short * "_")` rather than the original's
  `occursin`, so a short name can never cross-match another variable.

## Blockers

None recorded.
