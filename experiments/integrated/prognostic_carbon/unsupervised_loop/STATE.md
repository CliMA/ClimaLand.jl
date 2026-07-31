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

- **status:** **done** (2026-07-31) — harness ported, full 20-site battery green
  under BOTH LAI models, baselines recorded and committed
- **ported from:** `origin/ar/derecho_loop` (`single_site.jl`, `run_battery.pbs`,
  `test_sites.csv`)
- **lives at:** `experiments/integrated/prognostic_carbon/harness/`
  (`site_driver.jl`, `run_battery.pbs`, `test_sites.csv`) — deliberately NOT
  `experiments/integrated/era5/single_site.jl`, which is this branch's
  hard-coded desert experiment and is left untouched
- **PBS job ids:** 6967513 (4-site smoke test, PASSED 4/4, Exit_status=0, 7 min);
  6967717 (20-site, prescribed LAI) and 6967718 (20-site, prognostic Zhou LAI),
  both running
- **baseline recorded (GPP / LAI / Ra / Rh per site):** **yes, both** —
  prescribed 20/20 (job 6967717, `harness/baseline_prescribed_lai.tsv`) and
  prognostic 20/20 (job 6967718, `harness/baseline_prognostic_lai.tsv`)

### Stage-0 baseline, PRESCRIBED MODIS LAI (job 6967717, 20/20 PASS)

2000-09-01 → 2002-09-01, first year excluded. Fluxes g C m⁻² day⁻¹.
Committed in full at `harness/baseline_prescribed_lai.tsv`. NPP/GPP = 1 − Ra/GPP.

| site | biome | GPP | Ra | Rh | LAI | NPP/GPP |
|---|---|---|---|---|---|---|
| congo_basin | tropical rainforest | 8.71 | 4.01 | 1.69 | 5.16 | 0.540 |
| amazon_central | tropical rainforest | 8.16 | 3.81 | 1.87 | 4.79 | 0.534 |
| borneo | tropical rainforest | 7.44 | 3.58 | 0.92 | 4.09 | 0.519 |
| central_europe | temperate deciduous | 4.94 | 2.38 | 0.57 | 1.83 | 0.519 |
| n_australia_savanna | tropical savanna C4 | 4.42 | 2.20 | 1.12 | 1.17 | 0.502 |
| california_vaira | mediterranean | 3.83 | 2.00 | 0.58 | 1.74 | 0.478 |
| cerrado_brazil | tropical savanna C4 | 3.77 | 2.08 | 1.76 | 1.60 | 0.447 |
| ozark_us | temperate deciduous | 3.25 | 1.84 | 0.56 | 1.12 | 0.433 |
| pampas_argentina | temperate grassland | 3.10 | 1.80 | 0.97 | 1.03 | 0.419 |
| fennoscandia | boreal forest | 2.83 | 1.66 | 0.60 | 1.32 | 0.413 |
| iberia | mediterranean | 2.58 | 1.58 | 0.76 | 0.74 | 0.386 |
| us_great_plains | temperate grassland C4 | 2.08 | 1.45 | 0.52 | 0.66 | 0.306 |
| central_siberia | boreal forest | 1.61 | 1.27 | 0.32 | 1.28 | 0.207 |
| ne_china | temperate deciduous | 1.53 | 1.26 | 0.53 | 0.66 | 0.174 |
| canada_boreal | boreal forest | 1.37 | 1.20 | 0.62 | 0.71 | 0.126 |
| sahel | tropical savanna C4 | 1.20 | 1.09 | 0.55 | 0.37 | 0.089 |
| mojave_sw_us | desert | 1.04 | 1.07 | 0.17 | 0.16 | **−0.024** |
| alaska_north_slope | tundra | 0.98 | 1.07 | 0.00 | 0.40 | **−0.093** |
| sahara | desert | 0.00 | 0.81 | 0.07 | 0.00 | n/a |
| arabian | desert | 0.00 | 0.81 | 0.11 | 0.00 | n/a |

### Stage-0 baseline, PROGNOSTIC Zhou LAI (job 6967718, 20/20 PASS)

Full table at `harness/baseline_prognostic_lai.tsv`. Same window and units.
Compared with prescribed MODIS LAI above:

| site | LAI pres | LAI prog | GPP pres | GPP prog | NPP/GPP prog |
|---|---|---|---|---|---|
| congo_basin | 5.16 | 4.63 | 8.71 | 8.53 | 0.539 |
| amazon_central | 4.79 | 4.30 | 8.16 | 8.00 | 0.532 |
| borneo | 4.09 | 3.98 | 7.44 | 7.47 | 0.519 |
| central_europe | 1.83 | 2.23 | 4.94 | 5.15 | 0.521 |
| n_australia_savanna | 1.17 | **4.52** | 4.42 | 7.83 | 0.514 |
| iberia | 0.74 | **2.03** | 2.58 | 4.44 | 0.500 |
| cerrado_brazil | 1.60 | 3.05 | 3.77 | 5.17 | 0.495 |
| pampas_argentina | 1.03 | 2.44 | 3.10 | 4.70 | 0.494 |
| ozark_us | 1.12 | 2.52 | 3.25 | 4.69 | 0.493 |
| us_great_plains | 0.66 | **2.31** | 2.08 | 3.96 | 0.463 |
| ne_china | 0.66 | 1.66 | 1.53 | 3.16 | 0.425 |
| california_vaira | 1.74 | **0.82** | 3.83 | 2.91 | 0.417 |
| fennoscandia | 1.32 | 1.23 | 2.84 | 2.63 | 0.391 |
| canada_boreal | 0.71 | 1.24 | 1.37 | 1.99 | 0.293 |
| mojave_sw_us | 0.16 | 0.22 | 1.04 | 1.62 | 0.212 |
| central_siberia | 1.29 | 0.72 | 1.61 | 1.17 | **0.028** |
| alaska_north_slope | 0.40 | **0.14** | 0.98 | 0.49 | **−0.915** |
| sahel | 0.37 | **0.00** | 1.20 | **0.00** | n/a |
| sahara / arabian | 0.00 | 0.00 | 0.00 | 0.00 | n/a |

The two LAI models disagree substantially — prognostic runs much higher at
grassland/savanna/temperate sites and lower at `california_vaira`,
`central_siberia`, `alaska_north_slope`, and collapses `sahel` to bare ground.
That is the LAI model's own behaviour and **out of scope under rule 1**; it is
recorded because the carbon model must work under both, and the pools will
differ correspondingly. Do not "fix" it.

Under prognostic LAI the respiration problem persists at the same cold/arid
sites and is worse at `alaska_north_slope` (NPP/GPP −0.915, GPP halved to 0.49
while Ra stays 0.93), with `central_siberia` marginal at 0.028. The middle of
the distribution is healthier than under prescribed LAI (13 sites in 0.4–0.55).

**The current model already fails the stage-2 acceptance criterion.** That
criterion is "Ra neither collapses to zero nor exceeds GPP in the annual mean at
any site": Ra exceeds GPP at `mojave_sw_us` and `alaska_north_slope`, and is
positive with GPP = 0 at both true deserts. Four of 20 sites, all cold or arid.
Of the forest/grassland sites stage 2 must land in 0.3–0.6, three boreal/cold
ones miss low (`central_siberia` 0.21, `ne_china` 0.17, `canada_boreal` 0.13).

### Smoke-test baseline (job 6967513, 1 yr from 2000-09-01, no spinup excluded)

Fluxes in g C m⁻² day⁻¹, annual mean.

| site | biome | GPP | Ra | Rh | LAI |
|---|---|---|---|---|---|
| amazon_central | tropical rainforest | 8.83 | 4.08 | 1.60 | 4.78 |
| ozark_us | temperate deciduous | 3.75 | 2.04 | 0.49 | 1.13 |
| us_great_plains | temperate grassland C4 | 2.67 | 1.69 | 0.45 | 0.72 |
| sahara | desert | 0.00 | **1.00** | 0.06 | 0.00 |

Implied NPP/GPP = 1 − Ra/GPP: 0.54 amazon, 0.46 ozark, 0.37 great plains —
all inside the 0.3–0.6 band stage 2 must land in, so the *current* JULES Ra is
not wildly off where there is vegetation. The desert is the problem (below).
- **notes:** driver is ENV-parametrized (`SITE_LON/SITE_LAT/SITE_NAME/START/
  STOP/DT/SITE_OUTDIR/LAI_MODE/SPINUP_YEARS`). `LAI_MODE=prescribed` (MODIS,
  default, stages 1–3) or `prognostic` (Zhou, stage 4+). Writes
  `carbon_metrics.txt` per site; `run_battery.pbs` collects
  `baseline_summary.tsv`. GPP/Ra/Rh reported in g C m⁻² day⁻¹.

## Stage 1 — carbon pools in `biomass.jl`

- **status:** not_started — **unblocked, this is the next work**
- **gate satisfied:** stage 0 done; both baselines recorded, so rule 1 ("GPP and
  LAI unchanged") is now checkable against
  `harness/baseline_prescribed_lai.tsv` and `harness/baseline_prognostic_lai.tsv`
- **rule-1 check procedure:** rerun the battery under each LAI mode with the
  carbon model present and diff GPP and LAI per site against those two files.
  Any difference beyond round-off is a stage-1 bug by definition.
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
| 6967513 | 0 | 4-site smoke test of the ported harness, 1 yr, prescribed LAI | 2026-07-31 | **F, exit 0** | PASS 4/4 in 7 min; baseline table above | `.../battery_6967513.desched1/` |
| 6967717 | 0 | 20-site baseline, 2 yr, **prescribed** MODIS LAI, 1 yr spinup excluded | 2026-07-31 | **done** | PASS 20/20; table above; committed as `harness/baseline_prescribed_lai.tsv` | `.../battery_pc_base_pres_6967717.desched1/` |
| 6967718 | 0 | 20-site baseline, 2 yr, **prognostic** Zhou LAI, 1 yr spinup excluded | 2026-07-31 | **done** | PASS 20/20; table above; committed as `harness/baseline_prognostic_lai.tsv` | `.../battery_pc_base_prog_6967718.desched1/` |

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
- **2026-07-31, ROOT CAUSE: `autotrophic_respiration_Q10 = 1.0`, so maintenance
  respiration has no temperature response at all.** `Rpm = f_T·(R_leaf +
  Rd_ref·RAI + Rd_ref·μs·SAI)` with `f_T = Q10^((T−T_ref)/10)`, and the shipped
  default `Q10` is **1.0** — its own toml description says "fixed at 1.0 (not
  calibrated), keeping the maintenance baseline seasonally flat". So `f_T ≡ 1`.
  Measured: across `alaska_north_slope` (262 K), `mojave_sw_us` (281 K) and
  `sahel` (300 K), a 38 K span over which `Q10 = 2` would give a 14× range, Ra
  is 1.068 / 1.066 / 1.094 g C m⁻² day⁻¹ — constant to ±3%. At the two zero-LAI
  deserts Ra = `Rd_ref` × 1.0000085, i.e. exactly `Rd_ref` with `f_T = 1` and
  `RAI + μs·SAI = 1`.
  **Consequence for stage 2:** MODEL.md §2.1's `Rm` reuses this same `Q10`.
  Written assuming the temperature response is "already there", stage 2 would
  ship a temperature-insensitive `Rm` that reads correctly and is inert in
  practice. MODEL.md §8 claimed "2.0 (existing)" and has been corrected;
  setting `Q10 = 2.0` is a deliberate stage-2 change to report, not a silent fix.
- **2026-07-31, never call `qstat` inside a monitor or any wrapping script.**
  A monitor that tested `! qstat <jid>` for liveness reported `JOB_GONE` for
  6967718 while a bare `qstat` showed it Running at 33 min. `qstat` is in the
  sandbox's excludedCommands; wrapping it forces it back into the sandbox where
  it fails, and a failed `qstat` is indistinguishable from a finished job.
  Monitors must key on files the job writes; do liveness checks with bare
  `qstat` at the start of an iteration.
- **2026-07-31, the desert respires carbon it never fixed — quantified.** At
  `sahara`, GPP = 0.00 and LAI = 0.00, yet **Ra = 1.00 g C m⁻² day⁻¹**
  (≈ 365 g C m⁻² yr⁻¹). This is the JULES-style `AutotrophicRespirationModel`
  respiring *prescribed constant* area indices (`Rd_ref·RAI`,
  `Rd_ref·μs·SAI`, MODEL.md §2.1): SAI and RAI do not go to zero in a desert,
  so the plant pays maintenance forever with no photosynthesis and no pool to
  draw down. The carbon budget is open by construction. This is the strongest
  quantitative case for stage 2 (respire the pools, not the area indices), and
  it is also a prediction to check: after stage 2, Sahara Ra should fall to
  ≈ 0 because the pools there equilibrate near zero.
- **2026-07-31, per-job battery configs.** `run_battery.pbs` now sources
  `$PC_SCRATCH/battery_${PBS_JOBNAME}.conf` (falling back to `battery.conf`),
  and `RUNROOT` includes the job name. One shared conf file could not support
  two concurrent batteries — the second submission would silently run the
  first's configuration. Submit with `qsub -N pc_<name>`.
- **2026-07-31, PBS output file naming.** With `-j oe` and `-o <directory>`,
  Derecho writes `<jobid>.OU` (e.g. `6967513.desched1.OU`) into that directory,
  NOT `<jobname>.o<jobid>`. A monitor armed on the latter never fires and looks
  exactly like a still-running job. Watch `$PC_SCRATCH/<jobid>.desched1.OU`.
- **2026-07-31, `qstat` right after `qsub`.** A job id returned by `qsub` can be
  briefly absent from both `qstat -u` and `qstat -x -f <jid>` ("Unknown Job
  Id"). It is a propagation delay, not a failed submission — do not resubmit.
  Re-check before concluding anything.
- **2026-07-31, diagnostic name matching.** `output_short_name` is
  `"<short>_<schedule>_<suffix>"` (e.g. `gpp_1d_average`). The driver anchors
  its lookup with `startswith(name, short * "_")` rather than the original's
  `occursin`, so a short name can never cross-match another variable.

## Reporting — PR #1834 comment hygiene

- **Pinned status comment id `5145422196`.** This is the "opening comment": a
  live dashboard of stage status, running jobs and blockers. EDIT it every
  iteration (`gh pr comment --edit-last` will not find it — use
  `gh api -X PATCH repos/CliMA/ClimaLand.jl/issues/comments/5145422196 -f body=@file`).
  Never delete it, and never let it go stale.
- Iteration reports are posted as new comments below it. Keep at most ~10:
  before posting, list comments and delete the oldest surplus ones with
  `gh api -X DELETE repos/CliMA/ClimaLand.jl/issues/comments/<id>`, never
  touching `5145422196`.

## Blockers

None recorded.
