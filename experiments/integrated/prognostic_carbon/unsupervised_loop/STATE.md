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

- **status:** not_started
- **ported from:** `origin/ar/derecho_loop` (`single_site.jl`, `run_battery.pbs`,
  `test_sites.csv`)
- **PBS job ids:** —
- **baseline recorded (GPP / LAI / Ra / Rh per site):** no
- **notes:** —

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
| — | — | — | — | — | — | — |

## Decisions & findings

None recorded.

## Blockers

None recorded.
