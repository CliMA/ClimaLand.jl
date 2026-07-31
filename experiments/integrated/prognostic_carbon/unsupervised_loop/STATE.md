# Prognostic carbon — pipeline state

Single source of truth for where the work is. The loop reads this first every
iteration and commits it every iteration. See `loop_prompt.md` for what each
stage means and how it is gated, and `../MODEL.md` for the physics.

Status values: `not_started` | `in_progress` | `done` | `failed` | `blocked`

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
