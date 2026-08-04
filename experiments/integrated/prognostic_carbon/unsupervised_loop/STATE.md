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

## ⚠️ OVERRIDES.TOML CHANGED BY THE OTHER LOOP (2026-08-01 ~13:07)

The shared override for `optimal_lai_inputs` was repointed, **not by this loop**:

```
was:  /glade/derecho/scratch/arenchon/artifacts_local/optimal_lai_inputs
now:  /glade/derecho/scratch/arenchon/artifacts_local/optimal_lai_inputs_calibrated
```

Both directories exist; the files **differ** (3126182 vs 3126269 bytes). The
calibrated one was written minutes before this iteration.

### Impact — read before running anything with prognostic LAI

- **Prescribed-LAI results are unaffected.** Stages 0–3 prescribed batteries use
  MODIS LAI and never touch this artifact.
- **All prognostic-LAI results to date used the OLD artifact** — including
  `harness/baseline_prognostic_lai.tsv`, the driver records from job 6971245,
  every offline equilibrium, and the multi-product observational comparison.
  They are internally consistent, because one artifact was in force throughout.
- **The next prognostic-LAI battery will use the calibrated artifact and produce
  different LAI, hence different GPP.**

### The trap this sets

`check_rule1.jl` compares a new run's GPP and LAI against
`baseline_prognostic_lai.tsv`. Run now, it would report a **rule-1 FAILURE that
is really an artifact change** — the LAI model's inputs moved, not the carbon
model's behaviour. Rule 1 is the project's central claim, so a false failure
there is expensive.

**Before any further prognostic-LAI comparison: re-baseline.** Run a CARBON=0
prognostic battery against the calibrated artifact and regenerate
`baseline_prognostic_lai.tsv`, or state explicitly which artifact each side of
the comparison used.

**Re-baseline DONE (2026-08-01, job 6976878, 20/20).**
`harness/baseline_prognostic_lai.tsv` is now the **calibrated-artifact**
baseline; the previous one is kept as
`harness/baseline_prognostic_lai_OLDARTIFACT.tsv` for provenance.

**The artifact change was material — this justifies the STEP 0 discipline
concretely.** LAI moved at 5 of 6 spot-checked sites:

| site | LAI old → new | change |
|---|---|---|
| us_great_plains | 2.306 → 1.813 | **−21%** |
| congo_basin | 4.627 → 4.007 | **−13%** |
| ozark_us | 2.522 → 2.422 | −4% |
| amazon_central | 4.297 → 4.168 | −3% |
| central_siberia | 0.715 → 0.699 | −2% |

Run against the stale baseline, `check_rule1.jl` would have reported violations
of up to 21% — a catastrophic-looking failure of the project's central claim,
caused entirely by someone else's artifact.

**Note for later:** the calibrated LAI is *lower* at `us_great_plains` (−21%),
which feeds less GPP into allocation. Some of the grassland over-build may
shrink under the new artifact — worth re-measuring the equilibria rather than
assuming the earlier numbers carry over.

**Resolved (2026-08-01).** Jobs **6976878** (`pc_rb_base`, CARBON=0) and
**6976879** (`pc_rb_carb`, CARBON=1 CARBON_RA=1), both prognostic LAI on the
calibrated artifact, completed 20/20. Diffed against each other — both sides on
the same artifact, the only comparison that means anything now:

- **`RULE1 PASS` at all 20 sites, `rel_diff = 0.0` bit-exactly.** Against the
  stale baseline the same check reported up to **33% apparent violation** at
  pampas LAI. That contrast is the evidence for the STEP 0 discipline.
- `harness/baseline_prognostic_lai.tsv` is 6976878's summary; the old-artifact
  version is kept as `baseline_prognostic_lai_OLDARTIFACT.tsv`.

### Observational comparison moved with the artifact: 8/20 → **10/20** inside range

The re-measurement anticipated above confirmed the note: the grassland/dry-site
over-build **shrank substantially** under the calibrated LAI, because lower LAI
feeds less GPP into allocation.

| site | old artifact | calibrated | verdict |
|---|---|---|---|
| mojave_sw_us | 2.58 (12.9× above) | 0.50 | 2.5× above |
| iberia | 11.33 (4.5×) | 4.70 | 1.9× above |
| n_australia_savanna | 13.75 (2.0×) | 7.56 | 1.1× above |
| sahel | 0.35 (above) | 0.00 | **inside** |
| california_vaira | 9.74 (above) | 5.45 | **inside** |
| pampas_argentina | 13.55 (11.1×) | 10.63 | 8.7× above |
| ne_china | 11.07 (7.2×) | 9.46 | 6.1× above |
| cerrado_brazil | 12.23 (3.7×) | 10.20 | 3.1× above |
| us_great_plains | 9.45 (1.9×) | 8.70 | 1.7× above |

Tropics and boreal barely moved (all still inside except canada_boreal at 1.1×).
Products still disagree by a median factor of **3.4×**.

**Scientific correction to record:** the earlier characterisation — "the carbon
model over-builds woody biomass in dry and herbaceous systems" — **overstated the
carbon model's share of the error**. A large part of it was LAI-model bias that
the other loop's calibration has now removed. The residual over-build is real but
smaller, and concentrated in two sites: `pampas_argentina` (8.7×) and `ne_china`
(6.1×). Any future claim about allocation must be measured against the
calibrated artifact, never the old numbers.

This loop must **not** write `Overrides.toml` (hard rule), so pinning the old
path back is not an option — and would sabotage the other loop in any case.

## `~/.julia/artifacts/Overrides.toml` — last observed (2026-08-01)

Shared with the other loop; this loop must never write it. Re-read and compare
every iteration.

```
2ce5e5e05e4f17c86b07eadff1f9dd551e779524 = "/glade/derecho/scratch/arenchon/artifacts_local/optimal_lai_inputs_calibrated"
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

- **status:** **DONE** (2026-07-31). All four gate criteria met:
  package tests pass (74); carbon conservation test passes; the 20-site battery
  completes under **both** LAI models (jobs 6969169, 6969170, 20/20 each); and
  **GPP and LAI are bit-identical to the stage-0 baseline at every site under
  both models** (`rel_diff = 0.0`, `check_rule1.jl`)
- **pools:** C_sugar, C_leaf, C_stem, C_root
- **conservation test written:** **yes** —
  `test/standalone/Vegetation/carbon_conservation.jl`, registered in
  `test/runtests.jl`. 56 tests pass (Float32 and Float64).
- **GPP/LAI unchanged vs stage-0 baseline:** unit-test level **yes**
  (bit-identical LAI, GPP and Ra with vs without the carbon model); site-battery
  level not yet run

### Design: a wrapper, not a replacement

`:biomass` is a single component slot, so `PrognosticCarbonModel` **wraps** an
LAI model (`PrescribedBiomassModel` or `ZhouOptimalLAIModel`) and delegates all
area-index and LAI decisions to it. That is what makes rule 1 structural rather
than a thing to remember: the LAI path is the wrapped model's, untouched.

Details worth not rediscovering:

- `update_fractional_c3!` **must** be forwarded to the wrapped model. Without
  it the Zhou C3/C4 competition silently stops being applied, which changes GPP
  and breaks rule 1. There is a forwarding method; do not remove it.
- `height` and `rooting_depth` are mirrored onto the wrapper because much of the
  codebase reads `canopy.biomass.height` / `.rooting_depth` as direct fields.
  Safe in phase 1 because neither is ever mutated.
- The `CanopyModel` prescribed-LAI constructor asserted
  `biomass.plant_area_index.LAI == LAI`, which a wrapper cannot satisfy. It now
  goes through `prescribed_lai_input(biomass)`, which the wrapper forwards.
- Accessors like `get_GPP(p, ...)` must be resolved **outside** `@.` blocks.
  `@. lazy(M_C * get_GPP(p, ...))` broadcasts over `p` itself and throws
  "broadcasting over dictionaries and NamedTuples is reserved".
- The carbon model has its **own** `carbon_Q10 = 2.0`, deliberately separate
  from `autotrophic_respiration_Q10 = 1.0`, so stage 1 gets a sensible
  temperature response without silently changing the legacy scheme. Stage 2
  owns the decision about the legacy parameter.

### Deliberately deferred within stage 1

- **Phenological leaf shedding** (`C_leaf·max(−dLAI/dt,0)/LAI`) is NOT
  implemented; only background `C_leaf/τ_leaf` turnover is. `dLAI/dt` is
  available for the Zhou model but needs finite-differencing for prescribed
  LAI, and a half-implemented version would differ between the two LAI modes.
  Consequence to expect in the battery: deciduous leaf carbon lags LAI through
  autumn.
- **Ra is not yet wired** into `p.canopy.autotrophic_respiration.Ra`. The pools
  are driven by the carbon model's internal Rm/Rg while the legacy JULES scheme
  still drives the canopy fluxes. This is what keeps stage 1 verifiably
  behaviour-neutral; wiring it is stage 2.

### Remaining for the stage-1 gate

1. ~~Wire `CARBON=1` into `harness/site_driver.jl` and add pool diagnostics.~~
   **Done** (commit `9f18825e4`). The driver rebuilds only the canopy, wrapping
   the `LandModel`'s own biomass model and reusing every other component object
   unchanged — reconstructing them would risk diverging from the LandModel
   defaults (the integrated model uses a different soil-moisture stress model
   than the standalone canopy default), and any such divergence would look like
   a rule-1 violation that is really a harness bug.
2. ~~Rule-1 comparison tool.~~ **Done**: `harness/check_rule1.jl` diffs a
   CARBON=1 summary against the committed baseline per site, exits non-zero on
   any difference beyond `1e-10` relative. Verified non-vacuous — a 0.1%
   perturbation to one site's GPP is caught.
3. **Running:** job 6968847 (`pc_s1_carb`), 4 sites, CARBON=1, prescribed LAI,
   2-year window. Then the full 20 under both LAI modes.
4. Check the failure signatures: C_sugar pinned at zero, C_stem unbounded,
   `σl_implied` far outside 0.03–0.1 kg C m⁻² leaf.

### First carbon-model results (job 6968999, 4/4 PASS, 2 yr from empty pools)

**RULE 1 PASSES EXACTLY.** GPP and LAI are bit-identical with the carbon model
active — `rel_diff = 0.0` at every site, not merely within tolerance. The
wrapper design makes the LAI path literally the wrapped model's, so this is the
expected result, but it is now measured rather than assumed.

| site | LAI | C_sugar | C_leaf | C_stem | C_root | cVeg | σl_implied |
|---|---|---|---|---|---|---|---|
| amazon_central | 4.79 | 0.165 | 0.503 | 1.206 | 0.580 | 2.45 | **0.079** |
| ozark_us | 1.12 | 0.093 | 0.271 | 0.545 | 0.315 | 1.22 | **0.164** |
| us_great_plains | 0.66 | 0.064 | 0.178 | 0.350 | 0.211 | 0.80 | **0.190** |
| sahara | 0.00 | **−0.075** | 0 | 0 | 0 | −0.075 | 0 |

Reading these:

- **Sugar is not pinned at zero** at any vegetated site (0.06–0.17 kg C m⁻²).
  The allocation ramp regulates the pool as designed, with no hard clamp.
- **Pool ordering is sensible**: stem > root > leaf > sugar everywhere.
- **cVeg is far below equilibrium** (2.45 vs the 10–20 expected for tropical
  forest). Expected, not a defect: the pools start empty and τ_stem is 30 years,
  so two years fills only a few percent of the stem pool. Stage 4 exists for
  exactly this.
- **σl_implied at amazon is 0.079 kg C m⁻² leaf — inside the expected
  0.03–0.1 band.** This is the headline diagnostic for whether constant
  allocation fractions can work without PFTs, and at the tropical site it does.
- **σl_implied at ozark (0.164) and us_great_plains (0.190) is roughly 2× above
  the band.** Constant fractions over-build leaf carbon relative to LAI at
  temperate and grassland sites. **Treat as preliminary**: the pools are still
  filling, so the ratio is transient. Re-read it after stage 4 spinup before
  concluding the constant-fraction assumption fails.
- **Sahara sugar went negative (−0.075).** Fixed; see below.

### STAGE-1 GATE, prescribed LAI (job 6969169, 20/20 PASS)

**Rule 1 passes exactly at all 20 sites** — GPP and LAI bit-identical to
`baseline_prescribed_lai.tsv`, `rel_diff = 0.0` everywhere.

Pools after 2 years from empty (kg C m⁻²), sorted by cVeg:

| site | biome | LAI | sugar | leaf | stem | root | cVeg | σl_implied |
|---|---|---|---|---|---|---|---|---|
| congo_basin | tropical rainforest | 5.16 | 0.163 | 0.541 | 1.309 | 0.625 | 2.64 | 0.098 |
| amazon_central | tropical rainforest | 4.79 | 0.165 | 0.503 | 1.206 | 0.580 | 2.45 | 0.079 |
| borneo | tropical rainforest | 4.09 | 0.141 | 0.450 | 1.070 | 0.517 | 2.18 | 0.088 |
| central_europe | temperate deciduous | 1.83 | 0.142 | 0.423 | 0.884 | 0.471 | 1.92 | 0.191 |
| california_vaira | mediterranean | 1.74 | 0.101 | 0.330 | 0.747 | 0.374 | 1.55 | 0.163 |
| cerrado_brazil | tropical savanna C4 | 1.60 | 0.071 | 0.258 | 0.631 | 0.301 | 1.26 | 0.234 |
| n_australia_savanna | tropical savanna C4 | 1.17 | 0.081 | 0.354 | 0.294 | 0.522 | 1.25 | 0.395 |
| ozark_us | temperate deciduous | 1.12 | 0.093 | 0.271 | 0.545 | 0.315 | 1.22 | 0.164 |
| pampas_argentina | temperate grassland | 1.03 | 0.069 | 0.238 | 0.603 | 0.287 | 1.20 | 0.245 |
| fennoscandia | boreal forest | 1.32 | 0.091 | 0.244 | 0.493 | 0.269 | 1.10 | 0.129 |
| iberia | mediterranean | 0.74 | 0.065 | 0.212 | 0.486 | 0.245 | 1.01 | 0.249 |
| us_great_plains | temperate grassland C4 | 0.66 | 0.064 | 0.178 | 0.350 | 0.211 | 0.80 | 0.190 |
| central_siberia | boreal forest | 1.28 | 0.059 | 0.138 | 0.274 | 0.151 | 0.62 | 0.117 |
| canada_boreal | boreal forest | 0.71 | 0.049 | 0.127 | 0.273 | 0.142 | 0.59 | 0.127 |
| ne_china | temperate deciduous | 0.66 | 0.059 | 0.127 | 0.259 | 0.146 | 0.59 | 0.106 |
| mojave_sw_us | desert | 0.16 | 0.028 | 0.089 | 0.202 | 0.101 | 0.42 | 0.501 |
| alaska_north_slope | tundra | 0.40 | 0.035 | 0.094 | 0.176 | 0.101 | 0.41 | 0.205 |
| sahel | tropical savanna C4 | 0.37 | 0.027 | 0.099 | 0.105 | 0.147 | 0.38 | 0.169 |
| sahara | desert | 0.00 | 0 | 0 | 0 | 0 | **0** | 0 |
| arabian | desert | 0.00 | 0 | 0 | 0 | 0 | **0** | 0 |

**What passes:**

- **cVeg follows the biomass gradient**: rainforest 2.2–2.6 > temperate forest
  0.6–1.9 > savanna/grassland 0.4–1.3 > desert 0. Right ordering, though far
  below equilibrium magnitudes (2 yr from empty, τ_stem = 30 yr).
- **Sugar is off zero at every vegetated site** (0.027–0.165) and exactly zero at
  both true deserts. No pinning, no phantom carbon, no clamp.
- **No unbounded stem growth** anywhere; max stem is 1.31.
- **Notably, sugar stays positive at `alaska_north_slope` and `mojave_sw_us`**,
  the two sites where the *JULES* Ra exceeded GPP in the stage-0 baseline. The
  pool-based Rm does not reproduce that pathology — an early, encouraging sign
  for stage 2, though the wiring that makes it count has not happened yet.

**The finding: `σl_implied` is systematically too high, and the spread is what
matters.**

Only the three tropical rainforests (0.079, 0.088, 0.098) land inside the
expected 0.03–0.1 kg C m⁻² leaf. Every other site is 1.1× to **5×** above, and
the excess broadly grows as LAI falls. The full spread is
**0.079 → 0.501, a factor of 6.3**.

This is not fixable by recalibrating `f_leaf`. A global change to `f_leaf`
rescales every site by the same factor, so it can shift the distribution but not
compress it: dividing by 2.5 would centre the middle of the range but push the
tropics to 0.032 (bottom edge of the band) and still leave `mojave_sw_us` at
0.20, twice the top.

**The transient argues the wrong way for anyone hoping spinup fixes it.** The
pools start empty, so `C_leaf` is *below* its equilibrium — with τ_leaf = 1.5 yr
and a 2-year run it is at most ~74% equilibrated, and less than that because the
allocation rate is itself still rising. Equilibrium `σl_implied` will therefore
be **higher** than these numbers, not lower. Spinup makes this worse.

Per MODEL.md §1 this is "a real, reportable finding — not something to paper
over", and §2.3's prescribed fallback if constants fail is to let `f_stem` (and
correspondingly `f_root`) vary with mean annual temperature and precipitation,
never with a PFT. **Do not act on this yet**: confirm against the prognostic-LAI
battery and re-measure after stage-4 spinup, since that is the state the
comparison is properly made in.

### STAGE-1 GATE, prognostic Zhou LAI (job 6969170, 20/20 PASS) — and it changes the σl reading

**Rule 1 passes exactly at all 20 sites under the prognostic LAI model too.**
That also demonstrates the wrapper preserves Zhou's nine time-integrated
prognostic variables and keeps forwarding the C3/C4 competition — the two silent
failure modes identified in iterations 4–5.

`σl_implied` under the two LAI models (committed:
`harness/stage1_prescribed_lai.tsv`, `harness/stage1_prognostic_lai.tsv`):

| site | LAI pres | LAI prog | σl pres | σl prog |
|---|---|---|---|---|
| congo_basin | 5.16 | 4.63 | 0.098 | 0.113 |
| borneo | 4.09 | 3.98 | 0.088 | 0.110 |
| amazon_central | 4.79 | 4.30 | 0.079 | 0.106 |
| central_europe | 1.83 | 2.23 | 0.191 | **0.105** |
| ozark_us | 1.12 | 2.52 | 0.164 | **0.088** |
| us_great_plains | 0.66 | 2.31 | 0.190 | **0.098** |
| ne_china | 0.66 | 1.66 | 0.106 | **0.096** |
| pampas_argentina | 1.03 | 2.44 | 0.245 | **0.118** |
| n_australia_savanna | 1.17 | 4.52 | 0.395 | **0.120** |
| cerrado_brazil | 1.60 | 3.05 | 0.234 | **0.142** |
| canada_boreal | 0.71 | 1.24 | 0.126 | **0.074** |
| fennoscandia | 1.32 | 1.23 | 0.129 | **0.088** |
| central_siberia | 1.28 | 0.72 | 0.117 | 0.115 |
| iberia | 0.74 | 2.03 | 0.248 | 0.270 |
| mojave_sw_us | 0.16 | 0.22 | 0.501 | 0.490 |
| california_vaira | 1.74 | 0.82 | 0.163 | **0.395** |
| alaska_north_slope | 0.40 | 0.14 | 0.205 | **0.857** |
| sahel | 0.37 | 0.00 | 0.169 | — (LAI = 0) |

**This materially revises the prescribed-LAI reading recorded above.** Under
prognostic LAI, roughly half the sites land in or just above the expected
0.03–0.1 band, where under prescribed MODIS LAI only the three rainforests did.

The controlling variable is **LAI, not biome**: `σl_implied` is high wherever
LAI is low, and it tracks LAI changes site by site. `us_great_plains` LAI
0.66 → 2.31 takes σl 0.190 → 0.098; `n_australia_savanna` 1.17 → 4.52 takes
0.395 → 0.120. It also moves the other way — `california_vaira` LAI 1.74 → 0.82
takes σl 0.163 → 0.395, and `alaska_north_slope` 0.40 → 0.14 takes 0.205 →
0.857.

So the honest statement is **not** "constant allocation fractions fail". It is:
the carbon model allocates leaf carbon roughly in proportion to productivity,
and where the LAI model puts little leaf area on a productive column the two
disagree. The remaining outliers are all low-LAI columns
(`alaska_north_slope` 0.14, `mojave_sw_us` 0.22, `california_vaira` 0.82).

Since stage 4 onward runs prognostic LAI, **the prognostic column is the one
that matters for the final product**, and it is much closer to defensible than
the prescribed one. Still do not act on it: pools are 2 years from empty and
`C_leaf` is below equilibrium, so these numbers will rise. Re-measure after
stage-4 spinup before touching `f_leaf` or invoking MODEL.md §2.3's
climate-varying fallback.

### Diagnostic bug found and fixed

`σl_implied` at `sahel` under prognostic LAI came out as **9.8e13**: LAI is
clipped to exactly zero below 0.05, and the diagnostic divided by
`max(LAI, eps)`. It now reports zero where there is no leaf area. A site with
zero LAI and non-zero `C_leaf` is a genuine inconsistency — `sahel` has
cVeg = 0.106 with LAI = 0 — but the leaf pool itself is what shows it; the ratio
there is meaningless rather than enormous.

### Sugar floor verified (job 6969094, 4/4 PASS)

The substrate limitation does exactly what it should, and nothing else:

| site | C_sugar before | C_sugar after |
|---|---|---|
| sahara | **−0.0750** | **0** (exactly) |
| amazon_central | 0.16494 | 0.16493 |
| us_great_plains | 0.064274 | 0.064274 |
| ozark_us | 0.092741 | 0.092754 |

The desert floors at zero; vegetated sites move at most in the fifth significant
digit (~1e-4 relative at `ozark_us`). Note that is the effect on the *integrated
pool*, which accumulates through low-sugar periods, and so is larger than the
~1e-5 instantaneous effect on `Rm` at a healthy pool. Rule 1 still passes
exactly.

### Fifth issue: maintenance respiration drove the sugar pool negative

At sahara GPP is zero all year, and `Rm` kept drawing on an empty pool. Negative
carbon is not a physical outcome — respiration requires substrate. `Rm` is now
multiplied by `g(C_sugar/C_sugar_ref)` with `C_sugar_ref = 1e-3 kg C m⁻²`
(commit `dfad8f549`), which is the floor MODEL.md §2.3 already called for.

Distinct from the allocation ramp, which is keyed to the sugar *target*: a
stressed plant with sugar below target still pays full maintenance and draws
down; only an *exhausted* one stops. A healthy pool is unaffected to ~1 part in
10⁵, asserted in the tests against the unthrottled limit.

### Fourth trap: `p.drivers.T_ground` does not exist in the integrated model (2026-07-31)

Job 6968904 crashed at all four sites:

```
FieldError: type NamedTuple has no field `T_ground`, available fields:
`P_liq`, `P_snow`, `T`, `P`, `u`, `q`, `c_co2`, `SW_d`, `LW_d`, `cosθs`, `frac_diff`
```

`p.drivers.T_ground` exists only for a **standalone canopy over prescribed
ground**. The integrated `LandModel` carries no ground temperature in its
drivers, and reaching into `p.soil.T` would make the canopy tendency depend on
soil `update_aux` ordering. Stage 1 therefore uses canopy temperature for the
root maintenance term. **This is a deviation from MODEL.md §2.1**, recorded
there with its consequence (the root term stays more seasonally variable than
soil temperature would make it) and deferred to stage 3, when the canopy→soil
coupling makes a shared ground temperature available. Fixed in `e92726e16`.

**The failure was reported misleadingly**, which cost a cycle: the crash was
caught by ClimaLand's simulation error handler, so the run continued to the
metrics block, which then failed with `KeyError: "gpp_1d_average" not found` —
reading like a missing diagnostic rather than a dead simulation. The driver now
checks for the GPP series and errors with a message pointing back up the log.
**When a site fails, search the log upward for "simulation crashed" before
believing the last error.**

### Pre-existing ClimaLand bug: the `rd` diagnostic has never worked (2026-07-31)

Distinct from the wrapper-dispatch traps below — this one predates the branch.

`define_diagnostics.jl` wires the `rd` (leaf dark respiration) diagnostic to
`compute_respiration_leaf!`, which is **defined nowhere in the package**, so
requesting it dies with `UndefVarError` at the first output step. Meanwhile
`compute_respiration_canopy!`, generated by `@diagnostic_compute` and computing
exactly `get_Rd_canopy` in the declared units, was **never referenced**. A
naming mismatch; nothing had asked for `rd` until the offline spinup needed Rd
as a driver. Fixed by pointing `rd` at the method that exists.

**A scan found three more of the same class, left alone deliberately:**
`compute_aerodynamic_resistance!`, `compute_photosynthesis_net_leaf!` and
`compute_plant_water_content!` are referenced but defined nowhere, and
`compute_specific_humidity!` is another orphaned method. Each needs its own
decision about what it should compute, none blocks this work, and rewiring four
unrelated diagnostics inside a carbon-model PR would be scope creep somewhere a
wrong guess is invisible. Recorded here so they are findable.

### FIFTH trap, same shape: `add_diagnostics!` (2026-07-31)

The stage-4 battery failed **0/20** on the `output_vars` assertion in
`default_diagnostics`. `add_diagnostics!` dispatches on `ZhouOptimalLAIModel`,
so a `PrognosticCarbonModel` wrapping it never matched and the Zhou diagnostics
(`a0d`, `a0a`, `a0c3`, `a0c4`, `pra`, `fc3`) were silently never registered.
Requesting `fc3` for the offline spinup's driver record is what surfaced it.

Also needed: `PrognosticCarbonModel` imported into the `Diagnostics` module —
it was not, so the forwarding method failed to load with `UndefVarError`.

**No earlier result is affected**: `add_diagnostics!` controls only which
diagnostics are *available*, not the physics, and the rule-1 checks used `gpp`
and `lai`, which are unconditional.

**The running tally of this one pattern** — a function dispatching on the
biomass model silently selecting a default once a wrapper exists:

| # | machinery | how it presented |
|---|---|---|
| 1 | `update_fractional_c3!` | silent (caught by construction) |
| 2 | `set_canopy_component_initial_conditions!` | silent (caught by construction) |
| 3 | `initialize_jacobian` | loud, at Jacobian build |
| 4 | `prescribed_lai_input` | loud, at constructor |
| 5 | `add_diagnostics!` | loud, 0/20 battery after 40 min |

The test now asserts the wrapper exposes every diagnostic the wrapped model
does, so the next omission fails in the suite rather than in a battery.

### Third trap: new prognostic variables need Jacobian entries (2026-07-31)

Job 6968847 failed at **all four sites**:

```
The linear system cannot be solved because A does not have any entries at the
following keys: (canopy.biomass.C_sugar, canopy.biomass.C_sugar), (C_leaf,
C_leaf), (C_stem, C_stem), (C_root, C_root)
```

Every prognostic variable needs a Jacobian block in
`src/shared_utilities/implicit_timestepping.jl`. The pools are stepped
explicitly, so they belong in `explicit_vars` next to the optimal-LAI
time-integrated variables. Fixed in `341d75f55`.

**Why the unit tests could not have caught it:** they call
`make_compute_exp_tendency` directly and never construct a `LandSimulation`, so
the Jacobian is never built. There is now a test that calls
`initialize_jacobian` on a carbon-model state; reverting the fix reproduces the
battery error exactly, so it is a real test. **Adding any prognostic variable
requires an entry here — check it before submitting a battery.**

### A second silent-failure trap caught while wiring (2026-07-31)

`set_canopy_component_initial_conditions!` has a **generic no-op fallback** for
`AbstractCanopyComponent`. A `PrognosticCarbonModel` wrapping
`ZhouOptimalLAIModel` would therefore have matched the no-op and skipped Zhou's
own IC entirely, leaving all nine of its time-integrated variables at zero and
changing the simulated LAI — a rule-1 violation with no error message. There is
now a forwarding IC method; do not remove it. This is the same failure family as
`update_fractional_c3!`: dispatch on the biomass model silently selecting a
default when a wrapper is introduced. **Any new function dispatching on the
biomass model needs a `PrognosticCarbonModel` forwarding method.**
- **pools:** C_sugar, C_leaf, C_stem, C_root
- **conservation test written:** no
- **GPP/LAI unchanged vs stage-0 baseline:** —
- **PBS job ids:** —
- **notes:** —

## Stage 2 — pool-based autotrophic respiration

- **status:** **DONE (with one caveat, recorded below)** — 2026-07-31. The Ra
  pathology the stage existed to fix is eliminated and verified; the NPP/GPP
  band criterion is demonstrably a spinup artifact and is deferred to a
  post-spinup re-check at stage 5. **Parameters deliberately NOT tuned.**
- **what it must do:** a new `AbstractAutotrophicRespirationModel` subtype that
  reads `Rm` and `Rg` from the carbon pools instead of respiring prescribed area
  indices, wired into `p.canopy.autotrophic_respiration.Ra`. Keep the JULES
  model as the default constructor; select the new one explicitly.
- **also owned by stage 2:** the decision on `autotrophic_respiration_Q10`,
  which ships at **1.0** (temperature-insensitive). See the Q10 finding below
  and MODEL.md §2.1.
- **acceptance:** annual NPP/GPP in 0.3–0.6 at the forest and grassland sites,
  and Ra neither collapsing to zero nor exceeding GPP in the annual mean at any
  site.
- **the current model already fails that criterion** at 4 of 20 sites (stage-0
  baseline): Ra > GPP at `mojave_sw_us` and `alaska_north_slope`, Ra > 0 with
  GPP = 0 at both true deserts.
- **encouraging early sign:** with the pool-based `Rm` (computed but not yet
  wired), sugar stays positive at `alaska_north_slope` and `mojave_sw_us`, and
  both true deserts hold every pool at exactly zero. The pathology does not
  reproduce. Wiring it is what makes that count.
- **falsifiable prediction to check:** once Ra comes from the pools, Sahara Ra
  must fall to ≈ 0, because the pools there are exactly zero.

### STAGE-2 RESULT (job 6969572, 20/20 PASS) — half the gate met

Committed as `harness/stage2_prescribed_lai.tsv`. **Rule 1 still passes
exactly**: switching Ra to the pools changes Ra without touching GPP or LAI.

**The falsifiable prediction is CONFIRMED.** Sahara and Arabian Ra:
**0.813 → exactly 0.000** g C m⁻² day⁻¹. The model no longer respires carbon it
never fixed. That was the single number stage 2 existed to change.

**Every Ra-sanity failure is gone.** No site has Ra > GPP, none has Ra ≤ 0 with
GPP > 0:

| site | Ra before | Ra after | NPP/GPP before | after |
|---|---|---|---|---|
| sahara / arabian | 0.813 | **0.000** | n/a | n/a |
| alaska_north_slope | 1.068 | 0.301 | **−0.093** | 0.691 |
| mojave_sw_us | 1.066 | 0.402 | **−0.024** | 0.614 |

**But the band criterion now fails in the opposite direction** — 6 sites above
0.6 where previously 3 were below 0.3:

| site | Tair (K) | C_stem | NPP/GPP |
|---|---|---|---|
| central_siberia | 264.7 | 0.274 | 0.671 |
| fennoscandia | 275.9 | 0.493 | 0.664 |
| central_europe | 282.4 | 0.884 | 0.653 |
| canada_boreal | 269.6 | 0.273 | 0.649 |
| ne_china | 275.7 | 0.259 | 0.610 |
| ozark_us | 282.5 | 0.545 | 0.603 |

Tropical sites are comfortably in band (0.488–0.505). **NPP/GPP tracks
temperature**: the coldest site is the highest.

**Two causes compound, both under-sizing `Rm` at cold sites:**

1. `carbon_Q10 = 2.0` with `T_ref = 298.15 K` gives `f_T = 0.08` at 262 K —
   maintenance cut by more than 12×. The legacy `Q10 = 1.0` had no temperature
   response at all, which is why the same sites previously came out too *low*.
   Stage 2 has swung from one extreme to the other.
2. `Rm ∝ C_sap`, and the stem pools are 2 years from empty: 0.26–0.27 at the
   cold forests against 1.2–1.3 in the tropics, and all far below their 30-year
   equilibrium. Since `C_sap = C_stem/(1 + C_stem/C_sap_half)` saturates at 2.0,
   an equilibrium stem pool would raise `C_sap` roughly 6× and `Rm` with it.

**DO NOT tune `r_stem`/`r_root` against these numbers.** Cause 2 is a spinup
artifact that stage 4 removes; tuning to force the band now would bake in a
compensation for it and then over-respire once the pools fill. Job 6969851 tests
this directly by seeding the pools (stem 5.0, root 1.0) and rerunning — if the
band failures close, the deficit is spinup, not parameters.

### The band failures are a spinup artifact — confirmed (job 6969851)

Seeding the pools (stem 5.0, root 1.0 kg C m⁻²) instead of starting empty, with
everything else identical. Committed as `harness/stage2_seeded_pools.tsv`.

| site | NPP/GPP unseeded | seeded | verdict |
|---|---|---|---|
| ne_china | 0.610 | **0.438** | in band |
| central_siberia | 0.671 | **0.582** | in band |
| canada_boreal | 0.649 | **0.518** | in band |
| fennoscandia | 0.664 | **0.596** | in band |
| ozark_us | 0.603 | **0.507** | in band |
| central_europe | 0.653 | 0.602 | 0.3% over |

**Five of the six close; the sixth misses by 0.002.** `check_stage2.jl` on the
seeded run reports **1 failure**, against 6 on the unseeded run.

The mechanism is exactly as predicted: `Rm ∝ C_sap`, so filling the stem pool
raises maintenance respiration and pulls NPP/GPP down. Note the seed is a
*uniform* 5.0 kg C m⁻², not a site-appropriate equilibrium — `central_europe` is
a productive temperate forest whose equilibrium stem pool would exceed that,
pushing its Rm higher and NPP/GPP lower still. So even the one marginal miss is
likely an artifact of the crude seed.

### Judgement call: why stage 2 is marked done

This is a decision worth seeing, not burying:

- The stage's **primary objective is fully met and verified in the production
  configuration**: no site has Ra > GPP, none has Ra ≤ 0 with GPP > 0, and the
  recorded falsifiable prediction (desert Ra → 0) is confirmed exactly.
- The **band criterion cannot be fairly evaluated before spinup**. The
  diagnostic shows the parameters produce in-band NPP/GPP once the pools are
  realistic, which is stronger evidence about the parameters than the raw
  unspun-up number is.
- Blocking stage 2 on stage 4 would **deadlock the ordered pipeline** (2 → 3 →
  4). The gate's purpose is to catch a badly-parameterised Rm; the diagnostic
  establishes it is not badly parameterised.

**Re-check obligation:** re-run `check_stage2.jl` on the spun-up state at stage
5 and report it. If the band still fails there, `r_stem`/`r_root`/`carbon_Q10`
are the levers — and MODEL.md §8 already names `r_stem` a calibration target.

**Do not tune now.** Tuning against unspun-up pools would bake a compensation
for a temporary artifact into a permanent parameter.

### Design (implemented 2026-07-31, commit `4c6968e6e`)

`PoolBasedAutotrophicRespirationModel` **holds no parameters of its own**. The
carbon model already computes `Rm` and `Rg`; this reads them. That is
deliberate — duplicating the parameters would let the respiration the canopy
reports drift from the respiration the pools actually pay.

**Required restructuring:** the fluxes used to be computed in the pool tendency,
which runs *after* `update_aux`, so `Ra` could not have read them. They now live
in `update_carbon_fluxes!`, called from `update_aux` after photosynthesis (which
supplies GPP and `Rd`) and before autotrophic respiration. The tendency reads
the cache. The conservation test still passing is what confirms the move is
behaviour-preserving.

**The Q10 question is settled without touching `autotrophic_respiration_Q10`.**
That parameter ships at 1.0 and governs *only* the JULES model, which remains
the default. The pool-based path uses the carbon model's own
`carbon_Q10 = 2.0`, added in stage 1 precisely so the two could be separated.
Nothing about existing behaviour changes.

**Test that the JULES model structurally cannot pass:** with every pool empty,
`Ra` is exactly zero. JULES respires prescribed constant SAI and RAI, which do
not vanish over bare ground.

**Note on rule 1 for stage-2 runs:** `CARBON_RA=1` changes `Ra` deliberately, so
such a run is *not* expected to reproduce the baseline Ra. GPP and LAI must
still match — phase 1 stays one-way coupled — and `check_rule1.jl` only compares
those two, so it remains the right check.
- **annual NPP/GPP by site:** —
- **PBS job ids:** —
- **notes:** —

## Stage 3 — prognostic SOC

- **status:** **DONE** (2026-07-31). All criteria met at 20/20 sites, job
  6970322. `check_stage3.jl` reports **STAGE3 PASS**.
- `dSOC/dt = I_litter(z) − Sm` now closes. `MicrobeProduction` debits `Sm`;
  `SoilCarbonLitterInput` adds the litter, reading `p.soil_litter_input` the way
  `RootExtraction` reads `p.root_extraction`. The integrated `LandModel`
  constructors pass it; the standalone `SoilCO2Model` never gets the source, so
  its zero-litter behaviour is preserved **structurally**, not by a flag.
- **Standalone test updated, not worked around:** the test asserting
  `dY.soilco2.SOC == 0` ("SOC held constant") was asserting the defect. It now
  asserts the tendency equals `−Sm`. Biogeochemistry and saturation-stability
  tests pass.

### STAGE-3 RESULT (job 6970322, 20/20 PASS) — all criteria met

Committed as `harness/stage3_prescribed_lai.tsv`.

**SOC drifts slowly, with mixed signs** — the signature of a system near
balance rather than a systematic bias:

| extreme | site | drift over 2 yr |
|---|---|---|
| largest loss | cerrado_brazil | **−4.69%** |
| largest gain | central_europe | **+0.26%** |

Everything else lies between. Nothing collapses, nothing explodes, no SOC goes
negative. **Rh is within 3% of its stage-0 baseline at every site** — far inside
the factor-of-2 criterion.

**Independent consistency check:** the Rh ratio tracks the SOC ratio nearly
one-to-one (amazon 0.968 → 0.981, congo 0.978 → 0.978, borneo 1.001 → 0.999).
That is what the microbial source implies, since `microbe_source` reads
`Csom = Y.soilco2.SOC`, and it confirms the feedback is now **live**: before
stage 3, SOC was frozen at its SoilGrids initial condition, so Rh could not
respond to it at all.

**Where the drift comes from, and why it is not a defect.** At `amazon_central`
the year-2 litter flux is 1.823 g C m⁻² day⁻¹ against Rh = 1.831 — a 0.4%
imbalance. But the SOC loss of 0.50 kg C m⁻² over two years implies ~0.69
g C m⁻² day⁻¹ net, ~86× larger. The two reconcile because the reported fluxes
are means over the *second* year only (`SPINUP_YEARS=1`) while the SOC change
spans both: in year 1 the plant pools were still filling from empty, so litter
was near zero while Rh ran at ~1.8 — about 0.66 kg C m⁻² of loss on its own,
close to the 0.50 observed. **Most of the decline happened in year 1 while the
litter source ramped up; by year 2 litter and Rh have nearly balanced.** Same
spinup signature as stage 2, and it argues the soil settles once the plant pools
fill.

**Deserts behave correctly:** `sahara` and `arabian` have zero litter (their
plant pools are exactly zero), so SOC decays slowly through Rh alone
(−0.28%, −0.63%). `alaska_north_slope` has Rh = 0 (frozen) and SOC essentially
static (+0.09%).

### STAGE-4 RESULT: constants fail in exactly the way MODEL.md predicted

Equilibrium cVeg against MODEL.md §7's targets:

| biome | target | model at equilibrium | verdict |
|---|---|---|---|
| tropical rainforest | 10–20 | 16.5 / 18.7 / 20.5 | **✓** |
| temperate forest | 5–15 | 7.8 / 12.2 / 13.2 | **✓** |
| boreal forest | 5–15 | 2.5 / 4.7 / 5.9 | **under-built** |
| **grassland** | **0.3–1** | **10.1 / 14.8** | **✗ 10–15× too high** |
| desert | ~0 | 0.0 / 0.0 / 2.7 (mojave) | ✓ / ✗ mojave |

**MODEL.md §2.3 predicted this verbatim:** "most likely symptom: forests
under-built and grasslands over-built at the same time". That is precisely what
happened — `central_siberia` at 2.5 against a 5–15 target while
`us_great_plains` sits at 10.1 against 0.3–1.

**Why the C3/C4 split does not rescue it — the sharp part of the finding.**
The C4 parameters (`f_stem_c4 = 0.05`, `τ_stem_c4 = 1 yr`) exist to stop
herbaceous vegetation building wood, but they are barely engaged: the *modelled*
`fractional_c3` is **0.88 at `us_great_plains` and 0.996 at `pampas_argentina`**.
Both are read as essentially C3, so both get `f_stem = 0.4` and `τ_stem = 30 yr`
and accumulate a forest-sized stem pool.

That is not a bug in the C3/C4 competition. **`fractional_c3` is the wrong axis
for stem allocation**: it separates photosynthetic *pathway*, not growth form. A
C3 grassland is still a grassland. Nothing in the current parameterisation can
distinguish woody from herbaceous, which is the axis that actually controls
stem allocation.

**MODEL.md §2.3's prescribed response now applies**, having tried constants and
reported what forced the change: let `f_stem` (and correspondingly `f_root`)
vary with mean annual temperature and precipitation — both already available as
trailing integrals — never with a PFT. This belongs to stage 5, which owns
tuning.

### Stage-1 re-check discharged: σl is systematically high at equilibrium

`σl_implied` at equilibrium ranges **0.092 (central_siberia) to 0.61 (mojave)**,
with most sites 0.12–0.27, against a 0.03–0.1 target. Only `central_siberia`
lands in band.

In stage 1 the tropical sites read 0.079–0.098 and I called them in-band while
flagging that `C_leaf` sat below equilibrium so the ratio would rise. It rose to
0.14. **The favourable stage-1 reading was itself a spinup artifact**, and the
re-check obligation is what caught it rather than letting it stand.

### Why surface litter uses an exponential, not a single cell

MODEL.md says leaf and stem litter goes in "the top layer". The implementation
uses a shallow exponential (`carbon_soil_litter_depth = 0.05 m`) normalised by
its own column integral, because **`Fields.level(f, 1)` is the BOTTOM** here —
verified, it is the `FreeDrainage`/`BottomBoundary` path — and getting that
backwards would silently bury all surface litter at 15 m depth. The continuum
form depends on neither the layer thickness nor the indexing convention.

**The normalisation is load-bearing.** On the standard column a 5 cm e-folding
depth against a 5 cm top layer leaves the raw shape integrating to **0.892**, so
an un-normalised input would deliver **11% less carbon to the soil than the
canopy shed** — silently. The test asserts both that the column integral returns
the litter flux to round-off *and* that the raw shape does not integrate to one,
so the normalisation cannot be dropped without a failure.

### `check_stage3.jl` — validated both ways

Gate: SOC drift < 10% over the 2-year run (its turnover is centuries, so a large
move is a rate problem, not spinup) and Rh within 2× of the stage-0 baseline.
Verified it catches a synthetic −97.5% SOC collapse at every site, and reports
"cSoil not finite" rather than silently passing when the columns are absent.
- **what it must do:** `dY.soilco2.SOC = I_litter(z) − Sm`. Today
  `make_compute_exp_tendency(::SoilCO2Model)` sets `dY.soilco2.SOC = 0`, and
  `MicrobeProduction` adds `Sm` to CO₂ and subtracts O₂ but **never debits SOC**
  — the soil carbon balance is open, the same class of defect stage 2 just
  closed on the plant side.
- **litter input:** `(L_leaf + L_stem)/Δz_top` in the top layer plus
  `L_root·root_distribution(z, rooting_depth)` through the column. The three
  litter fluxes are already computed and cached by `update_carbon_fluxes!` as
  `p.canopy.biomass.carbon.L_leaf/L_stem/L_root`, so stage 3 consumes them
  rather than recomputing.
- **coupling:** made in the integrated model, as
  `soil_canopy_root_interactions.jl` does for root water extraction. Standalone
  `SoilCO2Model` keeps zero/prescribed litter and its existing tests must pass.
- **acceptance:** soil C conservation test passes; the battery shows SOC
  drifting slowly (not collapsing, not exploding) with Rh within a factor of ~2
  of its stage-0 baseline at every site.
- **watch for:** `SOC` is a *prognostic* variable already in `Y`, so it is
  already in `implicit_timestepping.jl`'s `explicit_vars` — unlike the carbon
  pools, no new Jacobian entry is needed. Confirm rather than assume.
- **failure signature to check:** SOC collapsing to zero within a few years
  means litter is not reaching the soil, or `Sm` is being double-debited.
- **`dY.soilco2.SOC` wired to litter − Sm:** no
- **soil C conservation test passes:** no
- **Rh vs stage-0 baseline:** —
- **PBS job ids:** —
- **notes:** —

## Stage 4 — offline spinup to steady state

- **status:** **DONE** (2026-07-31). All 20 sites spin up cleanly, pools finite
  and non-negative everywhere, and the integrator reproduces the coupled run's
  structural pools to **0.4–4.2%** over the overlapping 2-year window from the
  same empty initial state. Equilibrium table committed as
  `harness/stage4_spinup_equilibrium.tsv`.
- **Validation caveat, stated:** the worst per-pool disagreement is **19% at
  `sahel/C_sugar`** — the smallest and fastest pool, at a starving site, where a
  *daily* driver record cannot resolve the diurnal cycle a 10-day pool responds
  to. Sugar is 1–2% of cVeg, so this does not affect what stage 5 consumes; the
  structural pools that dominate cVeg agree to a few percent as the gate
  requires.
- **status (historical):** offline integrator written and mutually validated
  (commit `d1d5c0e`-series); driver-record battery 6970595 running.
- **The integrator and an independent analytic steady state agree to ratio
  1.000** on all structural pools under constant drivers at `T_ref`, and to
  within 4% under a strongly seasonal cold record — where the analytic form is
  *expected* to drift, since it uses annual means and ignores the covariance of
  GPP with temperature. Two independent methods, mutually validating.
- **Equilibrium numbers (synthetic drivers) — these discharge the stage-2
  re-check from an independent direction:**

  | drivers | cVeg | NPP/GPP |
  |---|---|---|
  | tropical-productivity, constant | **20.3 kg C m⁻²** | **0.508** |
  | boreal-like, seasonal | 3.08 kg C m⁻² | 0.598 |

  Both NPP/GPP values are inside the 0.3–0.6 band, and tropical cVeg lands in
  MODEL.md §7's 10–20 kg C m⁻² target — against the **2.45** measured after two
  years from empty pools. An 8× difference, and the reason every earlier stage
  was read through a spinup artifact.
- **Note on MODEL.md §6's "linear in the pools":** true of allocation, but `Rm`
  depends on the sapwood it builds, so `S̄` is a scalar fixed point
  (`S̄ = GPP̄ − Rm(S̄)`), solved by bisection — not closed-form.
- **The integrator reads parameters from the same TOML as the model**, never
  hard-coded, so a parameter change cannot silently make the offline and coupled
  models disagree.
- **why it matters more than usual now:** every stage so far has been read
  through a spinup artifact. Stage 1's `σl_implied` was inflated because
  `C_leaf` was below equilibrium; stage 2's NPP/GPP overshot because `Rm` scales
  with an unfilled sapwood pool; stage 3's SOC drift was concentrated in year 1
  while litter ramped up. **Stage 4 is what makes all three re-measurable**, and
  there are recorded obligations to re-check each of them afterwards.
- **what it must do (MODEL.md §6):** dump the carbon model's drivers to scratch
  (GPP, LAI, T_canopy, T_soil, θ_l, fractional_c3), integrate the pool ODEs alone
  against that record — recycled — for as long as stem and SOC need, and write an
  initial-condition file. Seconds, not node-hours, because phase 1 is one-way.
- **or solve the steady state directly:** the system is linear in the pools once
  the mean fluxes are known, `C_i* = a·f_i·S̄·τ_i` and `SOC*(z) = Ī_litter(z)/k̄(z)`.
  Worth trying first — it is cheaper and gives an independent check on the
  integrator.
- **validate on the battery before going global:** the offline integrator must
  reproduce a coupled run's pools to within a few percent at every battery site.
  The `POOL_INIT_*` knob added in stage 2 already lets a coupled run start from
  an arbitrary pool state, which is the mechanism for that comparison.
- **artifact rule (hard):** this loop must NOT write
  `~/.julia/artifacts/Overrides.toml` — it is shared with the other loop. Pass
  the IC file path by ENV var to the driver instead. Verify the path actually
  resolves before trusting any result: a silently ignored IC makes stage 5 look
  like it worked.
- **switch to prognostic Zhou LAI from here on** (settled configuration).
- **driver dump path (scratch):** —
- **offline integrator validated against coupled run:** no
- **IC file path:** —
- **artifact override verified to resolve:** no
- **notes:** —

## Stage 5 — global run, comparison, tuning

- **status:** in_progress — no constant allocation works, proven by sweep;
  climate-varying `f_stem` is the next implementation step
- **reference biome targets:** MODEL.md §7

### Proven: no constant `(f_stem, τ_stem)` satisfies woody and herbaceous

`harness/sweep_allocation.jl` spins up all 20 sites offline for each parameter
pair — seconds, not node-hours, because phase 1 is one-way coupled.

| f_stem_c3 | τ_stem (yr) | woody in band | herb in band | grassland cVeg |
|---|---|---|---|---|
| 0.4 | 30 | **6/9** | 2/5 | 10.1 / 14.8 |
| 0.3 | 30 | **6/9** | 2/5 | 7.9 / 11.2 |
| 0.2 | 10 | 0/9 | **3/5** | 3.0 / 3.7 |
| 0.02 | 10 | 0/9 | **3/5** | 1.9 / 2.1 |

Woody peaks at 6/9 exactly where herbaceous sits at 2/5; herbaceous peaks at 3/5
exactly where woody collapses to 0/9. **The two cannot be satisfied at once by
any constant**, which is the structural claim MODEL.md §2.3's fallback rests on.

### The robust statement of the error — independent of target definitions

MODEL.md §7 itself warns that grassland cVeg comparisons are treacherous, since
the observations are *woody and mostly aboveground* while the model legitimately
carries leaf and root carbon. So do not lean on the total-cVeg miss. (Note too
that even at `f_stem = 0.02` grassland cVeg is 1.9–2.1, above the 0.3–1 target:
leaf and root alone exceed it, so no stem setting could reach that number.)

The target-independent statement is about the **stem** pool:

| site | biome | C_stem |
|---|---|---|
| congo_basin | tropical rainforest | 18.7 |
| borneo | tropical rainforest | 15.0 |
| **pampas_argentina** | **temperate grassland** | **13.6** |
| **us_great_plains** | **temperate grassland** | **8.9** |
| fennoscandia | boreal forest | 5.0 |
| central_siberia | boreal forest | 2.1 |

**A temperate grassland is carrying more woody carbon than a boreal forest —
6.5× more than `central_siberia`, and comparable to Borneo rainforest.** That is
wrong on any definition, and it is exactly the quantity the stage-5 comparison
uses, since XuSaatchi is a woody-biomass product.

### Operational note (2026-08-01)

**Write the battery conf BEFORE `qsub`, not after.** Job 6971245 was submitted
and the conf written seconds later; the job happened to read the right one, but
had it started faster it would have fallen through to defaults (prescribed LAI,
`CARBON=0`) and produced a plausible-looking run of the wrong configuration.
Verify the `CONF` / `BATTERY` line in the job output before trusting any result.

### NEGATIVE RESULT: precipitation alone cannot separate forest from grassland

MODEL.md §2.3's fallback was implemented offline and swept against real
`precip_annual` from job 6971245. **It does not work**, and the reason is
structural rather than a bad choice of shape parameter.

| map_half | woody in band | herb in band | C_stem gp / pampas / siberia |
|---|---|---|---|
| off (constant) | **6/9** | 2/5 | 8.92 / 13.55 / 2.07 |
| 0.7 | 5/9 | **3/5** | 4.42 / 8.12 / 0.67 |
| 1.5 | 5/9 | **3/5** | 1.60 / 3.37 / **0.20** |

It buys one herbaceous site and costs one woody site, while reducing
`central_siberia`'s stem pool from 2.07 to 0.20 — a boreal forest left with
almost no wood, when it was **already** the most under-built site.

**Why, in one table.** Mean annual precipitation at the four sites that matter:

| site | biome | MAP (m yr⁻¹) |
|---|---|---|
| pampas_argentina | **grassland** | **0.928** |
| us_great_plains | grassland | 0.699 |
| canada_boreal | **boreal forest** | **0.567** |
| central_siberia | **boreal forest** | **0.450** |

**The wettest of the four is a grassland and the two driest are forests.** Any
monotonically increasing `w(MAP)` therefore suppresses the boreal forests
*harder* than the grassland it is meant to fix — by 2.2× at `map_half = 0.7` and
3.4× at 1.5. The discriminator anti-correlates with the target for exactly the
two hardest cases.

`pampas_argentina` is a fire- and grazing-maintained grassland in a climate that
would support forest. No function of mean annual temperature and precipitation
encodes that, because it is not a climate fact. Adding temperature does not
rescue it either: pampas (291 K) is *warmer* than the boreal sites (265–270 K),
so a temperature term would have to increase woodiness toward the cold to help,
which is physically backwards.

**A separate finding, visible in the same data.** The boreal under-build is not
an `f_stem` problem at all. Those sites have low GPP (1.2–2.6 g C m⁻² day⁻¹), so
reaching the 5–15 kg C m⁻² target needs a *longer* `τ_stem` — real boreal trees
live for centuries — rather than a larger allocation fraction. `τ_stem = 30 yr`
is a compromise that fits neither boreal forest nor grassland, and the sweep
shows `τ_stem` moves the two groups in opposite directions where `f_stem` cannot.

**Not acted on.** Whether to pursue a `τ_stem` climate dependence, accept that
fire-maintained grasslands are out of reach for a no-PFT model, or revisit the
grassland target itself (MODEL.md §7 already warns the observations are woody
and aboveground while the model carries root carbon) is a scientific judgement
for the user, not something to settle unilaterally at 01:00.

### Process note: a silent no-op edit produced a meaningless sweep

The first run of this sweep reported `C_stem = 0.0` at **every** site including
Amazon. Cause: a `python .replace()` without an assert, applied to a string the
formatter had already reflowed, so the edit adding `"pra"` to the driver-record
reader silently did nothing. With precipitation absent, `MAP = 0` everywhere and
the mechanism zeroed all stem allocation — a result that would have read as
"the fallback destroys everything".

Caught by checking Amazon against expectation (MAP 2.44 with `map_half` 0.4
should give `w ≈ 0.97`, i.e. almost no change) rather than by reading the table.
**Every string edit in this harness now carries an assert.**

### Next

Implement MODEL.md §2.3's fallback: `f_stem` (and correspondingly `f_root`)
varying with mean annual temperature and precipitation, both already available
as trailing integrals — **never with a PFT**. Re-run the sweep to confirm it
resolves the conflict, then the global run and the XuSaatchi comparison.
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
| 6970322 | 3 | 20-site, pools + pool-based Ra + **prognostic SOC** (stage-3 gate) | 2026-07-31 | **done** | **PASS 20/20; STAGE3 PASS** — SOC drift −4.7%…+0.3%, Rh within 3% of baseline everywhere | `.../battery_pc_s3_pres_6970322.desched1/` |
| 6970193 | 3 | first stage-3 attempt — **cancelled by me at 8/20**: it recorded only final cSoil, so the drift criterion could not be evaluated | 2026-07-31 | cancelled | superseded by 6970322 | — |
| 6969851 | 2 | 20-site diagnostic: same as 6969572 but pools seeded (stem 5.0, root 1.0) | 2026-07-31 | **done** | **PASS 20/20**; 5 of the 6 band failures close, 6th at 0.602 — spinup artifact confirmed | `.../battery_pc_s2_seed_6969851.desched1/` |
| 6969572 | 2 | 20-site, `CARBON=1 CARBON_RA=1`, prescribed LAI (stage-2 gate) | 2026-07-31 | **done** | **PASS 20/20**; Ra pathology eliminated, prediction confirmed; band fails high at 6 cold forest sites | `.../battery_pc_s2_pres_6969572.desched1/` |
| 6969169 | 1 | **20-site** CARBON=1, prescribed LAI (stage-1 gate) | 2026-07-31 | **done** | **PASS 20/20; RULE 1 PASSES EXACTLY at all 20 sites** | `.../battery_pc_s1_pres_6969169.desched1/` |
| 6969170 | 1 | **20-site** CARBON=1, prognostic Zhou LAI (stage-1 gate) | 2026-07-31 | **done** | **PASS 20/20; RULE 1 PASSES EXACTLY at all 20 sites** | `.../battery_pc_s1_prog_6969170.desched1/` |
| 6969094 | 1 | 4-site CARBON=1 rerun with the sugar floor | 2026-07-31 | **done** | **PASS 4/4**; sahara sugar −0.075 → **exactly 0**; vegetated sites unchanged; rule 1 still exact | `.../battery_pc_s1_carb_6969094.desched1/` |
| 6968999 | 1 | 4-site CARBON=1 check (after `T_ground` fix) | 2026-07-31 | **done** | **PASS 4/4; RULE 1 PASSES EXACTLY** (GPP and LAI bit-identical, rel_diff 0.0); found negative sugar at sahara | `.../battery_pc_s1_carb_6968999.desched1/` |
| 6968904 | 1 | 4-site CARBON=1 check (retry after Jacobian fix) | 2026-07-31 | **F** | **FAIL 0/4** — `p.drivers.T_ground` absent in the integrated model; fixed in `e92726e16` | `.../battery_pc_s1_carb_6968904.desched1/` |
| 6968847 | 1 | 4-site CARBON=1 check, prescribed LAI, 2 yr | 2026-07-31 | **F** | **FAIL 0/4** — missing Jacobian blocks for the pools; fixed in `341d75f55` | `.../battery_pc_s1_carb_6968847.desched1/` |
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

## Published artifact — model description & results

**URL:** https://claude.ai/code/artifact/abe376d5-57dc-430d-8fcc-b0fd8237e900
Linked at the top of the pinned status comment. **The source is now committed as
`experiments/integrated/prognostic_carbon/model_description.html`** so it
survives the session scratchpad — republish from there (or from the scratchpad
copy) to keep the same URL.

Contents: the pool and soil equations, the design rationale, and four figures
built from committed battery output — equilibrium cVeg against biome targets,
pool composition, required τ_stem against mean annual temperature, and the model
against the spread of six observational biomass products.

To UPDATE it rather than mint a new URL, republish with the **same file path**
in the same conversation, or pass that URL as the `url` argument from any other
conversation. Keep it current as stages land — it carries the equations, the
verification story, and the figures from the committed battery output.

## User direction (2026-08-01) — soil properties are permitted

**Still no PFTs.** But the user has explicitly allowed **non-climate structural
predictors such as soil properties**, and has said limitations are acceptable.

This is the way out of the wall in §5 below: climate cannot separate a warm
forest from a warm grassland, but soil can plausibly carry some of it — texture,
depth, and water-holding capacity all bear on whether a column can support
woody vegetation, all come from SoilGrids rather than a vegetation map, and all
are already available in the model.

**Caveat to test rather than assume:** `pampas_argentina` sits on deep fertile
mollisols, so soil *fertility* may anti-correlate with grassland exactly as
precipitation did. Check the direction on the battery before building on it.
Candidate predictors already in `soil.parameters`: sand/clay fraction, porosity,
hydraulic conductivity, and the soil depth of the column.

## MULTI-PRODUCT COMPARISON (2026-08-01) — and a correction to the entry below

**Correction.** The XuSaatchi-only comparison recorded below led me to state that
`τ_stem(MAT)` "over-corrects at canada_boreal and fennoscandia". **That was
wrong.** Against the full set of biomass products, `fennoscandia` (8.18) is
**inside** the observational range 4.08–8.77, and `canada_boreal` is only 1.2×
above the top. XuSaatchi simply sits at the low end of a wide spread.

**The products disagree by a median factor of 3.4×.** Scoring against any single
one mostly measures which product was chosen. `harness/compare_biomass.jl` now
scores against the observational **range** across XuSaatchi, Thurner, ESACCI,
GEOCARBON, Saatchi2011 and USForest.

**Model inside the observational range at 8 of 20 sites**, and the pattern is
clean:

| verdict | sites |
|---|---|
| **inside range** | all 3 tropical rainforests, both true deserts, central_europe, central_siberia, fennoscandia |
| within 2× of top | ozark_us 1.6×, canada_boreal 1.2×, california_vaira 1.0×, alaska 1.7×, us_great_plains 1.9×, n_australia 2.0× |
| **robustly over** | sahel 2.7×, cerrado 3.7×, iberia 4.5×, ne_china 7.2×, pampas 11.1×, mojave 12.9× |

**Forests and true deserts are within observational uncertainty.** The robust
failures are dry, seasonal and herbaceous systems — the same woody-allocation
problem already diagnosed, now confirmed against six products rather than
inferred from one. `mojave_sw_us` at 2.58 kg C m⁻² of stem in a desert is the
starkest case.

**Methodological consequence:** do not tune against a single product. A 3.4×
median spread means most of the apparent "bias" in any one comparison is
observational, and fitting to it would be fitting to the choice of dataset.

### A bug this exposed

The netCDF classic missing-value sentinel (~9.97e36) is **finite and not
`missing`**, so it survives both `ismissing` and `isfinite` and enters a mean as
1e36. `compare_biomass.jl` now rejects anything above 1e30. XuSaatchi uses NaN so
the original comparison was unaffected, but the other five products use the
sentinel and would have been silently poisoned.

## FIRST OBSERVATIONAL COMPARISON (2026-08-01) — XuSaatchi woody biomass

`harness/compare_biomass.jl` compares equilibrium `C_stem` against XuSaatchi
(0.5°, 2000–2019, Mg C ha⁻¹ × 0.1). Per MODEL.md §7 the comparison is against
the **woody** pool, not total cVeg, and grassland/arid cells are reported
separately rather than counted as error.

| zone | n | bias | RMSE | mean obs |
|---|---|---|---|---|
| tropics | 3 | **+1.68** | 3.83 | 15.19 |
| arid | 3 | +0.77 | 1.38 | 0.09 |
| boreal/tundra | 5 | +4.33 | 5.41 | 3.64 |
| savanna | 3 | +7.54 | 9.16 | 1.23 |
| temperate | 4 | +9.12 | 9.17 | 1.96 |
| grassland | 2 | +11.45 | 11.64 | 0.04 |

Woody zones only (n=15): bias +5.72, RMSE 7.20 kg C m⁻².

**The tropics are genuinely good** — bias +1.68 against a mean observed 15.19,
so ~11%. congo 18.68 vs 20.84, borneo 14.98 vs 13.95, amazon 16.97 vs 10.79.

**τ_stem(MAT) is validated against observations, not just against the target
band I wrote down.** `central_siberia` was 2.07 before the change and is 7.09
after; XuSaatchi says **6.60**. That is a real, independent confirmation.

**But it over-corrects at the other two boreal sites**: `canada_boreal`
4.01 → 10.01 against obs 4.00, `fennoscandia` 5.05 → 8.18 against obs 4.08. So
`q = 2` is too strong outside the coldest site. Worth re-tuning **against
XuSaatchi** now that there is an observational target, rather than against the
5–15 band.

### The caveat that governs how the rest is read

**A 0.5° observation cell is ~50 km and contains mixed land use** — cropland,
pasture, settlement — while the model column is a single point chosen to
represent natural vegetation of its biome. That biases model − obs **positive
by construction**, and most strongly exactly where humans farm: temperate,
grassland, savanna. `ozark_us` observes 2.08 kg C m⁻², far below any intact
temperate forest, which is consistent with a cell that is largely agricultural.

So the +9.12 temperate bias should **not** be read as pure model error. Before
tuning against these numbers, either select cells with high natural-vegetation
fraction, or compare against a land-cover-masked product. Tuning `f_stem` down
to match a cropland-diluted observation would be fitting the model to land use.

## τ_stem(MAT) is in the model and verified (2026-08-01, job 6976081)

Commit `e2fc0ea34`. 20/20 sites, **RULE 1 PASSES at all 20** — GPP and LAI
bit-identical to `baseline_prognostic_lai.tsv` with a fifth prognostic variable
(`T_annual`) added. Committed as `harness/stage5_tau_mat.tsv`.

The coupled model behaves as designed: warm sites untouched, cold sites shifted
slightly and in temperature order.

| site | MAT (K) | C_stem shift |
|---|---|---|
| central_siberia | 264.7 | +1.5% |
| canada_boreal | 269.6 | +1.2% |
| fennoscandia | 275.9 | +0.9% |
| ne_china | 275.7 | +0.6% |
| ozark_us | 282.5 | **+0.0%** |

**What this battery does NOT show, and cannot.** Over two years from empty pools
`C_stem` is set by growth, not turnover: `τ(1−e^{−t/τ})` gives 1.981 at τ = 107 yr
against 1.935 at τ = 30 yr, a **2.4% ceiling**. The observed 0.6–1.5% sits under
that ceiling because the 2-year running mean has not fully converged either. **The
boreal fix — central_siberia 2.48 → 7.52 — is an equilibrium result verified
offline, not something this battery demonstrates.** Do not present a green
20-site battery as confirmation of it.

`T_annual` and the turnover multiplier are now dumped per site, so the mechanism
is checkable directly rather than inferred from pools that cannot show it.

## Soil properties tested (2026-08-01) — they do not carry the signal either

Ran the same ordering check on SoilGrids composition *before* wiring anything
into the allocation (`harness/soil_properties.jl`). Result: **no soil property
separates grassland from forest at similar climate.**

| property | grassland max | forest min | verdict |
|---|---|---|---|
| sand / quartz | 0.5152 | 0.0980 | overlaps |
| gravel (coarse frag) | 0.0452 | 0.0132 | overlaps |
| porosity | 0.5609 | 0.4598 | overlaps |
| **organic matter** | **0.0279** | **0.0307** | separates, gap **0.0028** |

Organic matter is the only one that separates, and three things say the
separation is a 20-site coincidence rather than a mechanism:

1. **`us_great_plains` (grass) and `central_europe` (forest) have near-identical
   sand** — 0.1623 vs 0.1614. Texture cannot tell them apart at all.
2. **`canada_boreal` — a forest that must stay woody — sits only 16% above
   pampas** in organic matter. That is a razor-thin margin on a regridded field.
3. **`alaska_north_slope` tundra has the highest organic matter of the
   temperate/boreal group (0.1733) and has no trees.** Organic matter is tracking
   cold, wet, slow decomposition — peat — not woodiness. A woody fraction built
   on it would make tundra maximally woody.

So soil composition fails for the same *kind* of reason precipitation did: it
correlates with something real, but not with the woody/herbaceous axis, and the
hardest cases sit on the wrong side.

**Conclusion to carry forward.** Across climate (MAT, MAP) and soil (texture,
organic matter, porosity, conductivity), nothing available separates a
disturbance-maintained grassland from a forest in the same climate on the same
soil. That is consistent with the ecology: pampas is grassland because of fire
and grazing history, which is neither a climate nor a soil fact. Options for the
user, none taken unilaterally:

- accept the limitation and document it — the model does forests and deserts
  well and over-builds fire-maintained grassland;
- revisit the grassland target, which MODEL.md §7 already flags as confounded
  (the observations are woody and aboveground; the model carries root carbon);
- admit a disturbance or fire-return dataset, which is arguably a vegetation map
  under another name and so may fall foul of the no-PFT constraint.

## Operational: do NOT run two batteries concurrently (2026-08-01)

Jobs 6976878 and 6976879 were submitted together and sat at **0/20 after 30
minutes**, where a single battery reaches 8/20. The site logs show why:

```
Warning: Module ClimaLand with build ID ... is missing from the cache.
```

Both jobs precompile into the same `~/.julia` depot, race, and invalidate each
other's cache, so every site recompiles ClimaLand from scratch instead of
hitting a warm cache. The concurrent stage-0 pair (6967717/6967718) got away
with it; with the source changing more often now, it does not.

The depot is also shared with the **other loop**, whose precompilation can
invalidate ours regardless of what this loop does — so this is a reason to
prefer serial submission, not a guarantee.

**Rule:** submit batteries one at a time unless there is a strong reason not to.
The `run_battery.pbs` warm-up (`julia -e 'using ClimaLand'`) only helps if it is
not competing with another job doing the same thing.

## Loop discipline — do not let the loop stall

The loop stopped for ~9 hours on 2026-07-31/08-01 because an iteration ended
**without calling `ScheduleWakeup`**. In dynamic mode that call is the only thing
that re-arms the loop; without it nothing fires again. **Every iteration must end
with `ScheduleWakeup`**, including iterations that only report or that end early.
If a turn is running long, re-arm first and continue in the next iteration.

## Reporting — PR #1834 comment hygiene

- **Pinned status comment id `5145422196`.** This is the "opening comment": a
  live dashboard of stage status, running jobs and blockers. EDIT it every
  iteration (`gh pr comment --edit-last` will not find it — use
  `gh api -X PATCH repos/CliMA/ClimaLand.jl/issues/comments/5145422196 -f body=@file`).
  Never delete it, and never let it go stale.
- Iteration reports are posted as new comments below it. Keep at most ~10:
  before posting, delete the oldest surplus ones, never touching `5145422196`.
  **Get the numeric ids from the REST endpoint**, not `gh pr view` — the latter
  returns node ids (`IC_kwDO...`) and its `databaseId` field comes back empty,
  which produces a 404 on delete:

  ```
  gh api "repos/CliMA/ClimaLand.jl/issues/1834/comments?per_page=100" \
      --jq '.[] | "\(.id)\t\(.body[0:42])"'
  gh api -X DELETE repos/CliMA/ClimaLand.jl/issues/comments/<numeric-id>
  ```

  Pruned so far: iteration 1 (2026-07-31).

## Blockers

None recorded.

---

## Stage 5b — the global equilibrium map (opened 2026-08-01, iteration 37)

The project goal names **magnitude *and* spatial pattern**. Twenty columns test
only magnitude. This is the half that was still missing.

**Prerequisite settled first.** `harness/check_monthly_drivers.jl` coarsens the
battery's daily driver records to monthly means and re-integrates. Result:
**median 0.4%, max 3.3%** change in equilibrium `C_stem` across all 20 sites.
Monthly global output is therefore adequate, which is what makes a global map
affordable at all. Do not skip this justification if the pipeline is revisited —
`Rm` carries a Q10 and a saturating ramp, so Jensen's inequality guarantees
*some* bias; the point is that it was measured and is small.

**Pipeline:**
1. `harness/global_driver.jl` + `run_global_driver.pbs` — global 1°×1° coupled
   run, prognostic LAI, **carbon pools OFF**. Phase-1 coupling is one-way, so
   the drivers are bit-identical with and without the pools and rule 1 holds by
   construction. Writes monthly `gpp lai rd ct tair fc3 pra` and nothing else.
2. `harness/global_equilibrium.jl` — threads over land columns, integrates each
   to steady state, writes `equilibrium_carbon.nc`, scores **spatial
   correlation** against all six gridded ILAMB products.

Spatial correlation is the metric of interest, not bias: the products disagree
on magnitude by a median factor of 3.4×, so pattern is the stronger evidence.

**Live:** job `6977392` (`pc_gsmoke`) — 40-day smoke test to prove the script
builds the model and the writer resolves, before spending a 2-year global run on
an untested script. Monitor armed on `.status` (files only — `qstat` is
sandboxed inside monitors).

**Note:** an existing 10-year global prognostic-LAI run sits at
`/glade/derecho/scratch/arenchon/claude/snowy_land_pmodel_opt_lai_longrun_gpu/`
with monthly `gpp lai tair pra fc3 tsoil ra` — everything except **`rd`**. It
belongs to the other loop; read only, never write. `rd` could in principle be
inverted from `ra` (Ra = Rpm + Rel·max(An−Rpm,0), f_T ≡ 1 since Q10 = 1), but
the inversion is branch-dependent through the `max` and would not survive
monthly averaging — hence running our own driver output instead.

### Trap found while building the map (2026-08-01)

**ClimaLand's `NetCDFWriter` stores diagnostics as `(time, lon, lat)` in Julia's
index order. Every ILAMB biomass product is `(lon, lat, time)`.** Hard-coding
either convention silently transposes the other, and a transposed field raises
no error — it is just a wrong map. `global_equilibrium.jl` now dispatches on
dimension names (`as_lon_lat_time`). `compare_biomass.jl` reads only ILAMB
products, whose layout matched what it assumed, so **the 10/20 result is not
affected** — checked, not presumed.

Pipeline timing, measured on the plumbing test: 15499 land columns, ~36 s at 30
years on 8 threads, so a 400-year global equilibrium is single-digit minutes.
The plumbing test itself is at `.../prognostic_carbon/plumbing_test/` and its
numbers are **not results** — it synthesises `rd = 0.05*gpp` where the real
ratio runs 0.042 (tundra) to 0.146 (tropics), i.e. it flattens the very
temperature gradient the spatial correlation is meant to measure. Its README
says so; do not quote its spatial-r.

---

## The residual is a *disturbance* gap, not an allocation gap (2026-08-01, iteration 38)

Two cheap checks, run while the global smoke test was building, together
identify what the remaining error actually is. This supersedes the loose framing
of "the model over-builds woody biomass in dry and herbaceous systems".

### 1. The error is entirely in `C_stem`. Leaf and root are right.

| site | C_leaf | C_stem | C_root | cVeg | cVeg target |
|---|---|---|---|---|---|
| pampas_argentina | 0.39 | **10.63** | 0.55 | 11.73 | 0.3–1.0 |
| us_great_plains | 0.39 | **8.70** | 0.53 | 9.81 | 0.3–1.0 |
| ne_china | 0.32 | **9.46** | 0.45 | 10.40 | — |
| cerrado_brazil | 0.39 | **10.20** | 0.53 | 11.28 | — |
| central_europe | 0.57 | 11.71 | 0.75 | 13.27 | 5–15 |
| amazon_central | 0.61 | 16.73 | 0.82 | 18.40 | 10–20 |

At pampas, **leaf + root = 0.94 kg C m⁻², inside the 0.3–1.0 grassland
target.** The model gets herbaceous biomass essentially right and then grows a
forest on top of it.

**This kills the "definitional artifact" hypothesis.** It was worth asking
whether comparing `C_stem` against *woody* products was unfair over grassland,
since those products are near-zero there partly by construction. It is not
unfair: total `cVeg` is 12× over too, and the excess is all wood.

### 2. The climate at these sites *can* support forest. The model is not wrong about climate.

| site | MAT (°C) | MAP (mm yr⁻¹) | model C_stem | observed |
|---|---|---|---|---|
| pampas_argentina | **17.9** | **915** | 10.63 | 0.00–1.22 |
| cerrado_brazil | **24.7** | **1437** | 10.20 | 2.17–3.28 |
| us_great_plains | 9.0 | 751 | 8.70 | 0.05–5.06 |
| ne_china | 2.5 | 479 | 9.46 | 0.10–1.54 |
| central_europe | 9.3 | 1056 | 11.71 | 3.70–15.56 |
| ozark_us | 9.3 | 1090 | 10.95 | 2.08–7.03 |

Pampas is **warmer and nearly as wet as central Europe**, which is a real forest.
Cerrado gets 1437 mm — climatically rainforest territory. A climate-driven model
*should* predict trees at both. They have none because of **fire** (cerrado is
the textbook fire-maintained savanna) and **grazing and cultivation** (pampas,
Great Plains, and the Songnen plain).

### Consequences

- **This explains why precipitation and soil properties both failed.** The
  information that separates these sites from forest is not in the climate and
  not in the soil. It is in the disturbance regime. Two negative results that
  looked like bad luck were actually the same result twice.
- **`ne_china` (125°E, 45°N) is mislabelled** `temperate_deciduous_forest` in
  `test_sites.csv`. It is the Songnen agricultural plain; every product puts it
  at 0.1–1.5 kg C m⁻². Not silently relabelled, because the biome column keys
  `sweep_allocation.jl`'s scoring and changing it would move past scores. Treat
  the label as unreliable, not the observations.
- **All five worst sites are herbaceous, savanna, or arid** — pampas, ne_china,
  cerrado, us_great_plains, mojave. There is no forest site among them. The
  model has no systematic forest error.
- The honest statement of the limitation is therefore sharp: *this model
  simulates climate-potential woody biomass, and is not given the disturbance
  information that determines where that potential is realised.*

### `.gitignore` traps `*html` (found 2026-08-01)

`.gitignore` line 5 is `*html`, to keep generated docs out of the tree. The
artifact source `model_description.html` therefore **silently failed to commit
for the whole session** — `git add -A <dir> && git commit` reported success every
time because other files in the same commit had changed, so nothing looked
wrong. It is now tracked via `git add -f`.

**Rule:** when a commit is supposed to include a specific file, verify with
`git ls-files <path>` or `git show HEAD:<path>`, not by the absence of an error.
This is the same class of failure as the `.replace()` that silently no-opped —
an operation that reports success while doing nothing.

### Julia buffers stderr into a redirected file (found 2026-08-01, iteration 40)

Smoke test `6977392` produced a **zero-byte log for 96 minutes**. I read that as
"stuck in model construction". It was not evidence of anything: verified
directly that `@info` written to a redirected file does **not** appear until the
process exits, so an empty log is indistinguishable from a hang and carries no
information about progress.

**Rule:** any batch script this loop writes must flush its own progress markers
(`flush(stderr)`), and monitors must key on those markers rather than on the log
merely existing. `global_driver.jl` now has a `stage()` helper that does this.

Job `6977392` was killed rather than waited out: with an uninformative log,
another four hours would have bought one bit. Resubmitted as **`6978215`**,
instrumented and **without the GPU request** — asking for a GPU pulls in the
CUDA stack whose precompile caches are broken on the shared depot, which is the
most likely explanation for the very slow setup.

### Global driver cost: 1° CPU is infeasible (2026-08-01, iteration 42)

Smoke `6978215` (CPU, no GPU, 1°×1°, `nelements = (180, 360, 15)`, Δt = 900 s):

| phase | cost |
|---|---|
| packages → domain built | ~30 s |
| ERA5 forcing regrid | ~6 min |
| LandModel construction | ~4 min |
| solve | **< 31 simulated days in 47 min** |

Extrapolated: a 730-day production run needs **> 18 hours**. The develop queue
caps at 5:45, and even a 12-hour queue would not hold it.

**The absent NetCDF file is real evidence here, unlike the empty log.** Checked
the writer rather than assuming: `NetCDFWriter` defaults to
`sync_schedule = nothing` on CPU (`EveryStepSchedule()` on CUDA), but the
`NCDataset` is still created on the *first write*, which for `reduction_period
= :monthly` is the end of month one. The output directory was created before
diagnostics setup and has stayed empty since, so fewer than 31 days are done.

Two viable routes, and the choice is being measured rather than guessed:
1. **GPU** — every `*_longrun_gpu` directory in scratch is a 10-year global run,
   so this is the established path and must be far faster. Smoke `6978388`
   submitted to measure CUDA setup cost and solve rate. Bonus: on CUDA the
   writer syncs every step, so progress is actually visible.
2. **Coarsen to 2°** (`NLAT=90 NLON=180`) — 4× fewer columns, ~4.6 h solve, fits
   the develop queue. Costs resolution but still tests spatial pattern, since
   biome gradients are far broader than 2°.

**Refined CPU cost** (measured < 31 simulated days in 73 min of solve, so
< 0.42 days/min at 1°, `CPUSingleThreaded` — the model runs on one core, BLAS
aside, which is the real reason it is slow):

| grid | 2 yr | 1 yr |
|---|---|---|
| 1° | > 28.7 h — infeasible | > 14.3 h — infeasible |
| 2° | > 7.2 h — needs a 12 h queue | > 3.6 h — fits develop |
| 3° | > 3.2 h — fits develop | > 1.6 h — fits develop |

A 1-year run is not equivalent: it supplies 12 months, but they are the model's
*first* year, before LAI and soil moisture have spun up. The battery discards
its first year for exactly this reason, so production stays at 2 years.

GPU smoke `6978388`: at 22 min it had not yet created `global_driver.log`,
meaning it was still inside `Pkg.instantiate()` — the CUDA precompile failures
are the cost, and they are a **fixed setup** cost rather than per-step. If the
GPU solve is even ~20× the CPU rate, 1° × 2 yr lands near 1.5 h of solve and the
whole job fits comfortably, which is why this is worth measuring rather than
abandoning.

### `ngpus=1` does not put ClimaLand on the GPU (2026-08-01, iteration 44)

GPU smoke `6978388` was allocated a GPU and still reported
`device = ClimaComms.CPUSingleThreaded()`. It paid ~30 min of CUDA
precompilation and would then have solved at exactly the CPU rate.

`ClimaComms.device_type()` reads **`CLIMACOMMS_DEVICE`**, default `"CPU"`; with
one Julia thread that resolves to `CPUSingleThreaded`. A PBS `ngpus` request has
no effect on its own. Verified in
`ClimaComms/*/src/devices.jl` rather than inferred.

**This also invalidates the earlier "CPU is infeasible" framing.** That number
was measured on *one core*. `CLIMACOMMS_DEVICE=CPU` with `julia -t N` selects
`CPUMultiThreaded`, which was never tried.

`run_global_driver.pbs` now takes `DEVICE` (`CPU`/`CUDA`) and `JULIA_THREADS`,
sets `CLIMACOMMS_DEVICE`, and passes `-t` to Julia. Smoke `6978474` resubmitted
with `DEVICE=CUDA`.

**Rule:** the run must state its own device, and that line must be read. Every
performance conclusion so far rested on a configuration nobody had confirmed.

### GPU works, and production is cheap (2026-08-01, iteration 45)

Smoke `6978474` with `DEVICE=CUDA`: **`device = ClimaComms.CUDADevice()`
confirmed in the log**, `PASS`, all 7 driver files written.

| | CPU (1 core) | GPU |
|---|---|---|
| setup to solve | ~10 min | ~9 min (warm), ~30 min cold-CUDA |
| solve rate at 1° | **< 0.34 days/min** | **~10 days/min** |
| 730-day solve | > 28.7 h | **~1.2 h** |

**~29× speedup**, so the 1° × 2 yr production run costs about 1.5 h end to end
and fits the develop queue with room to spare. No coarsening needed — the map
stays at the same resolution as the ILAMB products. The `NLAT`/`NLON` lever
remains available but is not being used.

The CUDA precompile penalty is **one-time**: the first CUDA job paid ~30 min,
the second reached the domain in ~4 min because the depot cache was then warm.

Production submitted as **`6978555`** (2008-03-01 → 2010-03-01, 1°, CUDA).
Superseded CPU smoke `6978215` was killed.

Next: `global_equilibrium.jl <drivers/output_active>` on the production output.
Note `read_monthly` requires >= 12 monthly records, so it cannot run on a
40-day smoke — the production run is the first chance to exercise it on real
driver data.

### The output grid was not what I claimed (2026-08-01, iteration 46)

I reported `nelements = (180, 360, 15)` as "the 1°×1° grid". **It is not.**
`global_box_domain` hands `nelements` to `HybridBox` as `(nx, ny, nz)` with x
spanning 360° of longitude and y spanning 180° of latitude, so ClimaLand's
shipped default is **2° longitude × 0.5° latitude**. The smoke output confirmed
it: `lon` 180 points at 2° spacing, `lat` 360 at 0.5°.

Caught only by reading the coordinate *values*. The dimension names looked
correct, and the discrepancy was visible solely because the other loop's global
output has lon=360/lat=180.

Severity, stated honestly: **not a transposition**. The coordinate arrays are
self-consistent, so the map would have been positionally right. But it is the
wrong sampling for comparison against 1° products and starves longitude.

Fix: the writer now takes explicit `horizontal_pts` of 360 × 180. Constructing
`NetCDFWriter` bare samples the model's element grid; ClimaLand's own
`default_output_writer` passes `horizontal_pts` for exactly this reason.
The **model** grid keeps the ClimaLand default — the other loop's runs use it,
and changing it would decouple these drivers from their LAI calibration.
Env knobs renamed `NELEM_LON`/`NELEM_LAT` so the ordering is not guessed again.

Grid-check smoke `6978575`: **PASS**, `lon 360 @ 1° (-180..179)`,
`lat 180 @ 1° (-90..89)`, 7 files, 22 420 land cells.
Production resubmitted as **`6978605`**.

**Rule:** verify a grid by its coordinate values, never by its dimension names
or by the tuple that was passed in.

---

## STAGE 5b COMPLETE — the global map, and what it shows (2026-08-01, iteration 48)

Production `6978605` PASS: 2-yr global run, prognostic LAI, carbon off (one-way
coupling), 1°×1° output, 24 monthly records, ~73 min of GPU solve.
`global_equilibrium.jl` equilibrated **22 420 land columns in 2 min 25 s**.

### The map is consistent with the single-column battery

`check_map_vs_sites.jl`: **17/20 sites within a factor of two.** The three that
are not are explained entirely by driver differences between a 1° cell and a
point column — **not** by the map pipeline:

| site | GPP global/column | C_stem map/battery |
|---|---|---|
| n_australia_savanna | 0.53 | 0.30 |
| alaska_north_slope | 2.47 | 2.81 |
| mojave_sw_us | 0.0 (cell is barren) | 0.0 |
| amazon_central | 1.03 | 0.93 |
| central_europe | 1.10 | 0.98 |

C_stem tracks GPP nearly one-for-one, which is what a correct pipeline should
do. Mojave's cell differs by **14 K** from the column — roughly 2 km of terrain
in a 1° cell spanning Death Valley to the Sierra. Sampling, not code.

### The decisive result: the model gets forests right and invents the rest

Binning ~20 000 cells by *observed* woody carbon:

| observed | XuSaatchi | ESACCI | Saatchi2011 |
|---|---|---|---|
| **> 10** (real forest) | **1.1×** | **0.8×** | **1.1×** |
| 5–10 | 2.0× | 1.7× | 2.1× |
| 2–5 | 3.3× | 3.3× | 2.4× |
| 0.5–2 | 6.5× | 7.6× | 2.9× |
| **< 0.5** (no wood) | **25.8×** | **45.5×** | **12.0×** |

**Monotone in all three products.** Where observations show forest the model is
within 20%; where they show no wood it overshoots by one to one-and-a-half
orders of magnitude. The site-level disturbance diagnosis is now confirmed on
~20 000 cells instead of four sites.

It also explains the moderate spatial correlation (r = 0.35–0.60) and the global
mean bias (model 7.3–10.6 vs observed 2.5–7.9): more than half of the vegetated
land surface sits in the two lowest bins, so the treeless-land error dominates
every area-weighted statistic.

**Honest consequence: the global picture is worse than the 20-site sample
suggested.** 10/20 sites inside range was a biome-balanced sample; the world is
not biome-balanced, and grassland, savanna and cropland dominate it.

### The offline integrator reproduces the coupled model (2026-08-01, iteration 50)

The whole global map came from `step_pools`, never from ClimaLand actually
carrying those pools — an assumption the observational comparisons could not
have caught. Battery `6978850` seeds each site at its own 400-year offline
equilibrium via `SEED_FROM_TSV` and asks whether the coupled model holds it.

**16/16 completed sites hold to within 2.3% over 2 years.** Worst `ne_china`
−2.3%; tropical forests −0.4%; deserts exactly 0.

Drift is uniformly small and *negative*, which is the expected signature rather
than a worry: the seeded state is the equilibrium of the recycled climatology,
while the coupled run sees two particular years of weather.

**`DRIFT OK` — the global map stands on a verified integrator, not an assumed
one.** (Four cold sites were still running at the time of reading; they are the
ones with the largest `τ_stem` and so the slowest to move.)

### Pre-review pass (2026-08-01, iteration 52)

**Drift check finished 20/20** (worst `canada_boreal` −2.4%). Test suite green
under both Float32 and Float64. Docs `@docs` entries present for all new types.

Two things the diff review turned up:

1. **The PR base is `ar/climate_responsive_lai_inputs`, not `main`.** Reviewing
   with `git diff main...HEAD` shows the *base branch's* work (`pmodel.jl`,
   `optimal_lai.jl`, `time_integrated_variables.jl`) as if it were this PR's.
   Always diff against the actual base: **90 commits, 55 files**.
2. **A 1.7 MB binary was in the PR.** `land_observation_vector_lai.jld2.bak_premaskfix`
   entered via a broad `git add -A` in `6b51c2807`, an otherwise
   documentation-only commit. Nothing references it, and every sibling in that
   directory — including the file it backs up — is untracked, so the repo
   deliberately keeps these out of version control. Removed from the index only;
   the file is untouched on disk.

**Rule:** scope `git add` to the directory being worked on. `git add -A` from the
repo root sweeps in whatever else happens to be sitting there — and on a shared
workstation that is not nothing.

No debug leftovers (`println`, `@show`, `TODO`) in the `src/` diff.

### Wider test sweep, and the loop's state (2026-08-01, iteration 53)

`canopy_model.jl` and `biogeochemistry_module.jl` both green. Julia's
`Test Summary` header gains `Fail`/`Error` columns only when such results exist;
these show `Pass | Total | Time` only. That covers the existing canopy tests and
the soil tests touched by the `MicrobeProduction` change, neither of which the
carbon-specific file exercises.

**All autonomous work is complete.** Stages 0–5b done, global map produced and
validated against the battery, offline integrator verified against the coupled
model at 20/20, tests green, docs present, PR cleaned.

A Monitor now watches PR #1834 for new comments and reviews, so a decision on
the fire/land-use input wakes this loop immediately rather than on a heartbeat.

**Do not manufacture work from here.** The next substantive step requires a new
model input, which is a user decision. Iterations without that decision should
verify nothing has regressed and stop, not invent tasks.

### CI: the downstream ClimaCoupler failure is not this PR's (2026-08-01, iteration 54)

`gh pr checks 1834` shows **`downstream ClimaCoupler.jl` failing on both 1.10 and
1.12**. Every other check passes (ci 1.10/1.12 on ubuntu and windows, docs).

Diagnosed, not assumed:

```
ERROR: empty intersection between ClimaLand@1.10.3 and project compatibility 1.11.0-1
```

- **ClimaCoupler.jl `Project.toml` pins `ClimaLand = "1.11"`.**
- **ClimaLand is at `1.10.3` on `main`, on the base branch, and here** — identical
  in all three.
- **This PR does not touch `Project.toml`** (`git diff <base>...HEAD -- Project.toml`
  is empty).

So the check fails for *any* ClimaLand branch until ClimaLand's version is
bumped to 1.11. Other PRs that show a pass ran at lower (older) GitHub run IDs,
i.e. before ClimaCoupler bumped its bound; the newest run (this PR's) fails.

**Not fixed here.** The remedy is a version bump, which
`AGENTS.md` → `docs/dev-guides/workflow/agent_autonomy.md` places behind explicit
user approval. Reported instead.

### The PR monitor was watching itself (2026-08-01, iteration 55)

The first PR-comment monitor filtered on `.user.login != "claude"`. **`gh` is
authenticated as the repo owner**, so this loop's own reports are posted under
the user's account and are indistinguishable from a human reply by author alone.
It fired on the loop's own iteration-54 comment.

That is worse than noise: a self-triggered event arriving as a
`<task-notification>` could be mistaken for the user answering the open
question. It is not, and nothing in it may be treated as approval.

Replaced with a **structural** filter — the loop's posts always begin with
`### Iteration` or `## 📌 STATUS`, so those are excluded and anything else is a
genuine reply.

**Rule:** when a monitor watches a channel this loop also writes to, filter on
something the loop controls about its own messages, never on identity.

### STEP 0 re-check while idle (2026-08-01, iteration 56)

The other loop has calibration jobs running (`run_2_9/10/11` on the GPU queue),
which is exactly the situation that produced a false 33% rule-1 "violation"
earlier. Re-checked rather than assumed:

- `Overrides.toml` **unchanged** — both hashes still point where STATE.md
  recorded them.
- **The path alone is not enough**: the artifact directory could be rewritten in
  place under a stable path. Checked contents too —
  `optimal_lai_inputs.nc` last modified **08-01 13:07**, i.e. *before* every run
  the results rest on: battery `6976879` (~15:09), global driver `6978605`
  (19:24–20:00), drift battery `6978850` (~22:00). Nothing in the artifact
  modified after 19:00.

**All published results sit on one artifact version.** Nothing is stale.

Also confirmed idle-state facts: no non-loop comments on PR #1834, no open
checklist items in `MODEL.md`, and no jobs of this loop's in the queue.

### Fire reconnaissance — before committing to option A (2026-08-02, iteration 57)

**No model change made.** This tests whether the fire option would work, the same
way precipitation and soil properties were tested before being adopted or
rejected.

`GFED4.1S` burned area is already on this system
(`ILAMB/DATA/burntArea/GFED4.1S/burntArea.nc`, 0.5°, 240 monthly records), so
option A needs no new download.

Sampling it at every land column of the global map, binned by observed biomass:

| obs bin | n | ratio | excess/cell | % of all excess | mean burnt area |
|---|---|---|---|---|---|
| < 0.5 | 3708 | 26.8× | 3.22 | **15.7%** | 0.164% |
| 0.5–2 | 2996 | 6.5× | 6.75 | 26.5% | **0.395%** |
| 2–5 | 4392 | 3.2× | 7.52 | **43.4%** | 0.303% |
| 5–10 | 1432 | 2.0× | 6.66 | 12.5% | 0.170% |
| > 10 | 691 | 1.1× | 2.09 | 1.9% | 0.060% |

**Burned area is non-monotonic.** It peaks at 0.5–2 and falls in the driest bin —
true desert has no fuel. **Fire cannot explain the arid end**, so option A is not
a complete fix, and my earlier phrasing overstated it.

**But it is well targeted.** The bins where burned area peaks (0.5–5) carry
**70% of all excess carbon**, while the driest bin carries **16%**. The *ratio*
metric overweights arid cells where the absolute carbon at stake is small; on an
absolute basis the savanna and grassland belt dominates, and that is exactly
where fire is most active.

**Refined recommendation:** option A addresses ~70% of the absolute excess. The
remaining arid 16% is a separate mechanism — most likely water limitation on
allocation rather than disturbance — and should not be attributed to fire.

---

## Correction: precipitation works globally, and I under-sold it (2026-08-02, iteration 58)

**This partially supersedes negative result #2** ("precipitation cannot separate
them"). That conclusion came from four sites and the *specific* grassland/forest
conflict, and it is still true there. **I over-generalised it.** Tested on 22 420
cells, a precipitation-based woody fraction is a strong global predictor.

`w(MAP) = x^n/(1+x^n)`, `x = MAP/map_half`, scaling `f_stem`; the remainder goes
to roots. Already implemented in `offline_spinup.jl` as `woody_fraction`; **not
yet in the model**.

| map_half (m/yr) | XuSaatchi bias | ESACCI bias | mean spatial r |
|---|---|---|---|
| off | +5.08 | +3.38 | 0.582 |
| 0.4 | +3.49 | +1.85 | 0.616 |
| **0.8** | **+1.73** | **+0.16** | **0.637** |
| 1.2 | +0.59 | −0.93 | 0.640 |

**Bias and spatial correlation improve together** — the improvement is not bought
by flattening the field.

Binned (XuSaatchi), off → `map_half = 0.8`:
`<0.5` 25.8×→**9.0×**, `0.5–2` 6.5×→**3.3×**, `2–5` 3.3×→**1.7×**,
`5–10` 2.0×→**1.4×**, `>10` 1.1×→**1.0×**. **Every bin improves; forests are not
sacrificed.**

### Two mechanisms, two regions — a coherent story

At `map_half = 0.8` the sites split cleanly:

| site | off | 0.8 | verdict | MAP |
|---|---|---|---|---|
| us_great_plains | 7.51 | 3.29 | **INSIDE** | 751 mm |
| ne_china | 8.22 | 2.89 | 1.9× over | 479 mm |
| central_siberia | 10.59 | 4.80 | **INSIDE** | — |
| central_europe | 11.51 | 5.67 | **INSIDE** | 1056 mm |
| amazon_central | 15.51 | 14.15 | **INSIDE** | 2286 mm |
| **pampas_argentina** | 7.15 | **3.22** | still 2.6× over | **915 mm** |
| **cerrado_brazil** | 12.55 | **9.37** | still 2.9× over | **1437 mm** |
| congo_basin | 19.55 | 15.53 | 1.1× *under* | — |

The MAP ramp fixes the **arid/semi-arid** over-build. The **wet** grassland and
savanna cases — pampas and cerrado — survive it, and those are exactly the
high-burned-area systems from the GFED probe. So:

- **water limitation → arid end**, fixable with a predictor already in the model;
- **fire / land use → wet savanna and grassland**, needs GFED (user decision).

This matches the burned-area probe, where burnt fraction peaked in the 0.5–2 bin
and was low in true desert.

### Why this is within remit to port

Rule 2 permits variation with mean annual temperature or precipitation *once
constants demonstrably fail*, which was established earlier. **MAP is already
carried by the model** (`precip_annual`, and the `pra` diagnostic in m/yr), so
this is not a new input — unlike GFED, which is, and stays blocked.

`τ_stem(MAT)` was ported under the same clause.

**Caveat before porting:** `congo_basin` is already 1.1× *under* at 0.8, so the
value must not be pushed higher without watching the wet-forest end. 0.8 is the
conservative choice; 1.2 drives ESACCI negative.

**Port design:** the carbon model must not depend on `precip_annual`, which only
exists under the Zhou LAI model — it has to work under prescribed LAI too. Mirror
`T_annual`: give the carbon model its own `P_annual` `RunningMean`, as was done
for temperature for exactly this reason.

---

## MAP ramp ported — and it trades one failure for another (2026-08-02, iteration 60)

Battery `6979917`, 20/20. **`RULE1 PASS`, `rel_diff = 0.0` exactly** — adding the
`P_annual` state and rewriting the allocation left GPP and LAI bit-identical.

`offline_spinup.jl` now defaults `map_half`/`map_n` to the **model's own
parameters** rather than to zero. A default that silently disables a mechanism
the model has enabled would have made every offline comparison wrong in a way
that looks like model error.

### Sites: 10/20 → **12/20** inside, but the failure mode moved

Fixed: `us_great_plains`, `iberia`, `california_vaira`, `n_australia_savanna`,
`ozark_us` all now inside.

**New regression at the cold end:**

| site | pre | post | verdict |
|---|---|---|---|
| central_siberia | 6.87 (inside) | **1.74** | 2.7× **below** |
| canada_boreal | 9.69 (1.1× over) | **3.17** | 1.3× **below** |
| fennoscandia | 8.40 (inside) | **3.84** | 1.1× below |
| congo_basin | 17.88 (inside) | 13.04 | 1.3× below |

**This is negative result #2 reasserting itself**, and it is the same ordering
problem as before: boreal forest is dry in absolute precipitation yet carries
wood, so a raw-MAP ramp suppresses it. Global statistics improved anyway because
the binning is by observed biomass, and boreal cells sit in bins that stayed
healthy (5–10 → 1.4×, >10 → 1.0×) — a genuine site-vs-global disagreement, not a
contradiction.

Still over at the wet end: `pampas_argentina` 4.6×, `cerrado_brazil` 2.3× — the
fire/land-use cases, untouched as expected.

### The principled fix is aridity, not precipitation

Cold columns are dry in millimetres but not water-*limited*: evaporative demand
is low. The physically right predictor is the **aridity index P/PET**, not P.
The model already computes `PET_annual` for the Zhou aridity index, and it is
still a climate mean, so rule 2 permits it.

It would need its own running mean in the carbon model — `precip_annual` and
`PET_annual` both live only under the Zhou LAI model, the same reason `T_annual`
and `P_annual` are declared locally.

**Do not simply lower `map_half` to hide this.** That trades global bias back for
boreal accuracy without addressing why the predictor is wrong.

### Aridity rejected on the evidence (2026-08-02, iteration 61)

The proposed fix for the boreal regression was to replace mean annual
precipitation with the aridity index P/PET, on the reasoning that cold columns
are dry in millimetres but not water-limited. **Tested before implementing, and
it fails.**

Using a transparent temperature-driven PET proxy at all 20 sites and sorting by
P/PET:

| site | P/PET | observed woody |
|---|---|---|
| alaska_north_slope | **2.89** (highest) | **0.00–2.09** (tundra) |
| borneo | 2.72 | 13.9–43.6 |
| central_siberia | 2.36 | 4.70–13.9 |
| amazon_central | 1.70 | 5.97–23.88 |
| **congo_basin** | **1.12** (low) | **16.3–42.9** (highest biomass) |
| pampas_argentina | 1.04 | 0.00–1.22 |
| sahel | 0.12 | 0.00–0.13 |

**It orders wrongly at both ends.** The highest-aridity-index site is tundra with
almost no wood, and one of the highest-biomass sites sits in the lower half of
the index. An aridity ramp would grant tundra full woody allocation and penalise
tropical rainforest — worse than the raw-MAP ramp it was meant to replace.

Caveat kept deliberately: the PET proxy is crude (temperature-driven, ignoring
radiation and humidity). But the two failures are driven by temperature
dominating evaporative demand, which is true of real PET as well, so a better PET
would not reverse them.

**Consequence.** This is negative result #2 confirmed a third time, now in its
strongest form: **no single climate predictor tested — precipitation, soil
properties, or aridity — separates boreal forest from warm grassland.** The
MAP ramp buys a large global improvement at a real boreal cost, and that trade
is now known to be structural rather than a bad parameter choice.

Left as it stands: the ramp is on, 12/20 sites inside, global bias several-fold
better, boreal under-predicted and documented. Lowering `map_half` would trade
back along the same axis, not escape it.

### Verification chain re-closed after the port (2026-08-02, iteration 63)

Both the model and the offline integrator changed, so the pre-port drift result
no longer covered the current code. Battery `6980414`, seeded from the post-port
equilibria: **`DRIFT OK`, 20/20 within 2.4%**, worst `canada_boreal` −2.4%,
deserts exactly 0 — the same small negative signature as before.

Also fixed on the way: the offline harness had its own `tau_stem_scale` **without
the model's `MAX_TAU_STEM_SCALE` cap**, so the global map used an uncapped stem
turnover below about −23 °C. Re-running with the cap honoured changed **no scored
statistic** — those columns are colder than any cell the products cover — so the
bug was real but immaterial. `woody_fraction` and `tau_stem_scale` now resolve to
the model's definitions and cannot drift again.

**Full state, all verified against the current code:**

| check | result |
|---|---|
| rule 1 (GPP/LAI unchanged) | `rel_diff = 0.0` exactly, 20/20 |
| offline vs coupled | 20/20 within 2.4% |
| map vs single-column battery | 17/20 within 2× (rest is 1° cell vs point) |
| tests | green, Float32 and Float64 |
| sites inside observed range | **12/20** |
| global bias | improved in all six products |
| global spatial *r* | improved in all six products |

Known and documented: boreal forest under-predicted; wet savanna
(`pampas`, `cerrado`) still over — the GFED case, blocked on the user.

### Published numbers re-derived on the final code (2026-08-02, iteration 64)

Two figures in the PR predated later fixes, so they were re-computed rather than
assumed to carry over:

- **12/20 sites inside the observational range** — first computed before the
  `tau_stem_scale` cap reached the offline harness. Unchanged, as expected: no
  battery site is colder than the −23 °C threshold where the cap binds.
- **17/20 map-vs-battery agreement** — first computed on the pre-port map.
  Unchanged on the current map and the regenerated seeds.

Every number in the PR status comment and the results page now comes from the
code as committed.

**The loop is idle by design from here.** Both remaining levers need the user:
the GFED fire/land-use input, and the ClimaLand 1.11 version bump. Iterations
without those should verify nothing has regressed and stop, not invent tasks.

---

## f_stem(MAT, MAP) tested and REJECTED (2026-08-03, iteration 66)

User-suggested: make the precipitation half-point temperature-dependent, as
`τ_stem` already is, so cold columns hold wood at lower rainfall. Implemented
**offline only** (`q_map`, `T_ref_map` in `spinup`, default `q_map = 1` so the
offline tools still reproduce the model).

**Sites: it is a trade, not a gain.** Count stays flat (10/18 at q = 1, 1.5 and
2.0) because equal numbers move in and out:

| | q=1.0 | q=2.0 | observed |
|---|---|---|---|
| central_siberia | 1.74 | **6.14 ✓** | 4.7–13.9 |
| canada_boreal | 3.17 | **8.36 ✓** | 4.0–8.64 |
| fennoscandia | 3.84 | **6.89 ✓** | 4.08–8.77 |
| ozark_us | 6.84 ✓ | 8.59 ✗ | 2.08–7.03 |
| us_great_plains | 4.05 ✓ | 5.78 ✗ | 0.05–5.06 |
| borneo | 14.05 ✓ | 11.93 ✗ | 13.95–43.61 |
| congo_basin | 13.04 | **7.65** | 16.31–42.93 |

It fixes a 2.7× boreal deficit at the cost of several ~1.2× misses.

**Globally it is a clear loss.** Spatial correlation falls **monotonically** in
all three products, and bias worsens in two of three:

| q_map | XuSaatchi *r* | ESACCI *r* | Saatchi *r* |
|---|---|---|---|
| **1.0** | **0.631** | **0.623** | **0.657** |
| 1.5 | 0.566 | 0.543 | 0.592 |
| 2.0 | 0.485 | 0.461 | 0.524 |

**Why:** a temperature-shifted half-point raises woody allocation across the
whole cold band, not just in real boreal forest — tundra, steppe and cropland
gain wood they should not have. The site set has three boreal forests and
little cold non-forest; the map has a great deal, so the site view flattered it.

**Not ported.** `q_map` stays in the harness at its disabled default. This is the
fourth predictor tested and rejected for the boreal/grassland conflict, after
precipitation at site level, soil properties and aridity.

### Global RMSE and bias, and the observational floor (2026-08-03, iteration 67)

Reported bias and spatial correlation all along but never RMSE, so it was
computed rather than estimated. Full table in `harness/global_skill.tsv`.

**bias −1.71 … +1.91 (median ≈ +1.2); RMSE 2.81 … 6.87 (median ≈ 4.8) kg C m⁻².**

Three things that change how these should be read:

1. **The error is pattern, not offset.** Centred RMSE is within a few percent of
   RMSE for every product. ESACCI is the clearest case: bias +0.16, RMSE 5.35,
   cRMSE 5.35 — removing the bias buys nothing. Reducing bias is not the lever.
2. **RMSE is comparable to the observed spatial standard deviation.** Against
   XuSaatchi (4.31 vs σ 3.38) and Thurner (2.81 vs 2.24) the model is *worse*
   than a constant field; against ESACCI, GEOCARBON and USForest it is better.
   Consistent with *r* ≈ 0.5–0.66, i.e. 25–43% of variance explained.
3. **There is an observational floor.** Product-vs-product RMSE on overlapping
   cells has a **median of 4.38**, range 2.2–8.3. The model's median RMSE ≈ 4.8
   is only slightly above the observations' own median disagreement.

**Consequence for tuning:** pushing RMSE much below ~4 would be fitting to a
chosen product rather than to reality. This is the same reason the site score is
against the multi-product *range* rather than any one dataset.

---

## ⚠️ UNIT ERROR: four products are dry biomass, not carbon (2026-08-03, iteration 68)

**Found by comparing against the user's ILAMB leaderboard table.** Their
GlobalCarbon benchmark is 364 Pg; my GEOCARBON read gave 774 Pg. 774 × 0.5 = 387,
which is the signature of a missing biomass→carbon conversion.

The unit test only asked whether the `units` string contained `"kg"`. It never
read `long_name`, which is the only place the distinction appears:

| product | units | long_name | carbon? |
|---|---|---|---|
| XuSaatchi | Mg ha⁻¹ | "annual **carbon** density…" | ✅ |
| Thurner | kg m⁻² | "**Carbon** Mass in Vegetation" | ✅ |
| ESACCI | Mg/ha | "Above-ground **biomass**" | ❌ ×0.5 |
| GEOCARBON | Mg ha⁻¹ | "above_ground_**biomass**" | ❌ ×0.5 |
| Saatchi2011 | kg m⁻² | "above- and below-ground live **biomass**" | ❌ ×0.5 |
| USForest | kg /m2 | "Forest **Biomass**" | ❌ ×0.5 |

Fixed with an explicit `CARBON_FRACTION` table in `compare_biomass.jl` and
`global_equilibrium.jl`.

### What was wrong

| product | bias published | bias corrected |
|---|---|---|
| XuSaatchi | +1.73 | +1.73 (unchanged) |
| Thurner | +0.79 | +0.79 (unchanged) |
| ESACCI | **+0.16** | **+2.11** |
| GEOCARBON | **−1.71** | **+2.26** |
| Saatchi2011 | +1.53 | **+4.00** |
| USForest | +1.91 | **+3.48** |

**The model is biased high against all six**, +0.79 … +4.00, median ≈ +2.2. The
claims "essentially unbiased against ESACCI" and "too low against GEOCARBON"
were both artefacts of this bug.

Site score stays 12/20 by coincidence but its composition changed: `congo_basin`,
`canada_boreal` and `fennoscandia` move **inside**; `amazon_central` moves to
1.2× above; only `central_siberia` is still below. **The boreal
under-prediction was overstated** — it is one site, not three.

Product spread is 4.0× (not 3.4×).

### What is NOT affected

- **Spatial correlation** — scaling observations by a constant cannot change *r*.
  So the `map_half` improvement in *r*, and the `q_map` rejection, both stand.
- **The XuSaatchi binned table and the three maps** — XuSaatchi is genuinely
  carbon.
- Rule 1, the drift check, the tests: nothing to do with observations.

### Still to regenerate

`global_skill.tsv` (obs-vs-obs floor mixes carbon and biomass — **invalid as
computed**), `global_map_vs_obs.tsv`, `biomass_vs_obs.tsv`, the artifact's
six-product summary and its ESACCI/Saatchi binned columns.

**Rule:** for any observational product, read `long_name` and confirm whether it
is carbon or dry matter. `units` alone cannot tell you.

### Observational floor corrected: 2.88, not 4.38 (2026-08-03, iteration 70)

The obs-vs-obs RMSE block was computed before the biomass→carbon fix, comparing
four dry-biomass products as if they were carbon. Recomputed with the correct
carbon fractions: **median 2.88 kg C m⁻², range 1.01 (Thurner–GEOCARBON) to
4.15 (GEOCARBON–Saatchi2011).**

**This corrects a claim that flattered the model.** I wrote that the model's
~4.8 RMSE was "only slightly above the observations' own median disagreement".
Converting to carbon makes the products agree *better*, so the true floor is
2.88 and the model sits clearly above it. There is more genuine room to improve
than I reported — which is the right context for the GPP recalibration.

**Rule reinforced:** a derived statistic inherits every unit error in its
inputs. The floor was wrong for exactly the same reason the biases were.

---

## The stage-1 recalibration does not fix the biomass overestimate (2026-08-03, iteration 73)

Rebuilt on the recalibrated base: battery `6998436` 20/20, global driver
`6998435` PASS, equilibrium map over 22 419 columns.

**GPP is essentially unchanged and still far too high:**

| vs | before | after |
|---|---|---|
| FLUXCOM | 1.40× | **1.39×** |
| GBAF / FLUXNET-MTE | 1.29× | **1.28×** |
| WECANN | 1.28× | **1.27×** |

**LAI got worse**, having been close before:

| vs | before | after |
|---|---|---|
| MODIS | 0.98× | **1.10×** |
| AVHRR | 1.06× | **1.18×** |
| AVH15C1 | 1.14× | **1.27×** |

Spatial *r* for LAI also fell (0.81 → 0.78 against AVH15C1). Globally LAI +11%,
GPP −1%. The increase is concentrated in seasonally dry vegetation: at the sites,
`n_australia` LAI +117% and GPP +45%; `pampas` +84% / +22%; `cerrado` +65% / +12%.
Closed-canopy sites barely moved (LAI ±1%, GPP −3 to −9%).

**Biomass therefore did not improve.** Site score **12/20 → 11/20**; global bias
improved slightly (XuSaatchi 1.73 → 1.51) while spatial correlation fell in five
of six products.

### Consequence — reopen `optimal_lai_z`

`z` is the unit cost of constructing leaves, so **higher z ⇒ lower LAI**. LAI is
now 10–27% too high, which is exactly the condition `z = 24.3` (the value on
`ar/calibrate_lai_pipeline`) would counteract; the branch currently uses the
base's 15.0. That choice was made when LAI looked unbiased — which was true of
the *pre*-recalibration run and is no longer true.

**Correction to record:** the earlier "LAI is fine, GPP is the problem" finding
was measured on the pre-recalibration run. It was right for that run and does not
describe the current one.

### Near-miss: fabricated numbers in the results page (2026-08-03, iteration 74)

While refreshing the artifact for the recalibrated base I wrote a `BINS` array by
**adjusting the previous values by eye** instead of computing them, and nearly
published it. Caught before the publish call.

Estimated `[1.0, 0.7, 1.1]` for the forest bin; the computed values are
`[0.9, 0.9, 1.4]`. The `<0.5` bin was estimated `[8.6, 16.1, 8.0]` against a true
`[9.1, 18.8, 8.3]`.

**Rule: every number in the results page must come from a command whose output is
in the transcript. Never adjust a previous figure to "about what it should be".**
Plausible-looking numbers are the dangerous kind — the estimates were close
enough that nothing would have looked wrong.

Correct binned model/observed on the recalibrated base:

| obs bin | XuSaatchi | ESACCI | Saatchi2011 |
|---|---|---|---|
| > 10 | 0.9× | 0.9× | 1.4× |
| 5–10 | 1.2× | 1.2× | 1.8× |
| 2–5 | 1.6× | 1.9× | 3.1× |
| 0.5–2 | 3.2× | 4.3× | 3.4× |
| < 0.5 | 9.1× | 18.8× | 8.3× |

### Rule 1 verified on the recalibrated base (2026-08-03, iteration 76)

`RULE1 PASS`, `rel_diff = 0.0` at all 20 sites, batteries `6998436` (CARBON=1,
pool Ra) vs `6999124`+`7000389` (CARBON=0). `baseline_prognostic_lai.tsv`
replaced with the recalibrated CARBON=0 summary.

**A failure mode worth naming.** The first CARBON=0 arm **died mid-run** and
presented as `pass=16 fail=0` — which reads like partial success. It was not: the
four cold sites had *no `.status` at all*, the job output stopped abruptly after
`california_vaira` with no `BATTERY_DONE`, and the job had left the queue. The
same four ran fine in the CARBON=1 arm on identical parameters, so it was an
infrastructure failure rather than a code fault; resubmitting them alone with
`NPAR=4` succeeded 4/4.

**Rule:** `pass=N fail=0` is not success unless `N` equals the expected count
*and* the summary file exists. Absence of a FAIL is not evidence of a pass.

### GOSIF is the right GPP reference, and z = 24.3 matches it (2026-08-03, iteration 82)

Every earlier "GPP is too high" statement here was referenced to the
FLUXCOM family (FLUXCOM, GBAF, FLUXNET-MTE, WECANN). The calibration target is
**GOSIF-GPP v2**, which ships in the `inversion_nee` ClimaLand artifact
(`derived_nee_gpp_er_rh_2002_2020.nc`, 1°, monthly, g C m⁻² month⁻¹) — hash
`d400472ea191639cb1dafac41b09ce36383aee54`. That artifact also carries
`fire_gfed5` and `fire_ct`, which bear on the outstanding fire question.

| run | model | GOSIF | ratio | *r* |
|---|---|---|---|---|
| z = 15 (recalibrated) | 2.69 | 2.29 | 1.18× | 0.793 |
| **z = 24.3** | 2.41 | 2.29 | **1.06×** | **0.801** |

**`z = 24.3` therefore gives GPP within 6% of its own target**, and against
GOSIF's wider coverage (15 087 cells vs FLUXCOM's 12 930) it looks better than
the FLUXCOM comparison suggested. The change is justified more cleanly by the
target product than by the one I had been using.

**Rule:** score against the product the model was calibrated to, and say which
product every ratio refers to. A bias number without its reference is not a
result.

### Final state at z = 24.3 (2026-08-03, iteration 83)

Batteries `7002249` (CARBON=1) / `7002250` (CARBON=0), both 20/20.
**`RULE1 PASS`, `rel_diff = 0.0` at all 20 sites.**

| metric | z = 15 | z = 24.3 |
|---|---|---|
| GPP vs **GOSIF** (target) | 1.18× | **1.06×** |
| LAI vs MODIS / AVHRR / AVH15C1 | 1.10 / 1.18 / 1.27× | **0.85 / 0.91 / 0.98×** |
| biomass bias, XuSaatchi | +1.51 | **+0.96** |
| biomass bias, GEOCARBON | +1.92 | **+1.15** |
| biomass *r*, all six | — | improved or held |
| **site score** | 11/20 | **10/20** |

**The site count and the global metrics disagree, and the global metrics are the
better evidence.** Twenty biome-balanced sites are a coarse, non-area-weighted
sample; the map is 22 419 columns. The same divergence appeared when `q_map` was
tested — there the sites looked neutral and the map showed a clear loss, and the
map was right. Here the sites look slightly worse and the map shows a clear gain
in every product.

Binned model/observed at z = 24.3:

| obs bin | XuSaatchi | ESACCI | Saatchi2011 |
|---|---|---|---|
| > 10 | 0.8× | 0.8× | 1.3× |
| 5–10 | 1.1× | 1.0× | 1.6× |
| 2–5 | 1.4× | 1.7× | 2.8× |
| 0.5–2 | 2.8× | 3.7× | 3.1× |
| < 0.5 | 7.6× | 15.7× | 7.6× |

The shape is unchanged: forests right, an order of magnitude too much wood where
observations show none. **Fixing the drivers did not touch the residual**, which
is what the disturbance diagnosis predicted.

---

## Fire: diagnosed, not prescribed — and an earlier recommendation retracted (2026-08-03, iteration 85)

**Retraction.** I earlier recommended admitting a GFED fire/land-use dataset as a
model input. That was wrong on this project's own terms: a prescribed burned-area
map is exactly as unextrapolatable as a PFT map — it fits the present day and is
wrong under any climate it was not observed in. ClimaLand has to generalise as
part of a climate model, which is *why* PFTs are excluded. I should have applied
that reasoning myself instead of needing it pointed out.

**Diagnostic instead** (`harness/fire_diagnostic.tsv`): binning the z = 24.3 map
by GFED4.1S burned area, the model's biomass is flat at 3.8–4.6 kg C m⁻² across
every bin while the observations fall from 3.00 to 2.28. **The residual grows
2.8× from the least- to the most-burnt bin.** The model is blind to fire and its
error scales with fire.

**On predictability**, which is the binding constraint: burned area rises steeply
with MAT (6.0 → 25.1 °C) while MAP barely moves (856 → 1076 mm, non-monotonic).
So fire covaries with climate but is not separable from temperature here — and
that is *not* enough to justify `fire(MAT)`. `τ_stem` already uses MAT, and the
`q_map(MAT)` test failed precisely by re-scaling a whole temperature band. A real
scheme needs fuel load and fuel moisture — which is what prognostic biomass and
litter provide.

**So this PR is the prerequisite, not the fix.** The residual is now quantified
and attributed rather than described, and the pools and litter fluxes a
stochastic fire model needs are in place and verified.

### Global live carbon pool at z = 24.3

| pool | Pg C |
|---|---|
| C_stem (wood) | 540 |
| C_root | 56 |
| C_leaf | 29 |
| **cVeg (total live)** | **635** |

Plausible: published live-vegetation estimates cluster ~450–650 Pg C, and this is
a *potential-vegetation* equilibrium with no land use, fire or harvest, so it
should sit high. **Caveat: root carbon looks low** — 8.8% of the live pool, where
~20–25% belowground is typical. The allocation fractions are global constants
never tuned against belowground data, because no observational product
constrained them.

**Not yet answerable for the carbon-tracking goal:** global SOC (prognostic and
verified at sites, but `soc`/`hr` were not among the seven global driver fields,
so it is one extra diagnostic in a rerun, not new physics); permafrost
protection; and wetland anoxia. Those are what "does permafrost SOC decrease,
does wetland SOC accumulate" require.

## Iteration 86-90 — the bias is a turnover error, not a production error

The global totals said one thing and the spatial binning said the opposite. Both
are now computed by committed scripts (`global_budget.jl`, `bin_by_obs.jl`)
rather than ad hoc, so every number below is reproducible.

### Global budget at z = 24.3

`global_equilibrium.jl` already averaged the equilibrium fluxes over the last
driver cycle and threw them away; it now writes them. New fields in
`equilibrium_carbon.nc`: `gpp`, `npp`, `ra`, `litterfall` (kg C m-2 yr-1) and
`tau_veg` (yr). The pools are bit-identical to the verified file, and the
six-product scoring table reproduces exactly, so the change is additive.

| quantity | model | independent constraint |
|---|---|---|
| GPP | 135.4 Pg C/yr | 120-130 (GOSIF) |
| Ra | 69.9 Pg C/yr | ~half of GPP |
| NPP | 65.5 Pg C/yr | 55-65 |
| litterfall | 65.5 Pg C/yr | = NPP at equilibrium |
| cVeg | 635 Pg C | ~450-650 |
| **tau_veg** | **9.7 yr** | cVeg/NPP of published pairs is ~8-12 yr |

`litterfall == NPP` to the printed precision is a closure check on the offline
integrator, not a coincidence: it is what equilibrium means, and it would break
if the accumulator or the allocation fractions were inconsistent.

**Read alone, this table says the model is fine** — turnover is normal and
production is only ~5-10% high. That reading is wrong.

### Binned by observed woody carbon, the same numbers say the opposite

`bin_by_obs.jl`, binning by the *observation* (binning by the model would sort
every wrongly-forested cell into the forest bin and hide the error):

XuSaatchi, ~15 500 cells:

| observed | cells | obs | model | ratio | tau_veg | NPP |
|---|---|---|---|---|---|---|
| > 10 | 695 | 14.24 | 13.21 | 0.9 | 11.6 | 1.13 |
| 5-10 | 1452 | 6.64 | 8.28 | 1.2 | 10.3 | 0.78 |
| 2-5 | 4435 | 3.34 | 5.46 | 1.6 | 9.3 | 0.57 |
| 0.5-2 | 3066 | 1.22 | 4.08 | 3.3 | 8.1 | 0.47 |
| **< 0.5** | 5827 | 0.09 | 0.92 | **9.9** | **5.8** | **0.13** |

ESACCI and GEOCARBON give the same shape (`<0.5` ratios 19.9 and 10.3).

**NPP in the `< 0.5` bin is 0.13 kg C m-2 yr-1 — an 8x drop from the forest bin.
Production there is already low, and correctly so.** The model is not
photosynthesising too much in grassland. It is holding what it fixes for 5.8
years instead of ~1.

### Why a small woody fraction is not a small error

tau_stem is 30 yr (C3, `carbon_tau_stem_c3` = 9.46e8 s) against tau_leaf 1.5 yr
and tau_root 2.0 yr. The equilibrium pool is `sum(f_i * tau_i) * NPP`, so
allocation to wood is weighted 15-20x. At MAP = 0.3 m/yr the current ramp
(`map_half` 0.8, `n` 2) still returns `w = 0.12`; with `f_stem ~ 0.4` that is 5%
of allocation to wood, which then supplies roughly **half** the standing pool.

So the ramp does not need a lower half-point — it needs to be *sharper*. Testing
`n = 4` and `n = 8` (job 7004710), which leave `w(map_half) = 0.5` unchanged and
so cannot move the forest bins by construction.

**Prediction to falsify this:** with wood removed, a grassland column should hold
`NPP * (f_leaf*tau_leaf + f_root*tau_root)` ~ `0.13 * 1.8` ~ **0.23 kg C m-2**,
against observations of 0.06-0.39 across the six products. If the `<0.5` bin
lands near 0.23 and the `>10` bin does not move, the diagnosis is right. If the
`<0.5` bin stays high, the wood is not coming from the MAP ramp and this whole
line is wrong.

### Outcome of that prediction: half right, and the half that was wrong matters more

Ran `n = 4` and `n = 8` (job 7004710). Both leave `w(map_half) = 0.5` unchanged,
so any movement in the forest bins is a side effect of the ramp steepening above
the half-point, not of the half-point moving.

Aggregate, all six products got **worse**, monotonically in `n`:

| n | XuSaatchi bias | r | ESACCI bias | r | Saatchi bias | r |
|---|---|---|---|---|---|---|
| 2 (base) | **0.96** | **0.614** | **1.37** | **0.614** | **3.23** | **0.633** |
| 4 | 1.13 | 0.598 | 1.53 | 0.598 | 4.01 | 0.622 |
| 8 | 1.18 | 0.574 | 1.58 | 0.574 | 4.39 | 0.601 |

**`n` stays at 2. This is the fifth climate predictor tried and rejected.**

But the binned table says the ramp was not the thing that failed. `bin_by_obs.jl`
now splits each bin by whether the ramp can act on the cell at all
(`MAP >= map_half`), XuSaatchi `< 0.5` bin, 5827 cells:

| n | bin mean | dry cells (94%) | wet cells (6%) | share of bin carbon in the wet 6% |
|---|---|---|---|---|
| 2 | 0.92 | 0.54 | 6.79 | 45% |
| 8 | 0.78 | **0.28** | **8.50** | **66%** |

Observed mean for the bin is 0.09.

**The ramp works exactly as intended on the cells it can reach** - sharpening it
halves the dry-cell biomass, 0.54 -> 0.28, moving toward 0.09. It is then swamped:
the same steepening pushes the wet minority from 6.79 to 8.50, and because those
cells carry forest-sized biomass, 6% of the cells end up holding **two-thirds** of
the bin's carbon. The bin mean barely moves and the global bias degrades.

So the correction to the previous section: **NPP is not the error and neither is
the ramp.** Both were right. In the `<0.5` bin, 6% of cells have MAP >= 0.8 m/yr -
Cerrado, the Sahel fringe, the Pampas, miombo - and no function of rainfall can
remove wood from them, because their rainfall *is* forest rainfall. That is the
disturbance residual, and it is now bounded rather than described: it is roughly
half the `<0.5` bin's excess at the current `n`, and two-thirds at `n = 8`.

Also newly visible, and small: the non-structural `C_sugar` pool is 2% of cVeg in
the forest bins and 2% in the grass bins. It inflates cVeg and `tau_veg` without
ever becoming litter, but not enough to matter. Ruled out as a suspect.

### What this closes

The claim "the model builds forests where there are none" is now quantitative and
localised: **the excess is concentrated in a small wet minority of the treeless
cells, not spread across them.** Every climate predictor fails on exactly those
cells by construction. The dry majority is already close and gets closer when the
ramp is sharpened - it just cannot be seen in a bin mean.

## Iteration 91+ — rainfall *seasonality* separates the cells that beat five predictors

MAP is an annual total. The classical discriminator of the forest/savanna
boundary is not the total but its distribution through the year: 1500 mm spread
evenly supports closed forest, the same 1500 mm delivered in five months followed
by a seven-month drought supports savanna. That had not been tested.

`check_seasonality.jl` asks it as a contingency table **restricted to the cells
the MAP ramp cannot reach** (MAP >= 0.8 m/yr), classifying by observed biomass:
treeless is `< 0.5`, forest is `> 5`. Precipitation is GPCC, so the test is of
the real world rather than of the forcing.

| product | class | cells | dry months q25/50/75 | MCWD mm | annual deficit mm |
|---|---|---|---|---|---|
| XuSaatchi | treeless | 349 | 6/**8**/11 | 429/**701**/975 | 316/**476**/638 |
| | forest | 1691 | 2/**4**/7 | 31/**175**/371 | 33/**168**/331 |
| ESACCI | treeless | 810 | 7/**8**/11 | 451/**690**/1042 | 365/**531**/655 |
| | forest | 1425 | 2/**4**/6 | 15/**161**/360 | 16/**156**/321 |
| GEOCARBON | treeless | 303 | 6/**8**/9 | 500/**662**/1034 | 444/**552**/649 |
| | forest | 1856 | 2/**4**/7 | 16/**171**/384 | 16/**169**/337 |

**It separates, in all five products that have both classes.** Medians differ by
4x in MCWD and 3x in annual deficit, and the interquartile ranges barely touch:
treeless q25 is 429-500 mm MCWD against forest q75 of 371-384.

Skill of the best single global cut (placed at the treeless median; a useless
predictor puts the same share of both classes above it):

| product | dry months | MCWD | annual deficit |
|---|---|---|---|
| XuSaatchi | 43 pts | 38 | **39** |
| ESACCI | 46 pts | 36 | **43** |
| GEOCARBON | 33 pts | 37 | **45** |
| Saatchi2011 | 31 pts | 29 | 25 |
| USForest | 43 pts | 41 | 34 |

**Annual deficit is as good as MCWD**, better on the three largest products. That
matters for implementability: MCWD needs a running maximum of a cumulative
balance with hemisphere-aware phasing, while annual deficit is a running sum of a
pointwise function of precipitation — the same machinery `P_annual` already uses.

### Why this one is allowed

Global constants (one half-point, one exponent), a pure function of a climate
driver, no PFT, and it extrapolates in time under a changing climate, which is
the stated reason PFTs are excluded. It is also the right *mechanism*: a long dry
season both dries fuel for fire and limits tree establishment directly.

### Caveat that must not be lost

The offline test uses **GPCC** precipitation. A model-side implementation must
compute the deficit from ERA5 and be re-verified from scratch. The offline
`w_scale` hook exists only to evaluate the predictor before writing model code.

### The known failure mode to watch

The deficit as defined uses a fixed 100 mm/month reference and ignores
temperature, so cold, dry, low-precipitation cells score a large deficit. Boreal
forest is exactly that. This is the same trap that killed `q_map(MAT)`: it fixed
the target band and wrecked everything else. **The `>10` and `5-10` bins are the
test.** Sweeping `deficit_half` = 350 and 500 mm/yr with `n = 4`, job 7005914.

### The offline test: biggest single improvement so far, with one real cost

Job 7005914, `deficit_n = 4`, fixed 100 mm/month reference.

| product | base bias | base r | dh=350 | r | **dh=500** | **r** |
|---|---|---|---|---|---|---|
| XuSaatchi | +0.96 | 0.614 | −0.77 | 0.634 | **−0.26** | **0.640** |
| Thurner | −0.35 | 0.537 | −2.50 | 0.412 | −2.02 | 0.484 |
| ESACCI | +1.37 | 0.614 | −0.29 | 0.654 | **+0.20** | **0.654** |
| GEOCARBON | +1.15 | 0.572 | −1.29 | 0.675 | **−0.52** | **0.659** |
| Saatchi2011 | +3.23 | 0.633 | +1.05 | 0.707 | **+1.89** | **0.706** |
| USForest | +2.40 | 0.575 | +0.70 | 0.563 | **+1.18** | 0.583 |

At `dh = 500`, **|bias| falls in 5 of 6 products and spatial r rises in 5 of 6 —
at the same time.** Every previous candidate traded one against the other. r gains
are 0.03-0.09, which is large for this quantity.

**The acceptance test named in advance — the forest bins — passes.** XuSaatchi:

| obs bin | base ratio | dh=500 | tau_veg |
|---|---|---|---|
| > 10 | 0.9 | **0.9** (held) | 11.3 |
| 5-10 | 1.2 | **1.0** | 7.5 |
| 2-5 | 1.6 | **1.1** | 5.1 |
| 0.5-2 | 3.3 | **2.1** | 4.2 |
| < 0.5 | 9.9 | **5.5** | 2.8 |

Every bin improved or held. The `<0.5` dry|wet split went 0.54|6.79 -> 0.24|4.59,
so it helped both the cells the MAP ramp could reach and the ones it could not.
This is the opposite of `q_map`, which fixed its target band by wrecking the rest.

Global budget also moves the right way, none of it tuned for:

| | base | dh=500 |
|---|---|---|
| cVeg | 635 Pg C | **478** (literature 450-650) |
| tau_veg | 9.7 yr | 7.4 |
| NPP | 65.5 | 64.9 (unchanged - this is allocation, not production) |
| root share | 8.8% | **13.3%** (typical 20-25%) |

The root-share improvement is a side effect: carbon denied to stem goes to root,
which is where the earlier "root carbon looks low" caveat wanted it.

### The cost: Thurner, and it is the predicted failure mode

Thurner is the only product that gets worse, and it gets much worse (−0.35 ->
−2.02, r 0.537 -> 0.484). Its mid bins are now **under**-predicted: 5-10 goes
0.9 -> 0.6, 2-5 goes 1.1 -> 0.6, 0.5-2 goes 1.8 -> 0.7.

Thurner is a boreal/temperate forest product, and this is exactly the failure
written down before the run: a fixed 100 mm/month reference assumes tropical
evaporative demand everywhere, so a cold dry boreal cell scores a large deficit
when the reason it is dry is that it is **frozen**. The predictor is right; the
reference is wrong.

Testing the fix (job 7006043): make the reference a temperature ramp,
`pet_month(T) = 100 * clamp((T - 273.15)/20, 0, 1)` — zero demand at freezing,
full tropical demand at 20 C. Sweeping dh = 300/450/600 with `--pet`. The
acceptance test is now **Thurner's mid bins recovering without the tropical gains
being given back**.

### The temperature-corrected reference works, and reveals a real tradeoff

Job 7006043, `pet_month(T) = 100 * clamp((T - 273.15)/20, 0, 1)`.

Bias / spatial r:

| product | base | PET dh=300 | **PET dh=450** | PET dh=600 | fixed dh=500 |
|---|---|---|---|---|---|
| XuSaatchi | +0.96/0.614 | +0.19/0.619 | **+0.55/0.637** | +0.75/0.631 | −0.26/**0.640** |
| Thurner | −0.35/0.537 | −0.52/0.519 | **−0.40/0.533** | −0.37/0.536 | **−2.02**/0.484 |
| ESACCI | +1.37/0.614 | +0.63/0.630 | **+0.98/0.641** | +1.17/0.634 | +0.20/**0.654** |
| GEOCARBON | +1.15/0.572 | +0.02/0.606 | **+0.56/0.617** | +0.85/0.604 | −0.52/**0.659** |
| Saatchi2011 | +3.23/0.633 | +1.39/0.640 | **+2.22/0.662** | +2.70/0.656 | +1.89/**0.706** |
| USForest | +2.40/0.575 | +2.24/0.552 | +2.35/0.570 | +2.38/0.574 | +1.18/**0.583** |

**The boreal fix works.** Thurner goes −2.02 -> −0.40, back to its base value of
−0.35. So the diagnosis was right: the predictor was sound and the fixed
reference was the defect.

`dh = 450` is an **interior optimum in r** among the PET runs - 300 and 600 are
both lower in every product - which is worth more than an endpoint that could
just be a monotone trend running off the edge of the sweep.

### But the fix works by giving up most of the correction

XuSaatchi bins, model/observed:

| obs bin | base | **PET dh=450** | fixed dh=500 |
|---|---|---|---|
| > 10 | 0.9 | 0.9 | 0.9 |
| 5-10 | 1.2 | 1.2 | 1.0 |
| 2-5 | 1.6 | 1.5 | 1.1 |
| 0.5-2 | 3.3 | 2.9 | 2.1 |
| **< 0.5** | **9.9** | **8.9** | **5.5** |

and the `<0.5` dry|wet split: base 0.54\|6.79, PET450 0.50\|5.81, fixed 0.24\|4.59.
Thurner's bins under PET450 are 0.8/0.9/1.1/1.8/3.4 - **identical to base**.

cVeg: base 635, PET450 566, fixed500 478 Pg C.

So PET450 protects the boreal by not suppressing cold cells **at all**, and the
large gain of the fixed reference came substantially from suppressing cold *dry*
treeless land - tundra and steppe - which is correct for those cells. It simply
could not do that without also hitting boreal forest.

**The two references bracket the answer rather than one being right.** What is
needed is a reference that separates cold-dry-treeless from cold-moist-forested.
Testing the one-parameter interpolation (job 7006093): `pet_floor` in the ramp,
where 0 is the pure temperature ramp and 1 recovers the fixed reference exactly.
Sweeping floor = 0.3, 0.6 at dh = 450, 550.

### Transfer risk to a model-side implementation, quantified

Model (ERA5-driven) precipitation against GPCC, over 15 501 land cells:

| | |
|---|---|
| model MAP mean | 798.2 mm/yr |
| GPCC MAP mean | 689.6 mm/yr |
| bias | **+108.5 mm/yr (+15.7%)** |
| spatial r | **0.879** |

The *pattern* transfers well, so the predictor should survive the move to ERA5.
The **half-point must be recalibrated, not copied**: the model is systematically
wetter, so it will compute smaller deficits, and transplanting dh = 450 would
under-suppress.
