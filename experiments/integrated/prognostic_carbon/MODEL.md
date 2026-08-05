# Prognostic carbon pools — model specification

A simple live/dead carbon model for ClimaLand. It lives in
`src/standalone/Vegetation/biomass.jl`, alongside `PrescribedBiomassModel` and
`ZhouOptimalLAIModel`.

**Design rule that everything else follows from: the carbon model does not
change GPP or LAI.** Both are already modelled (P-model GPP; Zhou optimal LAI or
prescribed MODIS LAI) and are treated here as *inputs*. The carbon model is a
follower: it consumes GPP, tracks pools, and returns respiration and litter.
The one deliberate exception is autotrophic respiration, which is reworked (§2.1).

**No PFTs.** Where a PFT would normally set a number (allocation fractions,
turnover times), use a global constant — optionally split C3/C4 — and only if
that proves insufficient, let it vary continuously with mean annual temperature
and mean annual precipitation. Never with a plant functional type map.

---

## 1. Pools

All in kg C m⁻² of ground. Prognostic, in `Y.canopy.biomass`.

| symbol | pool | in | out |
|---|---|---|---|
| `C_sugar` | labile / non-structural C | GPP | maintenance respiration, structure building |
| `C_leaf` | leaf structure | `a·f_leaf·S` | leaf turnover + phenological shedding |
| `C_stem` | stem + branch structure | `a·f_stem·S` | turnover |
| `C_root` | root structure | `a·f_root·S` | turnover |
| `SOC(z)` | soil organic carbon, kg C m⁻³ | leaf/stem/root litter | microbial decomposition |

`SOC` already exists as `Y.soilco2.SOC` but is currently frozen at its
SoilGrids initial condition — §4 turns it on.

```
cVeg  = C_sugar + C_leaf + C_stem + C_root
cSoil = ∫ SOC dz
```

`cVeg` replaces the current `compute_vegetation_carbon!`, which is the purely
diagnostic `σl·LAI + ηsl·h·SAI` from prescribed area indices.

### The leaf pool / LAI consistency check

`C_leaf` is prognostic here, but LAI is set by a separate model. The two are
physically linked by the specific leaf density `σl` (kg C m⁻² leaf). We do NOT
impose `C_leaf = σl·LAI` — that would make the leaf allocation fraction
meaningless. Instead we output the *implied* specific leaf density

```
σl_implied = C_leaf / LAI
```

as a diagnostic. It should land near the prescribed `σl` (~0.03–0.1
kg C m⁻² leaf, thin herbaceous leaves at the low end). If it does not, the
allocation fractions and the LAI model disagree, and that is a real, reportable
finding — not something to paper over. It is also the cleanest single number to
calibrate `f_leaf` against.

DECIDED: constant `f_leaf` with emergent SLA, as above. Do not slave `C_leaf` to
`σl·LAI`, and do not add a nudge without asking.

## 2. Fluxes

GPP is `get_GPP(p, canopy.photosynthesis)` in mol CO₂ m⁻² s⁻¹; multiply by
`M_C = 0.012011` kg mol⁻¹ to get kg C m⁻² s⁻¹.

```
dC_sugar/dt = GPP − Rm − S
dC_leaf/dt  = a·f_leaf·S − C_leaf/τ_leaf − shed
dC_stem/dt  = a·f_stem·S − C_stem/τ_stem
dC_root/dt  = a·f_root·S − C_root/τ_root
```

`S` is the total sugar drawn for structure building, `a ≈ 0.7` the construction
efficiency, `(1−a)·S` the growth respiration, and
`f_leaf + f_stem + f_root = 1`.

### 2.1 Maintenance respiration `Rm`

This replaces the JULES-style term in `AutotrophicRespirationModel`, which today
respires *prescribed constant* area indices (`Rd_ref·RAI`, `Rd_ref·μs·SAI`).
With real pools we respire the pools:

```
Rm = f_T(T_canopy)·(M_C·Rd_canopy + r_stem·C_sap) + f_T(T_soil)·r_root·C_root
C_sap = C_stem / (1 + C_stem/C_sap_half)     # living (sapwood) fraction
f_T(T) = Q10^((T − T_ref)/10)
```

- The leaf term stays `Rd_canopy`, the leaf dark respiration the photosynthesis
  scheme already computes. Reusing it avoids double-counting `Rd` (it is already
  inside `An = GPP − Rd`) and keeps leaf respiration consistent with GPP.
- The sapwood saturation is what keeps a 20 kg C m⁻² tropical forest from
  respiring itself to death: living wood does not scale linearly with total wood.
- The root term should use soil temperature, not canopy temperature. **Not yet
  implemented as specified (2026-07-31):** there is no ground temperature the
  canopy can read in both configurations — `p.drivers.T_ground` exists only for
  a standalone canopy over prescribed ground, and reaching into the integrated
  model's `p.soil.T` would make the canopy tendency depend on soil `update_aux`
  ordering. Stage 1 therefore uses canopy temperature for both terms, which
  leaves the root term more seasonally variable than it should be. Revisit at
  stage 3, when the canopy→soil coupling is built and a shared ground
  temperature becomes available.

**`Q10` ships as 1.0, not 2.0 — stage 2 must set it deliberately.**
`autotrophic_respiration_Q10` in `toml/default_parameters.toml` is `1.0`, whose
description states it is "fixed at 1.0 (not calibrated), keeping the maintenance
baseline seasonally flat". So `f_T = 1^x ≡ 1`: the existing maintenance
respiration has *no* temperature response at all. Measured in the stage-0
baseline (job 6967717): across `alaska_north_slope` (262 K), `mojave_sw_us`
(281 K) and `sahel` (300 K) — a 38 K span over which `Q10 = 2` would give a 14×
range — Ra is 1.068, 1.066 and 1.094 g C m⁻² day⁻¹, i.e. constant to ±3%. At the
two zero-LAI deserts Ra equals `Rd_ref` to five decimal places.

Since `Rm` above reuses this same `Q10`, writing §2.1 while assuming the
temperature response is "already there" would produce a temperature-insensitive
`Rm` that looks correct in the code and is inert in practice. Setting `Q10` to
2.0 is a change to respiration behaviour, so it belongs to stage 2 and must be
reported as such, not slipped in.

### 2.2 Growth respiration `Rg`

`Rg = (1−a)·S`, and `Ra = Rm + Rg` — the same
`p.canopy.autotrophic_respiration.Ra` the rest of the model already consumes.
The current `Rg = Rel·max(An − Rpm, 0)` disappears: growth respiration becomes
an actual cost of actual growth, not a fixed fraction of net photosynthesis.

### 2.3 Allocation `S` — and why sugar never hits zero

The allocation *fractions* `f_leaf, f_stem, f_root` are global constants
(optionally one set for C3 and one for C4, blended by the modelled
`fractional_c3`). The allocation *rate* is what regulates the sugar pool.

Treat the sugar pool as a buffer with a target size set by the standing live
biomass — real NSC is roughly 5–20% of live biomass, and it is drawn down each
spring and refilled each autumn without ever emptying:

```
C_sugar_target = c_nsc·(C_leaf + C_sap + C_root)
x = C_sugar / C_sugar_target
S = S_max·g(x),   S_max = C_sugar/τ_alloc
```

with `g` a smooth sigmoid-like ramp: `g(x) = x^n/(1 + x^n)` (n ≈ 2–4), so

- sugar well above target ⇒ `g → 1` ⇒ growth runs at `C_sugar/τ_alloc`, drawing
  the pool back down;
- sugar well below target ⇒ `g → 0` ⇒ growth stops, and only maintenance
  respiration draws on the pool.

That is a proportional controller: it makes sugar oscillate seasonally around
its target rather than saturating at zero or growing without bound, and it needs
no hard clamp. `Rm` is *not* switched off by `g` — a plant that cannot pay
maintenance really does draw its reserves down. Genuine carbon starvation (sugar
→ 0 and staying there) is then a physical outcome the model can express, and a
signal worth reporting, not a numerical failure. The only clamp is a floor at
`C_sugar ≥ 0` in the respiration term, for the explicit step.

**How that floor is implemented (added 2026-07-31, after the stage-1 battery).**
Reserves draw down to zero, not below: respiration requires substrate. `Rm` is
multiplied by the same ramp applied to an absolute scale,

```
Rm ← g(C_sugar/C_sugar_ref)·Rm,     C_sugar_ref = 1e-3 kg C m⁻²
```

so maintenance stops as the pool empties. This is distinct from `g(x)` in the
allocation rate, which is keyed to the *target*: a stressed plant with sugar
below target still pays full maintenance and draws down: only an *exhausted*
plant stops. `C_sugar_ref` is small enough that a healthy pool
(0.05–0.2 kg C m⁻²) is unaffected to about one part in 10⁵.

This was not hypothetical. The first stage-1 battery reached
**C_sugar = −0.075 kg C m⁻²** at the Sahara point, where GPP is zero all year
and maintenance kept drawing on an empty pool.

If a global constant `f` set proves insufficient — most likely symptom: forests
under-built and grasslands over-built at the same time — the fallback is to make
`f_stem` (and correspondingly `f_root`) vary with mean annual temperature and
mean annual precipitation, both of which are already available as trailing
integrals (§3). Try constants first, and report what forced the change.

### 2.4 Turnover and litter

```
L_leaf = C_leaf/τ_leaf  +  C_leaf·max(−dLAI/dt, 0)/LAI      → soil surface layer
L_stem = C_stem/τ_stem                                       → soil surface layer
L_root = C_root/τ_root       → distributed by root_distribution(z, rooting_depth)
```

The second leaf term is phenological shedding: a *relative* LAI loss rate
applied to the leaf pool. Without it a deciduous canopy would drop its LAI in
autumn while its leaf carbon sat there. With it, `C_leaf` inherits LAI's
seasonality without being rigidly slaved to it. `dLAI/dt` comes from the LAI
model's own tendency (`dY.canopy.biomass.LAI` for `ZhouOptimalLAIModel`);
for prescribed LAI it must be finite-differenced or read from the
`TimeVaryingInput` — see the open questions.

The background term `C_leaf/τ_leaf` matters independently: an evergreen tropical
canopy has near-constant LAI but a ~1.5 yr leaf lifespan, and without it there
would be no leaf litter in the tropics at all.

## 3. Climate proxies, if and when needed

Start with global constants (C3/C4 split where it helps). If turnover times or
allocation fractions must vary, use trailing climate, never PFTs:

- mean annual temperature — a `RunningMean` of air temperature;
- mean annual precipitation — `precip_annual` already exists in
  `ZhouOptimalLAIModel`;
- growing-season length and aridity — `growing_days`, `PET_annual`, also already
  there;
- C4 fraction — `1 − fractional_c3`, modelled dynamically on this branch.

Reuse the existing `TimeIntegratedVariable`s rather than duplicating them (see
commit `05773316d` on duplicate time-integrated variables); declare new ones
only for the prescribed-LAI case, where the Zhou model is absent.

## 4. Soil organic carbon

`Y.soilco2.SOC` is prognostic but pinned: `make_compute_exp_tendency(::SoilCO2Model)`
sets `dY.soilco2.SOC = 0`, and `MicrobeProduction` adds `Sm` to CO₂ and
subtracts O₂ but never debits SOC. Turning SOC on means:

```
dSOC/dt = I_litter(z) − Sm
```

- `Sm` (kg C m⁻³ s⁻¹) is the existing `microbe_source`. Debiting it closes a
  carbon balance that is currently open.
- `I_litter(z) = (L_leaf + L_stem)/Δz_top` in the top layer, plus
  `L_root·root_distribution(z, rooting_depth)` through the column.
- In standalone `SoilCO2Model` (no canopy) litter input is zero or prescribed;
  the coupling is made in the integrated model, the same way root water
  extraction is coupled in `soil_canopy_root_interactions.jl`.

**SOC does not equilibrate on a runnable timescale** — τ is centuries to
millennia. Plan for analytical spinup (§6), not brute force.

## 5. Coupling direction — one-way first (DECIDED)

Phase 1 is one-way: `SAI`, `RAI` and canopy height stay prescribed.

With `SAI`, `RAI` and canopy height left prescribed, the carbon model is a pure
follower: nothing it computes feeds back into LAI, GPP, or the water and energy
budgets. That is what makes §6 cheap, and it makes rule 1 verifiable — a run
with the carbon model must reproduce the same GPP and LAI as one without it,
bit for bit apart from the respiration change.

Phase 2, only after phase 1 validates, can close the loop by deriving the area
indices and height from the pools:

```
h   = h_max(climate)·(1 − exp(−c_h·C_stem))
SAI = C_stem / (ηsl·h)          # inverts the existing cStem = ηsl·h·SAI
RAI = C_root / σ_root
```

Scientifically this is the right destination — constant SAI/RAI/height are the
weakest assumption in the current canopy — but it feeds back on radiation,
roughness and root water uptake, hence on GPP. It is a separate, explicit step,
gated on the user's approval.

## 6. Reaching steady state

Because of §5, equilibrium can be found **offline and cheaply**:

1. Run the coupled model a few years, saving the drivers the carbon model needs:
   GPP, LAI, T_canopy, T_soil, θ_l, and the C3/C4 fraction.
2. Integrate the pool ODEs alone against that record, recycled, for as many
   centuries as stem and SOC need. Seconds, not node-hours.
3. Or solve the steady state directly — the system is linear in the pools once
   the mean fluxes are known: `C_i* = a·f_i·S̄·τ_i`, `SOC*(z) = Ī_litter(z)/k̄(z)`.
4. Write the result as an initial-condition file and restart the coupled run from
   it — the same pattern the LAI pipeline used to rebuild `optimal_lai_inputs.nc`.
   Verify the artifact override actually resolves before trusting it.

Brute-force spinup inside the coupled model is the fallback, not the plan.

**Forcing to recycle (DECIDED):** the earliest available ERA5 decade, looped.
Note the constraint discovered while checking: the ClimaLand ERA5 forcing
artifact (`forty_yrs_era5_land_forcing_data`) starts in **1979**, not 1940 — so
the earliest decade we can actually recycle is **1979–1989**. That is the least
climate-changed forcing available, not a pre-industrial one. Atmospheric CO₂
should be held at that era's value (~340 ppm) rather than present-day, or the
"steady state" bakes in modern CO₂ fertilisation; check how `c_co2` is set before
relying on it.

## 7. Validation

Target: global `cVeg` of the right magnitude and pattern — roughly 10–20 kg C m⁻²
in tropical forest, 5–15 boreal/temperate forest, 0.3–1 grassland, ~0 desert.
Secondary checks: `cSoil` against HWSD and against SoilGrids (where SOC started);
NPP/GPP ≈ 0.4–0.5 in the global annual mean; `σl_implied` (§1) against the
prescribed `σl`; root:shoot ratio, a well-observed constraint that falls straight
out of `f_root`/`τ_root`.

### Reference data — already on Derecho (verified 2026-07-31)

The `ilamb_data` ClimaLand artifact holds only ET, GPP and two radiation fields —
**no biomass**. The full ILAMB reference tree is staged at NCAR, though:

```
/glade/campaign/cesm/community/lmwg/diag/ILAMB/DATA/
  biomass/XuSaatchi2021/XuSaatchi.nc      0.5°, 2000–2019, Mg ha⁻¹ (×0.1 → kg C m⁻²)
                                          "annual carbon density of global live
                                          woody vegetation" — the primary target
  biomass/GEOCARBON/biomass.nc            pan-tropical + boreal AGB
  biomass/ESACCI/biomass.nc               0.5°, 2010 epoch
  biomass/{Thurner,Saatchi2011,Tropical,GLOBAL.CARBON,NBCD2000,USForest}/
  soilc/HWSD_M/soilc_0.5x0.5.nc           cSoilAbove1m, kg m⁻²
```

Read them directly from that path (read-only, no download, no sandbox network
needed). **Caveat that matters for the comparison:** these products are *woody*
and mostly *aboveground* biomass. Compare the woody part of the model
(`C_stem`, plus coarse root if separable) to XuSaatchi, not total `cVeg`, and
treat grassland cells — where the observation is near zero by construction but
the model legitimately carries leaf and root carbon — separately rather than as
model error.

## 8. Parameters (first cut)

| parameter | meaning | first guess |
|---|---|---|
| `a` | construction efficiency | 0.7 |
| `f_leaf, f_stem, f_root` | allocation fractions (C3) | 0.3 / 0.4 / 0.3 |
| `f_leaf, f_stem, f_root` | allocation fractions (C4) | 0.4 / 0.05 / 0.55 |
| `τ_leaf` | leaf lifespan | 1.5 yr |
| `τ_stem` | stem turnover | 30 yr (C3), 1 yr (C4) |
| `τ_root` | root turnover | 2 yr |
| `r_stem` | sapwood specific maintenance respiration | to fit |
| `r_root` | root specific maintenance respiration | to fit |
| `C_sap_half` | sapwood saturation scale | ~2 kg C m⁻² |
| `c_nsc` | target sugar as fraction of live biomass | 0.1 |
| `τ_alloc` | allocation timescale | ~10 d |
| `n` | allocation ramp sharpness | 3 |
| `Q10`, `T_ref` | temperature sensitivity | **2.0** (must be set), 298.15 K |

Realistically 4–6 of these set global biomass: `τ_stem`, `f_stem`, `f_root`,
`r_stem`, and `a`. Those are the calibration targets; the rest should be fixed
at literature values and left alone.

---

## 9. What was tried and rejected

Written after the fact. The specification above says what the model is; this says
what it is *not*, and why — so the reasoning survives without reading the
iteration log in `unsupervised_loop/STATE.md`.

### The central problem

The model builds forests correctly and builds them where there are none. Binned
by observed woody carbon over ~20 000 cells, model/observed is ~0.8–1.3x where
observations show real forest and 8–16x where they show none, monotone across
three independent products.

Localised further: the excess is **not** spread across treeless land. At the
current ramp, 6% of treeless cells — the ones whose rainfall is *forest*
rainfall — carry 45% of that bin's modelled carbon. Cerrado, the Sahel fringe,
the Pampas, miombo. They are treeless because of fire, grazing and cultivation.

### Predictors tested and rejected

1. **Site-level precipitation.** Cannot separate boreal forest from warm
   grassland. (Globally a MAP ramp *does* help — that was a correction to an
   earlier claim, and the ramp is in the model.)
2. **Soil properties.** Only organic matter separates, and tundra peat has the
   highest OM of all.
3. **Aridity, P/PET.** Orders wrongly: Alaskan tundra scores most humid, Congo
   in the dry half.
4. **`f_stem(MAT, MAP)`** via a temperature-dependent half-point. Fixes the
   boreal and monotonically *degrades* spatial correlation in all three
   products, because it re-scales a whole temperature band.
5. **Sharpening the MAP ramp exponent** (`n` = 4, 8). Works exactly as intended
   on the cells it can reach — dry-cell biomass halves toward the observations —
   and is then swamped: the same steepening pushes the wet minority further up
   until 6% of cells hold two-thirds of the bin's carbon.
6. **A prescribed burned-area map.** Withdrawn on principle rather than on
   performance: a GFED map is as unextrapolatable as a PFT map, and the point of
   this model is to run in a climate that has not been observed.

### What worked

**Rainfall seasonality** — the annual water deficit, not the annual total. Two
cells with the same MAP differ by an order of magnitude in it, and that is what
separates wet savanna from wet forest. Validated as a contingency table *before*
any model change: medians differ 3–4x with barely overlapping interquartile
ranges, in all five products carrying both classes.

The half-point `carbon_deficit_half_woody = 0.55 m` has since survived a GPP
recalibration, an optimal-LAI recalibration, and a merge that replaced the
aridity index — three independent perturbations of the drivers it was calibrated
against.

### Corrections made along the way

Recorded because each invalidated a published number:

- **Four of the six biomass products report dry matter, not carbon** — visible
  only in `long_name`, not `units`. The model is biased high against *all* six;
  the previously reported near-zero ESACCI bias and negative GEOCARBON bias were
  artefacts.
- **A two-month phase error** in the driver climatology: the monthly NetCDFs
  carry no reference date, so record 1 is the run's first month, not January.
  Harmless for equilibria (a uniform shift of a recycled cycle) but wrong the
  moment a model month meets an observational one.
- **Three unit errors of the same family** — `pra` is an annual total not a rate,
  `gpp` and `hr` are mol CO2 not kg C. The last gave a global Rh fifty times any
  physical value. Conversions are now read from the file's `units` attribute.
- **"Improves every product" was withdrawn**: 6/6 became 5/6 after the phase fix,
  on a margin smaller than the error that produced it.
