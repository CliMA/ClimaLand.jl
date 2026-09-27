# Surface flux investigation notes (session 32a00791)

## Harness
- `scratch/site_run.jl SITE TAG` (run from a checkout, `--project=.buildkite`) -> `scratch/series/SITE_TAG_model.csv`, `SITE_obs.csv`
  - `CL_TOML=file` parameter overrides; `CL_PMODEL=1` = global-default P-model + PModelConductance + PiecewiseMoistureStress
- `scratch/site_eval.jl TAG...` metrics; `scratch/site_diurnal.jl SITE TAG...` May–Sep diurnal cycle
- `scratch/site_budget.jl SITE TAG` model diurnal energy budget; `scratch/obs_budget.jl SITE UTC_OFF YEAR` raw FLUXNET budget
- `scratch/gpp_cmp.jl SITE TAG` midday monthly GPP
- `scratch/canopy_night.jl`, `canopy_night_cap.jl` offline SurfaceFluxes replication of big-leaf canopy flux
- `scratch/snow_run.jl SITE TAG` SnowMIP metrics (SWE, Tsnow)
- Model timestamps may be offset 15 min (US-MOz); evaluator floors to 30 min.
- US-Ha1 has no LE/H/SWU/LWU obs in comparison data (only gpp).

## Commits (branch ts/Tsfc-fix, pushed to origin/ts/Tsfc-fix)
- `8b02be249` Solve for soil skin temperature with half-cell and litter thermal resistance
- `500f007f4` Include porous-media tortuosity factor in `soil_conductance`
- `9360ca397` Wire soil skin temperature solve into `SoilCanopyModel` and `SoilSnowModel`
- `84b9ffb55` Fix off-diagonal soil Jacobian block `∂(ρe_int residual)/∂ϑ_l`
- `165161092` Soil evaporation (SL14 DSL, adsorption guard) and CLM5-style under-canopy exchange
- `669519dd1` Limit the stable Monin-Obukhov stability parameter (CLM `zetamaxstable`)
- `ad74ddac9` Use plant area index (`LAI + SAI`) for canopy radiation and heat exchange
- `1c347feb4` Apply P-model moisture stress instantaneously instead of via acclimation
- `a7c87a332` Use the Tuzet leaf-water-potential moisture stress by default in `LandModel` and `SoilCanopyModel`
- `39b130ca1` Add CLM5-style canopy interception of rain and wet-canopy evaporation (opt-in)
- `bd1b2a238` Solve snow surface temperature before snow-soil ground heat flux in `SoilSnowModel`
- `2fcf19cf3` Shut off soil vapor conductance smoothly at residual saturation and use `θ_l` in implicit aux
- `66fcc6b22` Enable CLM5 canopy interception by default in `LandModel` and guard low-`PAI` canopies

## Results (Farquhar/Medlyn site config, vs main)
| site | var | main bias/rmse/r | v4 bias/rmse/r |
|---|---|---|---|
| US-MOz | shf | -51/137/0.23 | -45/121/0.52 |
| US-MOz | lhf | 21/84/0.90 | 14/76/0.90 |
| US-NR1 | lhf | -42/82/0.52 | -40/77/0.66 |
| US-Var | shf | 5.1/40/0.955 | 1.9/36/0.958 |
| US-Var | lhf | 7.3/32/0.91 | 4.7/30/0.92 |
| US-Var | lwu | -11.7/40 | -5.9/43 (amp 1.19->1.39) |

## Finding 1: Farquhar/Medlyn site config GPP 2.4x obs at US-MOz
- => canopy LE too high, Bowen ~0. With P-model (global default) midday GPP matches obs May–Sep
  (e.g. Jun 2.02e-5 vs 2.16e-5). Iterate with CL_PMODEL=1 from here on.
- P-model at US-Var (Mediterranean grass) overestimates LE strongly (bias +41, amp 2.1) – moisture
  stress / P-model issue, not surface-energy physics; note for later.

## Finding 2: nocturnal runaway decoupling (MOST supercritical collapse)
- Offline: U=2 m/s, Tc-Ta=-3 K -> ζ=55, u*=0.016, H≈-0.7 W/m²; U=1 -> ζ=100 always.
- Model night Tc-Ta ≈ -5.5 K, cshf ≈ -8 (obs total H ≈ -19), LWu ~13-17 W/m² too low.
- Fix: StableLimitedUniversalFunctionParams (ζ ≤ zeta_max_stable in flux-profile relations; CLM zetamaxstable).
- Reviewer (3f4b555a): approve-with-changes; do sensitivity; check snow; consider biomass heat storage; c_h forwarded.

| run (P-model) | MOz shf | MOz lhf | MOz lwu | NR1 shf | NR1 lhf | NR1 lwu |
|---|---|---|---|---|---|---|
| pm_main | -23/101/0.68 | -11/79/0.88 | 2.3/15.2/0.970 | -15/97 | -33/69 | -6.3/20.2 |
| pm_v4 (ζ∞) | -13/89/0.74 | -22/81/0.91 | 2.8/22.3/0.941 | -21/98 | -31/65 | -2.9/26.4 |
| pm_z2 | -21/89/0.75 | -19/80/0.91 | 7.4/18.0/0.966 | -28/96 | -29/64 | 3.4/19.7 |
| pm_z05 | -26/91/0.74 | -17/80/0.90 | 10.9/16.9/0.979 | -32/96 | -26/64 | 6.4/18.1 |

- MOz night (May–Sep): obs H -19, NETRAD -45..-51, LWout 406–432.
  z05: H -35..-47, rn -47..-52, LWu 411–434, Tc-Ta -1.6..-1.9. z2: H -22..-30, rn -38..-42, LWu 403–427, Tc-Ta -3.2..-3.8.

## Remaining problems (US-MOz daytime)
- Soil absorbs ~142 W/m² net radiation at noon; soil G ≈ 70 W/m² (obs G ≈ 7-10); daytime LWu 471 vs obs 443-450.
- Radiative transfer (SW two-stream, LW emissivity) uses LAI only; SAI (=1 at MOz) ignored (CLM uses LAI+SAI).
  Candidate: plant-area-index RT with stem optics, stems in canopy heat capacity & heat exchange, APAR scaled by LAI/PAI.
- Daytime model LE 270 vs obs LE_CORR 423 (LE_F_MDS 272); H 186 vs 157 (H_F_MDS 101).

## Commit `669519dd1`: `zeta_max_stable = 0.5` (chosen; SnowMIP `cdp` SWE RMSE `0.080 -> 0.039`)

## Finding 3: US-MOz winter SHF badly underestimated (obs Nov–Apr 35–87 W/m², model -18..25)
- Leafless forest has no radiative canopy (SAI ignored). => PAI change (`ad74ddac9`):
  RT/emissivity/heat capacity/sensible heat use `LAI + SAI` with CLM stem optics; leaf `APAR = par.abs * LAI / PAI`;
  transpiration uses `LAI`. Global `SAI = 0` => no global change where `SAI = 0`.

## Finding 4 (CRITICAL): plant water non-conservation with P-model
- βm enters P-model only via acclimated Vcmax25/Jmax25 (τ = 1/pmodel_α ≈ 36 d); instantaneous gs
  does not respond to water stress. Transpiration not supply-limited.
- US-Var (Piecewise): trans 2–4 mm/d Jul–Sep, root uptake 0, ϑ_l < 0 unbounded, ψ clamped -2041 m,
  Weibull K_plant = 0 permanently. βm 0.4–0.5 all summer (θ_low = θ_r = 0). LE +100 W/m² vs obs.
- US-NR1 (P-model + Tuzet, pmt_z05): frozen soil in winter -> ψ falls -> K_plant -> 0 permanently;
  summer clhf 2–8 (vs 23–31 with Piecewise); lhf rmse 64 -> 75.
- Fix (`1c347feb4`): acclimate unstressed capacities, apply βm instantaneously to GPP/Rd/An/gs (Stocker 2020).
  Reviewer e8839695 consulted.

## Commits `ad74ddac9` (PAI) and `1c347feb4` (instantaneous βm)
| run | MOz shf | MOz lhf | NR1 shf | NR1 lhf | Var shf | Var lhf |
|---|---|---|---|---|---|---|
| pm_z05 | -26/91 | -17/80 | -32/96 | -26/64 | -31/82 | 42/101 (r .53) |
| pm_pai | -10/80 | -28/90 | -9/75 | -29/69 | – | – |
| pm_ib (Piecewise) | -9/81 | -28/91 | -8/75 | -30/69 | -26/72 | 35/87 |
| pmt_ib (Tuzet) | -12/79 | -25/86 | -2/77 | -37/79 | -1/38 | 4/32 (r .90) |
- pmt_ib: plant water conserved at all sites (root uptake = transpiration). Piecewise still leaks at Var/NR1.
- NR1 low LE: site K_sat_plant = 5e-9 (14x below default) and MODIS LAI 1.4 (in situ ~4) limit supply.
- MOz MODIS LAI winter 0.4–0.6 (already stems/understory?) ; summer 3.8.

## Commit `a7c87a332`: LandModel/SoilCanopyModel default soil moisture stress -> Tuzet
- US-Ha1 (GPP obs only): Jun-Jul midday GPP 2.15/2.06e-5 (Piecewise) -> 2.30/2.25e-5 (Tuzet), obs 3.2/3.1e-5;
  ψ ≥ -56 m, root uptake = transpiration. full_land (412), soil_canopy_lsm, conservation tests pass.

## Commits `39b130ca1` and `66fcc6b22`: CLM5 canopy interception (liquid, default in `LandModel`)
- Reviewer f5dc9ae4: design matches CLM5 FracWet/CanopyFluxes. Fixed: lake throughfall (water/energy leak),
  f_wet evaluated after interception (CLM order), check in keyword ctor. Accepted/documented: rain enthalpy
  ρe_l(T_air)(I-D) into canopy (root-uptake convention), dew switch discontinuity in ∂q/∂T, Δt must match.
| run | MOz shf | MOz lhf (mrmse) | NR1 shf | NR1 lhf | Var shf | Var lhf |
|---|---|---|---|---|---|---|
| pmt_ib | -12/79 | -25/86 (35.6) | -2/77 | -37/79 | -1.1/38.0 | 4.0/32.2 |
| pmti (f_wet_max .05) | -15/78 | -22/85 (32.3) | -3/76 | -36/79 | -1.5/37.6 | 4.2/32.4 |
| pmti1 (f_wet_max 1) | -17/79 | -20/87 (31.1), r .897 | | | | |
- MOz wet-canopy evaporation ≈ 5–6 W/m² annual (~70 mm/yr, ~8% of ET); literature 10–20% of P for
  deciduous forest -> CLM5 default is on the low side but f_wet_max=1 degrades rmse/r.

## MOz summer diagnosis (May–Sep, pmti)
- Midday: model H 228 / LE 256 / Tc-Ta 2.9 K / LWu 476; obs H_F_MDS 90–100 (corr 139–157), LE_F_MDS 272–294
  (corr 423–457), LWout 450. Bowen: model 0.9, obs 0.3 (both corrected and uncorrected).
- Site Farquhar+Medlyn config (fmi): summer LE matches (Jun 195 vs 187) but GPP ~2x obs, April LE 143 vs 63,
  SHF bias -51. P-model: GPP -23% (Jun midday), LE -35%. => remaining summer LE deficit is a photosynthesis /
  stomatal-conductance calibration issue (P-model), not surface exchange.
- Leaf-off (Feb–Apr, Oct–Nov): night model H -45..-65 vs obs H_F_MDS -15..-18 (corr -26..-30); obs night
  closure gap ~35 W/m² -> part of the "winter SHF deficit" is obs night-flux underestimation.

## Global 2-year low-res GPU runs (`lr_main`, `lr_branch`, `lr_int`, `lr_favail`)
- `lr_main` (`main` @ `0d2fa0ef9`), `lr_branch` (`ts/Tsfc-fix` @ `39b130ca1`, `NoInterception`), `lr_int`
  (`CLM5Interception` default in `LandModel`), and `lr_favail` (all fixes on `ts/Tsfc-fix`, including low-`PAI`
  interception guard and smooth residual-saturation soil vapor conductance shutoff) run on GPU (`180x360x15`,
  `Δt = 900 s`, 2008-03-01 to 2010-03-01).
- In `lr_favail`, **0 cells develop NaNs anywhere on Earth over all 24 months** (identical valid-cell mask to `main`).
- Full-domain global skill (`lr_main` -> `lr_favail`, year-2 monthly climatology, W/m²):
  - **LHF vs ERA5**:
    - Global (`-60..90°`): bias `+3.35 -> +0.54` (6x reduction), `rmse_ann` `18.98 -> 19.82`
    - Tropics (`-23..23°`): bias `+6.39 -> +1.86`
    - N. midlatitudes (`23..50°N`): bias `+2.17 -> +0.66`, `rmse_ann` `16.37 -> 15.82`, `rmse_mon` `28.01 -> 27.59`, `rmse_seas` `22.73 -> 22.60`
    - Boreal (`50..90°N`): `rmse_ann` `9.51 -> 8.14`, `rmse_mon` `20.22 -> 20.16`
    - S. midlatitudes (`-60..-23°`): bias `+6.06 -> +2.15`
  - **LHF vs CLASS**:
    - Global: bias `+12.78 -> +9.85`, `rmse_ann` `21.20 -> 20.78`, `rmse_mon` `30.44 -> 30.15`, `rmse_seas` `21.84 -> 21.84`
    - N. midlatitudes: bias `+12.09 -> +10.37`, `rmse_ann` `17.68 -> 16.37`, `rmse_mon` `29.05 -> 27.76`, `rmse_seas` `23.05 -> 22.42`
    - Boreal: bias `+7.40 -> +5.73`, `rmse_ann` `12.33 -> 9.71`, `rmse_mon` `23.35 -> 22.30`
    - S. midlatitudes: bias `+18.22 -> +14.14`, `rmse_ann` `26.18 -> 24.85`, `rmse_mon` `34.11 -> 33.21`
  - **LHF vs FLUXCOM**:
    - Global: bias `+1.49 -> -1.62`, `rmse_ann` `21.60 -> 23.03`
    - N. midlatitudes: bias `+4.13 -> +2.40`, `rmse_ann` `16.73 -> 16.32`, `rmse_mon` `29.00 -> 28.66`, `rmse_seas` `23.69 -> 23.57`
    - Boreal: bias `+3.22 -> +1.58`, `rmse_ann` `10.22 -> 7.71`, `rmse_mon` `18.44 -> 17.73`
    - S. midlatitudes: bias `+13.73 -> +9.80`, `rmse_ann` `24.22 -> 23.49`, `rmse_mon` `32.90 -> 32.50`
  - **SHF vs FLUXCOM / CLASS / ERA5**:
    - Seasonal-cycle RMSE (`rmse_seas`) improves across all three references globally and in extratropics:
      - FLUXCOM `rmse_seas`: global `21.29 -> 20.59` (`nmid` `21.63 -> 20.45`, `bor` `22.53 -> 20.86`, `smid` `20.95 -> 20.43`)
      - CLASS `rmse_seas`: global `19.36 -> 18.87` (`nmid` `17.74 -> 16.44`, `bor` `19.65 -> 18.53`; `bor` `rmse_ann` `10.91 -> 9.85`, `rmse_mon` `22.48 -> 20.99`)
      - ERA5 `rmse_seas`: global `21.14 -> 20.64` (`trop` `19.66 -> 19.45`, `nmid` `20.98 -> 20.20`, `bor` `22.65 -> 21.87`, `smid` `23.57 -> 23.55`; `bor` `rmse_ann` `13.04 -> 12.38`, `rmse_mon` `26.14 -> 25.13`; `smid` bias `+1.84 -> +0.14`)
  - **LWu vs ERA5**:
    - Global: bias `-6.45 -> -1.58` (4x reduction), `rmse_ann` `28.62 -> 27.71`, `rmse_mon` `39.14 -> 38.83`
    - N. midlatitudes: bias `-21.27 -> -15.91`, `rmse_ann` `41.36 -> 37.89`, `rmse_mon` `52.00 -> 49.48`
    - Boreal: bias `-11.16 -> -7.59`, `rmse_ann` `17.70 -> 15.61`, `rmse_mon` `36.78 -> 36.25`
  - **SWu vs ERA5**:
    - Global: bias `+0.69 -> +0.26`, `rmse_ann` `12.27 -> 12.12`, `rmse_mon` `22.29 -> 22.18`, `rmse_seas` `18.60 -> 18.57`

## Finding 5: Two localized NaN mechanisms in global runs (diagnosed and fixed in `2fcf19cf3` and `66fcc6b22`)
1. **High-latitude canopy interception NaNs (`lr_int` vs `lr_branch`)**:
   - In `lr_int`, ~85 high-latitude grid cells (`55–75°N`) developed `NaN` in `Y.canopy.energy.T` when seasonal
     MODIS `PAI` dropped toward `0` (`0 < PAI < 0.05`).
   - Root cause: `canopy_energy.jl` shuts off canopy sensible heat flux (`cshf = 0`) for `PAI < 0.05` where canopy
     heat capacity `ac_canopy = cw * Ω * PAI` is tiny (`-> 0`), whereas `CLM5Interception` remained active down to
     `PAI = 0`. Dividing the interception enthalpy flux `F_int = ρ_l e_l(T_air) (I - D)` or dew latent heat by
     `ac_canopy` caused explicit Euler temperature spikes (`> 550 K`) in `Y.canopy.energy.T`.
   - Fix (`src/standalone/Vegetation/interception.jl`, `src/integrated/land.jl`): Added `canopy_is_active(PAI) = PAI >= 0.05`
     and gated `intercepted_liq`, storage capacity in `drip`, `wetted_fraction`, `condensing` (dew-to-storage), and
     `interception_energy_flux` below `PAI = 0.05` (storing any residual store when `PAI` crosses below `0.05` as
     immediate drip). Verified in 2-year global GPU run (`lr_favail`): 0 high-latitude interception NaNs.
2. **Sahara/Sahel desert soil evaporation NaNs (`lr_branch`)**:
   - Reproduced in single column (`scratch/column_nan.jl`, `CL_LON=13 CL_LAT=17`, both `Float32` and `Float64`):
     top-cell `Y.soil.ϑ_l` drained from `0.08` through `θ_r = 0.0527` down to `-0.25`, causing top-cell heat
     capacity to go negative and `T_soil` to blow up (`379 K -> -91.6 K -> DomainError`).
   - Root cause:
     (a) With Swenson & Lawrence (2014) `DSL = d_ds ((αS_c - S_l)/(αS_c))^p`, `DSL` saturates at `d_ds = 0.015 m`
         as `S_l -> 0` (unlike old `((αS_c - S_l)/S_l)^p` which diverged to `∞`), so `g_soil ≈ 2.4e-4 m/s` stays finite.
     (b) `effective_saturation` clips `ϑ_l_safe = max(ϑ_l, θ_r + sqrt(eps(FT)))`. In `Float32` (`sqrt(eps) = 3.45e-4`),
         `S_min ≈ 1e-3`, where van Genuchten `ψ ≈ -33 m` and Kelvin `hr = exp(g ψ M_w / (R T)) ≈ 0.998` (even in
         `Float64`, `hr ≈ 0.14–0.993` at `ϑ_l = θ_r`). Thus neither `g_soil` nor `hr` shuts off evaporation when
         mobile water `ϑ_l - θ_r` is exhausted!
     (c) `make_update_implicit_aux` (`energy_hydrology.jl:442`) used `min(ν - θ_i, ϑ_l)` instead of
         `volumetric_liquid_fraction(ϑ_l, ν - θ_i, θ_r)`, allowing negative `ϑ_l` to produce negative heat capacity.
   - Fix (reviewed by physics subagent `5f7c530e`, committed in `2fcf19cf3`):
     (a) `make_update_implicit_aux` now updates and uses `p.soil.θ_l = volumetric_liquid_fraction(Y.soil.ϑ_l, ν - Y.soil.θ_i, θ_r)`,
         consistent with `make_update_aux` and `compute_jacobian!`.
     (b) `soil_surface_vapor_conductance!` computes unclipped `S_l_sfc = clamp((min(θ_l_sfc, ϑ_l_top) - θ_r_sfc) / (ν_sfc - θ_r_sfc), 0, 1)`
         from `Y.soil.ϑ_l`, and `soil_conductance` modulates vapor conductance by the smooth $C^1$ liquid availability
         factor `f_avail = S_l^2 / (S_l^2 + S_0^2)` (`S_0 = 0.01`) so `g_soil -> 0` smoothly as `ϑ_l -> θ_r`.
   - Verified in desert column (`CL_LON=13 CL_LAT=17 CL_FT=Float32`), FLUXNET sites (`US-Var` LHF bias/RMSE
     `+4.20/32.4 -> +4.12/31.88 W/m²`, SHF RMSE `37.6 -> 37.46 W/m²`), and 2-year global GPU run (`lr_favail`,
     0 new NaNs globally over 24 months).

## Finding 6: Surface temperature wiring audit across configurations (committed in `bd1b2a238`)
- Soil skin temperature (`p.soil.T_sfc`) is live and wired into standalone `EnergyHydrology`, `SoilSnowModel`,
  `SoilCanopyModel`, and `LandModel`.
- Found and fixed ordering bug in `SoilSnowModel` (`src/integrated/soil_snow_model.jl:186-210`):
  `make_update_boundary_fluxes` previously called `update_soil_snow_ground_heat_flux!` (which uses `p.snow.T_sfc` in
  `snow_T_bottom`) before `Snow.update_surf_temp!`. Now `Snow.update_surf_temp!` is called before
  `update_soil_snow_ground_heat_flux!`, matching `LandModel`.

## Finding 7: Data-adaptive litter thermal resistance (`r_litter * (LAI + SAI)`, amended into `8b02be249` and `9360ca397`)
- **Problem with previous formulation (`r_litter * (1 - exp(-(LAI + SAI)))`)**:
  - Used the 2D horizontal optical cover fraction `1 - exp(-PAI)`, which saturates near `1` for `PAI ≳ 2`.
  - Consequently, a dense deciduous forest (`US-MOz`, summer `PAI = 4.83`) received the same litter thermal resistance
    (`0.099 m² K/W` at default `r_litter = 0.1`) as a sparse canopy, leaving summer midday ground heat flux at
    `US-MOz` too large (`G ≈ 58 W/m²` vs observed `7–10 W/m²`), while `US-Var` (`PAI ≈ 1.6–2.2`) required a
    site-specific hardcoded override `r_litter = 0.2` in `US-Var.jl` and `vaira_paper.jl` to reach `r_litter_eff ≈ 0.16–0.18 m² K/W`.
- **Principled fix (reviewed by physics subagent `4c8d43a9`)**:
  - Vertical heat conduction through stacked litter/thatch layers is a 1D series thermal resistance
    `r_litter_eff = d_litter / κ_litter`, where litter thickness `d_litter` (or litter area index `L_litter = f_litter * PAI`)
    scales linearly with plant area index `PAI = LAI + SAI`.
  - Replaced `r_litter * (1 - exp(-(LAI + SAI)))` with `r_litter * (LAI + SAI)` in `LandModel` and `SoilCanopyModel`,
    with a single global calibratable parameter `litter_thermal_resistance = 0.08` (`m² K/W` per unit `PAI`) in
    `default_parameters.toml` and `uncalibrated_parameters.toml`, and removed the site-specific `r_litter` overrides
    from `US-Var.jl` and `vaira_paper.jl`.
- **End-to-end FLUXNET site evaluation (P-model global default config, 0 site overrides for `r_litter`)**:

| Site | Var | `pm_main` (bias / RMSE / `r` / amp / mRMSE) | `favail` (bias / RMSE / `r` / amp / mRMSE) | `rlitter_pai` (bias / RMSE / `r` / amp / mRMSE) |
|---|---|---|---|---|
| **US-MOz** | `shf` | `-23.1 / 101.0 / 0.676 / 0.49 / 33.1` | `-15.3 / 77.6 / 0.814 / 1.04 / 30.2` | `-14.5 / 78.7 / 0.812 / 1.09 / 30.4` |
| **US-MOz** | `lhf` | `-10.9 / 79.3 / 0.879 / 0.92 / 23.0` | `-21.6 / 84.9 / 0.911 / 0.69 / 32.3` | `-22.1 / 83.1 / 0.912 / 0.73 / 33.1` |
| **US-MOz** | `swu` | `+4.7 / 20.5 / 0.902 / 1.32 / 10.7` | `-2.1 / 8.1 / 0.968 / 0.96 / 7.3` | `-2.1 / 8.1 / 0.968 / 0.96 / 7.3` |
| **US-MOz** | `lwu` | `+2.3 / 15.2 / 0.970 / 1.31 / 3.8` | `+10.6 / 15.2 / 0.986 / 1.29 / 10.6` | `+10.5 / 15.7 / 0.984 / 1.33 / 10.5` |
| **US-NR1** | `shf` | `-15.4 / 97.3 / 0.851 / 0.73 / 24.2` | `-2.9 / 75.7 / 0.916 / 1.10 / 18.1` | `-2.8 / 76.5 / 0.915 / 1.12 / 18.1` |
| **US-NR1** | `lhf` | `-32.6 / 69.3 / 0.678 / 0.55 / 33.4` | `-35.7 / 78.5 / 0.552 / 0.33 / 38.5` | `-35.8 / 78.1 / 0.564 / 0.34 / 38.6` |
| **US-NR1** | `swu` | `+26.0 / 61.8 / 0.844 / 2.44 / 30.3` | `+2.9 / 15.0 / 0.919 / 1.18 / 4.8` | `+2.9 / 15.1 / 0.919 / 1.18 / 4.8` |
| **US-NR1** | `lwu` | `-6.3 / 20.2 / 0.937 / 1.71 / 7.5` | `+7.1 / 13.0 / 0.982 / 1.47 / 7.4` | `+7.0 / 13.4 / 0.980 / 1.49 / 7.4` |
| **US-Var** | `shf` | `-24.8 / 72.3 / 0.880 / 0.67 / 34.5` | `-1.6 / 37.5 / 0.955 / 0.96 / 14.3` | `-1.2 / 37.2 / 0.957 / 1.01 / 14.1` |
| **US-Var** | `lhf` | `+42.3 / 86.9 / 0.650 / 2.02 / 57.5` | `+4.1 / 31.9 / 0.904 / 1.10 / 14.4` | `+4.2 / 32.3 / 0.908 / 1.14 / 14.2` |
| **US-Var** | `swu` | `-7.9 / 21.2 / 0.971 / 0.89 / 11.7` | `-7.9 / 21.2 / 0.971 / 0.89 / 11.7` | `-7.9 / 21.2 / 0.971 / 0.89 / 11.7` |
| **US-Var** | `lwu` | `-16.9 / 40.8 / 0.801 / 1.03 / 26.8` | `-1.9 / 40.8 / 0.826 / 1.33 / 25.8` | `-2.3 / 41.5 / 0.826 / 1.36 / 25.7` |

