# CalLMIP Phase 1a DK-Sor — Results Summary

**Status:** To be filled after calibration run completes.

---

## Calibration performance (1997-2012)

| Variable | Prior RMSE | Prior bias | Posterior RMSE | Posterior bias | Units |
|---|---|---|---|---|---|
| NEE | TBD | TBD | TBD | TBD | gC/m²/d |
| Qle | TBD | TBD | TBD | TBD | W/m² |
| Qh  | TBD | TBD | TBD | TBD | W/m² |

## Validation performance (2013)

| Variable | Posterior RMSE | Posterior bias | Units |
|---|---|---|---|
| NEE | TBD | TBD | gC/m²/d |
| Qle | TBD | TBD | W/m² |
| Qh  | TBD | TBD | W/m² |

---

## Parameter shifts (prior mean → posterior mean)

| Parameter | Prior mean | Posterior mean | Posterior std |
|---|---|---|---|
| pmodel_α | TBD | TBD | TBD |
| pmodel_β_c3 | TBD | TBD | TBD |
| soilCO2_reference_rate | 6.0e-8 | TBD | TBD |
| O2_michaelis_constant | TBD | TBD | TBD |
| … | | | |

---

## Emulator validation

- Training set size: N_ens × N_outputs (33 members × (16×365×3) = see `emulate_sample.jl`)
- Holdout correlation r: TBD (target > 0.6)
- GP kernel: Matern-5/2

---

## Posterior uncertainty coverage

- NEE 90% coverage (obs in [p5, p95]): TBD%  (target ≥ 70%)
- Qle 90% coverage: TBD%

---

## Known issues

### Winter NEE deficit (DAMM)
The DAMM model overestimates winter soil respiration at DK-Sor.
DJF mean NEE — prior: TBD, posterior: TBD, obs: TBD gC/m²/d.
This is a known ClimaLand.jl v1.8.x limitation; the lognormal prior on
`soilCO2_reference_rate` partially mitigates the issue compared to the
old bounded-Gaussian prior.

### Not-available diagnostics
- `evspsblsoi` (bare soil evaporation): Not diagnosed in ClimaLand.jl main → NaN in NetCDF
- `cLiveBioAbove` (aboveground biomass): Not diagnosed → NaN in NetCDF
