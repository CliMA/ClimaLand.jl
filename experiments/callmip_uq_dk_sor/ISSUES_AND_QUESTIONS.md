# DK-Sor CalLMIP Phase 1a: Issues and Open Questions

Branch: `rb/callmip_phase1a`  
Site: DK-Sor (Danish beech forest, FLUXNET)  
Pipeline: calibration (EKI/UKI) → GP emulator (CES) → MCMC posterior → 50-member forward ensemble → CalLMIP NetCDF

---

## Pipeline overview

```
generate_observations.jl          ← build y_obs (NEE, Qle, Qh), noise_cov, save observations.jld2
         ↓
run_calibration.jl  (SLURM, 25 tasks, ~8h)
         ↓  TransformUnscented Kalman Inversion, 10 iterations, 25 sigma-points
         ↓  each iteration: forward_model() × 25 in parallel → observation_map() → EKP update
         ↓  saves iteration_NNN/eki_file.jld2 at each step
emulate_sample.jl  (SLURM, 1 node, ~4h)
         ↓  loads EKP, trains GP on iteration pairs, pCN-MH MCMC → posterior samples
         ↓  saves emulator, chain, posterior .jld2
run_posterior_ensemble.jl          ← 50 draws from posterior → forward runs → daily_diagnostics
analyze_posterior_ensemble.jl      ← seasonality plots
run_callmip_simulations.jl         ← CalLMIP-standard prior+posterior runs
write_callmip_netcdf.jl            ← produce submission NetCDF
```

**Model**: ClimaLand full land model — EnergyHydrology soil, CanopyModel (PModel photosynthesis, PlantHydraulics, TwoStream radiation), SnowModel, SoilCO₂ (DAMM biogeochemistry). DT = 900 s, 1-year spinup (2003), 10-year calibration window (2004–2013), 419 valid observation days (wind-filtered).

**Parameters calibrated** (12): 9 canopy/conductance (`moisture_stress_c`, `pmodel_cstar`, `pmodel_β`, `leaf_Cd`, `canopy_z_0m_coeff`, `canopy_z_0b_coeff`, `canopy_d_coeff`, `canopy_K_lw`, `canopy_emissivity`) + 3 DAMM soil-CO₂ (`soilCO2_pre_exponential_factor`, `michaelis_constant`, `O2_michaelis_constant`).

---

## Bugs fixed in this branch

### Bug 1 — RAI = 17.5 instead of 1.5 *(critical)*
**Files**: `model_interface.jl`, `callmip_model_interface.jl`, `run_prior_mean.jl`  
**Effect**: Root Area Index 10× too large → `Ra_root = Rd × μr × RAI` was ≈10× normal → ecosystem respiration ≫ GPP for all physically plausible parameters → NEE ≈ 0 always. EKI could not improve NEE across 10 iterations (RMSE oscillated 10³–10⁸, never converging). The calibrated "posterior" was just as bad as the prior for NEE.  
**Fix**: `RAI = FT(1.5)`.

### Bug 2 — LAI product: MODIS instead of Copernicus *(moderate)*
**File**: `model_interface.jl`  
**Effect**: Switched `"LAI"` → `"LAI_alternative"` in a previous commit (957be61a4). Copernicus peak LAI ≈ 2.97 in April vs. MODIS 1.36 — a 2× difference at the time of maximum canopy activity.  
**Fix**: Reverted to `"LAI"` (Copernicus Global Land Service), matching `run_dk_sor_default.jl`.

### Bug 3 — 14-param prior instead of 12 *(minor, was fixed)*
**Files**: `run_calibration.jl`, `emulate_sample.jl`, `run_callmip_simulations.jl`  
**Effect**: Two autotrophic-respiration (AR) parameters were included that are not calibrated by Alexis's reference run. Caused dimension mismatch when loading stale `eki_file.jld2`.  
**Fix**: Aligned to Alexis's exact 12-param prior.

### Bug 4 — DT = 450 s instead of 900 s *(minor)*
**File**: `run_calibration.jl`  
**Effect**: Calibration ran 2× slower than necessary. DT = 900 matches the reference simulation.  
**Fix**: `DT = Float64(900)`.

---

## Current blocker: MCMC in `emulate_sample.jl` fails to explore posterior

### Symptom
`optimize_stepsize` (Robbins–Monro scheduler) drives pCN step β → 10⁻⁹⁵ after 500 iterations with acceptance rate stuck at 0.0005 (≈0%). Chain of 100,000 steps produces zero variance — all samples are identical. Two SLURM jobs (62189075, 62318217) timed out at 4h without producing a valid posterior.

### Root cause hypothesis
The observation vector has **n = 1257 entries** (3 variables × 419 days), each with individual noise σ (0.341 gC/m²/d for NEE, 6.9 W/m² for Qle, 8.45 W/m² for Qh). The posterior log-likelihood is:

```
log p(y | θ) = -½ ∑ᵢ (yᵢ − Gᵢ(θ))² / σᵢ²
```

With 419 NEE terms each contributing ~(1/σ_NEE)² ≈ 8.6 to the precision, any O(1) change in unconstrained θ moves the predicted NEE by enough to cause rejection. The pCN proposal is `θ' = √(1−β²)·θ + β·ξ`, so even β = 10⁻³ shifts the GP prediction by enough to make log-likelihood drop drastically. The optimizer cannot find a regime where 23% acceptance is achievable because the posterior in unconstrained space is far sharper than any β > 0 can resolve.

Additionally, `noise_learn = false` in the GP means the emulator variance is not added to the likelihood — the MCMC sees only the observation noise, not GP uncertainty. For a nearly-converged EKI posterior (tight ensemble), the GP predictions extrapolate outside training data for most proposed θ, producing massive log-likelihood swings.

### What was tried
1. `init_stepsize` 0.1 → 1e-3, `max_iter` 40 → 500: still diverges to β = 10⁻⁹⁵
2. Fallback to `init_stepsize/10 = 1e-4`: same result

---

## Open questions for CliMA Calibrate / CES experts

### Q1. High-dimensional observations and MCMC
Our observation vector has 1257 entries (n_days × 3 variables). In CES examples, the observation dimension is typically O(10–100). **Is there a recommended approach for O(1000)-dimensional observation vectors?** Options we are aware of but uncertain about:
- Reduce to monthly/seasonal means before MCMC (sacrifice temporal information)
- Use only a random subset of days for MCMC (e.g. 50 of 419)
- Let CES itself reduce via its output PCA — but SVD truncation retained 787/1257 components at 95% variance, which is still large
- Use a different sufficient-statistic formulation for the likelihood

### Q2. `noise_learn` setting for pCN-MH
We use `GaussianProcess(gppackage; noise_learn = false)`. Should `noise_learn = true` be used so the GP's predictive variance enters the MCMC likelihood (making it less brittle)? Does the `obs_noise_cov` passed to `Emulator()` need to be modified when switching `noise_learn`?

### Q3. Is UKI final-ensemble spread a valid posterior substitute?
After 10 TransformUnscented iterations, the EKP maintains an approximate Gaussian posterior with mean `u_mean` and covariance `u_cov`. For publication-quality uncertainty quantification, is it acceptable to draw posterior samples directly from `N(u_mean_final, u_cov_final)` (in unconstrained space) without running MCMC through the GP emulator? Or does that throw away the non-Gaussian structure CES is trying to capture?

### Q4. GP emulator quality assessment
GP MSE: train = 51.4, test = 59.6 (in SVD-reduced output space). We have 150 training columns (6 iterations × 25 members) and 787 retained output components. **How should we interpret these MSE values?** Is there a rule of thumb for when they indicate the emulator is reliable enough to trust the MCMC likelihood? The train ≈ test MSE suggests low overfitting, but the absolute magnitude is hard to assess without knowing the scale of the output.

### Q5. Recommended workflow for UKI + CES
The EKP is a **TransformUnscented** (UKI) process, not standard EKI. `get_training_points(ekp, i)` returns `(u_{i-1}, G_i)` pairs. Are there known incompatibilities or caveats when using CES (which was designed for EKI) with UKI outputs? Specifically:
- The sigma-point ensemble at each iteration has different statistical structure than EKI particles (deterministic sigma grid, not random draws). Does this affect GP training?
- After 10 UKI iterations the ensemble collapses toward the optimum. Is training on the early (diffuse) iterations and testing on the late (tight) ones the right split?

### Q6. pCN vs. other samplers
Is `pCNMHSampling()` appropriate for 12-dimensional parameter space with a peaked posterior? CES also provides `RWMHSampling()` and `BarkerSampling()`. For a peaked, high-observation-dimensional problem, would `BarkerSampling()` (which adapts to gradient information) be significantly more efficient? What sampler does the CES team recommend for calibrated land-model posteriors?

### Q7. Number of training iterations
We use iterations 1–6 for training (7 for test) out of 10 completed UKI iterations. The concern is that the first few iterations have sigma points spread uniformly in prior space (far from the posterior), while the last few are tightly clustered. **What is the recommended train/test split, and does it matter whether early or late iterations are in the training set?**

### Q8. Observation noise covariance specification
The noise_cov is diagonal with constant σ per variable: σ_NEE = 0.341 gC/m²/d, σ_Qle = 6.9 W/m², σ_Qh = 8.45 W/m². These come from `generate_observations.jl`. **Are there known issues with using a diagonal (uncorrelated) noise covariance for temporally correlated flux observations in CES?** Should we inflate the noise or use a nugget term to improve MCMC mixing?

---

## What is working

| Stage | Status | Notes |
|---|---|---|
| Model runs (single forward sim) | ✅ | NEE goes negative in summer (carbon sink), physically sane |
| `generate_observations.jl` | ✅ | 419 valid days, 1257-entry y_obs |
| `run_calibration.jl` (EKI/UKI) | ✅ | 10 iter, 25 members, ~8h on HPC, converges for Qle/Qh |
| GP emulator training | ✅ | MSE train≈test, emulator cached |
| MCMC sampling | ❌ | pCN-MH acceptance → 0, β → 10⁻⁹⁵ |
| Downstream (posterior ensemble, CalLMIP NetCDF) | ⏳ | Blocked on MCMC |

## NEE calibration quality concern

Even with bugs fixed, `calibration_diagnostics.png` shows NEE RMSE oscillates wildly across iterations (10³–10⁸ range) instead of converging monotonically. Qle and Qh converge cleanly. This suggests the NEE signal is very hard for the model to fit within the 12-parameter space, possibly because:
- SOC initial conditions / DAMM parameters interact non-linearly with the canopy carbon flux in ways that produce degenerate sigma points for some parameter combinations  
- The `michaelis_constant` and `O2_michaelis_constant` priors may be numerically sensitive  
- The 419-day observation window may create a loss landscape with many local minima for NEE that confuse UKI sigma-point updates

---

## Files in this branch (new / modified)

### Modified source files
- `src/standalone/Soil/Biogeochemistry/Biogeochemistry.jl` — added SOC consumption tendency (`dY.soilco2.SOC -= Sm`) so carbon is conserved in DAMM
- `src/simulations/Simulations.jl` — relaxed `I <: SciMLBase.DEIntegrator` type constraint to `I` (fixes LandSimulation construction error with newer SciML)

### Calibration experiment (`experiments/calibrate_dk_sor/`)
- `model_interface.jl` — **LAI: Copernicus**, **RAI: 1.5**, docstring fix
- `run_calibration.jl` — DT=900, 12-param Alexis priors, correct output dir
- `run_prior_mean.jl` — RAI: 1.5
- `slurm_calibration.sh` — ntasks=25, removed `Pkg.update()` (saves ~3.5h)
- `generate_observations.jl` — (existing, unchanged)
- `plot_calibration_check.jl` — *(new)* diagnostic seasonality plot
- `run_alexis_test.jl` — *(new)* test run with Alexis's 12-param EKI posterior
- `slurm_alexis_test.sh` — *(new)* SLURM script for the above

### UQ experiment (`experiments/callmip_uq_dk_sor/`)
- `emulate_sample.jl` — 12-param priors, `find_latest_ekp_path`, N_train=8, init_stepsize=1e-3, max_iter=500
- `callmip_model_interface.jl` — RAI: 1.5
- `run_callmip_simulations.jl` — 12-param Alexis priors
- `run_posterior_ensemble.jl` — (existing)

### Submission scripts
- `submit_calibrate_then_ces.sh` — *(new)* submits calibration + CES as a `--dependency=afterok` chain
