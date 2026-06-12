# CES Phase 2 — Known Issues & Proposed Improvements
### Input for Oliver Dunbar (OLLie) · ClimaLand / CalLMIP · April 2026

---

## Background

These issues were exposed during the **CalLMIP Phase 1a DK-Sor** CES run
(10 EKI iterations, 33 members, GP emulator → 100 000-step pCN-MH MCMC → 50-member
posterior ensemble). The EKI calibration itself performed well (NEE RMSE −38%). The CES
emulator–sampler step produced a pathological posterior for several key parameters.
Everything is documented in full in `CALLMIP_Phase1a_Report.md` (Steps 3a–3f).

Relevant output files from this run:
```
experiments/callmip_uq_dk_sor/output_posterior_uq/
  ├── emulator_validation_scatter_its1to8.png   ← GP predicted vs actual
  ├── emulator_validation_rmse_its1to8.png      ← per-dimension RMSE
  ├── posterior_marginals_its1to8.png           ← 16-param MCMC posteriors
  └── posterior_its1to8.jld2                   ← 16 × 100 000 posterior samples
```

---

## Issue 1 — GP Coverage Gap (EKI Ensemble Collapse Before Emulation)

### What happened

EKI converged tightly to `pmodel_β ≈ 19` by iteration 9. All 260 GP training points
(33 members × 8 used iterations, `N_train = 8`) cluster near that optimum and are sparse
everywhere else. When the MCMC chain explores `pmodel_β ≈ 38` (near the prior mean), the
GP has no training data there and incorrectly extrapolates a comparably low-cost surface.
The chain settles near the prior mean → near-zero GPP year-round → flat NEE posterior median.

### Evidence

From `posterior_marginals_its1to8.png`:

| Parameter | EKI optimal | MCMC posterior mean | Verdict |
|---|---:|---:|---|
| `pmodel_β` | 19.1 | **38.4** | Bimodal — GP blind to prior region |
| `moisture_stress_c` | ~2.7 | ~0.4 | MCMC opposite direction from EKI |
| `canopy_z_0m_coeff` | 0.15 | 0.02 | 7× lower than EKI optimum |
| `soilCO2_pre_exponential_factor` | ~16 500 | ~75 000 | 4–5× higher than EKI optimum |
| `leaf_Cd` | ~0.05 | ~0.15 | 2–4× higher than EKI optimum |

This is **not** a MCMC sampling bug. The pCN-MH chain is mixing — it is just sampling
from a wrong (emulator-extrapolated) posterior.

### CES settings used

```julia
N_train       = 8        # EKI iterations used for GP training
chain_length  = 100_000  # pCN-MH steps
init_stepsize = 0.1      # pCN step size
ensemble_size = 33
n_iterations  = 10
```

### Proposed fixes

1. **More EKI iterations / larger ensemble before CES:** 15–20 iterations or 50–66 members
   would produce a denser, less collapsed training cloud around the true optimum.
2. **Iterative CES:** Run CES, check if MCMC mean ≈ EKI optimal. If not, seed a second
   EKI pass with MCMC samples and retrain the GP until convergence.
3. **Use all EKI iterations for training:** The current `N_train = 8` skips iteration 9
   (the most-converged). Including all iterations may not help coverage but is low-cost.
4. **GP prior regularisation:** Add a stronger length-scale prior so the GP does not
   extrapolate flat cost into unexplored regions — make the emulator "sceptical" outside
   the training hull.

---

## Issue 2 — High PCA Dimensions Poorly Emulated

### What happened

The GP is trained in PCA-decorrelated output space (1257 dimensions for daily NEE + Qle + Qh
over 1997–2012). Dimensions >800 show test RMSE up to 25 — 5–10× higher than training RMSE
— indicating the GP is overfitting noise-dominated principal components.

These high-index PCA dimensions encode fine-scale temporal variation: exactly where seasonal
NEE shape information (spring onset, summer peak, fall transition) is concentrated. The
emulator is least accurate precisely where it matters most.

### Evidence

`emulator_validation_rmse_its1to8.png`:
- Dims 0–600: test RMSE ≈ 1–5 (acceptable)
- Dims 600–800: increasing scatter
- Dims 800–1257: spikes to RMSE = 15–25

### Proposed fixes

1. **Truncate PCA** at a variance threshold (e.g. 95% or 99%) rather than a fixed dimension.
   For DK-Sor, 95% variance is likely captured by ~400–600 dims; emulating the rest adds
   noise without information.
2. **Weight PCA dimensions by explained variance** in the GP likelihood — discount high-index
   dims so they contribute less to parameter inference.
3. **Separate GP per flux:** Train three GPs (NEE, Qle, Qh) each with their own PCA space
   rather than one joint GP over all 3 × N_days output dimensions. This may reduce the
   effective output dimension significantly.

---

## Issue 3 — No MCMC Convergence Diagnostics

### What happened

`emulate_sample.jl` runs 100 000 pCN-MH steps and returns the full chain with no
convergence check. There is no R̂ (Gelman–Rubin) statistic, no effective sample size
(ESS), and no automatic trace plot. A badly mixed chain (e.g. stuck at prior mean due to
coverage gap Issue 1 above) is indistinguishable from a converged one in the current
outputs.

In this run the chain ran for 100 000 steps but `pmodel_β` stayed in the bimodal region
near 38 — without a trace plot this would have gone undetected.

### Proposed fixes

1. **Run multiple independent chains** (3–4) with random initialisations; compute
   R̂ per parameter; fail loudly if any R̂ > 1.1.
2. **Report ESS per parameter** (target > 1 000 effective samples per dimension).
3. **Auto-generate trace plots** for the top-5 highest-variance parameters after every
   CES run — these are the cheapest convergence check.
4. **Burn-in fraction:** Currently all 100 000 samples are used. Discard the first 20–50%
   as burn-in to reduce sensitivity to initialisation.

---

## Issue 4 — Structural ClimaLand Model Limitations (Not CES)

These biases survive calibration and are independent of CES quality. They belong to the
ClimaLand model development discussion, not to CES improvement.

| Bias | Observed (FLUXNET) | EKI-optimal model | Root cause |
|---|---|---|---|
| Winter respiration (Oct–Mar) | +1.5 to +2.4 gC m⁻² d⁻¹ | ≈ 0 | Cold-temperature heterotrophic resp. not represented |
| Energy flux constraint | — | Qle/Qh RMSE unchanged (+1% / −3%) | 16 biochem/canopy params decouple from energy partition |

*Raised for completeness — resolution requires model physics changes, not CES tuning.*

---

## Summary Table

| Issue | Severity | Effort to fix | Owner |
|---|---|---|---|
| 1. GP coverage gap | 🔴 High — corrupts posterior | Medium (more EKI iters or iterative CES) | OLLie / CES |
| 2. High PCA dims overfitted | 🟡 Medium — degrades emulator accuracy | Low (truncate PCA threshold) | OLLie / CES |
| 3. No MCMC diagnostics | 🟡 Medium — silent failures undetectable | Low (add R̂ / ESS reporting) | OLLie / CES |
| 4. Winter resp. / energy flux | 🟠 Science gap | High (new model physics) | ClimaLand team |
