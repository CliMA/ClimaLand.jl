# ClimaLand — CalLMIP Phase 1a Model and Calibration Description
*Section for the CalLMIP Phase 1 protocol paper (Section 12)*
*Site: DK-Sor · April 2026*

---

## ClimaLand

### i. Model description

ClimaLand (v1.6.0, branch `rb/callmip_phase1a`) is a modular, process-based land surface
model developed by the CliMA project. It couples a two-stream radiative transfer canopy
model (TSRT; [Braghiere et al., 2021](https://doi.org/10.1029/2020MS002354)), a prognostic
multi-layer soil hydrology and heat model based on Richards' equation, a two-leaf
photosynthesis model (P-model; [Stocker et al., 2020](https://doi.org/10.5194/gmd-13-1545-2020)),
autotrophic respiration parameterized via leaf dark respiration scaled by nitrogen ratios,
heterotrophic soil CO₂ efflux via the DAMM (Dual Arrhenius and Michaelis-Menten) model
[Davidson et al., 2012], a prognostic snow model, and the Medlyn stomatal conductance model.
The model is implemented in Julia using ClimaCore.jl (v0.14.46) for spatial discretization.
Parameter calibration uses ClimaCalibrate.jl (v0.1.4) and EnsembleKalmanProcesses.jl (v2.5.0).

**This is not a CMIP6 or CMIP7-FT configuration.** ClimaLand is a new model.

### ii. Process switches

N cycling, fire, dynamic vegetation, demography, and DOC/POC are **not** included in this
configuration. The canopy structure is prescribed (fixed single PFT, prescribed LAI — see
point xv). Soil organic carbon (SOC) is initialized via a prescribed exponential profile
but has no dynamic decomposition cascade (no litter pools, no turnover); SOC is therefore
fixed throughout the simulation. The active prognostic components are:
canopy, soil hydro-thermodynamics, snow, and soil CO₂ (DAMM heterotrophic respiration).

### iii. N cycling

Not included. N deposition and N fertilisation data were not used.

### iv. Use of site metadata

**PFT:** DK-Sor is classified as a temperate deciduous broadleaf forest (DBF).
ClimaLand uses a single-PFT, single-column configuration. PFT-specific parameter defaults
(stomatal slope, P-model electron transport scaling, canopy optical properties) are drawn
from ClimaLand's DBF parameter set, which is subsequently modified by the 16 calibrated
parameters.

**Soil depth and texture:** The soil column extends from the surface to −9 m, discretized
into 20 layers with a stretched grid (finer spacing near the surface). Soil texture is
specified via fixed volumetric fractions: total porosity ν = 0.45, residual moisture
θ_r = 0.07, organic matter fraction ν_ss,om = 0.03, quartz fraction ν_ss,quartz = 0.47,
gravel fraction ν_ss,gravel = 0.12. These values were set based on site metadata and are
held fixed (not calibrated).

**% cover:** A single PFT covering 100% of the column area is assumed, consistent with
the monoculture beech (*Fagus sylvatica*) stand at DK-Sor.

### v. Parameter estimation method

Calibration uses **Ensemble Kalman Inversion (EKI)** with the **Transform Unscented Kalman
Inversion** (TUKI) variant [Huang et al., 2022] as implemented in ClimaCalibrate.jl v0.1.4.
TUKI constructs a deterministic sigma-point ensemble from the prior, propagates it through
the forward model, and applies an EKI update with prior regularisation
(`impose_prior = true`). The ensemble size is **N = 33** (2×16 + 1 sigma points for 16
parameters). **10 EKI iterations** were run; each iteration is one forward model evaluation
per ensemble member followed by an EKI update. A random seed (MersenneTwister(1234)) was
fixed for reproducibility.

**Calibration model time step:** Δt = 900 s (for the EKI forward model evaluations).
The CalLMIP prior/posterior simulations were subsequently re-run at Δt = 450 s.

**Observation vector** **y** has length 3n, constructed by concatenating daily-mean
NEE (gC m⁻² d⁻¹), Qle (W m⁻²), and Qh (W m⁻²) over n valid calibration days
(2004–2013, see day-selection criteria below).

**Day selection:** Days with daily mean wind speed ≥ 5 m s⁻¹ were excluded
(high-wind conditions reduce representativeness of the single-column model). Days with
missing values in any of the three target variables were also excluded.

**Observation error covariance Γ** is diagonal. The variance for each variable is the
time-mean squared measurement uncertainty from the FLUXNET2015 file
(σ²_NEE = mean(σ²_uc,NEE,t), and similarly for Qle and Qh). The same scalar variance
is applied to all valid days.

**Cost function** (implicit in the EKI update) minimises:
$$\mathcal{J}(\boldsymbol{\theta}) = \|\mathbf{y} - \mathcal{G}(\boldsymbol{\theta})\|^2_{\mathbf{\Gamma}^{-1}} + \|\boldsymbol{\theta} - \boldsymbol{\mu}_0\|^2_{\mathbf{C}_0^{-1}}$$
where G(θ) is the ClimaLand forward model output, μ₀ and C₀ are the prior mean and
covariance. All three data streams are assimilated **simultaneously** at each iteration.

**Posterior distribution** is obtained via Calibrate–Emulate–Sample (CES
[Cleary et al., 2021]):

1. A **Gaussian Process** (GP) emulator with a squared-exponential ARD kernel
   (GaussianProcesses.jl) is trained on 264 model evaluations from EKI iterations 1–8
   (8 iterations × 33 members; iteration 0 withheld to avoid prior-draw outliers).
   Iteration 9 and a randomly withheld 20% of training points are used as a test set.
   Input space is reduced via PCA retaining 99% of variance; output space is reduced
   via PCA retaining 95% of variance; a separate GP is fitted per retained output
   dimension. Extreme outlier columns (blow-up members) are pruned before training.
2. **Preconditioned Crank–Nicolson Metropolis–Hastings** (pCN-MH) MCMC [Cotter et al.,
   2013] with 100 000 steps uses the GP emulator as a cheap surrogate for the
   likelihood. Step size is optimised automatically (2 000 trial steps) to target
   ~23% acceptance.
3. After thinning and removal of 7 physically implausible members (daily NEE >
   50 gC m⁻² d⁻¹), a **43-member** posterior ensemble is obtained.

The **EKI-optimal parameter set** (best ensemble member at iteration 9, selected by
minimum cost) is used as the point-estimate posterior for the CalLMIP submission.

> **Note on CalLMIP Cal/Val temporal split:** The EKI calibration training window was
> 2004–2013 (inclusive). The CalLMIP submission splits the model output at year 2012,
> designating 2003–2012 as the "Cal" period and 2013 as the "Val" (temporal validation)
> period. Consequently, the CalLMIP validation year (2013) was included in the EKI
> training data and does not constitute a true hold-out for the EKI calibration itself.
> This is a limitation of the current Phase 1a setup and will be corrected for Phase 1b.

### vi. Parameter selection

Parameters were selected by expert elicitation based on (a) known sensitivity of
NEE/Qle/Qh to the parameter, and (b) known prior biases at this site (e.g. the default
P-model β = 51 produced summer GPP approximately twice the FLUXNET observed value;
the default canopy displacement height coefficient d_coeff = 0.007 gives d ≈ 0.17 m
for a 25 m canopy — an order of magnitude too low for a tall forest). No formal
sensitivity analysis was performed. **16 parameters** were selected out of O(100) total
ClimaLand parameters, covering: canopy photosynthesis and conductance (9), DAMM
heterotrophic respiration kinetics (3: pre-exponential factor, two Michaelis constants),
autotrophic-respiration-related parameters (3: Arrhenius Eₐ, root N ratio, stem N ratio),
and canopy heat capacity (1).

> **Note on parameter grouping:** `soilCO2_activation_energy` (DAMM Arrhenius Eₐ) is
> listed under heterotrophic respiration in §vii but was added specifically to address
> the winter CO₂ source deficit alongside the two autotrophic N ratio parameters. The
> user-supplied description groups it with autotrophic respiration (yielding 3 DAMM + 3
> autotrophic-related); both groupings are valid. The total of 16 is unambiguous.

### vii. Parameter table

| Short name | Long name | Process | Prior mean | Prior σ | Bounds |
|---|---|---|---:|---:|---|
| `moisture_stress_c` | Moisture stress slope | Stomatal conductance | 0.50 | 0.30 | [0.01, 5.0] |
| `pmodel_cstar` | P-model c* scaling | Photosynthesis | 0.43 | 0.15 | [0.05, 2.0] |
| `pmodel_β` | P-model water cost ratio | Photosynthesis | 20.0 | 8.0 | [2.0, 80.0] |
| `leaf_Cd` | Leaf drag coefficient | Canopy aerodynamics | 0.10 | 0.05 | [0.005, 1.0] |
| `canopy_z_0m_coeff` | Roughness length coeff. (z₀ₘ/h) | Canopy aerodynamics | 0.10 | 0.04 | [0.01, 0.25] |
| `canopy_z_0b_coeff` | Scalar roughness coefficient | Canopy aerodynamics | 0.001 | 0.0005 | [1×10⁻⁵, 0.01] |
| `canopy_d_coeff` | Displacement height coeff. (d/h) | Canopy aerodynamics | 0.65 | 0.12 | [0.30, 0.92] |
| `canopy_K_lw` | LW extinction coefficient | Canopy radiative transfer | 0.85 | 0.25 | [0.1, 2.0] |
| `canopy_emissivity` | Canopy thermal emissivity | Canopy radiative transfer | 0.97 | 0.02 | [0.9, 1.0] |
| `soilCO2_pre_exponential_factor` | DAMM pre-exponential factor α | Heterotrophic respiration | 25 000 | 10 000 | [1 000, 200 000] |
| `michaelis_constant` | DAMM substrate Michaelis constant Kₘ | Heterotrophic respiration | 0.010 | 0.005 | [1×10⁻⁴, 0.1] |
| `O2_michaelis_constant` | DAMM O₂ Michaelis constant | Heterotrophic respiration | 0.010 | 0.005 | [1×10⁻⁴, 0.1] |
| `soilCO2_activation_energy` | DAMM Arrhenius activation energy Eₐ (J mol⁻¹) | Heterotrophic respiration | 61 000 | 15 000 | [30 000, 100 000] |
| `root_leaf_nitrogen_ratio` | Root-to-leaf maintenance resp. ratio μᵣ | Autotrophic respiration | 1.0 | 0.5 | [0.1, 5.0] |
| `stem_leaf_nitrogen_ratio` | Stem-to-leaf maintenance resp. ratio μₛ | Autotrophic respiration | 0.10 | 0.05 | [0.01, 0.5] |
| `ac_canopy` | Canopy heat capacity (J m⁻² K⁻¹) | Canopy energy balance | 2 500 | 1 500 | [500, 10 000] |

All priors are constrained Gaussians mapped to physical bounds via a logit-normal transform
as implemented in EnsembleKalmanProcesses.jl (`constrained_gaussian`).

### viii. Spin-up method

A one-year meteorological spin-up is applied: ClimaLand is initialised on 2003-01-01 with
prescribed initial conditions (see §ix) and forced with FLUXNET gap-filled half-hourly
meteorology through 2003-12-31. The 2003 output is discarded; the calibration observation
window begins 2004-01-01. The same spin-up procedure is applied for every EKI ensemble
member and every CES posterior member. No multi-year or equilibrium spin-up is performed;
ClimaLand does not currently support accelerated steady-state methods (e.g., matrix
spin-up). The 1-year spin-up is a pragmatic choice based on the timescale for soil
moisture to approach a forced quasi-equilibrium under realistic meteorology.

> **⚠️ Deviation from CalLMIP Protocol §8.** The protocol requires cycling the forcing
> until carbon stock equilibrium is reached. ClimaLand Phase 1a uses only a 1-year
> spin-up. Soil carbon is prescribed as a fixed exponential profile
> (SOC(z) = 0.5 + 14.5 × exp(z / τ_soc) kg C m⁻³, τ_soc = 1/ln(30) m)
> and does not evolve dynamically. This is a known limitation to be addressed in Phase 1b.

### ix. Initial conditions

All initial conditions are **prescribed directly** at the start of the simulation
(2003-01-01) and are **not** included in the parameter calibration:

- **Soil moisture:** 95% saturation — ϑ_l = θ_r + 0.95 × (ν − θ_r)
- **Soil temperature:** interpolated from the atmospheric forcing driver at t = 0
- **Snow:** zero (S = S_l = U = 0)
- **Canopy temperature:** from atmospheric forcing at t = 0
- **Canopy hydraulics:** initialized at saturation (ϑ_l = ν for all stem/leaf nodes)
- **Soil CO₂:** 412 ppm; O₂ fraction 0.21; SOC from prescribed exponential profile

Initial conditions are identical for the prior and posterior simulations.

### x. Multi-datastream calibration approach

NEE, Qle, and Qh are assimilated **simultaneously** in a single observation vector
**y** = [NEE₁, …, NEEₙ, Qle₁, …, Qleₙ, Qh₁, …, Qhₙ]ᵀ at each EKI iteration.
The cost function weights each variable by its respective time-mean squared measurement
uncertainty. No stepwise or sequential treatment is applied.

### xi. Operational parameter set

ClimaLand does not yet have an operational multi-site calibration system. For Phase 1a,
a single-site, single-PFT parameter set is produced for DK-Sor (Temperate DBF). The
EKI-optimal parameters from iteration 9 constitute the posterior point estimate submitted
to CalLMIP. Extension to a globally operational PFT-specific parameter set is planned for
Phase 1b using multi-site calibration.

### xii. PFT coverage

A single PFT (Temperate DBF) is calibrated for DK-Sor. The site is a monoculture beech
(*Fagus sylvatica*) forest; 100% DBF cover is assumed for the single-column domain.

### xiii. Other uncertainty sources

Only **measurement uncertainty** (FLUXNET2015 daily σ_uc for NEE, Qle, Qh) is propagated
into the observation error covariance Γ. Model structural error and parameter–forcing
interaction uncertainty are not explicitly represented. The diagonal Γ assumes temporal
independence of daily observations. Initial condition, boundary, and driver uncertainties
are not propagated. The CES posterior ensemble provides an estimate of parametric
uncertainty only; the ensemble collapses toward the prior mean due to a known GP emulator
coverage failure described in `CES_Phase2_Issues_for_OLLie.md`.

### xiv. Optional experiments

Only the required **Experiment 1** (simultaneous NEE + Qle + Qh calibration) was
performed for Phase 1a.

### xv. LAI treatment

LAI is **prescribed** from the remotely-sensed time series provided in the PLUMBER2
forcing file, specifically the **Copernicus Global Land Service** LAI product
(`LAI_alternative` variable in the NetCDF forcing file). The time series is interpolated
to the model time step using a `TimeVaryingInput` piecewise-linear interpolant
(Δt = 900 s for EKI calibration runs; Δt = 450 s for the CalLMIP prior/posterior
simulations). Prognostic LAI is not used in this configuration.

---

*Last updated: April 2026*
