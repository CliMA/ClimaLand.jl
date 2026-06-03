# DK-Sor Uncertainty Quantification Pipeline

Full CalibrateEmulateSample (CES) UQ pipeline for the DK-Sor single-site
calibration, supporting [CalLMIP Phase 1a](https://github.com/callmip-org/Phase1).

## Pipeline overview

```
run_calibration.jl               (already done — produces EKP in output/)
        │
        ▼
emulate_sample.jl                (train GP emulator, run MCMC posterior)
        │
        ▼
run_posterior_ensemble.jl        (50 forward runs sampled from posterior)
        │
        ▼
analyze_posterior_ensemble.jl    (uncertainty bands on NEE/Qle/Qh)
        │
        ▼
run_callmip_simulations.jl       (prior + posterior ClimaLand runs with full
        │                          CalLMIP diagnostic output)
        │
        ▼
write_callmip_netcdf.jl          (produces CalLMIP-compliant NetCDF files)
```

## Required inputs

| File | Source |
|------|--------|
| `output/iteration_NNN/eki_file.jld2` | `run_calibration.jl` |
| `observations.jld2` | `generate_observations.jl` |
| `DK_Sor/DK-Sor_daily_aggregated_*_Flux.nc` | FLUXNET data |
| `DK_Sor/DK-Sor_1997-2014_FLUXNET2015_Met.nc` | FLUXNET data |

## Dependencies

The `emulate_sample.jl` and `analyze_posterior_ensemble.jl` scripts require
`CalibrateEmulateSample.jl`. Add it to the project:

```julia
# From experiments/callmip_uq_dk_sor/
julia --project=. -e 'using Pkg; Pkg.add("CalibrateEmulateSample"); Pkg.instantiate()'
```

Or add it to `.buildkite/Project.toml` as described in
[issue #1602](https://github.com/CliMA/ClimaLand.jl/issues/1602).

`run_posterior_ensemble.jl` uses the `.buildkite` project (same as calibration).

## Running on Caltech HPC

```bash
# Step 1 (already done): calibration
sbatch slurm_calibration.sh

# Step 2: emulate + sample (single node, ~4 h)
sbatch slurm_emulate_sample.sh

# Step 3: forward ensemble (~12 h, 50 members)
sbatch --dependency=afterok:<step2_jobid> slurm_posterior_ensemble.sh

# Step 4: analysis (local or interactive)
julia --project=experiments/callmip_uq_dk_sor \
      experiments/callmip_uq_dk_sor/analyze_posterior_ensemble.jl

# Step 5: CalLMIP simulations (prior + posterior, ~24 h)
sbatch --dependency=afterok:<step4_jobid> slurm_callmip_simulations.sh

# Step 6: write CalLMIP NetCDF (local or interactive, seconds)
julia --project=experiments/callmip_uq_dk_sor \
      experiments/callmip_uq_dk_sor/write_callmip_netcdf.jl
```

## Output files

```
experiments/callmip_uq_dk_sor/
├── output_posterior_uq/
│   ├── emulator_its1to7.jld2               (trained GP emulator)
│   ├── mcmc_and_chain_its1to7.jld2         (MCMC chain)
│   └── posterior_its1to7.jld2              (posterior + EKI optimum)
│       ├── constrained_posterior            (n_params × n_mcmc samples)
│       └── constrained_ekp_optimal         (EKI best-estimate)
├── output_posterior_ensemble/
│   ├── iteration_000/member_NNN/
│   │   ├── parameters.toml                 (sampled parameter values)
│   │   └── daily_diagnostics.jld2          (NEE, Qle, Qh timeseries)
│   ├── posterior_parameter_samples.jld2
│   └── posterior_ensemble_diagnostics.jld2
├── output_posterior_analysis/
│   ├── posterior_uncertainty_nee.png
│   ├── posterior_uncertainty_qle.png
│   ├── posterior_uncertainty_qh.png
│   ├── posterior_uncertainty_all_fluxes.png
│   └── posterior_uncertainty_bands.jld2    (percentile bands for CalLMIP)
├── output_callmip_sims/
│   └── iteration_000/
│       ├── member_001/callmip_diagnostics.jld2  (prior: all 13 ALMA vars)
│       └── member_002/callmip_diagnostics.jld2  (posterior: all 13 ALMA vars)
└── callmip_output/                              ← submit these 4 files
    ├── ClimaLand.CalLMIP1.0_Expt1_DK-Sor_Cal_Prior.nc
    ├── ClimaLand.CalLMIP1.0_Expt1_DK-Sor_Cal_Posterior.nc
    ├── ClimaLand.CalLMIP1.0_Expt1_DK-Sor_Val_Prior.nc
    └── ClimaLand.CalLMIP1.0_Expt1_DK-Sor_Val_Posterior.nc
```

## Calibrated parameters (14 total)

| # | Name | Physical meaning |
|---|------|-----------------|
| 1 | `moisture_stress_c` | Soil moisture stress threshold |
| 2 | `pmodel_cstar` | P-model cost function parameter |
| 3 | `pmodel_β` | P-model stomatal slope |
| 4 | `leaf_Cd` | Leaf drag coefficient |
| 5 | `canopy_z_0m_coeff` | Momentum roughness length coefficient |
| 6 | `canopy_z_0b_coeff` | Scalar roughness length coefficient |
| 7 | `canopy_d_coeff` | Displacement height coefficient |
| 8 | `canopy_K_lw` | Longwave extinction coefficient |
| 9 | `canopy_emissivity` | Canopy emissivity |
| 10 | `root_leaf_nitrogen_ratio` | Root-to-leaf nitrogen ratio |
| 11 | `stem_leaf_nitrogen_ratio` | Stem-to-leaf nitrogen ratio |
| 12 | `soilCO2_pre_exponential_factor` | DAMM pre-exponential factor |
| 13 | `michaelis_constant` | DAMM Michaelis–Menten constant (O₂) |
| 14 | `O2_michaelis_constant` | DAMM Michaelis–Menten constant (O₂) |

## Key settings to tune

`emulate_sample.jl`:
- `N_train` (default 7): number of EKI iterations used to train the GP.
  Increase if emulator test error is large; decrease if training data is limited.
- `retain_var_out` (default 0.95): adjusts output PCA truncation.
  Higher = more GP models = slower but more accurate emulator.
- `chain_length` (default 100 000): MCMC chain length.

`run_posterior_ensemble.jl`:
- `N_SAMPLES` (default 50): number of posterior draws to run through the land model.
  Increase for better-resolved uncertainty bands (update `--ntasks` in SLURM script).

## CalLMIP Phase 1 output variables

The 13 required ALMA variables and how they are derived from ClimaLand diagnostics:

| ALMA name | ClimaLand short | ClimaLand units | ALMA units | Conversion |
|-----------|----------------|----------------|-----------|------------|
| NEE | `nee` | mol CO₂ m⁻² s⁻¹ | kg C m⁻² s⁻¹ | × 12×10⁻³ |
| Qle | `lhf` | W m⁻² | W m⁻² | ×1 |
| Qh | `shf` | W m⁻² | W m⁻² | ×1 |
| GPP | `gpp` | mol CO₂ m⁻² s⁻¹ | kg C m⁻² s⁻¹ | × 12×10⁻³ |
| Reco | `er` | mol CO₂ m⁻² s⁻¹ | kg C m⁻² s⁻¹ | × 12×10⁻³ |
| TVeg | `trans` | kg m⁻² s⁻¹ | kg m⁻² s⁻¹ | ×1 |
| ESoil | `soillhf` | W m⁻² | kg m⁻² s⁻¹ | ÷ 2.5×10⁶ |
| Qg | `soilrn − soillhf − soilshf` | W m⁻² | W m⁻² | ×1 |
| AvgSurfT | `ct` (canopy T); fallback: top `tsoil` | K | K | ×1 |
| SoilMoist | integrate `swc` × ρ_w | m³ m⁻³ per layer | kg m⁻² | Σ(θ·dz·1000) |
| LAI | `lai` | m² m⁻² | m² m⁻² | ×1 |
| TotAbovBioMass | `cveg` | kg C m⁻² | kg C m⁻² | ×1 |
| TotSoilCarb | integrate `soc` | kg C m⁻³ per layer | kg C m⁻² | Σ(SOC·dz) |

---

## Known issues & fixes (2026-03-23)

### 1. `JLD2` fails to precompile — missing `ScopedValues`

**Symptom**: `emulate_sample.jl` (job `dk_sor_emulate_sample`) fails immediately with:

```
ERROR: LoadError: ArgumentError: Package ScopedValues [7e506255-...] is required
but does not seem to be installed: Run `Pkg.instantiate()` to install all recorded
dependencies.
```

**Root cause**: The `ScopedValues` v1.6.0 package directory in the Julia depot was
created as an empty stub (`~/miniconda3/share/julia/packages/ScopedValues/JYd8U/`)
with no source files inside, causing `Pkg.instantiate()` to silently skip it and
`JLD2` (which depends on it) to fail precompilation.

Note: `ScopedValues` is not a stdlib in Julia 1.12.5 — it lives in `Base.ScopedValues`
internally, but the registered package version (`1.6.0`) must still be installed in the
depot for older packages that `using ScopedValues` directly.

**Fix**:

```bash
# 1. Load the correct Julia version (must be 1.12.x matching the Manifest)
module use /groups/esm/modules
module load climacommon

# 2. Remove the empty broken package directory
rm -rf ~/miniconda3/share/julia/packages/ScopedValues/JYd8U

# 3. Re-instantiate to download ScopedValues 1.6.0 properly
cd /path/to/ClimaLand.jl
julia --project=experiments/callmip_uq_dk_sor -e 'using Pkg; Pkg.instantiate()'
```

**Safety note**: Only `JYd8U` (v1.6.0) was removed. The `z27HA` directory
(`ScopedValues` v1.5.0, used by Julia 1.10 environments) is unaffected.

---

### 2. `pkgdir(ClimaLand)` resolves to depot package, not local checkout

**Symptom**: After fixing issue #1, the job fails with:

```
ERROR: LoadError: IOError: mkdir(".../miniconda3/share/julia/packages/ClimaLand/K0wUS/
experiments/callmip_uq_dk_sor/output_posterior_uq"; ...): no such file or directory
```

**Root cause**: `emulate_sample.jl` used `pkgdir(ClimaLand)` to locate the repo root.
When the script runs with `--project=experiments/callmip_uq_dk_sor`, Julia resolves
`ClimaLand` against that project's Manifest, which lists it as a registered package in
the depot — not the local dev checkout. So `pkgdir(ClimaLand)` returns the depot path
instead of the working tree.

**Fix** (applied to `emulate_sample.jl`):

```julia
# Before:
const climaland_dir = pkgdir(ClimaLand)

# After:
const climaland_dir = abspath(joinpath(@__DIR__, "..", ".."))
```

`@__DIR__` always refers to the directory containing the source file at load time,
making the path independent of how `ClimaLand` is resolved in the project environment.

NEE sign convention: NEE > 0 = net carbon efflux to atmosphere (consistent with
FLUXNET and ALMA standards).  ER and GPP are both positive-definite.

The posterior files (`_Posterior.nc`) additionally contain 5th / 50th / 95th
percentile variables (`NEE_p05`, `NEE_p50`, `NEE_p95`, `Qle_p05` … `Qh_p95`)
derived from the 50-member parameter ensemble.  These require that
`analyze_posterior_ensemble.jl` has been run before `write_callmip_netcdf.jl`;
if not, the posterior point-estimate variables are still complete.

**Temporal split** (set by `CAL_END_YEAR` in `write_callmip_netcdf.jl`, default 2012):
- Cal files: 2004–2012  
- Val files: 2013 (temporal validation year)

**File naming** (CalLMIP Phase 1 Protocol v1.2, Section 7):
`<model>.<version>_Expt<no.>_<site>_<Cal/Val>_<Prior/Posterior>.nc`

### Issue 3 — GP hyperparameter optimisation diverges to Inf (2026-03-24)

---

### Issue 4 — `optimize_stepsize` fails to converge; MCMC step-size implications (2026-03-24)

**⚠️ Needs discussion with Ollie (CES developer)**

**Symptom:** Job fails after MCMC sampling completes (~17 min, chain at 100%) with:

```
ERROR: LoadError: "Failed to choose suitable stepsize in 20 iterations."
```

The 100k-step MCMC chain itself ran successfully; only the *step-size search*
preceding it timed out.

**Root cause:** `optimize_stepsize` runs an iterative binary search targeting a
Metropolis-Hastings acceptance rate of 15–35%. When the emulator likelihood
surface is poorly conditioned (e.g. weak identifiability in some directions after
SVD decorrelation), the acceptance rate oscillates and never stabilises within
`max_iter` evaluations. The function throws rather than returning a best guess.

**Temporary fix** (applied to `emulate_sample.jl`):

```julia
# Before:
new_step = optimize_stepsize(rng, mcmc; init_stepsize, N = 2000, discard_initial = 0)

# After:
new_step = try
    optimize_stepsize(rng, mcmc; init_stepsize, N = 2000, discard_initial = 0, max_iter = 40)
catch e
    @warn "optimize_stepsize did not converge ($e). Falling back to init_stepsize=$init_stepsize."
    init_stepsize
end
```

`max_iter` is doubled to 40 first. If still unconverged, the fallback
`init_stepsize = 0.1` is used. The chain will run but may mix poorly if 0.1 is
far from the optimal step size.

**Open question for Ollie:**
- Is the pCNMH step-size search expected to struggle when only 12 parameters are
  active (reduced from 14) with the current emulator?
- Should we pass a `target_acc` tuned for pCN (which can tolerate higher
  acceptance rates than vanilla RW-MH)?
- Is there a recommended fallback or a better initialisation of `init_stepsize`
  for this problem size?
- Would switching to `RWMHSampling` or `BarkerSampling` be more robust here?

---

### Issue 3 — GP hyperparameter optimisation diverges to Inf (2026-03-24)

**Symptom:** Job fails after ~6 min with:

```
AssertionError: isfinite(phi_c) && isfinite(dphi_c)
```

Originating in `LineSearches/src/hagerzhang.jl` inside `GaussianProcesses.optimize!`.
The `.out` log shows one or more output dimensions printing `Params: [Inf, Inf, ...]`
(e.g. dim 635) before the crash — the LBFGS optimiser escaped to ±∞ for
nearly-flat/unidentifiable output dimensions, causing the next dimension's line
search to receive non-finite inputs and assert.

**Root cause:** No bounds were placed on the log-length-scale parameters, so the
LBFGS optimizer can drive them to ±∞ for dimensions where the signal-to-noise
ratio is very low after SVD decorrelation.

**Fix** (applied to `emulate_sample.jl`):

Replace the bare `optimize_hyperparameters!(emulator)` call with a bounded
version:

```julia
# Before:
optimize_hyperparameters!(emulator)

# After:
# Bound log-length-scales to [-10, 10] to prevent LBFGS divergence.
n_gp_inputs = size(get_inputs(train_pairs), 1)
kb = [fill(-10.0, n_gp_inputs + 1), fill(10.0, n_gp_inputs + 1)]
optimize_hyperparameters!(emulator; kernbounds = kb)
```

The `+1` accounts for the log-variance parameter that SEArd appends after the
`n_gp_inputs` log-length-scale parameters. The bounds `[-10, 10]` correspond to
length scales of `e^{-10} ≈ 4.5×10^{-5}` to `e^{10} ≈ 22000` in normalised
input space, which is wide enough to be non-restrictive for well-identified
dimensions while preventing numerical blow-up for flat ones.

## References

- Rennon, O. et al., CalibrateEmulateSample.jl:
  https://github.com/CliMA/CalibrateEmulateSample.jl
- ClimaLand calibration experiment (DK-Sor):
  https://github.com/CliMA/ClimaLand.jl/pull/1618
- ClimaLand UQ issue:
  https://github.com/CliMA/ClimaLand.jl/issues/1602
- CalLMIP Phase 1:
  https://github.com/callmip-org/Phase1
