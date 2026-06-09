# CalLMIP Phase 1a — DK-Sor (Sorø, Denmark)

**Branch:** `rb/callmip_phase1a_v2`
**ClimaLand.jl commit:** `afa000ceb` (main)
**Protocol:** CalLMIP Phase 1a, Scenario 1 (default)
**Site:** DK-Sor — FLUXNET2015 DK-Sor, 55.4859°N 11.6446°E, broadleaf forest (DBF)

---

## Pipeline overview

```
Step 0a  generate_observations.jl       Build 16 yearly obs windows (1997-2012)
Step 0b  run_prior_check.jl             Sanity-check prior mean simulation (GATE)
Step 1   run_calibration.jl             EKI calibration (TransformUnscented, 10 iters)
Step 2   emulate_sample.jl              GP emulator → pCN-MCMC posterior (CES v1)
Step 3   run_callmip_simulations.jl     Prior + posterior 1997-2014 simulations
Step 4   run_posterior_ensemble.jl      50-member posterior uncertainty ensemble
Step 5   write_callmip_netcdf.jl        Write CalLMIP-compliant NetCDF output
Step 6   (post) analyze_posterior_ensemble.jl  Coverage statistics
         (post) evaluate_calibration.jl         RMSE/bias table
         (post) plot_calibration_check.jl        Convergence + seasonal plots
         (post) plot_posterior_distributions.jl  Prior vs posterior histograms
```

---

## Execution (HPC / Slurm)

> All scripts must be submitted from the ClimaLand.jl repo root.
> Never `cd` into the experiment dir before `sbatch` — `$SLURM_SUBMIT_DIR` won't point to the repo root.

```bash
cd /central/scratch/renatob/ClimaLand.jl

# 0a. Generate observations (quick, run interactively or as sbatch)
julia --project=experiments/callmip_phase1a_v2 \
      experiments/callmip_phase1a_v2/generate_observations.jl

# 0b. Prior check (GATE — must pass before proceeding)
sbatch experiments/callmip_phase1a_v2/slurm/slurm_prior_check.sh

# 1. EKI calibration (~12 h, 32 tasks, 8 threads each)
sbatch experiments/callmip_phase1a_v2/slurm/slurm_calibration.sh

# 2. Emulate + sample (~4 h)
sbatch experiments/callmip_phase1a_v2/slurm/slurm_emulate_sample.sh

# 3. CalLMIP prior + posterior simulations (~4 h)
sbatch experiments/callmip_phase1a_v2/slurm/slurm_callmip_simulations.sh

# 4. Posterior ensemble (50 members, ~8 h)
sbatch experiments/callmip_phase1a_v2/slurm/slurm_posterior_ensemble.sh

# 5+6. Postprocessing (NetCDF + plots, ~1 h)
sbatch experiments/callmip_phase1a_v2/slurm/slurm_postprocess.sh
```

### Environment variables (already in each Slurm script)
```bash
export JULIA_DEPOT_PATH=/central/scratch/renatob/julia_depot
export JULIA_PROJECT=/central/scratch/renatob/ClimaLand.jl/experiments/callmip_phase1a_v2
export JULIA_NUM_THREADS=8
```

---

## Output files

| Path | Description |
|---|---|
| `experiments/callmip_phase1a_v2/observations.jld2` | 16 yearly obs windows |
| `experiments/callmip_phase1a_v2/output_eki/` | EKI iteration checkpoints |
| `experiments/callmip_phase1a_v2/output_fixed_window/` | Fixed-window EKI re-run for GP |
| `experiments/callmip_phase1a_v2/output_ces/` | GP emulator + posterior MCMC samples |
| `experiments/callmip_phase1a_v2/output_callmip_sims/` | Prior + posterior CalLMIP simulations |
| `experiments/callmip_phase1a_v2/output_posterior_ensemble/` | 50-member ensemble |
| `experiments/callmip_phase1a_v2/callmip_output/*.nc` | **CalLMIP NetCDF submission files** |
| `experiments/callmip_phase1a_v2/figures/` | Diagnostic plots |

### CalLMIP NetCDF files
- `ClimaLand._Phase1a_Scen1_DK-Sor_Cal_Prior.nc`
- `ClimaLand._Phase1a_Scen1_DK-Sor_Cal_Posterior.nc`

Both files: 6574 days (1997-01-01–2014-12-31), dims `(time,lat,lon)` with lat=1, lon=1.

---

## Calibrated parameters (16)

| # | Name | Prior type | Bounds |
|---|---|---|---|
| 1 | `pmodel_α` | Gaussian | [0,1] |
| 2 | `pmodel_β_c3` | Gaussian | [1,500] |
| 3 | `pmodel_β_c4` | Gaussian | [1,500] |
| 4 | `pmodel_cstar` | Gaussian | [0,1000] |
| 5 | `moisture_stress_c` | Gaussian | [0,20] |
| 6 | `soilCO2_reference_rate` | Lognormal | [0,∞) |
| 7 | `soilCO2_activation_energy` | Gaussian | [0,∞) |
| 8 | `michaelis_constant` | Gaussian | [0,∞) |
| 9 | `O2_michaelis_constant` | Lognormal | [0,∞) |
| 10 | `root_leaf_nitrogen_ratio` | Gaussian | [0,1] |
| 11 | `relative_contribution_factor` | Gaussian | [0,1] |
| 12 | `leaf_Cd` | Gaussian | [0,1] |
| 13 | `canopy_z_0m_coeff` | Gaussian | [0,1] |
| 14 | `canopy_z_0b_coeff` | Gaussian | [0,1] |
| 15 | `canopy_d_coeff` | Gaussian | [0,1] |
| 16 | `canopy_K_lw` | Gaussian | [0,1] |

---

## Known issues and limitations

1. **Winter NEE (DAMM respiration)**: The DAMM model tends to overestimate winter
   soil respiration. The lognormal prior on `soilCO2_reference_rate` (log-mean 6×10⁻⁸
   kg C/m³/s) improves the prior, but the posterior will still show a DJF NEE bias.
   This is a known limitation of the current ClimaLand heterotrophic respiration
   parameterisation and is documented in `RESULTS.md`.

2. **Aboveground biomass (`cLiveBioAbove`) and bare-soil evaporation (`evspsblsoi`)**:
   Not available in the current ClimaLand.jl `main` diagnostic output. Written as NaN
   in the NetCDF with attribute `note = "not available in current ClimaLand.jl version"`.

3. **Posterior GP emulator**: Training uses a fixed 1997-2012 observation window
   (consistent across all ensemble members) rather than the minibatched windows used
   during EKI. This is by design to ensure the GP input/output matrices are
   dimensionally consistent.

4. **pmodel_β_c3 collapse**: Monitor the posterior mean of `pmodel_β_c3`. If it falls
   below ~30 the model's C3 Vcmax will be unrealistically low; `emulate_sample.jl`
   emits a warning in this case.

5. **Timestep**: Model runs at Δt = 900 s. The PLUMBER2 forcing is half-hourly but
   ClimaLand interpolates linearly between half-hourly values.
