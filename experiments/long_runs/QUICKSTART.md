# Quick Reference: Calibrated uSPAC Long Run

## 🚀 Submit Job to Derecho

```bash
cd /glade/derecho/scratch/reich/ClimaLand.jl

# 2-year test run (default, uses .buildkite project with GPU support)
qsub experiments/long_runs/PBS_calibrated_uspac_longrun.pbs

# Check status
qstat -u reich

# View output
tail -f longrun_output.txt
```

## 📊 For 19-Year Full Run

Edit `PBS_calibrated_uspac_longrun.pbs` line 18:
```bash
# Uncomment this line:
export LONGER_RUN=""
```

Then submit:
```bash
qsub experiments/long_runs/PBS_calibrated_uspac_longrun.pbs
```

## 📁 Output Location

```
calibrated_uspac_longrun_cpu/
├── global_diagnostics/     # NetCDF files
├── annual_timeseries_*.png # Time series plots
├── heatmap_*.png          # Spatial patterns
├── leaderboard_*.png      # Performance metrics
└── parameters.toml        # Parameter log
```

## ⚙️ Key Parameters (Auto-loaded)

```
βkx_base   = 1.916   βψx50_base = 1.371   βΠR_base   = 0.334
βkx_coord  = 0.745   βψx50_slope = -1.765  βΠR_slope  = -0.892
```

Source: `toml/default_parameters.toml` (iteration_002 calibration)

## 📝 What This Run Does

✅ Uses calibrated hydraulic trait coordination  
✅ Traits vary spatially with aridity (P/ET0)  
✅ Validated against TRY database (P50 diff: 0.07 MPa)  
✅ uSPAC conductance model (not standard PModel)  
✅ Global Aridity Index as covariate

## 🔍 Monitor Job

```bash
# Job status
qstat -u reich
qstat -f <JOBID>

# Live output
tail -f longrun_output.txt
tail -f longrun_error.txt

# Check results
ls -lh calibrated_uspac_longrun_*/
```

## 🆘 Quick Troubleshooting

| Problem | Solution |
|---------|----------|
| Job pending | Check queue: `qstat -Q` |
| CUDA not found | PBS uses `.buildkite` project (includes CUDA) |
| Out of memory | Edit PBS: `mem=150GB` |
| Timeout | Edit PBS: `walltime=36:00:00` |
| Aridity file not found | Run from ClimaLand.jl root |
| Parameter not in TOML | Check: `grep βkx_base toml/default_parameters.toml` |

## 📖 Full Documentation

- **Long run details**: `experiments/long_runs/CALIBRATED_LONGRUN_README.md`
- **Integration info**: `experiments/calibration/CALIBRATION_INTEGRATION_SUMMARY.md`  
- **Template**: `experiments/calibration/template_long_run_with_calibration.jl`
- **Reference**: `experiments/calibration/model_interface.jl`

## ⏱️ Expected Runtime

- **2-year run**: ~6-8 hours (1 GPU, 8 CPUs, 50GB RAM)
- **19-year run**: ~24 hours (1 GPU, 8 CPUs, 100GB RAM)

## 🎯 Validation Results

✅ **P50**: Model median -2.42 MPa vs TRY -2.49 MPa (diff: 0.07 MPa)  
✅ **Range**: -3.94 to -0.67 MPa (realistic global variation)  
✅ **Pattern**: Dry → more negative P50 ✓
