# Calibrated uSPAC Long Run

This directory contains a production long run script that uses the calibrated hydraulic trait coordination parameters from iteration_002.

## Files

- **`calibrated_uspac_longrun.jl`** - Main simulation script
- **`PBS_calibrated_uspac_longrun.pbs`** - PBS job submission script for Derecho
- **`CALIBRATED_LONGRUN_README.md`** - This file

## What's Different from Standard Long Runs?

Standard long run scripts (`snowy_land.jl`, etc.) use:
- `PModelConductance` or basic `PlantHydraulicsModel` with Weibull conductivity
- Default parameters from `toml/default_parameters.toml`

**This calibrated version uses:**
- ✅ `uSPACConductancePi` conductance model (trait coordination framework)
- ✅ Calibrated β parameters from iteration_002 (automatically loaded from default TOML)
- ✅ Global Aridity Index as spatial covariate (P/ET0)
- ✅ Validated against TRY Plant Trait Database (P50 median diff: 0.07 MPa)

## Quick Start

### Option 1: Submit to PBS Queue (Recommended for Production)

For a **2-year test run**:
```bash
cd /glade/derecho/scratch/reich/ClimaLand.jl
qsub experiments/long_runs/PBS_calibrated_uspac_longrun.pbs
```

For a **19-year full run**, edit the PBS script first:
```bash
# Uncomment line 18 in PBS_calibrated_uspac_longrun.pbs:
# export LONGER_RUN=""
qsub experiments/long_runs/PBS_calibrated_uspac_longrun.pbs
```

Check job status:
```bash
qstat -u $USER
```

### Option 2: Interactive Run (For Testing)

```bash
cd /glade/derecho/scratch/reich/ClimaLand.jl
module load julia/1.11.2
module load climacommon  # For GPU support

# 2-year test run (CPU)
julia --project=.buildkite experiments/long_runs/calibrated_uspac_longrun.jl

# 19-year full run (CPU)
LONGER_RUN="" julia --project=.buildkite experiments/long_runs/calibrated_uspac_longrun.jl
```

## Configuration

### Duration
- **Default (2 years)**: March 2008 - March 2010
- **LONGER_RUN (19 years)**: March 2000 - March 2019

Set via environment variable:
```bash
export LONGER_RUN=""  # Enable 19-year run
# OR leave unset for 2-year test
```

### Resolution
- Horizontal: 101 spectral elements
- Vertical: 15 spectral elements
- Timestep: 450 seconds

### Output
Results saved to: `calibrated_uspac_longrun_cpu/` or `calibrated_uspac_longrun_gpu/`
- Global diagnostics: NetCDF files
- Visualizations: Annual timeseries, heatmaps, leaderboard plots
- Parameter log: `parameters.toml`

## Calibrated Parameters Used

The following parameters are automatically loaded from `toml/default_parameters.toml`:

```toml
βkx_base   = 1.916   # Base hydraulic conductivity coordination
βkx_coord  = 0.745   # Coordination slope (kx vs P50)
βψx50_base = 1.371   # Base P50 water potential
βψx50_slope = -1.765 # P50 aridity dependence
βΠR_base   = 0.334   # Base osmotic potential
βΠR_slope  = -0.892  # Osmotic potential aridity dependence
```

These were calibrated against SMAP soil moisture data and validated against the TRY Plant Trait Database.

## Key Features

### 1. Trait Coordination
Plant hydraulic traits (P50, hydraulic conductivity, osmotic potential) vary spatially based on:
- **Aridity** (P/ET0): Drier climates → more negative P50, adjusted conductivity
- **Coordination**: Safety-efficiency tradeoffs built into parameter structure

### 2. Aridity Field
Loaded from: `experiments/calibration/aridity.jld2`
- Global Aridity Index (Trabucco & Zomer 2019)
- Regridded to model resolution
- Ocean points masked using topographic threshold

### 3. Trait Distribution (Optional)
Set `use_trait_distribution = true` in the script (line 211) to enable:
- Within-gridcell trait heterogeneity
- Climate-dependent variance
- Trait correlations (safety-efficiency tradeoffs)
- Gaussian-Hermite quadrature integration (n=3 default)

## Resource Requirements

### 2-Year Run
- Walltime: ~6-8 hours
- Memory: ~50 GB
- GPUs: 1 (optional, can run on CPU)
- CPUs: 4-8

### 19-Year Run  
- Walltime: ~24 hours (adjust PBS script if needed)
- Memory: ~100 GB
- GPUs: 1 (recommended)
- CPUs: 8

## Monitoring Progress

### Check PBS job
```bash
qstat -u $USER           # Job status
qstat -f <JOBID>         # Detailed info
cat longrun_output.txt   # Standard output
cat longrun_error.txt    # Error output
```

### Monitor output files
```bash
ls -lh calibrated_uspac_longrun_*/global_diagnostics/
tail -f longrun_output.txt  # Watch progress
```

## Expected Output

### Diagnostics (NetCDF)
- Soil moisture, temperature profiles
- Canopy fluxes (GPP, transpiration, sensible heat)
- Snow water equivalent, albedo
- Stomatal conductance, water potential

### Visualizations
- `annual_timeseries_*.png` - Time series of key variables
- `heatmap_*.png` - Spatial patterns at final timestep  
- `leaderboard_*.png` - Model performance metrics

## Validation

The calibrated parameters have been validated:

✅ **P50 (Water Potential at 50% Conductance Loss)**
- Model median: -2.42 MPa
- TRY database median: -2.49 MPa
- Difference: 0.07 MPa (2.8% error)

✅ **Global Patterns**
- Realistic spatial variation (wet → -0.67 MPa, dry → -3.94 MPa)
- Captures literature ranges (Choat et al. 2012, Liu et al. 2019)

See `experiments/calibration/CALIBRATION_INTEGRATION_SUMMARY.md` for details.

## Troubleshooting

### Aridity file not found
```
ERROR: Aridity file not found at .../aridity.jld2
```
**Solution**: Make sure you're running from the ClimaLand.jl root directory:
```bash
cd /glade/derecho/scratch/reich/ClimaLand.jl
julia --project=. experiments/long_runs/calibrated_uspac_longrun.jl
```

### Parameter not found in TOML
```
ERROR: KeyError: key "βkx_base" not found
```
**Solution**: Verify parameters are in `toml/default_parameters.toml`:
```bash
grep "βkx_base" toml/default_parameters.toml
```
They should have been added automatically. If missing, see:
`experiments/calibration/CALIBRATION_INTEGRATION_SUMMARY.md`

### CUDA not found error
```
ERROR: ArgumentError: Package CUDA not found in current path
```
**Solution**: The PBS script has been updated to use `--project=.buildkite` which includes CUDA. Make sure you're using the latest version of the PBS script.

### Out of memory
**Solution**: Request more memory in PBS script:
```bash
#PBS -l select=1:ncpus=8:ngpus=1:mem=150GB
```

### Job timeout
**Solution**: Increase walltime in PBS script:
```bash
#PBS -l walltime=36:00:00
```

## Next Steps

### Analysis
Use ClimaAnalysis tools to compare with observations:
- FLUXNET sites
- ILAMB benchmarks  
- Satellite products (MODIS GPP, SMAP soil moisture)

### Sensitivity Tests
Modify parameters in the script:
- `use_trait_distribution = true` - Enable trait heterogeneity
- `n_quad = 5` - Increase quadrature accuracy
- Change variance parameters (σ_P50_base, etc.)

### Longer Runs
For multi-decadal runs:
- Adjust walltime and resources in PBS script
- Consider checkpoint/restart capabilities
- Monitor disk space for diagnostics

## References

- **Calibration**: SMAP soil moisture L3 passive product (2015-2023)
- **Validation**: TRY Plant Trait Database v6.0
- **Aridity**: Global Aridity Index (Trabucco & Zomer 2019)
- **Forcing**: ERA5 reanalysis (hourly)
- **LAI**: MODIS MCD15A2H v006

## Support

For questions or issues:
1. Check `experiments/calibration/CALIBRATION_INTEGRATION_SUMMARY.md`
2. See `experiments/calibration/template_long_run_with_calibration.jl` for examples
3. Refer to `experiments/calibration/model_interface.jl` for implementation details
