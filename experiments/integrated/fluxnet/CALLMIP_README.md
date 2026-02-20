# CalMIP DK-Sor Site Simulation Guide

This directory contains the setup for running CalMIP (Calibration Land Model Intercomparison Project) simulations for the DK-Sor (Denmark Sorø) flux tower site in ClimaLand.jl.

## Files Created

### 1. Site Configuration
- **Location**: `ext/fluxnet_simulations/DK-Sor.jl`
- Contains site-specific parameters:
  - Location: 55.486°N, 11.6446°E, UTC+1
  - Vegetation: Deciduous Broadleaf Forest (Beech)
  - Soil parameters: porosity, hydraulic conductivity, van Genuchten parameters
  - Canopy parameters: height (30m), Vcmax25, stomatal conductance

### 2. Simulation Runner
- **File**: `experiments/integrated/fluxnet/callmip_dksor.jl`
- Runs full year simulation (2008)
- Timestep: 450 seconds (7.5 minutes)
- Saves 16 diagnostic variables at half-hourly resolution

### 3. Forcing Data
- **Location**: `/net/sampo/data1/renatob/callmip_forcing/DK-Sor.csv`
- **Source**: Converted from `DK-Sor_1997-2014_FLUXNET2015_Met.nc`
- **Variables**: Air temperature, VPD, pressure, precipitation, wind, radiation, CO2
- **Period**: 1997-2015 (315,552 half-hourly timesteps)
- **Conversion Script**: `convert_netcdf_fast.py` (Python)

### 4. Observation Data
- **Location**: `/net/sampo/data1/renatob/DK-Sor_daily_aggregated_1997-2013_FLUXNET2015_Flux.nc`
- **Variables**: LE (latent heat), H (sensible heat), NEE (net ecosystem exchange)
- **Usage**: For model validation and comparison

### 5. Analysis Tools

#### Comparison Script
- **File**: `compare_callmip_observations.py`
- **Purpose**: Compare simulation outputs with flux observations
- **Features**:
  - Reads HDF5 simulation outputs
  - Reads NetCDF flux observations
  - Calculates statistics (RMSE, bias, correlation)
  - Creates timeseries and scatter plots

#### Monitoring Script
- **File**: `monitor_callmip.sh`
- **Purpose**: Monitor simulation progress
- **Usage**: `./monitor_callmip.sh`

## Running the Simulation

### Quick Start
```bash
cd $HOME/ClimaLand.jl
julia --project=. experiments/integrated/fluxnet/callmip_dksor.jl
```

### Background Run (for long simulations)
```bash
cd $HOME/ClimaLand.jl
nohup julia --project=. experiments/integrated/fluxnet/callmip_dksor.jl > callmip_run.log 2>&1 &
```

### Monitor Progress
```bash
./experiments/integrated/fluxnet/monitor_callmip.sh
# or
tail -f callmip_run.log
```

## Output Structure

### Simulation Outputs
- **Directory**: `experiments/integrated/fluxnet/DK-Sor/callmip/out/`
- **Format**: HDF5 files (one per diagnostic variable)
- **Frequency**: Half-hourly (30-minute intervals)

### Diagnostic Variables (16 total)
Key outputs include:
- `lhf.hdf5` - Latent heat flux (W/m²)
- `shf.hdf5` - Sensible heat flux (W/m²)
- `gpp.hdf5` - Gross primary production (μmol/m²/s)
- `sw.hdf5` - Soil water content (m³/m³)
- `st.hdf5` - Soil temperature (K)
- Additional canopy and soil variables

## Comparing with Observations

Once the simulation completes, compare with observations:

```bash
python3 experiments/integrated/fluxnet/compare_callmip_observations.py \
  --flux-file /net/sampo/data1/renatob/DK-Sor_daily_aggregated_1997-2013_FLUXNET2015_Flux.nc \
  --output-dir experiments/integrated/fluxnet/DK-Sor/callmip/out/ \
  --start-date 2008-01-01 \
  --end-date 2009-01-01 \
  --plot-dir experiments/integrated/fluxnet/DK-Sor/callmip/plots/
```

This will generate:
- Timeseries comparison plots for LE, H, NEE
- Scatter plots with statistics (RMSE, bias, R)
- Console output with quantitative metrics

## Data Processing Pipeline

### NetCDF → CSV Conversion
The forcing data was converted using:
```bash
python3 experiments/integrated/fluxnet/convert_netcdf_fast.py \
  /net/sampo/data1/renatob/DK-Sor_1997-2014_FLUXNET2015_Met.nc \
  /net/sampo/data1/renatob/callmip_forcing/DK-Sor.csv
```

### Data Flow
```
NetCDF Forcing → CSV → ClimaLand Simulation → HDF5 Outputs
                                                    ↓
NetCDF Observations ← ← ← ← ← ← Comparison Script
```

## Key Modifications to ClimaLand

1. **Artifacts.jl**
   - Added `callmip_data_path()` function
   - Points to `/net/sampo/data1/renatob/callmip_forcing/`

2. **data_processing.jl**
   - Updated `read_fluxnet_data()` to recognize CalMIP sites
   - Checks for "DK-Sor" and routes to CalMIP data path

3. **initial_conditions.jl**
   - Updated `set_fluxnet_ic!()` to handle CalMIP sites
   - Routes CalMIP sites to appropriate data path

## Simulation Details

### Model Components
- **Soil**: Energy-Hydrology model with biogeochemistry
  - 20 layers, exponentially stretched
  - Depth: 0 to -10m
  - van Genuchten hydraulic model

- **Canopy**: Coupled photosynthesis and plant hydraulics
  - Two-stream radiative transfer
  - Farquhar photosynthesis model
  - Medlyn stomatal conductance
  - Plant hydraulics with Weibull vulnerability curve

- **CO2**: Soil CO2 production and transport

### Timestepping
- Integration: IMEX (Implicit-Explicit) scheme
- Timestep: 450 seconds (7.5 minutes)
- Adaptive stepping with Newton solver

### Boundary Conditions
- Top: Atmospheric forcing from CSV data
- Bottom: Zero flux (deep drainage allowed through runoff)

## Troubleshooting

### Common Issues

1. **Package not found errors**
   - Solution: Use `--project=.` flag (not `--project=experiments`)
   - Ensures proper package environment

2. **Data path errors**
   - Check that CSV file exists at: `/net/sampo/data1/renatob/callmip_forcing/DK-Sor.csv`
   - Verify file permissions

3. **Slow compilation**
   - First run takes 2-5 minutes to compile
   - Subsequent runs are much faster
   - Use `--compile=min` to reduce compilation time

4. **Memory issues**
   - Year-long simulation uses ~1.5-2GB RAM
   - Reduce diagnostic frequency if needed

### Checking Simulation Status
```bash
# Check if running
ps aux | grep "julia.*callmip_dksor"

# Check progress
tail -f callmip_run.log

# Check output files
ls -lh experiments/integrated/fluxnet/DK-Sor/callmip/out/
```

## Performance Expectations

- **Compilation time**: 2-5 minutes (first run)
- **Simulation time** (1 year, half-hourly output):
  - Expected: 5-15 minutes on typical workstation
  - Varies with CPU speed and available cores
- **Output size**: ~100-500 MB total (all diagnostics)

## Next Steps

1. ✅ Run full year simulation (2008)
2. ⏳ Wait for completion (~10-15 min)
3. 📊 Compare with observations
4. 📈 Generate validation plots
5. 📝 Calculate performance metrics
6. 🔄 Extend to multi-year runs if needed

## References

- CalMIP Protocol: https://github.com/callmip-org/Phase1
- FLUXNET2015: https://fluxnet.org/
- DK-Sor Site: Deciduous beech forest in Denmark
- ClimaLand.jl: https://github.com/CliMA/ClimaLand.jl

## Contact

For issues or questions about this CalMIP setup, see the ClimaLand.jl repository or CalMIP protocol documentation.
