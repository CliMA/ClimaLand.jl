# CalMIP Quick Start Guide

## Summary

This guide provides a quick reference for setting up CalMIP sites in ClimaLand.jl.

## Required Files for Each Site

| File Type | Location | Purpose |
|-----------|----------|---------|
| CSV forcing data | `data/callmip_sites/SITE-ID.csv` | Meteorological forcing |
| Site configuration | `ext/fluxnet_simulations/SITE-ID.jl` | Site parameters |
| Simulation runner | `experiments/integrated/fluxnet/callmip_SITE-ID.jl` | Run script |

## Climate Variables Required (CSV columns)

```
TIMESTAMP_START,TIMESTAMP_END,TA_F,VPD_F,PA_F,P_F,WS_F,LW_IN_F,SW_IN_F,CO2_F_MDS
```

| Variable | Description | Units |
|----------|-------------|-------|
| `TIMESTAMP_START` | Start time (YYYYMMDDHHMM) | - |
| `TIMESTAMP_END` | End time (YYYYMMDDHHMM) | - |
| `TA_F` | Air temperature | °C |
| `VPD_F` | Vapor pressure deficit | hPa |
| `PA_F` | Air pressure | kPa |
| `P_F` | Precipitation | mm |
| `WS_F` | Wind speed | m/s |
| `LW_IN_F` | Longwave radiation (down) | W/m² |
| `SW_IN_F` | Shortwave radiation (down) | W/m² |
| `CO2_F_MDS` | CO2 concentration | μmol/mol |

## Quick Setup Steps

### 0. Prepare CalMIP forcing data

CalMIP CSV data files are NOT stored in git. Download from CalMIP and set up access:

```bash
# Option 1: Use environment variable (recommended)
export CALLMIP_DATA_PATH="/path/to/your/callmip/data"
mkdir -p $CALLMIP_DATA_PATH

# Option 2: Use default local path
mkdir -p ../callmip_data

# Download and place your CSV files there
```

### 1. Update Artifacts.jl (if adding new sites)
```julia
# Add to: src/Artifacts.jl
function callmip_data_path(site_ID; context = nothing)
    @assert site_ID ∈ ("DK-Sor", "YOUR-SITE")  # ← Add here (only if adding new sites)
    # ...
end
```

### 2. Create site config file
```bash
# Copy template:
cp ext/fluxnet_simulations/DK-Sor.jl ext/fluxnet_simulations/YOUR-SITE.jl
# Edit YOUR-SITE.jl with your site's parameters
```

### 3. Update FluxnetSimulationsExt.jl
```julia
# Add to: ext/FluxnetSimulationsExt.jl
include("fluxnet_simulations/YOUR-SITE.jl")  # ← Add here
```

### 4. Create simulation runner
```bash
# Copy template:
cp experiments/integrated/fluxnet/callmip_dksor.jl \
   experiments/integrated/fluxnet/callmip_yoursite.jl
# Update site_ID and dates in callmip_yoursite.jl
```

### 5. Run simulation
```bash
# Make sure CALLMIP_DATA_PATH is set or data is in ../callmip_data/
export CALLMIP_DATA_PATH="/path/to/callmip/data"

julia --project=experiments experiments/integrated/fluxnet/callmip_yoursite.jl
```

## Site Information Needed

For the site configuration file, you need:

**Location:**
- Latitude (decimal degrees)
- Longitude (decimal degrees)
- Time zone offset from UTC (hours)

**Tower/Measurements:**
- Measurement height (m)

**Domain:**
- Soil depth (m, typically 10m or deeper)
- Vertical discretization (typically 20 elements)

**Soil Properties:**
- Porosity (ν)
- Saturated hydraulic conductivity (K_sat)
- Van Genuchten parameters (α, n)
- Residual water content (θ_r)
- Soil texture fractions

**Vegetation:**
- Canopy height (m)
- Stem height (m)
- Leaf height (m)
- Rooting depth (m)
- Plant functional type
- Typical LAI range

## Example Site: DK-Sor

```julia
# Site ID
site_ID = "DK-Sor"

# Location (Denmark)
lat = 55.486°N
long = 11.6446°E
time_offset = +1  # CET (UTC+1)

# Tower
atmos_h = 57 m

# Vegetation
- Type: Deciduous Broadleaf Forest (Beech)
- h_canopy = 30 m
- h_stem = 25 m
- h_leaf = 5 m
- rooting_depth = 1.0 m
```

## CalMIP Data Sources

1. **Protocol & Documentation:**
   - https://github.com/callmip-org/Phase1

2. **Forcing Data:**
   - Model Evaluation Portal: https://modelevaluation.org
   - CalMIP GitHub: https://github.com/callmip-org/Phase1/tree/main/Data

3. **Observation Data (for calibration):**
   - https://github.com/callmip-org/Phase1/tree/main/Data

## File Naming Convention

**Use hyphens in file names and site_ID strings:**
```julia
site_ID = "DK-Sor"  # ✓ Correct
site_ID = "DK_Sor"  # ✗ Wrong
```

**Use underscores in Val types:**
```julia
Val(:DK_Sor)  # ✓ Correct
Val(:DK-Sor)  # ✗ Wrong
```

## Typical Parameter Values

### For Deciduous Broadleaf Forests

```julia
# Soil
soil_ν = 0.45 - 0.55       # Porosity
soil_K_sat = 1e-7 - 1e-6   # Hydraulic conductivity (m/s)
soil_vg_n = 1.5 - 2.5      # Van Genuchten n
soil_vg_α = 0.02 - 0.08    # Van Genuchten α (1/m)

# Canopy
Vcmax25 = 5e-5 - 8e-5      # Max carboxylation rate (mol/m²/s)
g1 = 120 - 180             # Medlyn g1 parameter
h_canopy = 15 - 35         # Canopy height (m)
rooting_depth = 0.5 - 2.0  # Root depth (m)
```

### For Grasslands

```julia
# Soil (similar ranges)
soil_ν = 0.40 - 0.50
soil_K_sat = 5e-7 - 5e-6

# Canopy
Vcmax25 = 4e-5 - 7e-5
g1 = 80 - 120
h_canopy = 0.3 - 1.0       # Much shorter
rooting_depth = 0.3 - 1.0  # Shallower
```

## Checking Your Setup

Before running a full simulation, verify:

```bash
# 1. Check CSV file exists and has correct format
head data/callmip_sites/YOUR-SITE.csv

# 2. Check Julia can load the site
julia --project=experiments -e '
using ClimaLand
site_ID = "YOUR-SITE"
site_val = Symbol(replace(site_ID, "-" => "_"))
println("Site: ", site_val)
'

# 3. Do a short test run (1 day)
# Edit your runner script to set:
# stop_date = start_date + Day(1)
```

## Common Errors

### "Site not found in assertion"
→ Add your site to `callmip_data_path` in `Artifacts.jl`

### "Cannot read CSV file"
→ Check path: `data/callmip_sites/YOUR-SITE.csv`

### "Val type not defined"
→ Include your site file in `FluxnetSimulationsExt.jl`

### "Method not found for Val{:YOUR_SITE}"
→ Use underscores in Val: `Val{:YOUR_SITE}` not `Val{:YOUR-SITE}`

## Next Steps

1. ✓ Complete setup (follow steps above)
2. Run test simulation (1-7 days)
3. Compare with observations
4. Set up calibration workflow
5. Follow CalMIP protocol for full runs

## Additional Resources

- Full documentation: `docs/CALLMIP_SETUP.md`
- Example site: `ext/fluxnet_simulations/DK-Sor.jl`
- Example runner: `experiments/integrated/fluxnet/callmip_dksor.jl`
- CalMIP website: https://callmip-org.github.io
