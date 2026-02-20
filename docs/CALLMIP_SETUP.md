# CalMIP Site Setup for ClimaLand.jl

This directory contains the setup for running CalMIP (California Land Model Intercomparison Project) sites in ClimaLand.jl.

## Overview

The CalMIP project is a model intercomparison and calibration effort to improve land surface models. This setup provides a template for adding CalMIP sites to ClimaLand.

**Current Status:** Template setup complete for DK-Sor (Denmark Sorø) test site from CalMIP Phase 1a.

## CalMIP Resources

- **GitHub Repository:** https://github.com/callmip-org/Phase1
- **Website:** https://callmip-org.github.io
- **Protocol:** [CalMIP Phase 1 Protocol v1.1](https://github.com/callmip-org/Phase1/blob/main/Protocol/CalLMIP_Phase1_Protocol_v1.1.pdf)

## File Structure

```
ClimaLand.jl/
├── src/
│   └── Artifacts.jl                      # Updated with callmip_data_path function
├── ext/
│   ├── FluxnetSimulationsExt.jl         # Updated to include CalMIP sites
│   └── fluxnet_simulations/
│       └── DK-Sor.jl                    # Site-specific parameters
└── experiments/
    └── integrated/
        └── fluxnet/
            └── callmip_dksor.jl         # Simulation runner script
```

**Note on Data Storage:** CalMIP forcing data files (CSV) are NOT stored in the git repository. 
Like other ClimaLand sites (US-MOz, US-Ha1, etc.), these should be managed through the 
ClimaArtifacts system or stored locally. See "Data Location" section below.

## Adding a New CalMIP Site

Follow these steps to add a new CalMIP site:

### 1. Prepare Climate Forcing Data (CSV file)

**Important:** CalMIP forcing data files should NOT be committed to the git repository.

Create a CSV file with the following columns and units:

| Column Name      | Description                           | Units       |
|-----------------|---------------------------------------|-------------|
| TIMESTAMP_START | Start of averaging period (YYYYMMDDHHMM) | -      |
| TIMESTAMP_END   | End of averaging period (YYYYMMDDHHMM)   | -      |
| TA_F            | Air temperature                       | °C          |
| VPD_F           | Vapor pressure deficit                | hPa         |
| PA_F            | Air pressure                          | kPa         |
| P_F             | Accumulated precipitation             | mm          |
| WS_F            | Wind speed                            | m/s         |
| LW_IN_F         | Downwelling longwave radiation        | W/m²        |
| SW_IN_F         | Downwelling shortwave radiation       | W/m²        |
| CO2_F_MDS       | CO2 concentration                     | μmol/mol    |

#### Data Location Options:

**Option 1: Environment Variable (Recommended for local testing)**
```bash
# Set the path to your CalMIP data directory
export CALLMIP_DATA_PATH="/path/to/your/callmip/data"
# Place your CSV files there: /path/to/your/callmip/data/DK-Sor.csv
```

**Option 2: Default Local Path**
```bash
# ClimaLand will look in: ../callmip_data/ (relative to ClimaLand.jl)
mkdir ../callmip_data
# Place your CSV files there
```

**Option 3: Add to ClimaArtifacts (Recommended for production)**
For permanent integration, CalMIP data should be added to the ClimaArtifacts system
(similar to how Fluxnet sites like US-MOz are stored). This requires:
1. Uploading data to Caltech's artifact storage
2. Updating the Artifacts.toml file
3. Modifying `callmip_data_path()` in src/Artifacts.jl to use `@clima_artifact`

**Where to get CalMIP data:**
- CalMIP GitHub repository: https://github.com/callmip-org/Phase1/tree/main/Data
- Model evaluation data portal: https://modelevaluation.org
- Or obtain from Box link: https://caltech.box.com/shared/static/otrr2y0rgjct7hqhmq214nb8qjsvqj5p.gz

Add your site to the `callmip_data_path` function in `src/Artifacts.jl`:

```julia
function callmip_data_path(site_ID; context = nothing)
    @assert site_ID ∈ ("DK-Sor", "YOUR-SITE-ID")  # Add your site here
    
    climaland_dir = pkgdir(@__MODULE__)
    data_path = joinpath(climaland_dir, "data", "callmip_sites", "$(site_ID).csv")
    return data_path
end
```

### 3. Create Site Information File

Create a new file in `ext/fluxnet_simulations/` named after your site (e.g., `YOUR-SITE-ID.jl`).

Use `ext/fluxnet_simulations/DK-Sor.jl` as a template. You need to define:

```julia
# Replace DK_Sor with your site ID (underscores instead of hyphens)
function FluxnetSimulations.get_domain_info(FT, ::Val{:YOUR_SITE_ID}; ...)
function FluxnetSimulations.get_location(FT, ::Val{:YOUR_SITE_ID}; ...)
function FluxnetSimulations.get_fluxtower_height(FT, ::Val{:YOUR_SITE_ID}; ...)
function FluxnetSimulations.get_parameters(FT, ::Val{:YOUR_SITE_ID}; ...)
```

**Key information needed:**
- Latitude and longitude
- Time zone offset from UTC
- Tower/measurement height
- Soil depth and properties
- Vegetation type and parameters
- Root depth and canopy height

### 4. Update FluxnetSimulationsExt.jl

Add your site file to the includes in `ext/FluxnetSimulationsExt.jl`:

```julia
include("fluxnet_simulations/YOUR-SITE-ID.jl")
```

### 5. Create Simulation Runner

Create a simulation file in `experiments/integrated/fluxnet/` (e.g., `callmip_yoursite.jl`).

Use `experiments/integrated/fluxnet/callmip_dksor.jl` as a template.

**Key modifications needed:**
- Update `site_ID = "YOUR-SITE-ID"`
- Adjust simulation dates: `start_date` and `stop_date`
- Modify output directory path if needed

### 6. Run the Simulation

First, ensure your CalMIP forcing data is accessible (via environment variable or local path):

```bash
# Option 1: Set environment variable
export CALLMIP_DATA_PATH="/path/to/callmip/data"

# Then run
cd ClimaLand.jl
julia --project=experiments experiments/integrated/fluxnet/callmip_yoursite.jl
```

## CalMIP-Specific Considerations

### Output Variables

CalMIP requires specific output variables. The current setup includes:

**1D (surface) variables:**
- GPP (Gross Primary Productivity)
- NEE (Net Ecosystem Exchange) - calculated as ER - GPP
- ET (Evapotranspiration)
- Qle (Latent heat flux)
- Qh (Sensible heat flux)
- Rnet (Net radiation)

**2D (profile) variables:**
- Soil water content
- Soil temperature

Refer to the CalMIP protocol for the complete list of required outputs.

### Calibration

The CalMIP project focuses on model calibration. After setting up your site:

1. Run the model with default/prior parameters
2. Compare outputs with observations (from CalMIP data repository)
3. Use calibration framework (see `experiments/calibration/`)
4. Follow CalMIP protocol for calibration procedures

### Data Sources

- **Forcing data:** Available from CalMIP repository or modelevaluation.org
- **Observation data:** Available from CalMIP Phase 1 data directory
- **Site metadata:** See CalMIP protocol for soil texture, PFT cover, etc.

## Example: DK-Sor Test Site

The DK-Sor (Denmark Sorø) site is set up as an example:

- **Location:** 55.486°N, 11.6446°E
- **Vegetation:** Deciduous Broadleaf Forest (Beech)
- **Purpose:** CalMIP Phase 1a test calibration

**Files:**
- Data: `data/callmip_sites/DK-Sor.csv` (example template only)
- Site config: `ext/fluxnet_simulations/DK-Sor.jl`
- Runner: `experiments/integrated/fluxnet/callmip_dksor.jl`

## Data Format Notes

### Forcing Data Requirements

The CSV forcing data follows FLUXNET format conventions:
- Timestamps in format YYYYMMDDHHMM
- Half-hourly or hourly temporal resolution
- Gap-filled meteorological data preferred
- Quality flags can be included but are not required

### Initial Conditions

The simulation will use generic initial conditions unless you specify custom ones:
- Soil moisture: From equilibrium or prior runs
- Soil temperature: From equilibrium or prior runs
- Canopy state: LAI from MODIS or prescribed

## Troubleshooting

### Common Issues

1. **Missing data file:** Ensure your CSV file is in `data/callmip_sites/`
2. **Site ID mismatch:** Use hyphens in site_ID string, underscores in Val{:SITE_ID}
3. **Forcing data format:** Verify column names match exactly (case-sensitive)
4. **Time zone issues:** Check time_offset sign (positive for east of UTC)

### Data Processing

If you need to convert CalMIP NetCDF data to CSV format:

```julia
# Example conversion script (adapt as needed)
using NCDatasets
using DataFrames, CSV

# Load CalMIP NetCDF file
ds = NCDataset("path/to/callmip_forcing.nc")

# Extract variables and convert to DataFrame
df = DataFrame(
    TIMESTAMP_START = ds["time_start"][:],
    TIMESTAMP_END = ds["time_end"][:],
    TA_F = ds["Tair"][:],
    # ... add other variables
)

# Save as CSV
CSV.write("data/callmip_sites/YOUR-SITE-ID.csv", df)
```

## References

- CalMIP Protocol: https://github.com/callmip-org/Phase1
- ClimaLand Documentation: https://clima.github.io/ClimaLand.jl/dev/
- FLUXNET Data Format: https://fluxnet.org/data/fluxnet2015-dataset/

## Contact

For CalMIP-specific questions, visit: https://callmip-org.github.io

For ClimaLand integration issues, open an issue at: https://github.com/CliMA/ClimaLand.jl/issues
