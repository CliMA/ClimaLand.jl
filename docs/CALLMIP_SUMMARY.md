# CalMIP Site Setup - Summary

## What Was Created

A complete template setup for running CalMIP (California Land Model Intercomparison Project) sites in ClimaLand.jl, using DK-Sor (Denmark Sorø) as an example test site.

## Files Created/Modified

### 1. Documentation
- **`docs/CALLMIP_SETUP.md`** - Complete setup guide with detailed instructions
- **`docs/CALLMIP_QUICKSTART.md`** - Quick reference guide for adding new sites
- **`docs/CALLMIP_SUMMARY.md`** (this file) - Overview summary

### 2. Data Files
- **`data/callmip_sites/DK-Sor.csv`** - Example CSV template with forcing data format
  - Contains 24 rows of example meteorological data
  - Shows required column structure and units
  - **Note:** Replace with actual CalMIP data before running real simulations

### 3. Core Code Updates

#### src/Artifacts.jl
- Added `callmip_data_path()` function to handle CalMIP site data paths
- Supports local data directory until CalMIP data is added to artifacts

#### ext/fluxnet_simulations/DK-Sor.jl (NEW)
- Complete site configuration for DK-Sor test site
- Includes functions for:
  - `get_domain_info()` - Soil depth and discretization
  - `get_location()` - Lat/lon and time zone
  - `get_fluxtower_height()` - Measurement height
  - `get_parameters()` - Soil, canopy, and plant hydraulics parameters

#### ext/FluxnetSimulationsExt.jl
- Added include for DK-Sor.jl site configuration

#### ext/fluxnet_simulations/data_processing.jl
- Modified `read_fluxnet_data()` to handle both Fluxnet and CalMIP sites

### 4. Simulation Runner
- **`experiments/integrated/fluxnet/callmip_dksor.jl`** (NEW)
  - Complete simulation setup for DK-Sor site
  - Includes soil, canopy, and coupled model configuration
  - Sets up diagnostics and output
  - Based on ozark_pft.jl template

## DK-Sor Site Specifications

### Location
- **Site:** Denmark Sorø
- **Coordinates:** 55.486°N, 11.6446°E
- **Time Zone:** UTC+1 (CET)
- **Tower Height:** 57 m

### Vegetation
- **Type:** Deciduous Broadleaf Forest
- **Dominant Species:** Beech (Fagus sylvatica)
- **Canopy Height:** 30 m
- **Stem Height:** 25 m
- **Leaf Layer:** 5 m
- **Rooting Depth:** 1.0 m

### Soil Domain
- **Depth:** 10 m (zmin = -10, zmax = 0)
- **Elements:** 20
- **Porosity:** 0.50
- **Hydraulic Conductivity:** 5×10⁻⁷ m/s

## Required CSV Data Format

The forcing data CSV must contain these columns:

| Column | Description | Units |
|--------|-------------|-------|
| TIMESTAMP_START | Start of period (YYYYMMDDHHMM) | - |
| TIMESTAMP_END | End of period (YYYYMMDDHHMM) | - |
| TA_F | Air temperature | °C |
| VPD_F | Vapor pressure deficit | hPa |
| PA_F | Air pressure | kPa |
| P_F | Precipitation | mm |
| WS_F | Wind speed | m/s |
| LW_IN_F | Downwelling LW radiation | W/m² |
| SW_IN_F | Downwelling SW radiation | W/m² |
| CO2_F_MDS | CO2 concentration | μmol/mol |

## How to Add a New CalMIP Site

### Quick Steps:

1. **Create forcing data CSV**
   ```
   data/callmip_sites/YOUR-SITE.csv
   ```

2. **Update Artifacts.jl**
   - Add site to `callmip_data_path()` assertion

3. **Create site configuration**
   ```
   ext/fluxnet_simulations/YOUR-SITE.jl
   ```
   (Copy and modify DK-Sor.jl)

4. **Update FluxnetSimulationsExt.jl**
   - Add include statement

5. **Update data_processing.jl**
   - Add site to `callmip_sites` tuple

6. **Create simulation runner**
   ```
   experiments/integrated/fluxnet/callmip_yoursite.jl
   ```
   (Copy and modify callmip_dksor.jl)

### Detailed Instructions
See `docs/CALLMIP_SETUP.md` for complete step-by-step guide.

## Where to Get CalMIP Data

### Forcing Data
1. **CalMIP GitHub Repository:**
   - https://github.com/callmip-org/Phase1/tree/main/Data

2. **Model Evaluation Portal:**
   - https://modelevaluation.org

3. **Box Link (example):**
   - https://caltech.box.com/shared/static/otrr2y0rgjct7hqhmq214nb8qjsvqj5p.gz

### Observation Data (for calibration)
- Available in CalMIP Phase 1 Data directory
- Includes NEE, Qle, Qh with uncertainties
- See Phase 1a-test folder for DK-Sor example

### Site Metadata
- Refer to CalMIP Phase 1 Protocol (PDF)
- Includes soil texture, PFT cover percentages
- Available as NetCDF global attributes

## CalMIP Output Variables

The simulation is configured to output:

**Surface (1D) variables:**
- `gpp` - Gross Primary Productivity
- `er` - Ecosystem Respiration  
- `et` - Evapotranspiration
- `shf` - Sensible Heat Flux
- `lhf` - Latent Heat Flux
- `rn` - Net Radiation
- `swu`, `lwu` - Upwelling radiation
- `gs` - Stomatal Conductance
- `ct` - Canopy Temperature
- Others...

**Profile (2D) variables:**
- `swc` - Soil Water Content (by depth)
- `tsoil` - Soil Temperature (by depth)
- `si` - Snow/Ice

## Running the Example

```bash
# Navigate to ClimaLand directory
cd ClimaLand.jl

# Run the DK-Sor simulation
julia --project=experiments experiments/integrated/fluxnet/callmip_dksor.jl
```

**Note:** The example CSV contains only template data. Replace it with actual CalMIP forcing data before running real simulations.

## Output

Results are saved to:
```
experiments/integrated/fluxnet/DK-Sor/callmip/out/
```

Includes:
- Timeseries plots (GPP, fluxes, soil moisture, temperature)
- Diurnal cycle analysis
- Comparison with observations (if available)

## Next Steps

1. **Get Real Data:**
   - Download actual DK-Sor forcing data from CalMIP repository
   - Replace template CSV with real data

2. **Test Run:**
   - Run short simulation (1-7 days) to verify setup
   - Check outputs and diagnostics

3. **Full Simulation:**
   - Extend to full simulation period (months to years)
   - Generate outputs required by CalMIP protocol

4. **Calibration:**
   - Use ClimaLand calibration framework
   - Follow CalMIP Phase 1 protocol for calibration targets
   - See `experiments/calibration/` directory

5. **Add More Sites:**
   - Follow the same pattern for other CalMIP sites
   - CalMIP Phase 1 includes multiple forest, grassland, and cropland sites

## Important Notes

### Site ID Naming Convention
- Use **hyphens** in file names and strings: `"DK-Sor"`
- Use **underscores** in Val types: `Val(:DK_Sor)`
- This is handled automatically by `replace_hyphen()` function

### Data Sources Must Match
- Forcing data temporal resolution should match model timestep
- Typically half-hourly (1800s) or hourly
- Gap-fill missing data before running

### Initial Conditions
- Generic initial conditions are used by default
- For better results, spin up the model for 1-2 years
- Or use initial conditions from previous runs

### Calibration Targets
- CalMIP specifies which fluxes to calibrate against
- Typically: NEE, Qle, Qh with uncertainties provided
- See CalMIP protocol for full list

## References

### CalMIP Project
- **Website:** https://callmip-org.github.io
- **GitHub:** https://github.com/callmip-org/Phase1
- **Protocol:** CalMIP Phase 1 Protocol v1.1 (PDF in repository)

### ClimaLand
- **Documentation:** https://clima.github.io/ClimaLand.jl/dev/
- **GitHub:** https://github.com/CliMA/ClimaLand.jl
- **Issues:** https://github.com/CliMA/ClimaLand.jl/issues

### FLUXNET
- **Data Format:** https://fluxnet.org/data/fluxnet2015-dataset/
- Used as basis for CalMIP forcing data format

## Troubleshooting

### Common Issues

**"Site not found" error:**
- Ensure site is added to all necessary files (see checklist above)
- Check spelling consistency (DK-Sor vs DK_Sor)

**"Cannot find CSV file":**
- Verify file path: `data/callmip_sites/SITE-ID.csv`
- Check file permissions

**"Method not found for Val":**
- Site configuration file not included in FluxnetSimulationsExt.jl
- Val type uses underscore, not hyphen

**"Data format error":**
- Verify CSV has exactly the required columns
- Check for missing headers
- Ensure numeric data is not quoted

### Getting Help

- **CalMIP questions:** Visit https://callmip-org.github.io
- **ClimaLand questions:** Open issue at ClimaLand.jl GitHub
- **Documentation:** Read `docs/CALLMIP_SETUP.md` for detailed guide

## File Checklist

When adding a new site, ensure you've created/modified:

- [ ] CSV forcing data file
- [ ] Site configuration (.jl file)
- [ ] Artifacts.jl (add site to assertion)
- [ ] FluxnetSimulationsExt.jl (include site file)
- [ ] data_processing.jl (add to callmip_sites tuple)
- [ ] Simulation runner script
- [ ] Test with short simulation
- [ ] Verify outputs

## Example Site Comparison

| Aspect | DK-Sor (CalMIP) | US-MOz (Fluxnet) |
|--------|-----------------|------------------|
| Location | Denmark | Missouri, USA |
| Vegetation | Beech Forest | Oak Forest |
| Canopy Height | 30 m | 18.5 m |
| Time Zone | UTC+1 | UTC-6 |
| Data Source | CalMIP | AmeriFlux |
| Purpose | Calibration Test | Model Validation |

Both use the same underlying ClimaLand model structure and follow similar setup procedures.

## Version Information

- **ClimaLand.jl:** Compatible with current main branch
- **CalMIP Protocol:** Phase 1 Protocol v1.1
- **Date Created:** February 2026
- **Status:** Template setup complete, awaiting real forcing data

---

**For detailed instructions, see:** `docs/CALLMIP_SETUP.md`  
**For quick reference, see:** `docs/CALLMIP_QUICKSTART.md`
