# Soil Albedo Coefficient Calibration

This directory contains scripts used to calibrate the `CompositionBasedSoilAlbedo` coefficients
for ClimaLand's soil model.

## Background

CLM-style soil color maps can have significant biases over bare soil regions, particularly in deserts
where albedo is often overestimated by 20-50%. The `CompositionBasedSoilAlbedo` parameterization
addresses this by computing albedo directly from soil properties using a physics-based approach.

## Model Formulation

The dry albedo is computed via logistic regression:

```
η = η₀ + c_om × ν_ss_om + c_vgn × vg_n + c_cf × ν_ss_gravel
α_dry = α_min + (α_max - α_min) × σ(η)
```

where:
- `σ(η) = 1/(1 + exp(-η))` is the logistic (sigmoid) function
- `ν_ss_om` is soil organic matter volume fraction
- `vg_n` is the van Genuchten n parameter (pore size distribution index)
- `ν_ss_gravel` is coarse fragments volume fraction
- `α_min = 0.04`, `α_max = 0.60` are physical albedo bounds

The logistic transformation ensures bounded output with smooth gradients.

Moisture darkening is applied nonlinearly:
```
α = α_dry × (1 - f_wet × S_e^β)
```

where `S_e` is effective saturation, `β = 0.5` captures rapid initial darkening.

## Data Sources

1. **CERES bareground albedo**: Satellite-derived annual mean shortwave albedo
   - Source: CERES EBAF (Loeb et al., 2018)
   - Resolution: 1° × 1°

2. **SoilGrids**: Soil composition data
   - Organic matter content (`nu_ss_om`)
   - Coarse fragments (`nu_ss_cf`)
   - Source: SoilGrids 2.0 (Poggio et al., 2021)

3. **van Genuchten parameters**: Soil hydraulic properties
   - n parameter
   - Strong correlation with soil texture

## Calibration Regions

We focused on desert regions where bare soil dominates:

| Region | Latitude | Longitude | Rationale |
|--------|----------|-----------|-----------|
| Sahara | 15°N - 35°N | 15°W - 35°E | Largest hot desert |
| Arabian | 15°N - 30°N | 35°E - 60°E | Sandy/rocky mix |
| Australian | 20°S - 30°S | 120°E - 145°E | Variable textures |
| Gobi | 38°N - 48°N | 90°E - 115°E | Cold desert |
| Kalahari | 20°S - 28°S | 18°E - 28°E | Semi-arid |

Total: ~1,432 valid grid points used for calibration.

## Calibrated Coefficients

The coefficients below were fitted via logistic regression, then bias-corrected
based on global ClimaLand simulations. The raw CERES bareground albedo data
appears to have a positive bias, so the intercepts (η₀) were adjusted downward
by ~0.3 to reduce the overall upwelling shortwave bias in coupled simulations.

### PAR Band
```julia
η₀_PAR = -3.3      # Bias-corrected intercept (raw fit: -3.04)
c_om_PAR = -0.13   # Organic matter (darkening)
c_vgn_PAR = 1.24   # van Genuchten n (texture)
c_cf_PAR = 0.15    # Coarse fragments (rocks)
```

### NIR Band
```julia
η₀_NIR = -3.4      # Bias-corrected intercept (raw fit: -3.10)
c_om_NIR = -0.14
c_vgn_NIR = 1.28
c_cf_NIR = 0.16
```

## Performance Metrics

| Metric | Value |
|--------|-------|
| R² | 0.51 | 
| RMSE | 0.074 | 

### Predictor Correlations

| Predictor | Correlation | Physical Interpretation |
|-----------|-------------|------------------------|
| vg_n | +0.71 | **Strongest**: sandier soils are brighter |
| Coarse fragments | +0.35 | Rocky surfaces are brighter |
| Organic matter | -0.25 | Organic matter darkens soil |

## Files

- `fit_composition_albedo.jl`: Main fitting script
  - Loads CERES + SoilGrids data
  - Extracts desert region data
  - Fits logistic regression
  - Outputs calibrated coefficients

- `evaluate_albedo_models.jl`: Model evaluation script
  - Compares composition-based vs mean baseline
  - Regional performance analysis
  - Predictor correlation analysis

## Usage

```julia
# Run the fitting procedure
include("fit_composition_albedo.jl")
results = main()

# Run the evaluation
include("evaluate_albedo_models.jl")
evaluation = main()
```

Note: The scripts require ClimaLand artifacts to be available. Data files are
automatically downloaded via `ClimaLand.Artifacts` when the scripts are run.

## Physical Interpretation

The fitted coefficients are physically meaningful:

1. **van Genuchten n** (`c_vgn > 0`): Higher n values indicate sandier, coarser soils.
   Sandy soils (primarily quartz) are bright, while clay soils are darker.
   This is the strongest single predictor (r = 0.71).

2. **Organic matter** (`c_om < 0`): Organic matter (humus) is dark and strongly
   absorbs visible light. More OM → lower albedo.

3. **Coarse fragments** (`c_cf > 0`): Rocky/gravelly surfaces tend to be brighter
   than fine-grained soils. Important for distinguishing sandy vs rocky deserts.

## References

- Lobell, D. B., & Asner, G. P. (2002). Moisture effects on soil reflectance.
  Soil Science Society of America Journal, 66(3), 722-727.

- Matthias, A. D., et al. (2000). Surface roughness effects on soil albedo.
  Soil Science Society of America Journal, 64(3), 1035-1041.
