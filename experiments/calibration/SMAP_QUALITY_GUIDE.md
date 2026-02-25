# SMAP Quality Flag Usage Guide

## Overview

SMAP L3 soil moisture retrievals include quality flags (0-3 scale):
- **0**: Recommended quality - best retrievals
- **1**: Good quality - some correction needed, still reliable
- **2**: Marginal quality - use with caution
- **3**: Poor quality - not recommended

## Two Approaches to Quality Control

### 1. **Hard Filtering** (existing functionality, enhanced)
Filter out observations below a quality threshold.

```julia
# Only keep "recommended" quality (flag = 0)
lons, lats, sm = calculate_smap_temporal_mean(files; quality_flag_threshold=0)

# Keep recommended + good (flags 0-1)
lons, lats, sm = calculate_smap_temporal_mean(files; quality_flag_threshold=1)

# Keep all except poor (flags 0-2)
lons, lats, sm = calculate_smap_temporal_mean(files; quality_flag_threshold=2)
```

**Pros**: Simple, removes unreliable data  
**Cons**: Wastes marginal data that could contribute with lower weight

### 2. **Quality Weighting** (NEW)
Include all data but weight by quality in averaging.

```julia
# Calculate quality-weighted mean with inverse weighting
result = calculate_smap_weighted_temporal_mean(
    files;
    quality_flag_threshold=3,      # Include all quality levels
    weighting_scheme=:inverse       # Best gets 1.0, poor gets 0.1
)

# Access results
sm_mean = result.sm_mean           # Weighted mean soil moisture
sm_std = result.sm_std             # Standard deviation (unweighted)
total_weight = result.total_weight # Sum of weights per pixel
n_obs = result.n_obs              # Number of observations per pixel
```

**Pros**: Uses all available data, gives appropriate importance  
**Cons**: More complex, requires choosing weighting scheme

## Weighting Schemes

| Scheme | Flag 0 | Flag 1 | Flag 2 | Flag 3 | Best For |
|--------|--------|--------|--------|--------|----------|
| `:inverse` | 1.0 | 0.5 | 0.25 | 0.1 | General use, balanced |
| `:exponential` | 1.0 | 0.368 | 0.135 | 0.050 | Aggressive deweighting |
| `:quadratic` | 1.0 | 0.444 | 0.111 | 0.0 | Gentle then steep drop |
| `:binary` | 1.0 | 0.0 | 0.0 | 0.0 | Equivalent to threshold=0 |

Test with:
```julia
# See all weight mappings
for flag in 0:3
    w = quality_flag_to_weight(flag; scheme=:inverse)
    println("Flag $flag -> weight $w")
end
```

## Recommended Workflow

### For Calibration/Data Assimilation

**Option A: Conservative (recommended for initial work)**
```julia
# Filter out poor quality, weight the rest
result = calculate_smap_weighted_temporal_mean(
    files;
    quality_flag_threshold=2,      # Reject poor quality (flag=3)
    weighting_scheme=:inverse       # Weight remaining by quality
)

# Use weights in calibration cost function
obs_data = create_valid_observation_mask_weighted(
    result.lons, result.lats, result.sm_mean, result.total_weight,
    model_coords;
    min_weight=0.5  # Only use pixels with decent total weight
)

# In cost function:
weighted_residuals = obs_data.observation_weights .* (model - obs_data.observation_vector)
cost = sum(weighted_residuals.^2)
```

**Option B: Maximize data coverage**
```julia
# Include all data with aggressive weighting
result = calculate_smap_weighted_temporal_mean(
    files;
    quality_flag_threshold=3,       # Include all
    weighting_scheme=:exponential   # Strong preference for good quality
)
```

**Option C: Quality before quantity**
```julia
# Hard filter for best quality only
lons, lats, sm = calculate_smap_temporal_mean(files; quality_flag_threshold=0)
```

### For Spatial Analysis/Climatology

```julia
# More inclusive for coverage
result = calculate_smap_weighted_temporal_mean(
    files;
    quality_flag_threshold=2,
    weighting_scheme=:inverse
)

# Can also examine uncertainty
high_uncertainty = result.sm_std .> 0.1  # High variability pixels
low_weight = result.total_weight .< 5.0   # Low confidence pixels
```

## Inspecting Quality Flags

Load a single file to examine quality distribution:
```julia
lons, lats, sm, qual = load_smap_grid_with_quality("path/to/smap/file.h5")

# Count by quality flag
for flag in 0:3
    n = count(qual .== flag)
    println("Flag $flag: $n pixels")
end

# Statistics by quality
for flag in 0:3
    mask = (qual .== flag) .& (.!isnan.(sm))
    if count(mask) > 0
        sm_subset = sm[mask]
        println("Flag $flag: mean=$(mean(sm_subset)), std=$(std(sm_subset))")
    end
end
```

## Combining with Existing Workflows

### Point-based extraction (existing `load_smap_sm`)
Still uses hard filtering (quality_flag == 0).

### Grid-based calibration (NEW)
```julia
# 1. Calculate weighted temporal mean
smap_result = calculate_smap_weighted_temporal_mean(
    files; 
    quality_flag_threshold=2,
    weighting_scheme=:inverse
)

# 2. Create observation mask with weights
obs_data = create_valid_observation_mask_weighted(
    smap_result.lons, 
    smap_result.lats, 
    smap_result.sm_mean, 
    smap_result.total_weight,
    model_coords;
    min_weight=0.5  # Optional: filter low-confidence pixels
)

# 3. Use in calibration
# obs_data.observation_vector -> SMAP values
# obs_data.observation_weights -> quality-based weights
# obs_data.smap_pixel_coords -> model points to average per SMAP pixel
```

## Tips

1. **Start conservative**: Use `quality_flag_threshold=0` or `1` to see how much data you get
2. **Check coverage**: Print `n_valid_pixels` to see if you need to relax thresholds
3. **Examine weights**: Look at `result.total_weight` distribution - low values indicate few/poor observations
4. **Validate**: Compare results from different thresholds/schemes to ensure stability
5. **Document choice**: Record which threshold and weighting scheme you use for reproducibility

## Example: Full Pipeline

```julia
using Dates
include("smap_utils.jl")

# 1. Find files for your period
files = find_smap_files(
    DateTime(2020, 6, 1), 
    DateTime(2020, 8, 31);
    data_dir="/path/to/smap"
)

# 2. Calculate weighted temporal mean
smap_result = calculate_smap_weighted_temporal_mean(
    files;
    quality_flag_threshold=2,       # Reject poor quality
    weighting_scheme=:inverse        # Weight by quality
)

# 3. Create observation mask for your model
obs_data = create_valid_observation_mask_weighted(
    smap_result.lons,
    smap_result.lats,
    smap_result.sm_mean,
    smap_result.total_weight,
    model_coords;
    min_weight=1.0  # Require at least unit total weight
)

# 4. Use in calibration or validation
println("Found $(length(obs_data.observation_vector)) valid SMAP pixels")
println("Coverage: $(obs_data.coverage_stats["land_coverage_pct"])% of land")
println("Mean weight: $(obs_data.coverage_stats["mean_weight"])")

# 5. In your cost function
function cost_function(model_output)
    # Average model output to SMAP scale
    model_averaged = [mean(model_output[group]) for group in obs_data.smap_pixel_groups]
    
    # Weighted residuals
    residuals = model_averaged - obs_data.observation_vector
    weighted_residuals = obs_data.observation_weights .* residuals
    
    return sum(weighted_residuals.^2)
end
```
