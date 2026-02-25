#!/usr/bin/env julia
# Quick test script for SMAP quality flag weighting
# Modify the parameters below to experiment

using Dates
using Statistics
include("smap_utils.jl")

# =============================================================================
# CONFIGURE THESE PARAMETERS
# =============================================================================

# Date range
START_DATE = DateTime(2020, 6, 1)
STOP_DATE = DateTime(2020, 6, 30)

# Data directory (or set SMAP_DATA_PATH environment variable)
DATA_DIR = get(ENV, "SMAP_DATA_PATH", "/net/sampo/data1/smap/L3_SM_P_E")

# Quality control settings
QUALITY_THRESHOLD = 2        # 0=best only, 1=good+, 2=marginal+, 3=all
WEIGHTING_SCHEME = :inverse  # :inverse, :exponential, :quadratic, :binary
MIN_WEIGHT = 0.5            # Minimum weight to include pixel

# =============================================================================
# RUN ANALYSIS
# =============================================================================

println("="^80)
println("SMAP Quality Flag Test")
println("="^80)
println("\nSettings:")
println("  Date range: $START_DATE to $STOP_DATE")
println("  Quality threshold: $QUALITY_THRESHOLD (0=best, 3=all)")
println("  Weighting scheme: $WEIGHTING_SCHEME")
println("  Min weight: $MIN_WEIGHT")
println()

# Find files
println("Searching for SMAP files...")
files = find_smap_files(START_DATE, STOP_DATE; data_dir=DATA_DIR)

if isempty(files)
    println("\n❌ No files found!")
    println("Check:")
    println("  1. Data directory exists: $DATA_DIR")
    println("  2. Files exist for date range")
    println("  3. File naming: SMAP_L3_SM_P_E_YYYYMMDD_R*.h5")
    exit(1)
end

println("✓ Found $(length(files)) files")

# Compare: No weighting vs Quality weighting
println("\n" * "="^80)
println("Comparison: Standard vs Quality-Weighted Averaging")
println("="^80)

println("\n1. Standard averaging (binary threshold at $QUALITY_THRESHOLD):")
lons1, lats1, sm1 = calculate_smap_temporal_mean(
    files; 
    quality_flag_threshold=QUALITY_THRESHOLD
)
n_valid1 = count(.!isnan.(sm1))
println("   Valid pixels: $n_valid1 ($(round(100*n_valid1/length(sm1), digits=2))%)")
if n_valid1 > 0
    valid_sm1 = sm1[.!isnan.(sm1)]
    println("   Mean SM: $(round(mean(valid_sm1), digits=4)) m³/m³")
    println("   Std SM: $(round(std(valid_sm1), digits=4)) m³/m³")
end

println("\n2. Quality-weighted averaging (scheme=$WEIGHTING_SCHEME):")
result = calculate_smap_weighted_temporal_mean(
    files;
    quality_flag_threshold=QUALITY_THRESHOLD,
    weighting_scheme=WEIGHTING_SCHEME
)
n_valid2 = count(.!isnan.(result.sm_mean))
println("   Valid pixels: $n_valid2 ($(round(100*n_valid2/length(result.sm_mean), digits=2))%)")
if n_valid2 > 0
    valid_sm2 = result.sm_mean[.!isnan.(result.sm_mean)]
    valid_weights = result.total_weight[result.total_weight .> 0]
    valid_nobs = result.n_obs[result.n_obs .> 0]
    
    println("   Mean SM: $(round(mean(valid_sm2), digits=4)) m³/m³")
    println("   Std SM: $(round(std(valid_sm2), digits=4)) m³/m³")
    println("   Mean total weight: $(round(mean(valid_weights), digits=2))")
    println("   Mean # obs/pixel: $(round(mean(valid_nobs), digits=1))")
end

# Show weight distribution for first file
println("\n" * "="^80)
println("Quality Flag Analysis (first file)")
println("="^80)

lons, lats, sm, qual = load_smap_grid_with_quality(files[1])

println("\nQuality flag distribution:")
for flag in 0:3
    n = count(qual .== flag)
    pct = round(100 * n / length(qual), digits=2)
    weight = quality_flag_to_weight(flag; scheme=WEIGHTING_SCHEME)
    
    flag_name = ["Recommended", "Good", "Marginal", "Poor"][flag+1]
    println("  Flag $flag ($flag_name): $n pixels ($pct%), weight=$weight")
end

# Soil moisture statistics by quality
println("\nSoil moisture by quality flag:")
for flag in 0:3
    mask = (qual .== flag) .& (.!isnan.(sm))
    n = count(mask)
    
    if n > 0
        sm_subset = sm[mask]
        flag_name = ["Recommended", "Good", "Marginal", "Poor"][flag+1]
        println("  Flag $flag ($flag_name): n=$n, mean=$(round(mean(sm_subset), digits=4)) m³/m³")
    end
end

# Demonstrate different weighting schemes
println("\n" * "="^80)
println("Impact of Different Weighting Schemes")
println("="^80)

schemes = [:inverse, :exponential, :quadratic, :binary]
scheme_results = Dict()

for scheme in schemes
    result = calculate_smap_weighted_temporal_mean(
        files;
        quality_flag_threshold=3,  # Include all
        weighting_scheme=scheme
    )
    
    valid_sm = result.sm_mean[.!isnan.(result.sm_mean)]
    if !isempty(valid_sm)
        scheme_results[scheme] = (
            mean=mean(valid_sm),
            std=std(valid_sm),
            n_valid=length(valid_sm)
        )
    end
end

println("\n$(lpad("Scheme", 12)) | $(lpad("Mean SM", 10)) | $(lpad("Std SM", 10)) | $(lpad("N Valid", 10))")
println("-"^50)
for scheme in schemes
    if haskey(scheme_results, scheme)
        r = scheme_results[scheme]
        println("$(lpad(string(scheme), 12)) | $(lpad(round(r.mean, digits=4), 10)) | $(lpad(round(r.std, digits=4), 10)) | $(lpad(r.n_valid, 10))")
    end
end

println("\n" * "="^80)
println("Summary")
println("="^80)
println("""
✓ Quality weighting allows using more data while maintaining quality control
✓ Good quality observations get higher weight in averaging
✓ Different schemes affect mean estimates (though usually small differences)
✓ Choose based on your application:
  - Conservative: binary or inverse with low threshold
  - Maximize coverage: exponential or quadratic with high threshold
  
Next steps:
1. Adjust QUALITY_THRESHOLD, WEIGHTING_SCHEME, MIN_WEIGHT above
2. Integrate weighted data into your calibration pipeline
3. See SMAP_QUALITY_GUIDE.md for detailed usage
""")
