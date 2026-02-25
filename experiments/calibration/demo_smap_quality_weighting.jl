# Demonstration of SMAP Quality Flag Filtering and Weighting
# This script shows how to use quality flags to:
# 1. Filter out bad quality data
# 2. Weight observations by quality in temporal averaging

using Dates
include("smap_utils.jl")

# ==============================================================================
# Example 1: Compare different quality flag thresholds (filtering only)
# ==============================================================================
function example_filtering()
    println("\n" * "="^80)
    println("Example 1: Quality Flag Filtering (Binary)")
    println("="^80)
    
    # Setup date range and find files
    start_date = DateTime(2020, 6, 1)
    stop_date = DateTime(2020, 6, 30)
    data_dir = get(ENV, "SMAP_DATA_PATH", "/net/sampo/data1/smap/L3_SM_P_E")
    
    files = find_smap_files(start_date, stop_date; data_dir)
    
    if isempty(files)
        println("No SMAP files found. Check SMAP_DATA_PATH.")
        return
    end
    
    println("\nFound $(length(files)) files for June 2020")
    
    # Compare different quality thresholds
    thresholds = [0, 1, 2, 3]  # 0=best only, 3=include all
    
    for threshold in thresholds
        println("\n--- Quality Flag Threshold: $threshold ---")
        
        lons, lats, sm_mean = calculate_smap_temporal_mean(
            files; 
            quality_flag_threshold=threshold
        )
        
        n_valid = count(.!isnan.(sm_mean))
        coverage_pct = round(100 * n_valid / length(sm_mean), digits=2)
        
        println("Valid pixels: $n_valid / $(length(sm_mean)) ($coverage_pct%)")
        
        if n_valid > 0
            valid_sm = sm_mean[.!isnan.(sm_mean)]
            println("SM range: [$(round(minimum(valid_sm), digits=3)), $(round(maximum(valid_sm), digits=3))] m³/m³")
            println("SM mean: $(round(mean(valid_sm), digits=3)) m³/m³")
        end
    end
end

# ==============================================================================
# Example 2: Quality-weighted temporal averaging
# ==============================================================================
function example_weighting()
    println("\n" * "="^80)
    println("Example 2: Quality-Weighted Temporal Averaging")
    println("="^80)
    
    # Setup
    start_date = DateTime(2020, 6, 1)
    stop_date = DateTime(2020, 6, 30)
    data_dir = get(ENV, "SMAP_DATA_PATH", "/net/sampo/data1/smap/L3_SM_P_E")
    
    files = find_smap_files(start_date, stop_date; data_dir)
    
    if isempty(files)
        println("No SMAP files found. Check SMAP_DATA_PATH.")
        return
    end
    
    # Compare different weighting schemes
    schemes = [:inverse, :exponential, :quadratic, :binary]
    
    for scheme in schemes
        println("\n--- Weighting Scheme: $scheme ---")
        
        result = calculate_smap_weighted_temporal_mean(
            files;
            quality_flag_threshold=3,  # Include all quality levels
            weighting_scheme=scheme
        )
        
        n_valid = count(.!isnan.(result.sm_mean))
        coverage_pct = round(100 * n_valid / length(result.sm_mean), digits=2)
        
        println("Valid pixels: $n_valid ($coverage_pct%)")
        
        if n_valid > 0
            valid_sm = result.sm_mean[.!isnan.(result.sm_mean)]
            valid_std = result.sm_std[.!isnan.(result.sm_std)]
            valid_weights = result.total_weight[result.total_weight .> 0]
            
            println("SM mean: $(round(mean(valid_sm), digits=3)) ± $(round(mean(valid_std), digits=3)) m³/m³")
            println("Avg total weight per pixel: $(round(mean(valid_weights), digits=2))")
            println("Avg observations per pixel: $(round(mean(result.n_obs[result.n_obs .> 0]), digits=1))")
        end
    end
end

# ==============================================================================
# Example 3: View quality flag weights
# ==============================================================================
function example_weight_mapping()
    println("\n" * "="^80)
    println("Example 3: Quality Flag Weight Mappings")
    println("="^80)
    
    quality_flags = 0:3
    schemes = [:inverse, :exponential, :quadratic, :binary]
    
    println("\nQuality Flag | " * join(string.(schemes), " | "))
    println("-" * "="^79)
    
    for flag in quality_flags
        weights = [quality_flag_to_weight(flag; scheme=s) for s in schemes]
        weight_strs = [lpad(round(w, digits=3), 8) for w in weights]
        
        flag_desc = if flag == 0
            "0 (best)    "
        elseif flag == 1
            "1 (good)    "
        elseif flag == 2
            "2 (marginal)"
        else
            "3 (poor)    "
        end
        
        println("$flag_desc | " * join(weight_strs, " | "))
    end
    
    println("\nWeighting Scheme Descriptions:")
    println("  - inverse:     1.0, 0.5, 0.25, 0.1")
    println("  - exponential: exp(-flag) - smooth decay")
    println("  - quadratic:   ((3-flag)/3)² - gentle then steep")
    println("  - binary:      Only flag=0 gets weight")
end

# ==============================================================================
# Example 4: Load a single file with quality flags
# ==============================================================================
function example_single_file_quality()
    println("\n" * "="^80)
    println("Example 4: Inspecting Quality Flags from a Single File")
    println("="^80)
    
    # Find one file
    data_dir = get(ENV, "SMAP_DATA_PATH", "/net/sampo/data1/smap/L3_SM_P_E")
    start_date = DateTime(2020, 6, 15)
    files = find_smap_files(start_date, start_date; data_dir)
    
    if isempty(files)
        println("No SMAP file found for 2020-06-15. Check SMAP_DATA_PATH.")
        return
    end
    
    file = files[1]
    println("\nLoading: $(basename(file))")
    
    # Load with quality flags
    lons, lats, sm, qual = load_smap_grid_with_quality(file)
    
    println("\nData loaded:")
    println("  Total pixels: $(length(sm))")
    println("  Valid SM: $(count(.!isnan.(sm)))")
    println("  Valid quality flags: $(count(qual .>= 0))")
    
    # Quality flag distribution
    valid_qual = qual[qual .>= 0]
    
    if !isempty(valid_qual)
        println("\nQuality Flag Distribution:")
        for flag in 0:3
            n = count(valid_qual .== flag)
            pct = round(100 * n / length(valid_qual), digits=2)
            println("  Flag $flag: $n pixels ($pct%)")
        end
    end
    
    # Statistics by quality flag
    println("\nSoil Moisture by Quality Flag:")
    for flag in 0:3
        mask = (qual .== flag) .& (.!isnan.(sm))
        n = count(mask)
        
        if n > 0
            sm_subset = sm[mask]
            sm_mean = round(mean(sm_subset), digits=3)
            sm_std = round(std(sm_subset), digits=3)
            println("  Flag $flag: n=$n, mean=$sm_mean ± $sm_std m³/m³")
        else
            println("  Flag $flag: n=0")
        end
    end
end

# ==============================================================================
# Main execution
# ==============================================================================
function main()
    println("SMAP Quality Flag Demonstration")
    println("This script shows how to filter and weight SMAP data by quality flags")
    
    # Check if data directory exists
    data_dir = get(ENV, "SMAP_DATA_PATH", "/net/sampo/data1/smap/L3_SM_P_E")
    if !isdir(data_dir)
        println("\n⚠️  WARNING: SMAP data directory not found: $data_dir")
        println("Set SMAP_DATA_PATH environment variable or update the default path.")
        println("\nTo set: export SMAP_DATA_PATH=/path/to/smap/data")
        return
    end
    
    # Run examples
    example_weight_mapping()
    example_single_file_quality()
    example_filtering()
    example_weighting()
    
    println("\n" * "="^80)
    println("Demo Complete!")
    println("="^80)
    println("\nKey Takeaways:")
    println("1. Use quality_flag_threshold to hard-filter bad data")
    println("2. Use weighting_scheme to give higher importance to good quality")
    println("3. Combine both: filter worst data (threshold=2) + weight remaining (:inverse)")
    println("4. Quality weights can be used in calibration cost functions")
end

# Run the demo
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
