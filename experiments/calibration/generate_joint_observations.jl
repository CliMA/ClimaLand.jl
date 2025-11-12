using JLD2
using Dates
using Statistics

include(joinpath(@__DIR__, "api.jl"))

"""
    generate_joint_observations(;
        lhf_filepath,
        smap_filepath,
        output_filepath,
        lhf_weight = 0.5,
        sm_weight = 0.5,
        normalize_by_variance = true
    )

Combine LH flux and SMAP soil moisture observation vectors for joint calibration.

# Arguments
- `lhf_weight`: Relative importance of LH flux (0-1, default 0.5 for equal contribution)
- `sm_weight`: Relative importance of SMAP (0-1, default 0.5 for equal contribution)
- `normalize_by_variance`: If true, also scale by observation variance to normalize units

The weighting ensures that LH and SMAP contribute equally to the total cost function,
regardless of the number of observations per variable.
"""
function generate_joint_observations(;
    lhf_filepath::String = "experiments/calibration/land_observation_vector.jld2",
    smap_filepath::String = "experiments/calibration/land_observation_vector_SMAP.jld2",
    output_filepath::String = "experiments/calibration/land_observation_vector_joint.jld2",
    lhf_weight::Float64 = 0.5,
    sm_weight::Float64 = 0.5,
    normalize_by_variance::Bool = true
)
    @info "Generating joint observation vector" lhf_filepath smap_filepath lhf_weight sm_weight
    
    # Load LH flux observations
    if !isfile(lhf_filepath)
        error("LH observation file not found: $lhf_filepath. Run generate_observations.jl first.")
    end
    lhf_data = JLD2.load(lhf_filepath)
    lhf_vector = lhf_data["observation_vector"]
    
    # Load SMAP observations (with pixel groupings)
    if !isfile(smap_filepath)
        error("SMAP observation file not found: $smap_filepath. Run generate_smap_observations.jl first.")
    end
    sm_data = JLD2.load(smap_filepath)
    sm_vector = sm_data["observation_vector"]
    sm_metadata = sm_data["metadata"]
    
    # Verify both use same date ranges
    lhf_ranges = get(lhf_data, "sample_date_ranges", nothing)
    sm_ranges = get(sm_data, "sample_date_ranges", nothing)
    
    if lhf_ranges != sm_ranges
        @warn "Date ranges mismatch between LH and SMAP observations!" lhf_ranges sm_ranges
    end
    
    # Calculate observation counts and variances
    n_lhf = length(lhf_vector)
    n_sm = length(sm_vector)
    
    var_lhf = var(lhf_vector)
    var_sm = var(sm_vector)
    
    @info "Observation statistics" n_lhf n_sm var_lhf var_sm
    
    # Calculate scaling factors to equalize contribution to cost function
    # Goal: lhf_weight * sum((lhf_scaled)²) ≈ sm_weight * sum((sm_scaled)²)
    
    # Method 1: Scale by sqrt(n) to normalize total contribution
    scale_lhf = sqrt(lhf_weight / n_lhf)
    scale_sm = sqrt(sm_weight / n_sm)
    
    # Method 2: Also normalize by variance (puts both on similar numerical scale)
    if normalize_by_variance
        scale_lhf *= 1.0 / sqrt(var_lhf)
        scale_sm *= 1.0 / sqrt(var_sm)
    end
    
    @info "Scaling factors" scale_lhf scale_sm ratio=scale_sm/scale_lhf
    
    # Apply scaling (this effectively creates weighted observations)
    lhf_scaled = lhf_vector .* scale_lhf
    sm_scaled = sm_vector .* scale_sm
    
    # Concatenate: [LH observations, SMAP observations]
    observation_vector = vcat(lhf_scaled, sm_scaled)
    
    # Also save unscaled for reference
    observation_vector_unscaled = vcat(lhf_vector, sm_vector)
    
    @info "Joint observation vector created" n_lhf n_sm total=length(observation_vector)
    
    # Calculate expected cost contributions
    expected_lhf_contribution = n_lhf * scale_lhf^2 * var_lhf
    expected_sm_contribution = n_sm * scale_sm^2 * var_sm
    total_expected = expected_lhf_contribution + expected_sm_contribution
    
    # Print spatial resolution info
    println("\n===== Spatial Resolution =====")
    println("LH flux: ERA5 native resolution (~0.25°)")
    println("SMAP SM: 9 km EASE-Grid (~0.09° at equator)")
    for (i, meta) in enumerate(sm_metadata)
        if haskey(meta["coverage_stats"], "avg_model_points_per_smap_pixel")
            avg_pts = meta["coverage_stats"]["avg_model_points_per_smap_pixel"]
            println("  Period $i: averaging ~$(round(avg_pts, digits=1)) model points per SMAP pixel")
        end
    end
    
    # Print weighting info
    println("\n===== Observation Weighting =====")
    println("LH flux observations:   $n_lhf")
    println("SMAP observations:      $n_sm")
    println("Ratio (LH/SMAP):        $(round(n_lhf/n_sm, digits=2))")
    println("\nScaling factors:")
    println("  LH flux scale:        $(round(scale_lhf, digits=6))")
    println("  SMAP scale:           $(round(scale_sm, digits=6))")
    println("\nExpected cost contributions:")
    println("  LH flux:              $(round(100*expected_lhf_contribution/total_expected, digits=2))%")
    println("  SMAP:                 $(round(100*expected_sm_contribution/total_expected, digits=2))%")
    
    # Save combined file
    JLD2.jldsave(
        output_filepath;
        observation_vector,  # Scaled/weighted version (use this for calibration)
        observation_vector_unscaled,  # For diagnostics
        lhf_vector_scaled = lhf_scaled,
        sm_vector_scaled = sm_scaled,
        lhf_vector_unscaled = lhf_vector,
        sm_vector_unscaled = sm_vector,
        lhf_metadata = get(lhf_data, "metadata", Dict()),
        sm_metadata = sm_metadata,
        sample_date_ranges = lhf_ranges,
        nelements = lhf_data["nelements"],
        # Weighting metadata
        lhf_weight,
        sm_weight,
        scale_lhf,
        scale_sm,
        normalize_by_variance,
        created = now(),
        description = "Joint observation vector: ERA5 LH + SMAP SM (weighted: LH=$lhf_weight, SM=$sm_weight, var_norm=$normalize_by_variance)"
    )
    
    @info "Saved joint observation vector" filepath=output_filepath
    
    # Print summary
    println("\n===== Joint Observation Summary =====")
    println("Total observations:    $(length(observation_vector))")
    println("LH fraction:           $(round(100*n_lhf/length(observation_vector), digits=2))%")
    println("SMAP fraction:         $(round(100*n_sm/length(observation_vector), digits=2))%")
    
    if haskey(sm_data, "metadata")
        println("\n===== SMAP Coverage by Period =====")
        for (i, meta) in enumerate(sm_data["metadata"])
            println("Period $i ($(meta["start_date"]) to $(meta["stop_date"])): $(meta["n_observations"]) obs, $(meta["coverage_stats"]["land_coverage_pct"])% coverage")
        end
    end
    
    return observation_vector
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    generate_joint_observations(
        lhf_weight = 0.5,   # 50% contribution from LH
        sm_weight = 0.5,    # 50% contribution from SMAP
        normalize_by_variance = true
    )
end