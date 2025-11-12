using Dates
using JLD2
import ClimaComms
import ClimaCore
import ClimaLand

include(joinpath(@__DIR__, "smap_utils.jl"))
include(joinpath(@__DIR__, "api.jl"))

"""
    generate_smap_observation_vector(;
        sample_date_ranges,
        nelements,
        mask_threshold = 0.75,
        quality_flag_threshold = 0,
        output_filepath = "experiments/calibration/land_observation_vector_SMAP.jld2"
    )

Generate observation vector from SMAP soil moisture data.
Handles SMAP's sparse coverage by tracking valid_indices separately for each period.

Creates a JLD2 file containing:
- observation_vector: Concatenated SMAP observations (NaNs removed)
- metadata: Array of Dict, one per period, containing:
  - valid_indices: Model grid indices with SMAP coverage for this period
  - n_observations: Count of valid obs
  - coverage_stats: Coverage percentages
"""
function generate_smap_observation_vector(;
    sample_date_ranges::Vector{Tuple{String,String}},
    nelements::Tuple{Int,Int},
    mask_threshold::Float64 = 0.75,
    quality_flag_threshold::Int = 0,
    output_filepath::String = "experiments/calibration/land_observation_vector_SMAP.jld2"
)
    FT = Float64
    context = ClimaComms.context()
    
    @info "Generating SMAP observation vector" n_periods=length(sample_date_ranges) nelements quality_flag_threshold
    
    # Create model domain to get grid and land mask
    domain = ClimaLand.Domains.global_domain(FT; nelements, mask_threshold, context)
    surface_space = domain.space.surface
    
    # Extract surface land mask
    land_mask_3d = ClimaLand.Domains.landsea_mask(domain; threshold = mask_threshold)
    land_mask_field = ClimaCore.Fields.zeros(FT, surface_space)
    parent_3d = parent(land_mask_3d)
    parent_2d = parent(land_mask_field)
    
    if size(parent_3d, 3) > 0
        parent_2d .= parent_3d[:, :, 1, :]
    else
        parent_2d .= parent_3d
    end
    
    # Extract model coordinates
    model_coords = extract_model_coords_from_space(surface_space, land_mask_field)
    
    @info "Model grid setup" total_points=length(model_coords.lons) land_points=sum(model_coords.land_mask)
    
    # Process each date range **SEPARATELY** (SMAP coverage varies by period)
    all_obs_vectors = []
    all_metadata = []
    
    for (idx, (start_str, stop_str)) in enumerate(sample_date_ranges)
        start_date = DateTime(start_str)
        stop_date = DateTime(stop_str)
        
        @info "Processing SMAP period $idx/$(length(sample_date_ranges))" start_date stop_date
        
        # Find SMAP files for this period
        smap_files = find_smap_files(start_date, stop_date)
        
        if isempty(smap_files)
            @warn "No SMAP files found for period: $start_date to $stop_date, skipping"
            # Create empty entry to maintain period indexing
            push!(all_metadata, Dict(
                "start_date" => start_str,
                "stop_date" => stop_str,
                "n_observations" => 0,
                "smap_pixel_groups" => Vector{Vector{Int}}(),
                "smap_pixel_indices" => Int[],
                "coverage_stats" => Dict("land_coverage_pct" => 0.0)
            ))
            push!(all_obs_vectors, Float64[])
            continue
        end
        
        @info "Found SMAP files" n_files=length(smap_files)
        
        # Calculate temporal mean over this period
        lons_smap, lats_smap, sm_smap = calculate_smap_temporal_mean(
            smap_files;
            quality_flag_threshold
        )
        
        # Create valid observation mask FOR THIS PERIOD
        # **NOW RETURNS PIXEL GROUPINGS**
        valid_data = create_valid_observation_mask(
            lons_smap,
            lats_smap,
            sm_smap,
            model_coords
        )
        
        # Store period-specific results
        push!(all_obs_vectors, valid_data.observation_vector)
        push!(all_metadata, Dict(
            "start_date" => start_str,
            "stop_date" => stop_str,
            "n_observations" => length(valid_data.observation_vector),
            "smap_pixel_groups" => valid_data.smap_pixel_groups,  # **NEW: for averaging model output**
            "smap_pixel_indices" => valid_data.smap_pixel_indices,
            "coverage_stats" => valid_data.coverage_stats,
            "n_smap_files" => length(smap_files)
        ))
        
        @info "Period processed" period_idx=idx n_obs=length(valid_data.observation_vector) coverage_pct=valid_data.coverage_stats["land_coverage_pct"] avg_points_per_pixel=valid_data.coverage_stats["avg_model_points_per_smap_pixel"]
    end
    
    # Concatenate all periods
    observation_vector = vcat(all_obs_vectors...)
    
    @info "Total observations across all periods" n_obs=length(observation_vector) n_periods=length(sample_date_ranges)
    
    # Verify observation vector matches metadata
    expected_length = sum(m["n_observations"] for m in all_metadata)
    if length(observation_vector) != expected_length
        error("Observation vector length mismatch: got $(length(observation_vector)), expected $expected_length")
    end
    
    # Save to file
    JLD2.jldsave(
        output_filepath;
        observation_vector,
        metadata = all_metadata,  # **Array of Dict, one per period**
        nelements,
        sample_date_ranges,
        quality_flag_threshold,
        mask_threshold,
        created = now(),
        description = "SMAP L3 soil moisture observation vector with period-specific valid_indices"
    )
    
    @info "Saved SMAP observation vector" filepath=output_filepath total_obs=length(observation_vector) n_periods=length(all_metadata)
    
    # Print summary statistics
    for (i, meta) in enumerate(all_metadata)
        @info "Period $i summary" dates="$(meta["start_date"]) to $(meta["stop_date"])" n_obs=meta["n_observations"] coverage=meta["coverage_stats"]["land_coverage_pct"]
    end
    
    return observation_vector
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    # Use your actual calibration date ranges
    include(joinpath(@__DIR__, "run_uspac_calibration.jl"))  # To get generate_seasonal_date_ranges
    
    generate_smap_observation_vector(
        sample_date_ranges = generate_seasonal_date_ranges(2015, 2023),
        nelements = (101, 15),
        mask_threshold = 0.75,
        quality_flag_threshold = 0,
        output_filepath = "experiments/calibration/land_observation_vector_SMAP.jld2"
    )
end