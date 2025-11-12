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

Creates a JLD2 file containing:
- observation_vector: Concatenated valid SMAP observations (no NaNs)
- valid_indices: Indices to extract same pixels from model output
- valid_mask: Boolean mask for model grid
- metadata: Configuration and statistics
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
    
    # Get land mask - use the landsea_mask directly without trying to convert it
    # The mask is already on the subsurface space, we need to work with it differently
    land_mask_3d = ClimaLand.Domains.landsea_mask(domain; threshold = mask_threshold)
    
    # Extract surface values by taking the top layer
    # Get the surface space coordinates
    coords = ClimaCore.Fields.coordinate_field(surface_space)
    
    # Create a land mask field on the surface space
    # We'll extract this by interpolating from the 3D mask
    land_mask_field = ClimaCore.Fields.zeros(FT, surface_space)
    
    # Get parent arrays to work with
    # The 3D mask has the same horizontal structure, just with vertical levels
    # We can copy the top level values
    parent_3d = parent(land_mask_3d)
    parent_2d = parent(land_mask_field)
    
    # Copy top level of 3D mask to 2D surface mask
    # Assuming IJFH layout: [I, J, F, H] where we want the top of the vertical column
    if size(parent_3d, 3) > 0  # Check vertical dimension exists
        # Take the top level (index 1 in vertical direction)
        parent_2d .= parent_3d[:, :, 1, :]
    else
        # If no vertical dimension, the mask might already be 2D
        parent_2d .= parent_3d
    end
    
    # Extract model coordinates
    model_coords = extract_model_coords_from_space(surface_space, land_mask_field)
    
    @info "Model grid setup" total_points=length(model_coords.lons) land_points=sum(model_coords.land_mask)
    
    # Process each date range
    all_obs_vectors = []
    all_valid_indices = []
    all_metadata = []
    
    for (idx, (start_str, stop_str)) in enumerate(sample_date_ranges)
        start_date = DateTime(start_str)
        stop_date = DateTime(stop_str)
        
        @info "Processing SMAP period $idx/$(length(sample_date_ranges))" start_date stop_date
        
        # Find SMAP files
        smap_files = find_smap_files(start_date, stop_date)
        
        if isempty(smap_files)
            error("No SMAP files found for period: $start_date to $stop_date")
        end
        
        @info "Found SMAP files" n_files=length(smap_files)
        
        # Calculate temporal mean over the period
        lons_smap, lats_smap, sm_smap = calculate_smap_temporal_mean(
            smap_files;
            quality_flag_threshold
        )
        
        # Create valid observation mask
        valid_data = create_valid_observation_mask(
            lons_smap,
            lats_smap,
            sm_smap,
            model_coords
        )
        
        # Store results
        push!(all_obs_vectors, valid_data.observation_vector)
        push!(all_valid_indices, valid_data.valid_indices)
        push!(all_metadata, Dict(
            "start_date" => start_str,
            "stop_date" => stop_str,
            "n_observations" => length(valid_data.observation_vector),
            "coverage_stats" => valid_data.coverage_stats
        ))
        
        @info "Period processed" period_idx=idx n_obs=length(valid_data.observation_vector) coverage_pct=valid_data.coverage_stats["land_coverage_pct"]
    end
    
    # Concatenate all periods
    observation_vector = vcat(all_obs_vectors...)
    
    # For simplicity, use indices from first period
    # (assumes same spatial coverage across periods)
    # If coverage varies, you'll need to handle this differently
    valid_indices = all_valid_indices[1]
    
    @info "Total observations" n_obs=length(observation_vector) n_periods=length(sample_date_ranges)
    
    # Save to file
    JLD2.jldsave(
        output_filepath;
        observation_vector,
        valid_indices,
        valid_mask = model_coords.land_mask,  # Full mask for reference
        nelements,
        sample_date_ranges,
        quality_flag_threshold,
        mask_threshold,
        metadata = all_metadata,
        created = now()
    )
    
    @info "Saved SMAP observation vector" filepath=output_filepath size=length(observation_vector)
    
    return observation_vector
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    # Match your calibration configuration
    generate_smap_observation_vector(
        sample_date_ranges = [("2016-01-01", "2016-03-31")],  # Example: Q1 2016
        nelements = (101, 15),
        mask_threshold = 0.75,
        quality_flag_threshold = 0,
        output_filepath = "experiments/calibration/land_observation_vector_SMAP.jld2"
    )
end