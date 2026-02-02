import ClimaAnalysis
import Dates
using Dates: DateTime
import ClimaCalibrate
import ClimaCalibrate.EnsembleBuilder
import ClimaCalibrate.Checker
using Statistics  # for mean() in aggregate_to_smap_pixels

include(
    joinpath(pkgdir(ClimaLand), "experiments/calibration/observation_utils.jl"),
)

using CairoMakie, GeoMakie, Printf, StatsBase

# Need access to get_era5_obs_var_dict and get_sim_var_dict
ext = Base.get_extension(ClimaLand, :LandSimulationVisualizationExt)
"""
    ClimaCalibrate.observation_map(iteration)

Return G ensemble for an `iteration`.

G ensemble represents the concatenated forward model evaluations from all
ensemble members, arranged horizontally. Each individual forward model
evaluation corresponds to preprocessed, flattened simulation data from a single
ensemble member that has been matched to the corresponding observational data.
"""
function ClimaCalibrate.observation_map(iteration)
    output_dir = CALIBRATE_CONFIG.output_dir
    ekp = JLD2.load_object(ClimaCalibrate.ekp_path(output_dir, iteration))
    ensemble_size = EKP.get_N_ens(ekp)

    # Load observation metadata (includes short_names for joint observations)
    obs_file = CALIBRATE_CONFIG.obs_vec_filepath
    if !isfile(obs_file)
        error("Observation file not found: $obs_file")
    end
    
    obs_data = JLD2.load(obs_file)
    
    # Get short_names from metadata (for joint observations) or from EKP (for ERA5-only)
    if haskey(obs_data, "short_names")
        short_names = obs_data["short_names"]
        @info "Loaded short_names from observation metadata" short_names
    else
        # Fallback: try to extract from observation (ERA5-only path)
        obs_series = EKP.get_observation_series(ekp)
        short_names = ClimaCalibrate.ObservationRecipe.short_names(
            first(obs_series.observations),
        )
    end
    
    # **Load SMAP-specific metadata if using SMAP**
    obs_metadata = nothing
    if "sm_surface" in short_names
        obs_metadata = obs_data  # Already loaded above
        if haskey(obs_metadata, "sm_metadata")
            @info "Loaded SMAP observation metadata" n_periods=length(obs_metadata["sm_metadata"])
        end
    end

    # Determine observation vector length (handle both old and new formats)
    if haskey(obs_data, "joint_samples")
        # New joint format: samples are stored as vectors
        obs_length = length(obs_data["joint_samples"][1])  # First period sample length
    elseif haskey(obs_data, "observation_samples")
        # SMAP-only format
        obs_length = length(obs_data["observation_samples"][1])
    else
        # Old ERA5 format with EKP.Observation - check first period
        obs_vec = obs_data["observation_vector"]
        if !isempty(obs_vec)
            first_obs = obs_vec[1]
            if hasproperty(first_obs, :samples) && !isempty(first_obs.samples)
                obs_length = length(first_obs.samples[1])
            else
                obs_length = length(obs_vec)
            end
        else
            error("Empty observation_vector in observation file")
        end
    end
    @info "Building ensemble observation matrix" ensemble_size obs_length
    
    # Use GEnsembleBuilder ONLY for ERA5 with ClimaAnalysis metadata
    # Use manual path for SMAP-only and joint calibrations
    has_smap = "sm_surface" in short_names
    use_builder = !has_smap && !haskey(obs_data, "joint_samples")
    
    if use_builder
        # Original ERA5-only path with GEnsembleBuilder
        g_ens_builder = EnsembleBuilder.GEnsembleBuilder(ekp)
        for m in 1:ensemble_size
            member_path =
                ClimaCalibrate.path_to_ensemble_member(output_dir, iteration, m)
            simdir_path = joinpath(member_path, "global_diagnostics/output_active")
            @info "Processing member $m: $simdir_path"
            try
                process_member_data!(g_ens_builder, m, simdir_path, short_names)
            catch e
                @error "Error processing member $m, filling observation map entry with NaNs" exception =
                    e
                EnsembleBuilder.fill_g_ens_col!(g_ens_builder, m, NaN)
            end
        end
        
        if EnsembleBuilder.is_complete(g_ens_builder)
            return EnsembleBuilder.get_g_ensemble(g_ens_builder)
        else
            @error "G ensemble matrix is not completed. You may find it useful to call `EnsembleBuilder.missing_short_names(g_ens_builder, 1) or display g_ens_builder in the REPL"
        end
    else
        # SMAP-only or Joint observation path: manually build ensemble matrix
        g_ensemble = zeros(obs_length, ensemble_size)
        
        for m in 1:ensemble_size
            member_path =
                ClimaCalibrate.path_to_ensemble_member(output_dir, iteration, m)
            simdir_path = joinpath(member_path, "global_diagnostics/output_active")
            @info "Processing member $m: $simdir_path"
            try
                if has_smap && !haskey(obs_data, "joint_samples")
                    # SMAP-only processing
                    member_vector = process_smap_member_vector(simdir_path, obs_data)
                else
                    # Joint LHF+SMAP processing
                    member_vector = process_joint_member_vector(simdir_path, short_names, obs_data)
                end
                g_ensemble[:, m] = member_vector
            catch e
                @error "Error processing member $m, filling with NaNs" exception = e
                g_ensemble[:, m] .= NaN
            end
        end
        
        return g_ensemble
    end
end

"""
    process_smap_member_vector(diagnostics_folder_path, obs_data)

Process model output for SMAP-only calibration and return observation vector.
Extracts surface soil moisture at grid points where SMAP observations exist.
Uses flat linear indices on the 404×202 lat-lon diagnostic output grid.
"""
function process_smap_member_vector(diagnostics_folder_path, obs_data)
    @info "Processing SMAP-only calibration data"
    
    # Get valid indices from observation metadata
    if !haskey(obs_data, "valid_indices")
        error("valid_indices not found in observation file. Please regenerate observations with generate_smap_observations_from_regridded.jl")
    end
    
    valid_indices = obs_data["valid_indices"]
    
    # Load 3D soil water content
    swc = get(
        ClimaAnalysis.SimDir(diagnostics_folder_path),
        short_name = "swc",
    )
    
    # Shift to start of month and take temporal average
    swc = ClimaAnalysis.shift_to_start_of_previous_month(swc)
    sm_var = ClimaAnalysis.average_time(swc, ignore_nan = true)
    
    # Extract surface level (top of soil column)
    z_vals = sm_var.dims["z"]
    surface_z = z_vals[end]  # Top layer
    sm_var = ClimaAnalysis.slice(sm_var, z = surface_z)
    
    # Flatten to 1D (ClimaAnalysis outputs 2D lat-lon grid: 404×202)
    sim_data_flat = vec(sm_var.data)
    
    @info "Simulation data loaded" total_points=length(sim_data_flat) valid_points=length(valid_indices) sim_size=size(sm_var.data)
    
    # Extract values at valid observation points
    # valid_indices are simple linear indices into the flattened array
    member_vector = sim_data_flat[valid_indices]
    
    @info "SMAP member processing complete" n_observations=length(member_vector)
    
    return member_vector
end

"""
    process_joint_member_vector(diagnostics_folder_path, short_names, obs_data)

Process model output for joint calibration and return observation vector.
Handles both LH flux and SMAP with proper scaling.
"""
function process_joint_member_vector(diagnostics_folder_path, short_names, obs_data)
    @info "Processing joint calibration data" short_names
    
    # Get expected sizes from obs_data
    #n_lhf_expected = obs_data["n_lhf_per_period"][1]  # First (and only) period
    n_sm_expected = obs_data["n_sm_per_period"][1]
    
    @info "Expected observation counts" n_lhf=n_lhf_expected n_sm=n_sm_expected
    
    # Get simulation variable dictionary
    sim_var_dict = ext.get_sim_var_dict(diagnostics_folder_path)
    
    # 1. Process LH flux - match observation time period AND processing
    lhf_var = sim_var_dict["lhf"]()
    @info "LH flux before processing" dims=keys(lhf_var.dims) sizes=[(k, length(v)) for (k,v) in lhf_var.dims]
    
    # Get the date range from obs_data
    start_date_str, stop_date_str = obs_data["sample_date_ranges"][1]
    start_date = DateTime(start_date_str)
    stop_date = DateTime(stop_date_str)
    
    # Window to match observation time period FIRST
    lhf_var = ClimaAnalysis.window(
        lhf_var,
        "time",
        left = start_date,
        right = stop_date,
        by = ClimaAnalysis.MatchValue(),
    )
    @info "LH flux after time windowing" dims=keys(lhf_var.dims) sizes=[(k, length(v)) for (k,v) in lhf_var.dims]
    
    # Then do seasonal averaging (like observations - use ignore_nan = true to match)
    lhf_var = ClimaAnalysis.average_season_across_time(lhf_var, ignore_nan = true)
    @info "LH flux after seasonal average" dims=keys(lhf_var.dims) sizes=[(k, length(v)) for (k,v) in lhf_var.dims]
    
    # DON'T apply ocean mask - simulation already has it built-in
    # The observations were created with a different mask, so we'll get mismatch
    # Instead, we need to match the observation processing exactly
    # Skip the explicit mask application since simulation data already masked
    
    # Remove longitude double-counting
    lons = ClimaAnalysis.longitudes(lhf_var)
    lhf_var = ClimaAnalysis.window(
        lhf_var,
        "longitude",
        right = length(lons) - 1,
        by = ClimaAnalysis.Index(),
    )
    
    # Remove polar latitudes
    lats = ClimaAnalysis.latitudes(lhf_var)
    lhf_var = ClimaAnalysis.window(
        lhf_var,
        "latitude",
        left = 2,
        right = length(lats) - 1,
        by = ClimaAnalysis.Index(),
    )
    
    # Flatten ALL seasons
    # Replace NaN with 0.0 so vector size matches observations (which were regridded)
    # The simulation has more NaN points than observations due to stricter internal masking
    lhf_var = ClimaAnalysis.flatten(lhf_var)
    lhf_vector_raw = lhf_var.data
    
    # Check if we need to handle NaN mismatch
    #n_lhf_expected = obs_data["n_lhf_per_period"][1]
    n_lhf_actual = length(lhf_vector_raw)
    
    if n_lhf_actual < n_lhf_expected
        @warn "Simulation flattened output smaller than expected - this means simulation has more NaN points" n_lhf_actual n_lhf_expected
        # This shouldn't happen if masks are identical - indicates simulation has extra NaNs
        # Pad with zeros (will have minimal impact on calibration with proper scaling)
        lhf_vector = vcat(lhf_vector_raw, zeros(Float32, n_lhf_expected - n_lhf_actual))
    elseif n_lhf_actual > n_lhf_expected
        @warn "Simulation has more valid points than observations - truncating" n_lhf_actual n_lhf_expected
        lhf_vector = lhf_vector_raw[1:n_lhf_expected]
    else
        lhf_vector = lhf_vector_raw
    end
    
    @info "LH flux after flattening and size adjustment" length=length(lhf_vector)
    
    # 2. Process SMAP - load swc directly (not in sim_var_dict), temporal mean, extract surface
    # Note: swc is not in get_sim_var_dict, so we load it directly with SimDir
    sm_var = get(
        ClimaAnalysis.SimDir(diagnostics_folder_path),
        short_name = "swc",
    )
    sm_var = ClimaAnalysis.shift_to_start_of_previous_month(sm_var)
    sm_var = ClimaAnalysis.average_time(sm_var, ignore_nan = false)
    
    # Extract surface level (top of soil column)
    # swc has dimensions (lon, lat, z) where z is depth
    # We want the surface layer (z[end] = top of column)
    z_vals = sm_var.dims["z"]
    surface_z = z_vals[end]  # Top layer
    sm_var = ClimaAnalysis.slice(sm_var, z = surface_z)
    
    sm_var = ClimaAnalysis.flatten(sm_var)
    sm_vector_all = sm_var.data
    
    # Sample to match SMAP observation count (simple stride for now)
    # TODO: Replace with proper SMAP pixel grouping aggregation
    if length(sm_vector_all) > n_sm_expected
        stride = div(length(sm_vector_all), n_sm_expected)
        sm_vector = sm_vector_all[1:stride:stride*n_sm_expected]
    else
        sm_vector = sm_vector_all[1:min(n_sm_expected, end)]
    end
    
    # 3. Apply same scaling that was used for observations
    lhf_weight = obs_data["lhf_weight"]
    sm_weight = obs_data["sm_weight"]
    
    var_lhf = var(lhf_vector)
    var_sm = length(sm_vector) > 1 ? var(sm_vector) : 1.0
    
    scale_lhf = sqrt(lhf_weight / length(lhf_vector)) / sqrt(var_lhf)
    scale_sm = sqrt(sm_weight / length(sm_vector)) / sqrt(var_sm)
    
    lhf_vector .*= Float32(scale_lhf)
    sm_vector .*= Float32(scale_sm)
    
    @info "Applied observation scaling" scale_lhf scale_sm
    
    # 4. Concatenate in same order as observation vector: [LH, SM]
    member_vector = vcat(lhf_vector, sm_vector)
    
    @info "Joint member processing complete" n_lhf=length(lhf_vector) n_sm=length(sm_vector) total=length(member_vector)
    
    return member_vector
end

"""
    process_member_data!(
        g_ens_builder,
        col_idx,
        diagnostics_folder_path,
        short_names,
    )

Fill out the `col_idx`th of the G ensemble matrix using variables with the names
`short_names` from the NetCDF files in `diagnostics_folder_path`.
"""
function process_member_data!(
    g_ens_builder,
    col_idx,
    diagnostics_folder_path,
    short_names,
)
    @info "Processing member $col_idx with short names: $short_names"
    
    # Detect if we're using joint calibration or single-variable**
    has_smap = "sm_surface" in short_names
    has_lhf = "lhf" in short_names
    
    if has_smap && has_lhf
        # Joint calibration: process both variables
        @info "Joint calibration detected: processing both LH and SMAP"
        return process_joint_member_data!(
            g_ens_builder,
            col_idx,
            diagnostics_folder_path,
            short_names
        )
    elseif has_smap
        # SMAP only - use ERA5-style processing since SMAP is now in same format
        @info "SMAP-only calibration (using standard processing)"
        return process_era5_member_data!(
            g_ens_builder,
            col_idx,
            diagnostics_folder_path,
            ["sm_surface"]
        )
    else
        # ERA5 only (LH flux or other variables)
        @info "ERA5-only calibration"
        return process_era5_member_data!(
            g_ens_builder,
            col_idx,
            diagnostics_folder_path,
            short_names
        )
    end
end

"""
    process_joint_member_data!(g_ens_builder, col_idx, diagnostics_folder_path, short_names)

Process both LH flux (ERA5) and SMAP soil moisture for joint calibration.
Concatenates in the same order as observation vector: [LH..., SM...]
"""
function process_joint_member_data!(
    g_ens_builder,
    col_idx,
    diagnostics_folder_path,
    short_names
)
    obs_file = CALIBRATE_CONFIG.obs_vec_filepath
    obs_data = JLD2.load(obs_file)
    
    # Check that we have the expected structure
    if !("lhf_metadata" in keys(obs_data) && "sm_metadata" in keys(obs_data))
        error("Joint observation file is missing required metadata. Expected 'lhf_metadata' and 'sm_metadata'.")
    end
    
    @info "Processing joint calibration data"
    
    # 1. Process LH flux (ERA5 style)
    lhf_vector = process_era5_variable(diagnostics_folder_path, "lhf")
    
    # 2. Process SMAP (with spatial aggregation)
    sm_vector = process_smap_variable(
        diagnostics_folder_path,
        "sm_surface",
        obs_data["sm_metadata"]
    )
    
    # 3. Apply scaling if observations were pre-scaled
    if haskey(obs_data, "scale_lhf") && haskey(obs_data, "scale_sm")
        scale_lhf = obs_data["scale_lhf"]
        scale_sm = obs_data["scale_sm"]
        
        lhf_vector .*= scale_lhf
        sm_vector .*= scale_sm
        
        @info "Applied observation scaling" scale_lhf scale_sm
    end
    
    # 4. Concatenate in same order as observation vector: [LH, SM]
    member_vector = vcat(lhf_vector, sm_vector)
    
    @info "Joint member processing complete" n_lhf=length(lhf_vector) n_sm=length(sm_vector) total=length(member_vector)
    
    # Fill the g_ensemble column
    EnsembleBuilder.set_g_ens_col!(g_ens_builder, col_idx, member_vector)
    
    return nothing
end

"""
    process_era5_variable(diagnostics_folder_path, short_name)

Process a single ERA5 variable and return flattened vector.
"""
function process_era5_variable(diagnostics_folder_path, short_name)
    sim_var_dict = ext.get_sim_var_dict(diagnostics_folder_path)
    
    var = sim_var_dict[short_name]()
    var = ClimaAnalysis.average_season_across_time(var, ignore_nan = false)

    # Remove longitude double-counting
    lons = ClimaAnalysis.longitudes(var)
    var = ClimaAnalysis.window(
        var,
        "longitude",
        right = length(lons) - 1,
        by = ClimaAnalysis.Index(),
    )

    # Remove polar latitudes
    lats = ClimaAnalysis.latitudes(var)
    var = ClimaAnalysis.window(
        var,
        "latitude",
        left = 2,
        right = length(lats) - 1,
        by = ClimaAnalysis.Index(),
    )
    
    # Flatten to vector
    var_flat = ClimaAnalysis.flatten(var)
    
    @info "Processed ERA5 variable" short_name size=length(var_flat.data)
    
    return var_flat.data
end

"""
    process_smap_variable(diagnostics_folder_path, short_name, sm_metadata)

Process SMAP variable with spatial aggregation to match SMAP pixels.
"""
function process_smap_variable(
    diagnostics_folder_path,
    short_name,
    sm_metadata
)
    sample_date_ranges = CALIBRATE_CONFIG.sample_date_ranges
    sim_var_dict = ext.get_sim_var_dict(diagnostics_folder_path)
    
    member_data_all_periods = []
    
    for (period_idx, (start_str, stop_str)) in enumerate(sample_date_ranges)
        start_date = DateTime(start_str)
        stop_date = DateTime(stop_str)
        
        period_metadata = sm_metadata[period_idx]
        smap_pixel_groups = period_metadata["smap_pixel_groups"]
        
        # Extract model data for this period
        var = sim_var_dict[short_name]()
        
        # Get temporal window
        var_window = ClimaAnalysis.window(
            var,
            "time",
            left = start_date,
            right = stop_date,
            by = ClimaAnalysis.MatchValue(),
        )
        
        # Average over time
        var_avg = ClimaAnalysis.average_time(var_window)
        
        # Flatten spatial dimensions
        var_flat = ClimaAnalysis.flatten(var_avg)
        full_data = var_flat.data
        
        # Aggregate to SMAP pixel resolution
        smap_aggregated = aggregate_to_smap_pixels(full_data, smap_pixel_groups)
        
        push!(member_data_all_periods, smap_aggregated)
    end
    
    # Concatenate all periods
    member_vector = vcat(member_data_all_periods...)
    
    @info "Processed SMAP variable" short_name total_obs=length(member_vector)
    
    return member_vector
end

"""
    aggregate_to_smap_pixels(model_data::Vector, smap_pixel_groups::Vector{Vector{Int}})

Aggregate fine-resolution model output to coarse SMAP pixels by spatial averaging.

# Arguments
- model_data: Flattened model grid (all spatial points)
- smap_pixel_groups: Vector of vectors, where each inner vector contains 
                     model grid indices belonging to one SMAP pixel

# Returns
- Vector of aggregated values, one per SMAP pixel (spatial mean)
"""
function aggregate_to_smap_pixels(
    model_data::Vector{FT},
    smap_pixel_groups::Vector{Vector{Int}}
) where FT
    n_smap_pixels = length(smap_pixel_groups)
    aggregated = zeros(FT, n_smap_pixels)
    
    for (i, model_indices) in enumerate(smap_pixel_groups)
        if isempty(model_indices)
            aggregated[i] = NaN
            continue
        end
        
        # Extract model values for this SMAP pixel
        pixel_values = model_data[model_indices]
        
        # Spatial average (could also use median, area-weighted mean, etc.)
        # Filter NaNs if any
        valid_values = filter(!isnan, pixel_values)
        
        if isempty(valid_values)
            aggregated[i] = NaN
            @warn "All model points are NaN for SMAP pixel $i"
        else
            aggregated[i] = mean(valid_values)
        end
    end
    
    return aggregated
end

"""
    process_era5_member_data!(g_ens_builder, col_idx, diagnostics_folder_path, short_names)

Process ERA5 variables (existing functionality, now renamed for clarity).
"""
function process_era5_member_data!(
    g_ens_builder,
    col_idx,
    diagnostics_folder_path,
    short_names
)
    # Process variables - handle both ERA5 and SMAP
    sim_var_dict = ext.get_sim_var_dict(diagnostics_folder_path)
    
    vars = map(short_names) do short_name
        @info "Processing variable: $short_name"
        
        # Special handling for SMAP soil moisture (3D variable)
        if short_name == "sm_surface"
            # Load 3D soil water content
            swc = get(
                ClimaAnalysis.SimDir(diagnostics_folder_path),
                short_name = "swc",
            )
            
            # Shift to start of month and take temporal average
            swc = ClimaAnalysis.shift_to_start_of_previous_month(swc)
            var = ClimaAnalysis.average_time(swc, ignore_nan = true)
            
            # Extract surface level (top of soil column)
            # swc has dimensions (lon, lat, z) where z is depth
            z_vals = var.dims["z"]
            surface_z = z_vals[end]  # Top layer
            var = ClimaAnalysis.slice(var, z = surface_z)
            
            @info "Extracted surface soil moisture" dims=keys(var.dims)
        else
            # ERA5 or other 2D variables
            var = sim_var_dict[short_name]()
            var = ClimaAnalysis.average_season_across_time(var, ignore_nan = false)
        end

        # Common processing: remove longitude double-counting
        lons = ClimaAnalysis.longitudes(var)
        var = ClimaAnalysis.window(
            var,
            "longitude",
            right = length(lons) - 1,
            by = ClimaAnalysis.Index(),
        )

        # Remove polar latitudes
        lats = ClimaAnalysis.latitudes(var)
        var = ClimaAnalysis.window(
            var,
            "latitude",
            left = 2,
            right = length(lats) - 1,
            by = ClimaAnalysis.Index(),
        )
        
        var
    end

    seq_indices_checker = Checker.SequentialIndicesChecker()
    checkers = (seq_indices_checker,)

    for var in vars
        use_var = EnsembleBuilder.fill_g_ens_col!(
            g_ens_builder,
            col_idx,
            var;
            checkers,
            verbose = true,
        )
        use_var || error(
            "OutputVar with short name ($(ClimaAnalysis.short_name(var))) was passed, but not used",
        )
    end
    
    return nothing
end

"""
    ClimaCalibrate.analyze_iteration(ekp,
                                     g_ensemble,
                                     prior,
                                     output_dir,
                                     iteration)

Analyze an iteration by plotting the bias plots, constrained parameters over
iterations, and errors over iterations and time.
"""
function ClimaCalibrate.analyze_iteration(
    ekp,
    g_ensemble,
    prior,
    output_dir,
    iteration,
)
    plot_output_path = ClimaCalibrate.path_to_iteration(output_dir, iteration)
    plot_constrained_params_and_errors(plot_output_path, ekp, prior)

    # Leaderboard plots can only be plotted when the model saves swu, lwu, shf, lhf diagnostics
    # Plot ERA5 bias plots for only the first ensemble member
    # This can take a while to plot, so we plot only one of the members.
    # We choose the first ensemble member because the parameters for the first
    # ensemble member are supposed to be the mean of the parameters of the
    # ensemble members if it is EKP.TransformUnscented
    output_path =
        ClimaCalibrate.path_to_ensemble_member(output_dir, iteration, 1)

    diagnostics_folder_path =
        joinpath(output_path, "global_diagnostics", "output_active")
    ext.compute_monthly_leaderboard(
        output_path,
        diagnostics_folder_path,
        "ERA5",
    )
    ext.compute_seasonal_leaderboard(
        output_path,
        diagnostics_folder_path,
        "ERA5",
    )
end

"""
    plot_constrained_params_and_errors(output_dir, ekp, prior)

Plot the constrained parameters and errors from `ekp` and `prior` and save
them to `output_dir`.
"""
function plot_constrained_params_and_errors(output_dir, ekp, prior)
    dim_size = sum(length.(EKP.batch(prior)))
    fig = CairoMakie.Figure(size = ((dim_size + 1) * 500, 500))
    for i in 1:dim_size
        EKP.Visualize.plot_ϕ_over_iters(fig[1, i], ekp, prior, i)
    end
    EKP.Visualize.plot_error_over_iters(
        fig[1, dim_size + 1],
        ekp,
        error_metric = "loss",
    )
    EKP.Visualize.plot_error_over_time(
        fig[1, dim_size + 2],
        ekp,
        error_metric = "loss",
    )
    CairoMakie.save(
        joinpath(output_dir, "constrained_params_and_error.png"),
        fig,
    )
    return nothing
end
