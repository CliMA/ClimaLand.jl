import ClimaAnalysis
import Dates
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

    # Determine which observation is used by the short names
    # This assumes that observations do not differ in the variables that are
    # being calibrated
    obs_series = EKP.get_observation_series(ekp)
    short_names = ClimaCalibrate.ObservationRecipe.short_names(
        first(obs_series.observations),
    )
    
    # **NEW: Load SMAP observation metadata if using SMAP**
    obs_metadata = nothing
    if "sm_surface" in short_names
        obs_file = CALIBRATE_CONFIG.obs_vec_filepath
        if isfile(obs_file)
            obs_metadata = JLD2.load(obs_file)
            @info "Loaded SMAP observation metadata" valid_obs=length(obs_metadata["valid_indices"])
        else
            error("SMAP observation file not found: $obs_file")
        end
    end

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
    
    # **NEW: Detect if we're using joint calibration or single-variable**
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
        # SMAP only
        @info "SMAP-only calibration"
        obs_file = CALIBRATE_CONFIG.obs_vec_filepath
        obs_data = JLD2.load(obs_file)
        return process_smap_member_data!(
            g_ens_builder,
            col_idx,
            diagnostics_folder_path,
            ["sm_surface"],  # Only SM
            obs_data
        )
    else
        # ERA5 only (LH flux)
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
    # **EXISTING CODE** - just renamed from process_member_data
    era5_obs_vars = ext.get_era5_obs_var_dict()
    
    for short_name in short_names
        short_name in keys(era5_obs_vars) || error(
            "Variable $short_name does not appear in the observation dataset."
        )
    end

    sim_var_dict = ext.get_sim_var_dict(diagnostics_folder_path)
    vars = map(short_names) do short_name
        var = sim_var_dict[short_name]()
        var = ClimaAnalysis.average_season_across_time(var, ignore_nan = false)

        # To prevent double counting along the longitudes since -180 and 180
        # degrees are the same point
        lons = ClimaAnalysis.longitudes(var)
        var = ClimaAnalysis.window(
            var,
            "longitude",
            right = length(lons) - 1,
            by = ClimaAnalysis.Index(),
        )

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
