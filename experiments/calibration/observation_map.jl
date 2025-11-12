import ClimaAnalysis
import Dates
import ClimaCalibrate
import ClimaCalibrate.EnsembleBuilder
import ClimaCalibrate.Checker

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
    nelements = CALIBRATE_CONFIG.nelements
    
    @info "Short names: $short_names"
    
    # **NEW: Handle SMAP vs ERA5**
    if "sm_surface" in short_names && !isnothing(obs_metadata)
        @info "Processing SMAP soil moisture observations"
        return process_smap_member_data(
            diagnostics_folder_path,
            short_names,
            current_minibatch,
            obs_metadata
        )
    else
        # Existing ERA5 processing
        return process_era5_member_data(
            diagnostics_folder_path,
            short_names,
            current_minibatch,
            nelements
        )
    end
end

"""
    process_smap_member_data(diagnostics_folder_path, short_names, current_minibatch, obs_metadata)

Process model soil moisture and subset to match SMAP observation locations.
"""
function process_smap_member_data(
    diagnostics_folder_path,
    short_names,
    current_minibatch,
    obs_metadata
)
    valid_indices = obs_metadata["valid_indices"]
    sample_date_ranges = CALIBRATE_CONFIG.sample_date_ranges
    nelements = CALIBRATE_CONFIG.nelements
    
    sim_var_dict = ext.get_sim_var_dict(diagnostics_folder_path)
    
    vars = map(short_names) do short_name
        var = sim_var_dict[short_name]()
        var = ClimaAnalysis.average_season_across_time(var, ignore_nan = false)
        
        # Apply ocean mask
        ocean_mask = make_ocean_mask(nelements)
        var = ocean_mask(var)
        
        var
    end
    
    # Flatten and subset to SMAP observation locations
    flattened_data = map(current_minibatch) do idx
        start_date, stop_date = sample_date_ranges[idx]
        flat_data = map(vars) do var
            var = ClimaAnalysis.window(
                var,
                "time",
                left = start_date,
                right = stop_date,
                by = ClimaAnalysis.MatchValue(),
            )
            
            # Flatten the entire grid
            full_flat = ClimaAnalysis.flatten(var).data
            
            # **CRITICAL: Subset to match SMAP observation locations**
            smap_subset = full_flat[valid_indices]
            
            @info "Subsetted model output to SMAP locations" full_size=length(full_flat) subset_size=length(smap_subset)
            
            smap_subset
        end
        flat_data
    end
    
    member = vcat(vcat(flattened_data...)...)
    @info "Size of SMAP-matched member data: $(length(member))"
    
    return member
end

"""
    process_era5_member_data(diagnostics_folder_path, short_names, current_minibatch, nelements)

Process ERA5 variables (existing functionality).
"""
function process_era5_member_data(
    diagnostics_folder_path,
    short_names,
    current_minibatch,
    nelements
)
    # **EXISTING CODE** from your original process_member_data
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

    # This check should be used, because fill_g_ens_col! is not aware of the
    # meaning of the time dimension (e.g. seasonal averages vs monthly
    # averages). For example, without this check, if the simulation data contain
    # monthly averages and metadata track seasonal averages, then no error is
    # thrown, because all dates in metadata are in all the dates in var.
    seq_indices_checker = Checker.SequentialIndicesChecker()
    checkers = (seq_indices_checker,)

    # fill_g_ens_col! will remove the spinup and is mask-aware.
    # g_ens_builder contain the metadata from the observations, so
    # fill_g_ens_col! will only choose values over temporal and spatial
    # coordinates that exist in the observational data
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
