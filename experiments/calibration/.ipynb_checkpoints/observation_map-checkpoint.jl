import ClimaAnalysis
import Dates
import ClimaCalibrate

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
    current_minibatch = EKP.get_current_minibatch(ekp)
    obs = EKP.get_obs(ekp)
    single_obs_len = sum(length(obs))
    ensemble_size = EKP.get_N_ens(ekp)

    obs_series = EKP.get_observation_series(ekp)
    # Determine which observation is used by the short names
    short_names = ClimaCalibrate.ObservationRecipe.short_names(
        first(obs_series.observations),
    )

    G_ensemble = Array{Float64}(undef, single_obs_len, ensemble_size)
    for m in 1:ensemble_size
        member_path =
            ClimaCalibrate.path_to_ensemble_member(output_dir, iteration, m)
        simdir_path = joinpath(member_path, "global_diagnostics/output_active")
        @info "Processing member $m: $simdir_path"
        try
            G_ensemble[:, m] .=
                process_member_data(simdir_path, short_names, current_minibatch)

        catch e
            @error "Error processing member $m, filling observation map entry with NaNs" exception =
                e
            G_ensemble[:, m] .= NaN
        end
    end

    return G_ensemble
end

"""
    process_member_data(diagnostics_folder_path, short_names, current_minibatch)

Process the data of a single ensemble member and return a single column of the
G ensemble matrix.
"""
function process_member_data(
    diagnostics_folder_path,
    short_names,
    current_minibatch,
)
    sample_date_ranges = CALIBRATE_CONFIG.sample_date_ranges
    nelements = CALIBRATE_CONFIG.nelements
    @info "Short names: $short_names"
    era5_obs_vars = ext.get_era5_obs_var_dict()
    # for non-era5 data uncomment and integrate:
    # obs_dict =
    # if all(in(keys(ext.get_era5_obs_var_dict())), short_names)
    #     ext.get_era5_obs_var_dict()
    # elseif Base.hasproperty(ext, :get_smap_obs_var_dict) &&
    #        all(in(keys(ext.get_smap_obs_var_dict())), short_names)
    #     ext.get_smap_obs_var_dict()
    # else
    #     error("Short names not found in ERA5 or SMAP observation dicts")
    # end
    for short_name in short_names
        short_name in keys(era5_obs_vars) || error(
            "Variable $short_name does not appear in the observation dataset. Add the variable to get_era5_obs_var_dict",
        )
    end

    sim_var_dict = ext.get_sim_var_dict(diagnostics_folder_path)
    vars = map(short_names) do short_name
        var = sim_var_dict[short_name]()
        var = ClimaAnalysis.average_season_across_time(var, ignore_nan = false)

        ocean_mask = make_ocean_mask(nelements)
        var = ocean_mask(var)

        # Replace all NaNs on land with the mean value
        # This is needed to stop small NaN values from stopping the calibration
        nanmean_land_val = mean(!isnan, var.data)
        var = ClimaAnalysis.replace(
            val -> isnan(val) ? nanmean_land_val : val,
            var,
        )
        var = ocean_mask(var)

        # To prevent double counting along the longitudes since -180 and 180
        # degrees are the same point
        lons = ClimaAnalysis.longitudes(var)
        var = ClimaAnalysis.window(
            var,
            "longitude",
            right = length(lons) - 1,
            by = ClimaAnalysis.Index(),
        )

        # Exclude the poles
        lats = ClimaAnalysis.latitudes(var)
        var = ClimaAnalysis.window(
            var,
            "latitude",
            left = 2,
            right = length(lats) - 1,
            by = ClimaAnalysis.Index(),
        )
    end
    # Flatten and concatenate the data for each minibatch
    # Note that we implicitly remove spinup when windowing is done
    @info "Current minibatch is $current_minibatch"
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
            ClimaAnalysis.flatten(var).data
        end
        flat_data
    end
    member = vcat(vcat(flattened_data...)...)
    @info "Size of member is $(length(member))"
    return member
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
