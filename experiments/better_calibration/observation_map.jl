import ClimaAnalysis
import Dates
import ClimaCalibrate

include(
    joinpath(
        pkgdir(ClimaLand),
        "experiments/better_calibration/observation_utils.jl",
    ),
)
include(
    joinpath(
        pkgdir(ClimaLand),
        "ext/land_sim_vis/leaderboard/data_sources.jl",
    ),
)
include(
    joinpath(
        pkgdir(ClimaLand),
        "ext/land_sim_vis/leaderboard/leaderboard.jl",
    ),
)

"""
    ClimaCalibrate.observation_map(iteration)

Return G ensemble for an `iteration`.
"""
function ClimaCalibrate.observation_map(iteration)
    (; output_dir) = get_config()
    ekp = JLD2.load_object(ClimaCalibrate.ekp_path(output_dir, iteration))
    current_minibatch = EKP.get_current_minibatch(ekp)
    obs = EKP.get_obs(ekp)
    single_obs_len = sum(length(obs))
    @info "The length of the member is $single_obs_len"
    # @info single_member_len
    ensemble_size = EKP.get_N_ens(ekp)
    # Name of observation is determined by short names of OutputVars
    # TODO: I am not sure how I feel about using short names like this
    # TODO: Should add a function to be able to extract short names like this
    obs_series = EKP.get_observation_series(ekp)
    short_names = split(obs_series.observations[1].names[1], ";")

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

# TODO: It might make sense to add sample_date_ranges and nelements as
# parameters to this function
"""
    process_member_data(diagnostics_folder_path, short_names, current_minibatch)

Process the data of a single ensemble member.
"""
function process_member_data(
    diagnostics_folder_path,
    short_names,
    current_minibatch,
)
    (; sample_date_ranges, nelements) = get_config()
    @info "Short names: $short_names"
    for short_name in short_names
        short_name in ("lhf", "shf", "swu", "lwu") || error(
            "Do not know how to process member data with variable with the short name $short_name",
        )
    end

    sim_var_dict = get_sim_var_dict(diagnostics_folder_path)
    # TODO: This is annoying to keep track of how observational and simulational
    # data is processed differently even though they are mostly similar
    # Process each variable to be seasonal averages across time
    # TODO: Also, do not like how the preprocessing is done in completely
    # different files which make it difficult to follow due to code locality
    vars = map(short_names) do short_name
        var = sim_var_dict[short_name]()
        var =
            ClimaAnalysis.average_season_across_time(var, ignore_nan = false)

        ocean_mask = make_ocean_mask(nelements)
        var = ocean_mask(var)

        # Replace all NaNs on land with the mean value
        # This is needed to stop small NaN values from stopping the calibration
        # TODO: This is a stopgap. It might be better to leave them as NaNs. As a temporary
        # solution, we replace with 500.0
        nanmean_land_val = 500.0 # mean(filter(!isnan, var.data))
        var = ClimaAnalysis.replace(val -> isnan(val) ? nanmean_land_val : val, var)
        var = ocean_mask(var)

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
    end
    # Flatten and concatenate the data for each minibatch
    # Note that we implicitly remove spinup when windowing is done
    @info "Current minibatch is $current_minibatch"
    flattened_data = map(current_minibatch) do idx
        start_date, end_date = sample_date_ranges[idx]
        flat_data = map(vars) do var
            var = ClimaAnalysis.window(
                var,
                "time",
                left = start_date,
                right = end_date,
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

    # Only plot using the first ensemble member
    output_path =
        ClimaCalibrate.path_to_ensemble_member(output_dir, iteration, 1)

    # Plot ERA5 bias plots for only the first ensemble member
    # This can take a while to plot, so we plot only one of the members.
    # We choose the first ensemble member because the parameters are supposed to
    # be the mean of the parameters of the ensemble members if it is
    # EKP.TransformUnscented
    diagnostics_folder_path = joinpath(output_path, "global_diagnostics")
    try
        compute_monthly_leaderboard(
            output_path,
            diagnostics_folder_path,
            "ERA5",
        )
        compute_seasonal_leaderboard(
            output_path,
            diagnostics_folder_path,
            "ERA5",
        )
    catch e
        @error "Error in `analyze_iteration`" error = e
    end
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
    EKP.Visualize.plot_error_over_iters(fig[1, dim_size + 1], ekp)
    EKP.Visualize.plot_error_over_time(fig[1, dim_size + 2], ekp)
    CairoMakie.save(
        joinpath(output_dir, "constrained_params_and_error.png"),
        fig,
    )
    return nothing
end
