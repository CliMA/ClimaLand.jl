import ClimaAnalysis
import Dates
import ClimaCalibrate

include("observation_utils.jl")

function ClimaCalibrate.observation_map(iteration)
    ekp = JLD2.load_object(ClimaCalibrate.ekp_path(output_dir, iteration))
    current_minibatch = EKP.get_current_minibatch(ekp)
    obs = EKP.get_obs(ekp)
    single_obs_len = sum(length(obs))
    single_member_len = single_obs_len * length(current_minibatch)
    ensemble_size = EKP.get_N_ens(ekp)
    # Name of observation is determined by short names of OutputVars
    # TODO: I am not sure how I feel about using short names like this
    short_names = split(obs_series.observations[1].names[1], ";")

    G_ensemble = Array{Float64}(undef, single_member_len, ensemble_size)
    for m = 1:ensemble_size
        member_path = ClimaCalibrate.path_to_ensemble_member(output_dir, iteration, m)
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


function process_member_data(diagnostics_folder_path, short_names, current_minibatch)
    for short_name in short_names
        short_name in ("lhf", "shf", "swu", "lwu") || error(
            "Do not know how to process member data with variable with the short name of $short_name",
        )
    end

    sim_var_dict = get_sim_var_dict(diagnostics_folder_path)
    # TODO: This is annoying to keep track of how observational and simulational data is
    # processed differently even though they are mostly similar
    # Process each variable to be seasonal averages across time
    vars = map(short_names) do short_name
        sim_var = sim_var_dict[short_name]()

        # Throw away first SPINUP months
        first_date = first(ClimaAnalysis.dates(var))
        var = ClimaAnalysis.window(
            var,
            left = first_date + SPINUP,
            by = ClimaAnalysis.MatchValue(),
        )

        var = ClimaAnalysis.average_season_across_time(var, ignore_nan = true)
        ocean_mask = make_ocean_mask(NELEMENTS)
        var = ocean_mask(var)
        var = window(var, "longitude", right = length(lons) - 1, by = ClimaAnalysis.Index())
    end

    # TODO: Verify that current_minibatch is what I think it is
    # Flatten and concatenate the data for each minibatch
    flattened_data = map(current_minibatch) do idx
        start_date, end_date = SAMPLE_DATE_RANGES[idx]
        map(vars) do var
            var = ClimaAnalysis.window(
                var,
                "time",
                left = start_date,
                right = end_date,
                by = ClimaAnalysis.MatchValue(),
            )
            ClimaAnalysis.flatten(var).data
        end
    end
    return vcat(vcat(flattened_data...)...)
end

function ClimaCalibrate.analyze_iteration(ekp, g_ensemble, prior, output_dir, iteration)
    plot_output_path = ClimaCalibrate.path_to_iteration(output_dir, iteration)
    plot_constrained_params_and_errors(plot_output_path, ekp, prior)

    # Only plot using the first ensemble member
    output_path = ClimaCalibrate.path_to_ensemble_member(output_dir, iteration, 1)

    # TODO: This need to change and I don't know to what :(
    diagnostics_folder_path = joinpath(output_path, "model_config")
    try
        compute_leaderboard(output_path, diagnostics_folder_path, spinup_months)
    catch e
        @error "Error in `analyze_iteration`" error = e
    end
end

function plot_constrained_params_and_errors(output_dir, ekp, prior)
    dim_size = sum(length.(EKP.batch(prior)))
    fig = CairoMakie.Figure(size = ((dim_size + 1) * 500, 500))
    for i = 1:dim_size
        EKP.Visualize.plot_ϕ_over_iters(fig[1, i], ekp, prior, i)
    end
    EKP.Visualize.plot_error_over_iters(fig[1, dim_size+1], ekp)
    EKP.Visualize.plot_error_over_time(fig[1, dim_size+2], ekp)
    CairoMakie.save(joinpath(output_dir, "constrained_params_and_error.png"), fig)
    return nothing
end
