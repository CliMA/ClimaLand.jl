import ClimaAnalysis
import ClimaUtilities.ClimaArtifacts: @clima_artifact
import GeoMakie
import CairoMakie
import Dates

include("data_sources.jl")

"""
    compute_leaderboard(leaderboard_base_path, diagnostics_folder_path)

Plot the biases against observations and simulation data.

The paramter `leaderboard_base_path` is the path to save the leaderboards and bias plots and
`diagnostics_folder_path` is the path to the simulation data.
"""
function compute_leaderboard(leaderboard_base_path, diagnostics_folder_path)
    @info "Error against observations"

    # Set up dict for storing simulation and observational data after processing
    sim_obs_comparsion_dict = Dict()

    # Print dates for debugging
    _, var_func = first(sim_var_dict)
    var = var_func()
    output_dates =
        Dates.DateTime(var.attributes["start_date"]) .+
        Dates.Second.(ClimaAnalysis.times(var))
    @info "Working with dates:"
    @info output_dates

    short_names = keys(sim_var_dict)
    for short_name in short_names
        @info short_name
        # Simulation data
        sim_var = sim_var_dict[short_name]()

        # Observational data
        obs_var = obs_var_dict[short_name](sim_var.attributes["start_date"])

        # Remove first spin_up_months from simulation
        spin_up_months = 12
        spinup_cutoff = spin_up_months * 31 * 86400.0
        ClimaAnalysis.times(sim_var)[end] >= spinup_cutoff &&
            (sim_var = ClimaAnalysis.window(sim_var, "time", left = spinup_cutoff))

        obs_var = ClimaAnalysis.resampled_as(obs_var, sim_var)

        sim_obs_comparsion_dict[short_name] = (sim_var, obs_var)
    end

    months = Dates.monthname.(1:12)

    # Plot monthly comparsions
    for short_name in short_names
        fig = CairoMakie.Figure(size = (650 * 6, 450 * 2))
        sim_var, obs_var = sim_obs_comparsion_dict[short_name]
        mask_vis, mask_fn = mask_dict[short_name](sim_var, obs_var)
        times = ClimaAnalysis.times(sim_var) |> copy
        times = reshape(times, (2, 6))
        for ((indices, t), month) in zip(pairs(times), months)
            layout = fig[Tuple(indices)...] = CairoMakie.GridLayout()
            ClimaAnalysis.Visualize.plot_bias_on_globe!(
                layout,
                ClimaAnalysis.slice(sim_var, time = t),
                ClimaAnalysis.slice(obs_var, time = t),
                mask = mask_vis,
                cmap_extrema = compare_vars_biases_plot_extrema[short_name],
                apply_mask_to_var = mask_fn,
                more_kwargs = Dict(:mask => ClimaAnalysis.Utils.kwargs(color = :white)),
            )
            CairoMakie.Label(layout[0, 1], month, tellwidth = false, fontsize = 30)
        end
        CairoMakie.save(joinpath(leaderboard_base_path, "$short_name.png"), fig)
    end

    # Plot month (x-axis) and global bias and global RMSE (y-axis)
    fig = CairoMakie.Figure(size = (450 * length(keys(sim_var_dict)), 600))
    for (col, short_name) in enumerate(keys(sim_var_dict))
        sim_var, obs_var = sim_obs_comparsion_dict[short_name]
        _, mask_fn = mask_dict[short_name](sim_var, obs_var)
        # Compute globas bias and global rmse
        rmse_vec = []
        bias_vec = []
        times = ClimaAnalysis.times(sim_var) |> copy
        for t in times
            g_rmse = ClimaAnalysis.global_rmse(
                ClimaAnalysis.slice(sim_var, time = t),
                ClimaAnalysis.slice(obs_var, time = t),
                mask = mask_fn,
            )
            g_bias = ClimaAnalysis.global_bias(
                ClimaAnalysis.slice(sim_var, time = t),
                ClimaAnalysis.slice(obs_var, time = t),
                mask = mask_fn,
            )
            push!(rmse_vec, g_rmse)
            push!(bias_vec, g_bias)
        end
        ax_rmse = CairoMakie.Axis(
            fig[1, col],
            title = "Global RMSE for $short_name",
            xlabel = "Month",
            ylabel = "Global RMSE ($(ClimaAnalysis.units(sim_var)))",
            xticks = (1:12, Dates.monthname.(1:12) .|> x -> x[1:3]),
        )
        CairoMakie.lines!(ax_rmse, 1:12 |> collect, rmse_vec)
        CairoMakie.scatter!(ax_rmse, 1:12 |> collect, rmse_vec)

        ax_bias = CairoMakie.Axis(
            fig[2, col],
            title = "Global Bias for $short_name",
            xlabel = "Month",
            ylabel = "Global Bias ($(ClimaAnalysis.units(sim_var)))",
            xticks = (1:12, Dates.monthname.(1:12) .|> x -> x[1:3]),
        )
        CairoMakie.lines!(ax_bias, 1:12 |> collect, bias_vec)
        CairoMakie.scatter!(ax_bias, 1:12 |> collect, bias_vec)
    end
    CairoMakie.save(joinpath(leaderboard_base_path, "global_rmse_and_bias_graphs.png"), fig)
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) < 2
        error(
            "Usage: julia leaderboard.jl <Filepath to save leaderboard and bias plots> <Filepath to simulation data>",
        )
    end
    leaderboard_base_path = ARGS[begin]
    diagnostics_folder_path = ARGS[begin+1]
    compute_leaderboard(leaderboard_base_path, diagnostics_folder_path)
end
