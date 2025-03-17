import ClimaAnalysis
import GeoMakie
import CairoMakie
import Dates

include("data_sources.jl")

"""
    compute_leaderboard(leaderboard_base_path, diagnostics_folder_path)

Plot the biases and a leaderboard of various variables defined over longitude, latitude, and
time.

The parameter `leaderboard_base_path` is the path to save the leaderboards and bias plots,
and `diagnostics_folder_path` is the path to the simulation data.

Loading and preprocessing simulation data is done by `get_sim_var_dict`. Loading and
preprocessing observational data is done by `get_obs_var_dict`. The masks are normalizing
the global RMSE and bias are determined by `get_mask_dict`. The ranges of the bias plots are
determined by `get_compare_vars_biases_plot_extrema`. See the functions defined in
data_sources.jl.
"""
function compute_leaderboard(leaderboard_base_path, diagnostics_folder_path)
    sim_var_dict = get_sim_var_dict(diagnostics_folder_path)
    obs_var_dict = get_obs_var_dict()
    mask_dict = get_mask_dict()
    compare_vars_biases_plot_extrema = get_compare_vars_biases_plot_extrema()

    @info "Error against observations"

    # Set up dict for storing simulation and observational data after processing
    # and for storing the month we are interested in
    sim_obs_comparsion_dict = Dict()
    times_dict = Dict() # map short name to (times, dates)

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

        # Remove first spin_up_months months from simulation
        spin_up_months = 12
        spinup_cutoff = spin_up_months * 30 * 86400.0
        ClimaAnalysis.times(sim_var)[end] >= spinup_cutoff && (
            sim_var =
                ClimaAnalysis.window(sim_var, "time", left = spinup_cutoff)
        )

        # Determine which 12 months to compute RMSE on
        # For a very long run, there is a possibility of the final year of the observational
        # data to be earlier than the final year of the simulation data
        sim_times = ClimaAnalysis.times(sim_var)
        obs_times = ClimaAnalysis.times(obs_var)
        sim_time_idx = max(length(sim_times) - 11, 1)
        obs_time_idx = max(length(obs_times) - 11, 1)
        if obs_times[obs_time_idx] < sim_times[sim_time_idx]
            left_time = obs_times[obs_time_idx]
            right_time = obs_times[end]
            times = obs_times[obs_time_idx:end]
            dates =
                Dates.DateTime(obs_var.attributes["start_date"]) .+
                Dates.Second.(times)
        else
            left_time = sim_times[sim_time_idx]
            right_time = sim_times[end]
            times = sim_times[sim_time_idx:end]
            dates =
                Dates.DateTime(sim_var.attributes["start_date"]) .+
                Dates.Second.(times)
        end

        # Get 12 or less months of data
        sim_var = ClimaAnalysis.window(
            sim_var,
            "time",
            left = times[begin],
            right = times[end],
        )
        obs_var = ClimaAnalysis.window(
            obs_var,
            "time",
            left = times[begin],
            right = times[end],
        )

        obs_var = ClimaAnalysis.resampled_as(obs_var, sim_var)

        sim_obs_comparsion_dict[short_name] = (sim_var, obs_var)
        times_dict[short_name] = (times, dates)
    end

    # Plot monthly comparsions
    for short_name in short_names
        # Get the right number of months for plotting
        sim_var, obs_var = sim_obs_comparsion_dict[short_name]
        times, monthly_dates = times_dict[short_name]
        num_times = length(monthly_dates)
        months_and_years = (
            (Dates.monthname(date), Dates.year(date)) for date in monthly_dates
        )

        fig = CairoMakie.Figure(size = (650 * ceil(num_times / 2), 450 * 2))
        mask = mask_dict[short_name](sim_var, obs_var)
        times = vcat(
            times,
            Array{Union{Missing, eltype(times)}}(missing, 12 - num_times),
        )
        times = reshape(times, (2, 6))
        for ((indices, t), (month, year)) in zip(pairs(times), months_and_years)
            layout = fig[Tuple(indices)...] = CairoMakie.GridLayout()
            ClimaAnalysis.Visualize.plot_bias_on_globe!(
                layout,
                ClimaAnalysis.slice(sim_var, time = t),
                ClimaAnalysis.slice(obs_var, time = t),
                cmap_extrema = compare_vars_biases_plot_extrema[short_name],
                mask = mask,
            )
            CairoMakie.Label(
                layout[0, 1],
                month * " $year",
                tellwidth = false,
                fontsize = 30,
            )
        end
        CairoMakie.save(
            joinpath(leaderboard_base_path, "$(short_name)_bias_plot.png"),
            fig,
        )
    end

    # Plot month (x-axis) and global bias and global RMSE (y-axis)
    fig = CairoMakie.Figure(size = (450 * length(keys(sim_var_dict)), 600))
    for (col, short_name) in enumerate(keys(sim_var_dict))
        sim_var, obs_var = sim_obs_comparsion_dict[short_name]
        times, monthly_dates = times_dict[short_name]
        years = unique(Dates.year.(monthly_dates))
        num_times = length(times)
        mask = mask_dict[short_name](sim_var, obs_var)
        # Compute globas bias and global rmse
        rmse_vec = []
        bias_vec = []
        for t in times
            g_rmse = ClimaAnalysis.global_rmse(
                ClimaAnalysis.slice(sim_var, time = t),
                ClimaAnalysis.slice(obs_var, time = t),
                mask = mask,
            )
            g_bias = ClimaAnalysis.global_bias(
                ClimaAnalysis.slice(sim_var, time = t),
                ClimaAnalysis.slice(obs_var, time = t),
                mask = mask,
            )
            push!(rmse_vec, g_rmse)
            push!(bias_vec, g_bias)
        end
        ax_rmse = CairoMakie.Axis(
            fig[1, col],
            title = "Global RMSE for $short_name",
            xlabel = "Month and Year: $years",
            ylabel = "Global RMSE ($(ClimaAnalysis.units(sim_var)))",
            xticks = (
                1:num_times,
                Dates.monthname.(monthly_dates) .|> x -> x[1:3],
            ),
        )
        CairoMakie.lines!(ax_rmse, 1:num_times |> collect, rmse_vec)
        CairoMakie.scatter!(ax_rmse, 1:num_times |> collect, rmse_vec)

        ax_bias = CairoMakie.Axis(
            fig[2, col],
            title = "Global Bias for $short_name",
            xlabel = "Month and Year: $years",
            ylabel = "Global Bias ($(ClimaAnalysis.units(sim_var)))",
            xticks = (
                1:num_times,
                Dates.monthname.(monthly_dates) .|> x -> x[1:3],
            ),
        )
        CairoMakie.lines!(ax_bias, 1:num_times |> collect, bias_vec)
        CairoMakie.scatter!(ax_bias, 1:num_times |> collect, bias_vec)
    end
    CairoMakie.save(
        joinpath(leaderboard_base_path, "global_rmse_and_bias_graphs.png"),
        fig,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) != 2
        error(
            "Usage: julia leaderboard.jl <Filepath to save leaderboard and bias plots> <Filepath to simulation data>",
        )
    end
    leaderboard_base_path = ARGS[begin]
    diagnostics_folder_path = ARGS[begin + 1]
    compute_leaderboard(leaderboard_base_path, diagnostics_folder_path)
end
