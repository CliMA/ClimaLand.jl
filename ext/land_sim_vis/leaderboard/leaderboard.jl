"""
    compute_monthly_leaderboard(leaderboard_base_path,
                                diagnostics_folder_path,
                                data_source)

Plot the biases and a monthly leaderboard of various variables defined over longitude,
latitude, and time. The observational data is determined by `data_source` and can either be
`ILAMB` or `ERA5`.

The parameter `leaderboard_base_path` is the path to save the leaderboards and bias plots,
and `diagnostics_folder_path` is the path to the simulation data.

Loading and preprocessing simulation data is done by `get_sim_var_dict`. Loading and
preprocessing observational data is done by `get_obs_var_dict`. The masks are normalizing
the global RMSE and bias are determined by `get_mask_dict`. The ranges of the bias plots are
determined by `get_compare_vars_biases_plot_extrema`. See the functions defined in
data_sources.jl.
"""
function compute_monthly_leaderboard(
    leaderboard_base_path,
    diagnostics_folder_path,
    data_source,
)
    sim_var_dict = get_sim_var_dict(diagnostics_folder_path)
    obs_var_dict = get_obs_var_dict(data_source)
    mask_dict = get_mask_dict(data_source)

    compare_vars_biases_plot_extrema = get_compare_vars_biases_plot_extrema()
    short_names = intersect(keys(sim_var_dict), keys(obs_var_dict))
    issubset(short_names, keys(mask_dict)) ||
        error("Not all variables ($short_names) have a mask $(keys(mask_dict))")

    @info "Error against observations"

    # Set up dict for storing simulation and observational data after processing
    # and for storing the month we are interested in
    sim_obs_comparsion_dict = Dict()

    # Print dates for debugging
    _, var_func = first(sim_var_dict)
    var = var_func()
    output_dates =
        Dates.DateTime(var.attributes["start_date"]) .+
        Dates.Second.(ClimaAnalysis.times(var))
    @info "Working with dates:"
    @info output_dates

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

        # Get the first valid time and last valid time
        obs_times = ClimaAnalysis.times(obs_var)
        sim_times = ClimaAnalysis.times(sim_var)
        min_time = maximum(first.((obs_times, sim_times)))
        max_time = minimum(last.((obs_times, sim_times)))

        # Window OutputVars to restrain the times to those that are the same between
        # both OutputVars
        sim_var = ClimaAnalysis.window(
            sim_var,
            "time",
            left = min_time,
            right = max_time,
        )
        obs_var = ClimaAnalysis.window(
            obs_var,
            "time",
            left = min_time,
            right = max_time,
        )

        obs_var = ClimaAnalysis.shift_longitude(obs_var, -180.0, 180.0)
        obs_var = ClimaAnalysis.resampled_as(obs_var, sim_var)

        sim_obs_comparsion_dict[short_name] = (sim_var, obs_var)
    end

    # Plot monthly comparsions
    for short_name in short_names
        sim_var, obs_var = sim_obs_comparsion_dict[short_name]

        # Grab the last 12 months if possible for plotting
        times = ClimaAnalysis.times(sim_var)
        min_idx = max(1, length(times) - 11)
        times = times[min_idx:end]
        monthly_dates =
            Dates.DateTime(sim_var.attributes["start_date"]) .+
            Dates.Second.(times)
        num_times = length(times)
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
            joinpath(
                leaderboard_base_path,
                "$(data_source)_$(short_name)_bias_plot.png",
            ),
            fig,
        )
    end

    # Plot month (x-axis) and global bias and global RMSE (y-axis)
    fig = CairoMakie.Figure(size = (250 + 450 * length(short_names), 900))
    fig_rmse_bias = fig[1, 1] = CairoMakie.GridLayout()
    for (col, short_name) in enumerate(short_names)
        sim_var, obs_var = sim_obs_comparsion_dict[short_name]
        times = ClimaAnalysis.times(sim_var)
        mask = mask_dict[short_name](sim_var, obs_var)

        sim_vec = [
            ClimaAnalysis.weighted_average_lonlat(
                ClimaAnalysis.apply_oceanmask(
                    ClimaAnalysis.slice(sim_var, time = t),
                ),
            ).data[] for t in times
        ]
        rmse_vec = [
            ClimaAnalysis.global_rmse(
                ClimaAnalysis.slice(sim_var, time = t),
                ClimaAnalysis.slice(obs_var, time = t),
                mask = mask,
            ) for t in times
        ]
        bias_vec = [
            ClimaAnalysis.global_bias(
                ClimaAnalysis.slice(sim_var, time = t),
                ClimaAnalysis.slice(obs_var, time = t),
                mask = mask,
            ) for t in times
        ]

        ax_sim = CairoMakie.Axis(
            fig_rmse_bias[1, col],
            title = "Global sim lonlat averages for $short_name",
            xlabel = "Month",
            ylabel = "Global lonlat averages ($(ClimaAnalysis.units(sim_var)))",
            xticks = (1:12, Dates.monthabbr.(1:12)),
        )
        ax_rmse = CairoMakie.Axis(
            fig_rmse_bias[2, col],
            title = "Global RMSE for $short_name",
            xlabel = "Month",
            ylabel = "Global RMSE ($(ClimaAnalysis.units(sim_var)))",
            xticks = (1:12, Dates.monthabbr.(1:12)),
        )
        ax_bias = CairoMakie.Axis(
            fig_rmse_bias[3, col],
            title = "Global bias for $short_name",
            xlabel = "Month",
            ylabel = "Global bias ($(ClimaAnalysis.units(sim_var)))",
            xticks = (1:12, Dates.monthabbr.(1:12)),
        )

        # Partition months, rmse_vec, and bias_vec into years
        months =
            Dates.month.(
                Dates.DateTime(sim_var.attributes["start_date"]) .+
                Dates.Second.(times),
            )
        months_split, sim_vec_split, rmse_vec_split, bias_vec_split =
            partition_by_val(12, months, sim_vec, rmse_vec, bias_vec)

        # Plot each year with the earlier year being more transparent than the later years
        axes = (ax_sim, ax_rmse, ax_bias)
        num_years = length(months_split)
        for (curr_year, (months, sim_vec, rmse_vec, bias_vec)) in enumerate(
            zip(months_split, sim_vec_split, rmse_vec_split, bias_vec_split),
        )
            alpha = curr_year / num_years
            data_vecs = [sim_vec, rmse_vec, bias_vec]

            for (ax, data_vec) in zip(axes, data_vecs)
                CairoMakie.lines!(
                    ax,
                    months,
                    data_vec,
                    alpha = alpha,
                    color = :blue,
                )
            end
        end

        # Compute the average over each of the months
        num_months = 12
        num_years = length(sim_vec_split)
        average_per_months = (
            begin
                season_to_avg = compute_group_averages(months, data_vec)
                [get(season_to_avg, month, NaN) for month in 1:12]
            end for data_vec in (sim_vec, rmse_vec, bias_vec)
        )
        for (ax, data_vec) in zip(axes, average_per_months)
            CairoMakie.lines!(ax, 1:12, data_vec, color = :orange)
        end
    end

    # Add a legend for the meaning of the black line
    blue_line = CairoMakie.LineElement(color = :blue)
    orange_line = CairoMakie.LineElement(color = :orange)
    CairoMakie.Legend(
        fig[1, 2],
        [blue_line, orange_line],
        [
            "More transparent - earlier years\nLess transparent - later years",
            "Average over each month",
        ],
    )

    CairoMakie.save(
        joinpath(
            leaderboard_base_path,
            "$(data_source)_global_rmse_and_bias_graphs.png",
        ),
        fig,
    )
end

"""
    compute_seasonal_leaderboard(leaderboard_base_path,
                                 diagnostics_folder_path,
                                 data_source)

Plot the biases and a seasonal leaderboard of various variables defined over longitude,
latitude, and time. The observational data is determined by `data_source` and can either be
`ILAMB` or `ERA5`.

The parameter `leaderboard_base_path` is the path to save the leaderboards and bias plots,
and `diagnostics_folder_path` is the path to the simulation data.

Loading and preprocessing simulation data is done by `get_sim_var_dict`. Loading and
preprocessing observational data is done by `get_obs_var_dict`. The masks are normalizing
the global RMSE and bias are determined by `get_mask_dict`. The ranges of the bias plots are
determined by `get_compare_vars_biases_plot_extrema`. See the functions defined in
data_sources.jl.
"""
function compute_seasonal_leaderboard(
    leaderboard_base_path,
    diagnostics_folder_path,
    data_source,
)
    # Get everything we need from data_sources.jl
    sim_var_dict = get_sim_var_dict(diagnostics_folder_path)
    obs_var_dict = get_obs_var_dict(data_source)
    short_names = intersect(keys(sim_var_dict), keys(obs_var_dict))

    # Need to intialize mask function
    mask_dict = get_mask_dict(data_source)
    issubset(short_names, keys(mask_dict)) ||
        error("Not all variables ($short_names) have a mask $(keys(mask_dict))")
    # Store the mask functions after initialization
    mask_fn_dict = Dict()

    compare_vars_biases_plot_extrema = get_compare_vars_biases_plot_extrema()

    # Set up dict for storing simulation and observational data after processing
    # Map short name to Dict which maps season to tuple of OutputVars
    sim_obs_season_comparsion_dict = Dict()
    # Map short name to time series of time averages for each season
    sim_obs_time_avg_over_seasons_comparsion_dict = Dict()
    seasons = ["ANN", "MAM", "JJA", "SON", "DJF"]

    spin_up_months = 12
    for short_name in short_names
        @info short_name
        # Simulation data
        sim_var = sim_var_dict[short_name]()

        # Observational data
        obs_var = obs_var_dict[short_name](sim_var.attributes["start_date"])

        # Make masking function
        mask_fn_dict[short_name] = mask_dict[short_name](sim_var, obs_var)

        # Remove first spin_up_months from simulation if possible
        spinup_cutoff = spin_up_months * 31 * 86400.0
        ClimaAnalysis.times(sim_var)[end] >= spinup_cutoff && (
            sim_var =
                ClimaAnalysis.window(sim_var, "time", left = spinup_cutoff)
        )

        # Determine which times can be used
        sim_times = ClimaAnalysis.times(sim_var)
        obs_times = ClimaAnalysis.times(obs_var)
        min_time = maximum(first.((sim_times, obs_times)))
        max_time = minimum(last.((sim_times, obs_times)))

        # Window sim_var and obs_var
        sim_var = ClimaAnalysis.window(
            sim_var,
            "time",
            left = min_time,
            right = max_time,
        )
        obs_var = ClimaAnalysis.window(
            obs_var,
            "time",
            left = min_time,
            right = max_time,
        )

        # Resample
        obs_var = ClimaAnalysis.shift_longitude(obs_var, -180.0, 180.0)
        obs_var = ClimaAnalysis.resampled_as(obs_var, sim_var)
        sim_var_seasons = (sim_var, ClimaAnalysis.split_by_season(sim_var)...)
        obs_var_seasons = (obs_var, ClimaAnalysis.split_by_season(obs_var)...)

        # Take time average
        sim_var_seasons = (
            !isempty(sim_var) ? ClimaAnalysis.average_time(sim_var) : sim_var for sim_var in sim_var_seasons
        )
        obs_var_seasons = (
            !isempty(obs_var) ? ClimaAnalysis.average_time(obs_var) : obs_var for obs_var in obs_var_seasons
        )

        # Save observation and simulation data
        sim_obs_season_comparsion_dict[short_name] = Dict(
            season => (sim_var_s, obs_var_s) for
            (season, sim_var_s, obs_var_s) in
            zip(seasons, sim_var_seasons, obs_var_seasons)
        )

        # Compute time averages across seasons
        obs_var_seasons_over_time =
            ClimaAnalysis.split_by_season_across_time(obs_var)
        sim_var_seasons_over_time =
            ClimaAnalysis.split_by_season_across_time(sim_var)

        # Take time average of each season, it is reasonable to assume that
        # there is no missing months between the first and last points in time
        time_averages_sim_var =
            ClimaAnalysis.average_time.(sim_var_seasons_over_time)
        time_averages_obs_var =
            ClimaAnalysis.average_time.(obs_var_seasons_over_time)

        sim_obs_time_avg_over_seasons_comparsion_dict[short_name] =
            (time_averages_sim_var, time_averages_obs_var)
    end

    # Make global seasonal averages over all years (after removing spinup)
    # Rows correspond to short names
    # Cols correspond to seasons
    groups = ["ANN", "MAM", "JJA", "SON", "DJF"]
    fig_bias = CairoMakie.Figure(;
        size = (600 * length(groups), 400 * length(short_names)),
    )
    for (row_idx, short_name) in enumerate(short_names)
        CairoMakie.Label(
            fig_bias[row_idx, 0],
            short_name,
            tellheight = false,
            fontsize = 30,
        )
        for (col_idx, group) in enumerate(groups)
            sim_var, obs_var = sim_obs_season_comparsion_dict[short_name][group]
            isempty(sim_var) && break
            layout = fig_bias[row_idx, col_idx] = CairoMakie.GridLayout()
            sim_var.attributes["short_name"] = "mean $(ClimaAnalysis.short_name(sim_var))"
            ClimaAnalysis.Visualize.plot_bias_on_globe!(
                layout,
                sim_var,
                obs_var,
                cmap_extrema = compare_vars_biases_plot_extrema[short_name],
                mask = mask_fn_dict[short_name],
            )
        end
    end

    # Plot the labels for the short names
    for (col_idx, group) in enumerate(groups)
        CairoMakie.Label(
            fig_bias[0, col_idx],
            group,
            tellwidth = false,
            fontsize = 30,
        )
    end

    # Add a title
    titlelayout = CairoMakie.GridLayout(
        fig_bias[-1, 1:length(groups)],
        halign = :center,
        tellwidth = false,
    )
    CairoMakie.Label(
        titlelayout[1, 1],
        "Annual and seasonal biases over all years excluding spin up ($spin_up_months months)",
        fontsize = 40,
    )
    CairoMakie.save(
        joinpath(
            leaderboard_base_path,
            "$(data_source)_seasonal_avg_over_all_years.png",
        ),
        fig_bias,
    )


    # Add plot of time average for simulation and bias data excluding spin up
    # Rows correspond to short names
    # Cols correspond to "SIM" and "ANN"
    groups = ["SIM", "ANN"]
    fig_sim_ann = CairoMakie.Figure(;
        size = (600 * length(groups), 400 * length(short_names)),
    )
    for (row_idx, short_name) in enumerate(short_names)
        CairoMakie.Label(
            fig_sim_ann[row_idx, 0],
            short_name,
            tellheight = false,
            fontsize = 30,
        )
        for (col_idx, group) in enumerate(groups)
            if group == "SIM"
                sim_var, _ = sim_obs_season_comparsion_dict[short_name]["ANN"]
                isempty(sim_var) && break
                layout = fig_sim_ann[row_idx, col_idx] = CairoMakie.GridLayout()
                ClimaAnalysis.Visualize.contour2D_on_globe!(
                    layout,
                    sim_var,
                    mask = mask_fn_dict[short_name],
                )
            else
                sim_var, obs_var =
                    sim_obs_season_comparsion_dict[short_name][group]
                isempty(sim_var) && break
                layout = fig_sim_ann[row_idx, col_idx] = CairoMakie.GridLayout()
                ClimaAnalysis.Visualize.plot_bias_on_globe!(
                    layout,
                    sim_var,
                    obs_var,
                    cmap_extrema = compare_vars_biases_plot_extrema[short_name],
                    mask = mask_fn_dict[short_name],
                )
            end
        end
    end

    # Plot the labels for the short names
    for (col_idx, group) in enumerate(groups)
        CairoMakie.Label(
            fig_sim_ann[0, col_idx],
            group,
            tellwidth = false,
            fontsize = 30,
        )
    end

    # Add a title
    titlelayout = CairoMakie.GridLayout(
        fig_sim_ann[-1, 1:length(groups)],
        halign = :center,
        tellwidth = false,
    )
    CairoMakie.Label(
        titlelayout[1, 1],
        "Time average for simulation and bias data\nexcluding spin up ($spin_up_months months)",
        fontsize = 40,
    )

    CairoMakie.save(
        joinpath(
            leaderboard_base_path,
            "$(data_source)_sim_annual_time_avg_over_all_years.png",
        ),
        fig_sim_ann,
    )

    # Make plot with seasons on x-axis and RMSE and bias on the y-axis
    # Rows correspond to short names
    # Cols correspond to SIM and ANN
    # Set up figure to plot on
    fig = CairoMakie.Figure(size = (450 * length(short_names), 900))
    fig_rmse_bias = fig[1, 1] = CairoMakie.GridLayout()


    # Loop over sim_obs_season_over_time_comparsion_dict
    for (col, short_name) in enumerate(short_names)
        sim_vars, obs_vars =
            sim_obs_time_avg_over_seasons_comparsion_dict[short_name]
        mask = mask_fn_dict[short_name]

        # Get season and compute global bias and global rmse
        seasons = [sim_var.attributes["season"] for sim_var in sim_vars]
        sim_vec = [
            ClimaAnalysis.weighted_average_lonlat(
                ClimaAnalysis.apply_oceanmask(sim_var),
            ).data[] for sim_var in sim_vars
        ]
        rmse_vec = [
            ClimaAnalysis.global_rmse(sim_var, obs_var, mask = mask) for
            (sim_var, obs_var) in zip(sim_vars, obs_vars)
        ]
        bias_vec = [
            ClimaAnalysis.global_bias(sim_var, obs_var, mask = mask) for
            (sim_var, obs_var) in zip(sim_vars, obs_vars)
        ]

        # Partition by seasons
        # Map each season to a number for plotting
        season_to_num = Dict("MAM" => 1, "JJA" => 2, "SON" => 3, "DJF" => 4)
        seasons = [season_to_num[season] for season in seasons]
        seasons_split, sim_vec_split, bias_vec_split, rmse_vec_split =
            partition_by_val(4, seasons, sim_vec, bias_vec, rmse_vec)

        # Set up three axes
        ax_sim = CairoMakie.Axis(
            fig_rmse_bias[1, col],
            title = "Global sim lonlat averages for $short_name",
            xlabel = "Season",
            ylabel = "Global lonlat averages ($(ClimaAnalysis.units(first(sim_vars))))",
            xticks = (1:4, ["MAM", "JJA", "SON", "DJF"]),
        )
        ax_rmse = CairoMakie.Axis(
            fig_rmse_bias[2, col],
            title = "Global RMSE for $short_name",
            xlabel = "Season",
            ylabel = "Global RMSE ($(ClimaAnalysis.units(first(sim_vars))))",
            xticks = (1:4, ["MAM", "JJA", "SON", "DJF"]),
        )
        ax_bias = CairoMakie.Axis(
            fig_rmse_bias[3, col],
            title = "Global bias for $short_name",
            xlabel = "Season",
            ylabel = "Global bias ($(ClimaAnalysis.units(first(sim_vars))))",
            xticks = (1:4, ["MAM", "JJA", "SON", "DJF"]),
        )

        # Plot on axes
        axes = (ax_sim, ax_rmse, ax_bias)
        num_years = length(seasons_split)
        for (curr_year, (seasons, sim_vec, rmse_vec, bias_vec)) in enumerate(
            zip(seasons_split, sim_vec_split, rmse_vec_split, bias_vec_split),
        )
            alpha = curr_year / num_years
            data_vecs = [sim_vec, rmse_vec, bias_vec]

            for (ax, data_vec) in zip(axes, data_vecs)
                CairoMakie.lines!(
                    ax,
                    seasons,
                    data_vec,
                    alpha = alpha,
                    color = :blue,
                )
            end
        end

        # Compute the average over each of the seasons
        num_years = length(sim_vec_split)
        average_per_seasons = (
            begin
                season_to_avg = compute_group_averages(seasons, data_vec)
                [get(season_to_avg, season, NaN) for season in 1:4]
            end for data_vec in (sim_vec, rmse_vec, bias_vec)
        )

        for (ax, data_vec) in zip(axes, average_per_seasons)
            CairoMakie.lines!(ax, 1:4, data_vec, color = :orange)
        end
    end

    # Add a legend for the meaning of the black line
    blue_line = CairoMakie.LineElement(color = :blue)
    orange_line = CairoMakie.LineElement(color = :orange)
    CairoMakie.Legend(
        fig[1, 2],
        [blue_line, orange_line],
        [
            "More transparent - earlier years\nLess transparent - later years",
            "Average over each season",
        ],
    )

    CairoMakie.save(
        joinpath(
            leaderboard_base_path,
            "$(data_source)_seasonal_global_rmse_and_bias_graphs.png",
        ),
        fig,
    )
end

"""
    partition_by_val(val, x_vec, y_vecs...)

Partition `x_vec` and `y_vecs` into subarrays, so that the last value of each
subarray is `val` if possible.

Examples
=========

```jldoctest
julia> partition_by_val(4, [3,4,1,2,3,4,1,2], collect(1:8))
2-element Vector{Vector{Vector{Int64}}}:
 [[3, 4], [1, 2, 3, 4], [1, 2]]
 [[1, 2], [3, 4, 5, 6], [7, 8]]
```
"""
function partition_by_val(val, x_vec, y_vecs...)
    for vec in y_vecs
        length(x_vec) != length(vec) && error(
            "Length of $x_vec ($(length(x_vec))) and length of $vec ($(length(vec))) are not the same",
        )
    end
    start_and_end_indices = vcat(0, findall(==(val), x_vec))
    if !(length(x_vec) in start_and_end_indices)
        push!(start_and_end_indices, length(x_vec))
    end
    vecs = (x_vec, y_vecs...)
    ret_vecs = [Vector{eltype(vec)}[] for vec in vecs]
    for idx in eachindex(start_and_end_indices[1:(end - 1)])
        for (vec, ret_vec) in zip(vecs, ret_vecs)
            start_idx = start_and_end_indices[idx] + 1
            end_idx = start_and_end_indices[idx + 1]
            push!(ret_vec, vec[start_idx:end_idx])
        end
    end
    return ret_vecs
end

"""
    compute_group_averages(groups, vals)

Return a dictionary mapping values in `groups` to average of the corresponding
values in `vals`.

Examples
=========

```jldoctest
julia> compute_group_averages([1,2,3,4,1,2], [1,2,100,200,11,-2])
Dict{Int64, Float64} with 4 entries:
  4 => 200.0
  2 => 0.0
  3 => 100.0
  1 => 6.0
```
"""
function compute_group_averages(groups, vals)
    length(groups) != length(vals) &&
        error("Length of $groups and $vals are not the same")
    group_to_vals = Dict{eltype(groups), Vector{eltype(vals)}}()
    for (group, val) in zip(groups, vals)
        if group ∉ keys(group_to_vals)
            group_to_vals[group] = [val]
        else
            push!(group_to_vals[group], val)
        end
    end
    return Dict(group => mean(vals) for (group, vals) in group_to_vals)
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) != 3
        error(
            "Usage: julia leaderboard.jl <Filepath to save leaderboard and bias plots> <Filepath to simulation data> <\"ERA5\" or \"ILAMB\">",
        )
    end
    leaderboard_base_path = ARGS[begin]
    diagnostics_folder_path = ARGS[begin + 1]
    data_source = ARGS[begin + 2]

    # Error handling
    isdir(leaderboard_base_path) ||
        error("$leaderboard_base_path is not a directory")

    compute_monthly_leaderboard(
        leaderboard_base_path,
        diagnostics_folder_path,
        data_source,
    )
    compute_seasonal_leaderboard(
        leaderboard_base_path,
        diagnostics_folder_path,
        data_source,
    )
end
