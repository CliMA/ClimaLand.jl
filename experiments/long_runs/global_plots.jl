using ClimaUtilities.ClimaArtifacts
import ClimaDiagnostics
import ClimaAnalysis as CAN
import ClimaAnalysis.Visualize as viz
import ClimaUtilities
import ClimaUtilities.ClimaArtifacts: @clima_artifact
import ClimaParams as CP
using ClimaLand
using ClimaLand.Snow
using ClimaLand.Soil
using ClimaLand.Canopy
import ClimaLand
import ClimaLand.Parameters as LP
using CairoMakie
import GeoMakie
using Dates
using Statistics

# Load functions to access the data
include("data_global_plots.jl")

"""
    make_global_plots(
        data_dir;
        root_save_path = pwd(),
        plot_net_radiation = false,
    )

Make a set of global plots comparing simulation data to observational data.
The plots include contour plots of the annually averaged simulation and
observational data, as well as a bias plot.

By default, the plots use the second year of data from the simulation; this
can be changed by modifying the `year_ind` variable. The observational data
plotted is from 2008.

Also by default, the plots include the following variables:
- Latent heat flux (LE)
- Sensible heat flux (H)
- Upward longwave radiation (LWᵤ)
- Upward shortwave radiation (SWᵤ)

If `plot_net_radiation = true`, LWᵤ and SWᵤ are replaced with net
radiation fluxes (LWₙ, SWₙ).
This isn't done by default because the net radiation flux observation data
takes some time to compute, since it uses the downward radiative fluxes,
which come from a large hourly dataset.
"""
function make_global_plots(
    data_dir;
    root_save_path = pwd(),
    plot_net_radiation = false,
)
    # Get information about the variables to plot
    short_names, title_stubs, levels_dict = get_var_info(plot_net_radiation)

    # Load the data
    sim_var_dict = get_sim_var_dict(CAN.SimDir(data_dir))
    obs_var_dict = get_obs_var_dict()

    # Create figure for all plots
    num_cols = 3
    fig = CairoMakie.Figure(
        size = (500 * num_cols, 405 * length(short_names)),
        figure_padding = 50,
    )

    for (idx, short_name) in enumerate(short_names)
        # Plot figures in odd rows, colorbars in even rows
        fig_row = (idx - 1) * 2 + 1

        # Access simulation and observation data for this variable
        sim_var = CAN.apply_oceanmask(
            CAN.shift_to_start_of_previous_month(sim_var_dict[short_name]()),
        )
        sim_var_times = CAN.times(sim_var)
        obs_var = CAN.apply_oceanmask(
            obs_var_dict[short_name](sim_var.attributes["start_date"]),
        )

        # Check that the variables are 2D (+time)
        length(sim_var.dims) == 3 ||
            error("Can only plot 2D simulation variables")
        length(obs_var.dims) == 3 ||
            error("Can only plot 2D observational variables")

        title_stub = title_stubs[short_name]
        units_label = "(" * sim_var.attributes["units"] * ")"

        # Select the time window of data to plot
        # Use data from year 1 of the simulation
        year_ind = 1
        sim_var_window = CAN.window(
            sim_var,
            "time",
            left = (year_ind - 1) * 366 * 86400 + 30 * 86400, # 1 year left of year_ind, in seconds.
            right = year_ind * 366 * 86400, # 1 year right of year_ind, in seconds
        )

        obs_var_window = CAN.window(
            obs_var,
            "time",
            left = 0, # observation data starts at 2008
            right = 366 * 86400, # 1 year of observation data, in seconds
        )

        # Load the visualization mask
        mask = viz.oceanmask()

        ### GLOBAL HEATMAPS (plotted by default)
        # Average the data over this window of time
        sim_var_annual_average = CAN.average_time(sim_var_window)
        obs_var_annual_average = CAN.average_time(obs_var_window)

        # Get colorbar limits to share between heatmaps
        clims, levels, nlevels =
            get_clims(sim_var_annual_average, obs_var_annual_average)

        # Make colorbar tick labels
        ticklabels = map(x -> string(x), levels)
        ticks = (levels, ticklabels)

        # Choose colormap based on variable
        # Treat SHF differently because we include negative and positive values
        if short_name in ["shf"]
            # Make sure that 0 is at the center
            cmap = viz._constrained_cmap(
                Makie.cgrad(:PRGn_9).colors,
                clims[1],
                clims[2];
                categorical = true,
            )
        else
            cmap = Makie.cgrad(
                :Greens_4,
                nlevels - 1;
                categorical = true,
                rev = false,
            )
        end

        # Plot simulation data heatmap
        sim_title =
            fig_row == 1 ?
            CairoMakie.rich("ClimaLand, annually averaged", fontsize = 28) : "" # title of the figure
        p_loc_sim = [fig_row, 1] # column 1
        label = "$title_stub $units_label"
        plot_var_contour!(
            fig,
            sim_var_annual_average,
            p_loc_sim,
            sim_title,
            levels,
            cmap,
            mask,
            ticks,
            label,
        )

        # Plot observational data heatmap
        obs_title =
            fig_row == 1 ?
            CairoMakie.rich("ERA5, annually averaged", fontsize = 28) : "" # title of the figure
        p_loc_obs = [fig_row, 2] # column 2
        plot_var_contour!(
            fig,
            obs_var_annual_average,
            p_loc_obs,
            obs_title,
            levels,
            cmap,
            mask,
            ticks,
            label,
        )

        ### GLOBAL BIAS PLOTS
        # Compute the bias
        apply_mask = CAN.apply_oceanmask
        bias_var = CAN.bias(
            sim_var_annual_average,
            obs_var_annual_average,
            mask = apply_mask,
        )
        global_bias = round(bias_var.attributes["global_bias"], sigdigits = 3)
        rmse = round(
            CAN.global_rmse(
                sim_var_annual_average,
                obs_var_annual_average,
                mask = apply_mask,
            ),
            sigdigits = 3,
        )
        # units = CAN.units(bias_var)
        # bias_var.attributes["long_name"] *= " (RMSE: $rmse $units, Global bias: $global_bias $units)"

        # Set up colormap
        levels_bias = levels_dict["$(short_name)_bias"]
        clims_bias = (levels_bias[1], levels_bias[end])
        color_scheme = :vik
        min_level, max_level = clims_bias
        # Make sure that 0 is at the center
        cmap = Visualize._constrained_cmap(
            Makie.cgrad(color_scheme).colors,
            min_level,
            max_level;
            categorical = true,
        )

        # Set up plot parameters
        bias_title =
            fig_row == 1 ?
            CairoMakie.rich("ClimaLand vs ERA5 bias", fontsize = 28) : "" # title of the figure
        ticklabels = map(x -> string(Int(round(x))), levels_bias)#; digits = digits_to_round)), levels)
        ticks = (levels_bias, ticklabels)
        p_loc_bias = [fig_row, 3] # column 3

        # Plot the bias contour plot
        plot_var_contour!(
            fig,
            bias_var,
            p_loc_bias,
            bias_title,
            levels_bias,
            cmap,
            mask,
            ticks,
            label,
        )
    end
    # Save the figure
    save_name = joinpath(root_save_path, "global_plots_bias.png")
    CairoMakie.save(save_name, fig)
    @show save_name
    return nothing
end

function make_seasonal_plots(
    data_dir;
    root_save_path = pwd(),
    plot_net_radiation = false,
)
    # Get information about the variables to plot
    short_names, title_stubs, levels_dict = get_var_info(plot_net_radiation)

    # Load the data
    sim_var_dict = get_sim_var_dict(CAN.SimDir(data_dir))
    obs_var_dict = get_obs_var_dict()

    # create figure for all plots
    num_cols = 2
    num_rows = 2
    fig = CairoMakie.Figure(size = (550num_cols, 350num_rows))

    for (idx, short_name) in enumerate(short_names)
        # Determine row and column index for 2x2 grid plotting
        row_idx = mod(idx, 2) + 1 # odd var idx in even rows, offset by 1 for title
        col_idx = idx <= fld(length(short_names), 2) ? 1 : 2 # first half of vars in first column (floor division)

        # Access simulation and observation data for this variable
        sim_var = CAN.apply_oceanmask(
            CAN.shift_to_start_of_previous_month(sim_var_dict[short_name]()),
        )
        sim_var_times = CAN.times(sim_var)
        obs_var = CAN.apply_oceanmask(
            obs_var_dict[short_name](sim_var.attributes["start_date"]),
        )
        title_stub = title_stubs[short_name]
        units_label = "(" * sim_var.attributes["units"] * ")"
        label = "$title_stub $units_label"

        # Check that the times are aligned between simulation and observation data
        @assert all(
            sim_var.dims["time"][1:12] - obs_var.dims["time"][1:12] .== 0,
        )

        # At each point in time, compute the global average of the variable
        n_months = 12
        sim_var_global_average = zeros(n_months)
        obs_var_global_average = zeros(n_months)
        for i in 1:n_months
            # This assumes all variables are 2D (no z dimension)
            sim_slice_args = Dict(:time => sim_var_times[i])
            obs_slice_args = Dict(:time => sim_var_times[i])

            sim_var_sliced = CAN.slice(sim_var; sim_slice_args...)
            sim_var_masked = CAN.apply_oceanmask(sim_var_sliced)
            sim_var_global_average[i] =
                CAN.weighted_average_lonlat(sim_var_masked).data[1]

            obs_var_sliced = CAN.slice(obs_var; obs_slice_args...)
            obs_var_masked = CAN.apply_oceanmask(obs_var_sliced)
            obs_var_global_average[i] =
                CAN.weighted_average_lonlat(obs_var_masked).data[1]
        end

        month_names = [
            "Jan",
            "Feb",
            "Mar",
            "Apr",
            "May",
            "Jun",
            "Jul",
            "Aug",
            "Sep",
            "Oct",
            "Nov",
            "Dec",
        ]

        # Set up axis
        ax = Axis(
            fig[row_idx, col_idx], # last column
            ylabel = label,
            height = 250,
            xgridvisible = false,
            ygridvisible = false,
            xticks = (1:1:12, month_names),
        )

        # Plot simulation data
        CairoMakie.lines!(
            ax,
            sim_var_global_average,
            color = :blue,
            linewidth = 3,
            label = "ClimaLand",
        )

        # Plot observational data
        CairoMakie.scatter!(
            ax,
            obs_var_global_average,
            color = :orange,
            label = "ERA5",
        )

        # Add legend to top left plot
        row_idx == col_idx == 1 &&
            CairoMakie.axislegend(ax, position = :rt, labelsize = 18)

    end
    # Add title at the end so it spans all columns
    Label(
        fig[0, :],
        "ClimaLand vs ERA5 seasonal cycle",
        fontsize = 24,
        font = "TeX Gyre Heros Bold Makie",
    )

    save_name = joinpath(root_save_path, "seasonal_cycle.png")
    CairoMakie.save(save_name, fig)
    @show save_name
    return nothing
end

"""
    get_var_info(plot_net_radiation)

Get information about the variables to plot, including variable names,
titles to use, and bias colorbar levels. Only plot net radiation fluxes
if explicitly requested, since they take some time to compute.
"""
function get_var_info(plot_net_radiation)
    short_names =
        plot_net_radiation ? ["lhf", "shf", "lwn", "swn"] :
        ["lhf", "shf", "lwu", "swu"]

    title_stubs = Dict(
        "lhf" => "LE", #"Latent heat flux",
        "shf" => "H", #"Sensible heat flux",
        "lwu" => "LWᵤ", #"Upward longwave radiation",
        "swu" => "SWᵤ", #"Upward shortwave radiation",
        "lwn" => "LWₙ", #"Net longwave radiation",
        "swn" => "SWₙ", #"Net shortwave radiation",
    )

    # These have been estimated by visual inspection, and may need to be adjusted
    levels_dict = Dict(
        "lhf_bias" => collect(-50:10:30),
        "shf_bias" => collect(-30:10:50),
        "lwu_bias" => collect(-40:5:20),
        "swu_bias" => collect(-20:5:40),
        "lwn_bias" => collect(-40:5:40),
        "swn_bias" => collect(-40:5:40),
    )

    return short_names, title_stubs, levels_dict
end

# TODO this function can be replaced by `CAN.weighted_average_latlon(masked_var)` once a release is made
function compute_global_average(masked_var)
    latitude_name = ClimaAnalysis.latitude_name(masked_var)
    lon_name = ClimaAnalysis.longitude_name(masked_var)

    land_data = ClimaAnalysis.apply_oceanmask(masked_var).data
    lat_data = masked_var.dims[latitude_name]
    mask = .~isnan.(land_data)
    nlon = length(masked_var.dims[lon_name])
    resized_lat_data = transpose(repeat(lat_data, 1, nlon))

    return sum(land_data[mask] .* cosd.(resized_lat_data[mask])) /
           sum(cosd.(resized_lat_data[mask]))
end

"""
    get_clims(sim_var, obs_var)

Get colorbar limits to use for heatmaps of the simulation and observation data.
We want to use the same colorbar limits for both heatmaps, so that the colors are
consistent between the two. To do this, we combine the data from both and
use the 5th and 95th percentiles of the combined data as the limits.
Depending on the range of the data, we use a step size of 10 or 40 for the colorbar levels.

Returns:
    clims::Tuple{Int, Int} - The colorbar limits to use for the heatmaps
    levels::Vector{Int} - The levels to use for the colorbar
    nlevels::Int - The number of levels in the colorbar
"""
function get_clims(sim_var, obs_var)
    # This is a hacky way to remove NaNs in data, but it works for now
    sim_var_no_nans = deepcopy(sim_var.data)
    obs_var_no_nans = deepcopy(obs_var.data)
    sim_var_no_nans[isnan.(sim_var_no_nans)] .= 0
    obs_var_no_nans[isnan.(obs_var_no_nans)] .= 0

    # Combine the simulation and observation data
    combined_data = vcat(vec(sim_var_no_nans), vec(obs_var_no_nans))

    # Use the 5th and 95th percentiles as the colorbar limits
    min = Int(round(quantile(combined_data, 0.05); digits = -1))
    max = Int(round(quantile(combined_data, 0.95); digits = -1))

    # Choose the step size for colorbar levels based on the range of the data
    diff = max - min
    step_size = diff < 150 ? 10 : 40

    max = (diff % step_size != 0) ? max + step_size - (diff % step_size) : max
    clims = (min, max)
    levels = collect(min:step_size:max)
    nlevels = length(levels)

    return clims, levels, nlevels
end

"""
    plot_var_contour!(fig, var, p_loc, title, levels, cmap, mask, ticks, label)

Plot a heatmap of the variable `var` on the figure `fig` at position `p_loc`.
The heatmap is plotted with the given `levels` and `cmap`, and the `mask` is overlaid.
A colorbar is added to the figure with the given `ticks` and `label`.
"""
function plot_var_contour!(
    fig,
    var::CAN.OutputVar,
    p_loc,
    title,
    levels,
    cmap,
    mask,
    ticks,
    label,
)
    ax = GeoMakie.GeoAxis(
        fig[p_loc[1], p_loc[2]];
        title = title,
        xticklabelsvisible = false, # don't show lat labels
        yticklabelsvisible = false, # don't show lon labels
        xgridvisible = false, # don't show lat grid
        ygridvisible = false, # don't show lon grid
        height = 235,
        # ylabel = label, # plot variable label on y-axis for leftmost column (sim)
        # ylabelvisible = true,
        # ylabelpadding = 0,
    )
    lon = CAN.longitudes(var)
    lat = CAN.latitudes(var)
    plot = Makie.contourf!(
        ax,
        lon,
        lat,
        var.data;
        rasterize = true,
        levels = levels,
        colormap = cmap,
        extendhigh = :auto,
        extendlow = :auto,
    )
    Makie.poly!(ax, mask; color = :white) # plot mask
    Makie.lines!(ax, GeoMakie.coastlines(); color = :black) # plot coastlines

    # Add colorbar
    Makie.Colorbar(
        fig[p_loc[1] + 1, p_loc[2]], # place colorbar below plot
        plot,
        label = label,
        vertical = false, # horizontal colorbar
        flipaxis = false, # label underneath colorbar
        labelsize = 20,
        width = 400, # a little smaller
        tellwidth = false; # make colorbar width indep of plot width
        alignmode = Mixed(top = -20),
        ticks = ticks,
    )
    Makie.rowgap!(fig.layout, 0) # move colorbar closer to plot
end
