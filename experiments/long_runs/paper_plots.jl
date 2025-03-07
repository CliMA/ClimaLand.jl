using ClimaUtilities.ClimaArtifacts
import ClimaDiagnostics
import ClimaAnalysis
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

root_path = joinpath(pwd(), "snowy_land_longrun_gpu")
!isdir(root_path) && mkdir(root_path)
# outdir = "snowy_land_longrun_gpu/output_active" # on local
# Set outdir to wherever diagnostics are saved
# outdir = "snowy_land_longrun_gpu-3720-a_larger2_ksat_alpha/output_active" # on local
# outdir = "snowy_land_longrun_gpu-3761-alpha1/output_active" # on local
# outdir = "snowy_land_longrun_gpu-3771-new_default/output_active"
# outdir = "snowy_land_longrun_gpu-3761-new_default/output_active" # on local
outdir = "snowy_land_longrun_gpu-3777-alpha2/output_active" # on local
root_path = outdir

short_names = ["lhf", "shf", "lwu", "swu", "lwn", "swn"]
title_stubs = Dict(
    "lhf" => "Latent heat flux",
    "shf" => "Sensible heat flux",
    "lwu" => "Upward longwave radiation",
    "swu" => "Upward shortwave radiation",
    "lwn" => "Net longwave radiation",
    "swn" => "Net shortwave radiation",
)
# Define levels for contour colorbars
levels_dict = Dict(
    "lhf_bias" => collect(-50:10:30),
    "shf_bias" => collect(-30:10:50),
    "lwu_bias" => collect(-40:5:20),
    "swu_bias" => collect(-20:5:40),
    "lwn_bias" => collect(-40:5:40),
    "swn_bias" => collect(-40:5:40),
)

include("data_paper_plots.jl")

function compute_global_average(masked_var)
    latitude_name = ClimaAnalysis.latitude_name(masked_var)
    lon_name = ClimaAnalysis.longitude_name(masked_var)

    land_data = ClimaAnalysis.apply_oceanmask(masked_var).data
    lat_data = masked_var.dims[latitude_name]
    mask = .~isnan.(land_data)
    nlon = length(masked_var.dims[lon_name])µ
    resized_lat_data = transpose(repeat(lat_data, 1, nlon))

    return sum(land_data[mask] .* cosd.(resized_lat_data[mask])) /
           sum(cosd.(resized_lat_data[mask]))
end

function make_paper_figures(
    root_path,
    outdir,
    short_names,
    title_stubs;
    plot_bias = false,
    plot_seasonal = false,
)
    # Set up for comparison to data (copied from leaderboard.jl)
    # use sim_var and obs_var together for the seasonal plot because they already have the same units :)
    sim_var_dict = get_sim_var_dict(ClimaAnalysis.SimDir(outdir))
    obs_var_dict = get_obs_var_dict()

    # create figure for all plots
    num_cols = plot_bias || plot_seasonal ? 3 : 2
    fig = CairoMakie.Figure(
        size = (500num_cols, 405 * length(short_names)),
        figure_padding = 50,
    )

    for (idx, short_name) in enumerate(short_names)
        # Plot figures in odd rows, colorbars in even rows
        fig_row = (idx - 1) * 2 + 1

        title_stub = title_stubs[short_name]

        # Access simulation data in the time we want
        sim_var = sim_var_dict[short_name]()
        sim_var_times = ClimaAnalysis.times(sim_var)

        obs_var = obs_var_dict[short_name](sim_var.attributes["start_date"])



        sim_var_global_average = zeros(12)
        obs_var_global_average = zeros(12)
        for i in 1:12
            sim_slice_args =
                ClimaAnalysis.has_altitude(sim_var) ?
                Dict(:z => 1, :time => sim_var_times[i]) :
                Dict(:time => sim_var_times[i]) # if has altitude, take first layer
            obs_slice_args =
                ClimaAnalysis.has_altitude(sim_var) ?
                Dict(:z => 1, :time => sim_var_times[i]) :
                Dict(:time => sim_var_times[i]) # if has altitude, take first layer

            sim_var_sliced = ClimaAnalysis.slice(sim_var; sim_slice_args...)
            sim_var_masked = ClimaAnalysis.apply_oceanmask(sim_var_sliced)
            sim_var_global_average[i] = compute_global_average(sim_var_masked)

            obs_var_sliced = ClimaAnalysis.slice(obs_var; obs_slice_args...)
            obs_var_masked = ClimaAnalysis.apply_oceanmask(obs_var_sliced)
            obs_var_global_average[i] = compute_global_average(obs_var_masked)
        end

        i = 1 # use first year of simulation
        kwarg_z = ClimaAnalysis.has_altitude(sim_var) ? Dict(:z => 1) : Dict() # if has altitude, take first layer
        sim_var_sliced = ClimaAnalysis.slice(sim_var; kwarg_z...)
        sim_var_window = ClimaAnalysis.window(
            sim_var_sliced,
            "time",
            left = (i - 1) * 366 * 86400 + 30 * 86400, # 1 year left of year i, in seconds.
            right = i * 366 * 86400, # 1 year right of year i, in seconds
        )
        units_label = "(" * sim_var.attributes["units"] * ")"

        # Access observation data in the time we want
        obs_var_sliced = ClimaAnalysis.slice(obs_var; kwarg_z...)
        obs_var_window = ClimaAnalysis.window(
            obs_var_sliced,
            "time",
            left = 0, # observation data starts at 2008
            right = 366 * 86400, # 1 year of observation data, in seconds
        )

        ## GLOBAL HEATMAP
        # sim_var_annual_average below contains annually-averaged data in the provided window
        sim_var_annual_average = ClimaAnalysis.average_time(sim_var_window)
        # ClimaAnalysis.center_longitude!(sim_var_annual_average, 0)

        # obs_var_annual_average below contains annually-averaged data in the provided window
        obs_var_annual_average = ClimaAnalysis.average_time(obs_var_window)
        # ClimaAnalysis.center_longitude!(obs_var_annual_average, 0)

        t = sim_var_times[1]

        # Get colorbar limits to share between heatmaps
        sim_var_no_nans = deepcopy(sim_var_annual_average.data)
        sim_var_no_nans[isnan.(sim_var_no_nans)] .= 0 # TODO this is a hacky way to remove NaNs in data
        # sim_extrema = extrema(sim_var_no_nans)
        # obs_extrema = extrema(obs_var_annual_average.data)

        combined_data =
            vcat(vec(sim_var_no_nans), vec(obs_var_annual_average.data))
        min = Int(round(quantile(combined_data, 0.05); digits = -1))
        max = Int(round(quantile(combined_data, 0.95); digits = -1))
        diff = max - min
        step_size = diff < 150 ? 10 : 40

        # Increase max to be a multiple of step_size
        max =
            (diff % step_size != 0) ? max + step_size - (diff % step_size) : max
        clims = (min, max)
        levels = collect(min:step_size:max)
        nlevels = length(levels)

        # extrema = (min(sim_extrema[1], obs_extrema[1]), max(sim_extrema[2], obs_extrema[2]))
        # we MUST pass the same clims to the sim plot, obs plot, AND colorbars
        # levels = levels_dict[short_name]
        # nlevels = length(levels)
        # clims = (levels[1], levels[end])

        # make colorbar tick labels
        ticklabels = map(x -> string(x), levels)#; digits = digits_to_round)), levels)
        ticks = (levels, ticklabels)

        # Plot simulation data heatmap
        sim_title =
            fig_row == 1 ?
            CairoMakie.rich("ClimaLand, annually averaged", fontsize = 28) : "" # title of the figure

        # Treat SHF differently because we include negative and positive values
        if short_name in ["shf"]
            # Plot simulation data heatmap
            viz.plot_bias_on_globe!(
                fig,
                sim_var_annual_average,
                sim_var_annual_average, # obs_var_annual_average,
                false, # plot_bias
                levels,
                p_loc = (fig_row, 1), # plot in the first column
                plot_colorbar = true,
                mask = viz.oceanmask(),
                cmap_extrema = clims,
                more_kwargs = Dict(
                    :mask => ClimaAnalysis.Utils.kwargs(color = :white),
                    :plot => ClimaAnalysis.Utils.kwargs(rasterize = true),
                    :axis => ClimaAnalysis.Utils.kwargs(
                        title = sim_title,
                        xticklabelsvisible = false, # don't show lat labels
                        yticklabelsvisible = false, # don't show lon labels
                        xgridvisible = false, # don't show lat grid
                        ygridvisible = false, # don't show lon grid
                        height = 235,
                        ylabel = "$(sim_var.attributes["long_name"]) $units_label", # plot variable label on y-axis for leftmost column (sim)
                        ylabelvisible = true,
                        ylabelpadding = 0,
                    ),
                    :cb => ClimaAnalysis.Utils.kwargs(
                        vertical = false, # horizontal colorbar
                        label = "$(sim_var.attributes["long_name"]) $units_label",
                        labelsize = 20,
                        flipaxis = false, # label underneath colorbar
                        height = 15,
                        width = 400, # a little smaller
                        tellwidth = false, # make colorbar width indep of plot width
                        alignmode = Mixed(top = -20),
                    ),
                ),
            )
        else
            # Choose color scheme based on variable and clims
            colormap = Makie.cgrad(
                :Greens_4,
                nlevels - 1;
                categorical = true,
                rev = false,
            )
            colors = colormap.colors.colors
            contour = viz.contour2D_on_globe!(
                fig,
                sim_var_annual_average,
                # ylabel = "$(sim_var.attributes["long_name"]) $units_label", # plot variable label on y-axis for leftmost column (sim)
                p_loc = (fig_row, 1), # plot in the first column
                plot_colorbar = false,
                mask = viz.oceanmask(),
                more_kwargs = Dict(
                    :mask => ClimaAnalysis.Utils.kwargs(color = :white),
                    :plot => ClimaAnalysis.Utils.kwargs(
                        rasterize = true,
                        # colorrange = clims,
                        colormap = colormap,
                        levels = levels,
                        extendhigh = :auto,
                        extendlow = :auto,
                    ),
                    :axis => ClimaAnalysis.Utils.kwargs(
                        title = sim_title,
                        xticklabelsvisible = false, # don't show lat labels
                        yticklabelsvisible = false, # don't show lon labels
                        xgridvisible = false, # don't show lat grid
                        ygridvisible = false, # don't show lon grid
                        height = 235,
                        ylabel = "$(sim_var.attributes["long_name"]) $units_label", # plot variable label on y-axis for leftmost column (sim)
                        ylabelpadding = -5,
                        ylabelvisible = true,
                        # limits = clims,
                    ),
                ),
            )
            Makie.Colorbar(
                fig[fig_row + 1, 1],
                label = "$(sim_var.attributes["long_name"]) $units_label",
                labelsize = 20,
                vertical = false, # horizontal colorbar
                flipaxis = false, # label underneath colorbar
                height = 15,
                width = 400, # a little smaller
                tellwidth = false, # make colorbar width indep of plot width
                colorrange = clims, # use global min/max across sim and obs plots
                colormap = colormap,
                ticks = ticks,
                highclip = last(colors),
                lowclip = first(colors),
                alignmode = Mixed(top = -6),
            )
        end

        # Plot observational data heatmap
        obs_title =
            fig_row == 1 ?
            CairoMakie.rich("ERA5, annually averaged", fontsize = 28) : "" # title of the figure
        if short_name in ["shf"]
            viz.plot_bias_on_globe!(
                fig,
                obs_var_annual_average,
                obs_var_annual_average, # obs_var_annual_average,
                false, # plot_bias
                levels,
                p_loc = (fig_row, 2), # plot in the first column
                plot_colorbar = true,
                # colorbar_label = "$(sim_var.attributes["long_name"]) $units_label",
                mask = viz.oceanmask(),
                cmap_extrema = clims,
                more_kwargs = Dict(
                    :mask => ClimaAnalysis.Utils.kwargs(color = :white),
                    :plot => ClimaAnalysis.Utils.kwargs(rasterize = true),
                    :axis => ClimaAnalysis.Utils.kwargs(
                        title = obs_title,
                        xticklabelsvisible = false, # don't show lat labels
                        yticklabelsvisible = false, # don't show lon labels
                        xgridvisible = false, # don't show lat grid
                        ygridvisible = false, # don't show lon grid
                        height = 235,
                        ylabelvisible = false,
                    ),
                    :cb => ClimaAnalysis.Utils.kwargs(
                        vertical = false, # horizontal colorbar
                        label = "$(sim_var.attributes["long_name"]) $units_label",
                        labelsize = 20,
                        flipaxis = false, # label underneath colorbar
                        height = 15,
                        width = 400, # a little smaller
                        tellwidth = false, # make colorbar width indep of plot width
                        alignmode = Mixed(top = -20),
                    ),
                ),
            )
        else
            # Choose color scheme based on variable and clims
            colormap = Makie.cgrad(
                :Greens_4,
                nlevels - 1;
                categorical = true,
                rev = false,
            )
            colors = colormap.colors.colors

            viz.contour2D_on_globe!(
                fig,
                obs_var_annual_average,
                p_loc = (fig_row, 2), # plot in the second column
                plot_colorbar = false,
                mask = viz.oceanmask(),
                more_kwargs = Dict(
                    :mask => ClimaAnalysis.Utils.kwargs(color = :white),
                    :plot => ClimaAnalysis.Utils.kwargs(
                        rasterize = true,
                        # colorrange = clims,
                        colormap = colormap,
                        levels = levels,
                        extendhigh = :auto,
                        extendlow = :auto,
                    ),
                    :axis => ClimaAnalysis.Utils.kwargs(
                        title = obs_title,
                        xticklabelsvisible = false, # don't show lat labels
                        yticklabelsvisible = false, # don't show lon labels
                        xgridvisible = false, # don't show lat grid
                        ygridvisible = false, # don't show lon grid
                        height = 235,
                        ylabelvisible = false,
                    ),
                ),
            )
            Makie.Colorbar(
                fig[fig_row + 1, 2],
                label = "$(sim_var.attributes["long_name"]) $units_label",
                labelsize = 20,
                vertical = false, # horizontal colorbar
                flipaxis = false, # label underneath colorbar
                height = 15,
                width = 400, # a little smaller
                tellwidth = false, # make colorbar width indep of plot width
                colorrange = clims, # use global min/max across sim and obs plots
                colormap = colormap,
                ticks = ticks,
                highclip = last(colors),
                lowclip = first(colors),
                alignmode = Mixed(top = -6),
            )
        end
        # Makie.colgap!(fig.layout, 0, fig_row) # reduce gap between plot/colorbar pairs
        # CairoMakie.save(joinpath(root_path, "$(short_name)_$(t)-annual_avg-obs.pdf"), fig_global_obs)

        @assert !(plot_bias && plot_seasonal) # only one of these can be true
        ## BIAS PLOT
        if plot_bias
            levels_bias = levels_dict["$(short_name)_bias"]
            clims_bias = (levels_bias[1], levels_bias[end])

            bias_title =
                fig_row == 1 ?
                CairoMakie.rich("ClimaLand vs ERA5 bias", fontsize = 28) : "" # title of the figure
            bias_colorbar_label = "$(sim_var.attributes["long_name"]) $units_label"
            viz.plot_bias_on_globe!(
                fig,
                sim_var_annual_average,
                obs_var_annual_average,
                true, # plot_bias
                levels_bias,
                p_loc = (fig_row, 3), # plot in the third column
                cmap_extrema = clims_bias,
                mask = viz.oceanmask(),
                more_kwargs = Dict(
                    :mask => ClimaAnalysis.Utils.kwargs(color = :white),
                    :plot => ClimaAnalysis.Utils.kwargs(
                        rasterize = true,
                        # colorrange = clims,
                    ),
                    :axis => ClimaAnalysis.Utils.kwargs(
                        title = bias_title,
                        xticklabelsvisible = false, # don't show lat labels
                        yticklabelsvisible = false, # don't show lon labels
                        xgridvisible = false, # don't show lat grid
                        ygridvisible = false, # don't show lon grid
                        height = 235,
                        ylabelvisible = false,
                    ),
                    :cb => ClimaAnalysis.Utils.kwargs(
                        vertical = false, # horizontal colorbar
                        label = bias_colorbar_label,
                        labelsize = 20,
                        flipaxis = false, # label underneath colorbar
                        height = 15,
                        width = 400, # a little smaller
                        tellwidth = false, # make colorbar width indep of plot width
                        alignmode = Mixed(top = -20),
                    ),
                ),
            )
        end

        ## SEASONAL CYCLE
        if plot_seasonal
            # data_sources.jl has observational data for "gpp", "lwu", and "et" only - maybe separate short_names loop for this
            # Simulation data

            # var_global_average below is a vector of vector, one for each year of simulation, containing monthly global average of var.
            # i represent a year, from 1 to last year
            # for more details on the ClimaAnalysis functions, see ClimaAnalysis docs.

            # ~only compute seasonal cycle for last year so we skip spinup~
            # compute seasonal cycle for second to last year so we skip spinup AND have data for dec after off-by-one correction (shift_to_start_of_previous_month)
            sim_var_global_average =
                ClimaAnalysis.average_lon(
                    ClimaAnalysis.weighted_average_lat(
                        ClimaAnalysis.apply_oceanmask(sim_var_window),
                    ),
                ).data

            # fig_seasonal_cycle = CairoMakie.Figure(size = (600, 400))
            seasonal_title =
                fig_row == 1 ?
                CairoMakie.rich(
                    "ClimaLand vs ERA5 seasonal cycle",
                    fontsize = 28,
                ) : "" # title of the figure

            ax = Axis(
                fig[fig_row, 3],
                title = seasonal_title,
                ylabel = title_stub * " $units_label",
                height = 250,
                xgridvisible = false,
                ygridvisible = false,
                xticks = (
                    1:1:12,
                    [
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
                        "Dev",
                    ],
                ),
            )
            # [
            # plot model output
            CairoMakie.lines!(
                ax,
                sim_var_global_average,
                color = :blue,
                linewidth = 3,
                label = "ClimaLand",
            )

            # Add comparison to observational data (copied from leaderboard.jl)
            # Observational data

            # var_global_average below is a vector of vector, one for each year of simulation, containing monthly global average of var.
            # i represent a year, from 1 to last year
            # for more details on the ClimaAnalysis functions, see ClimaAnalysis docs.
            obs_var_global_average =
                ClimaAnalysis.average_lon(
                    ClimaAnalysis.weighted_average_lat(
                        ClimaAnalysis.apply_oceanmask(obs_var_window),
                    ),
                ).data

            CairoMakie.scatter!(
                ax,
                obs_var_global_average,
                color = :orange,
                label = "ERA5",
            )
            CairoMakie.axislegend(ax, position = :rt)
        end


    end
    save_name = joinpath(root_path, "combined_figures.pdf")
    save_name =
        plot_bias ? joinpath(root_path, "combined_figures_bias.pdf") : save_name
    save_name =
        plot_seasonal ? joinpath(root_path, "combined_figures_seasonal.pdf") :
        save_name
    CairoMakie.save(save_name, fig)
    @show save_name
    return nothing
end

make_paper_figures(
    root_path,
    outdir,
    short_names,
    title_stubs,
    plot_bias = true,
    plot_seasonal = false,
)
