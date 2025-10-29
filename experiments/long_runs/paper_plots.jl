###########################################################################
## Global plots comparing long run outputs to ERA5 data
# This script generates global heatmaps of ClimaLand outputs,
# ERA5 data, and their difference (bias) for 4 variables:
# latent heat flux (lhf or LE), sensible heat flux (shf or H),
# net longwave radiation (LWn), and net shortwave radiation (SWn).
#
## Usage:
# Run this script after running the snowy_land_pmodel_longrun_gpu experiment.
# The buildkite run number of the long run experiment is specified at the
# top of this script as `longrun_number`. The experiment output is expected
# to be located in the default output directory on clima, but this can be
# changed by modifying the `outdir` variable.
# To generate the plots once you have the long run data, simply run:
# `julia --project=.buildkite experiments/long_runs/paper_plots.jl`
#
# Note that this script requires a specific branch of ClimaAnalysis,
# which can be found under the tag `climaland-v1-paper-plots-tag`
# in that repository.
###########################################################################

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
using GeoMakie
using Dates
using Statistics

# Specify the longrun buildkite number
longrun_number = 4281

root_path = joinpath(pwd(), "snowy_land_longrun_gpu-$longrun_number")
!isdir(root_path) && mkdir(root_path)

outdir = "snowy_land_longrun_gpu-$longrun_number" # local
outdir = "/scratch/clima/slurm-buildkite/climaland-long-runs/$longrun_number/climaland-long-runs/snowy_land_pmodel_longrun_gpu/global_diagnostics/output_active" # on clima

short_names = ["lhf", "shf", "lwn", "swn"]#, "lwu", "swu"]
title_stubs = Dict(
    "lhf" => "LE",#"Latent heat flux",
    "shf" => "H",#"Sensible heat flux",
    "lwu" => "LWᵤ",#"Upward longwave radiation",
    "swu" => "SWᵤ",#"Upward shortwave radiation",
    "lwn" => "LWₙ",#"Net longwave radiation",
    "swn" => "SWₙ",#"Net shortwave radiation",
)
# Define levels for contour colorbars
levels_dict = Dict(
    "lhf_bias" => collect(-30:5:10),
    "shf_bias" => collect(-10:5:30),
    "lwu_bias" => collect(-40:5:20),
    "swu_bias" => collect(-20:5:40),
    "lwn_bias" => collect(-20:5:20),
    "swn_bias" => collect(-20:5:20),
)

include("data_paper_plots.jl")


function make_paper_figures(
    root_path,
    outdir,
    short_names,
    title_stubs;
    plot_bias = false,
)
    # We want to run comparison for years 2017-2020
    # Simulation was run from 2000-03-01 to 2020-03-01.
    # Output is saved to the first day of the next month, then shifted back to the previous month
    #  so for the last output date (2020-03-01 shifted to 2020-03-01) the data is for February 2020.
    #  so we compare to data from 2017-03-01 to 2020-02-01 (inclusive)
    comparison_start_date = DateTime("2017-03-01")
    comparison_end_date = DateTime("2020-02-01")

    # Set up for comparison to data (copied from leaderboard.jl)
    # use sim_var and obs_var together for the seasonal plot because they already have the same units :)
    sim_var_dict = get_sim_var_dict(ClimaAnalysis.SimDir(outdir))
    obs_var_dict = get_obs_var_dict(comparison_start_date, comparison_end_date)

    # create figure for all plots
    num_cols = plot_bias ? 3 : 2
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
        @assert DateTime(sim_var.attributes["start_date"]) <
                comparison_start_date

        obs_var = obs_var_dict[short_name](sim_var.attributes["start_date"])
        obs_var_times = ClimaAnalysis.times(obs_var)
        @assert DateTime(obs_var.attributes["start_date"]) <
                comparison_start_date

        sim_var_window = ClimaAnalysis.window(
            sim_var,
            "time",
            left = comparison_start_date,
            right = comparison_end_date,
        )
        units_label = "(" * sim_var.attributes["units"] * ")"

        # Access observation data in the time we want
        obs_var_window = ClimaAnalysis.window(
            obs_var,
            "time",
            left = comparison_start_date,
            right = comparison_end_date,
        )

        # Shift longitude from [0, 360] to [-180, 180]
        obs_var_window =
            ClimaAnalysis.shift_longitude(obs_var_window, -180.0, 180.0)
        # Resample obs_var to sim_var grid
        obs_var_window =
            ClimaAnalysis.resampled_as(obs_var_window, sim_var_window)

        start_date = DateTime(sim_var.attributes["start_date"])
        @assert start_date + Second(ClimaAnalysis.times(sim_var_window)[1]) ==
                comparison_start_date
        @assert start_date + Second(ClimaAnalysis.times(obs_var_window)[1]) ==
                comparison_start_date
        @assert start_date + Second(ClimaAnalysis.times(sim_var_window)[end]) ==
                comparison_end_date
        @assert start_date + Second(ClimaAnalysis.times(obs_var_window)[end]) ==
                comparison_end_date


        ## GLOBAL HEATMAP
        # sim_var_annual_average below contains annually-averaged data in the provided window
        sim_var_annual_average = ClimaAnalysis.average_time(sim_var_window)

        # obs_var_annual_average below contains annually-averaged data in the provided window
        obs_var_annual_average = ClimaAnalysis.average_time(obs_var_window)

        # Display RMSE and bias values between sim and obs
        rmse = ClimaAnalysis.global_rmse(
            sim_var_annual_average,
            obs_var_annual_average,
            mask = ClimaAnalysis.apply_oceanmask,
        )

        @show "$short_name RMSE: $rmse "

        bias = ClimaAnalysis.global_bias(
            sim_var_annual_average,
            obs_var_annual_average,
            mask = ClimaAnalysis.apply_oceanmask,
        )
        @show "$short_name Bias: $bias"


        t = sim_var_times[1]

        # Get colorbar limits to share between heatmaps
        sim_var_no_nans = deepcopy(sim_var_annual_average.data)
        sim_var_no_nans[isnan.(sim_var_no_nans)] .= 0 # TODO this is a hacky way to remove NaNs in data

        combined_data =
            vcat(vec(sim_var_no_nans), vec(obs_var_annual_average.data))
        min = Int(round(quantile(combined_data, 0.02); digits = -1))
        max = Int(round(quantile(combined_data, 0.98); digits = -1))
        diff = max - min
        step_size = diff < 150 ? 10 : 20

        # Increase max to be a multiple of step_size
        max =
            (diff % step_size != 0) ? max + step_size - (diff % step_size) : max
        clims = (min, max)
        levels = collect(min:step_size:max)
        nlevels = length(levels)

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
                        ylabel = "$title_stub $units_label", # plot variable label on y-axis for leftmost column (sim)
                        ylabelvisible = true,
                        ylabelpadding = 0,
                    ),
                    :cb => ClimaAnalysis.Utils.kwargs(
                        vertical = false, # horizontal colorbar
                        label = "$title_stub $units_label",
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
                # ylabel = "$title_stub $units_label", # plot variable label on y-axis for leftmost column (sim)
                p_loc = (fig_row, 1), # plot in the first column
                plot_colorbar = false,
                mask = viz.oceanmask(),
                more_kwargs = Dict(
                    :mask => ClimaAnalysis.Utils.kwargs(color = :white),
                    :plot => ClimaAnalysis.Utils.kwargs(
                        rasterize = true,
                        # colorrange = clims,
                        colormap = colormap,
                        # levels = levels, TODO add these back
                        # extendhigh = :auto,
                        # extendlow = :auto,
                    ),
                    :axis => ClimaAnalysis.Utils.kwargs(
                        title = sim_title,
                        xticklabelsvisible = false, # don't show lat labels
                        yticklabelsvisible = false, # don't show lon labels
                        xgridvisible = false, # don't show lat grid
                        ygridvisible = false, # don't show lon grid
                        height = 235,
                        ylabel = "$title_stub $units_label", # plot variable label on y-axis for leftmost column (sim)
                        ylabelpadding = -5,
                        ylabelvisible = true,
                        # limits = clims,
                    ),
                ),
            )
            Makie.Colorbar(
                fig[fig_row + 1, 1],
                label = "$title_stub $units_label",
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
                # colorbar_label = "$title_stub $units_label",
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
                        label = "$title_stub $units_label",
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
                        # levels = levels, TODO add these back
                        # extendhigh = :auto,
                        # extendlow = :auto,
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
                label = "$title_stub $units_label",
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

        ## BIAS PLOT
        if plot_bias
            levels_bias = levels_dict["$(short_name)_bias"]
            clims_bias = (levels_bias[1], levels_bias[end])

            bias_title =
                fig_row == 1 ?
                CairoMakie.rich("ClimaLand vs ERA5 bias", fontsize = 28) : "" # title of the figure
            bias_colorbar_label = "$title_stub $units_label"
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

    end
    save_name = joinpath(root_path, "combined_figures.png")
    save_name =
        plot_bias ?
        joinpath(root_path, "combined_figures_bias_2017-2020-shifted.png") :
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
)
