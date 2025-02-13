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

root_path = joinpath(pwd(), "snowy_land_longrun_gpu")
!isdir(root_path) && mkdir(root_path)
# 3057 (prev) -> 3310 (newer) -> 3363 (sc = 0)
# outdir = "/scratch/clima/slurm-buildkite/climaland-long-runs/3057/climaland-long-runs/snowy_land_longrun_gpu/global_diagnostics/output_active/" # on clima
# outdir = "/scratch/clima/slurm-buildkite/climaland-long-runs/3310/climaland-long-runs/snowy_land_longrun_gpu/global_diagnostics/output_active/" # on clima
# outdir = "/scratch/clima/slurm-buildkite/climaland-long-runs/3363/climaland-long-runs/snowy_land_longrun_gpu/global_diagnostics/output_active/" # on clima
outdir = "snowy_land_longrun_gpu/output_active" # on local
# simdir = ClimaAnalysis.SimDir(outdir)

short_names = ["lhf", "shf", "lwu", "swu"]
# units_labels = Dict(
#     "lhf" => "(W/m²)",
#     "shf" => "(W/m²)",
#     "lwu" => "(W/m²)",
#     "swu" => "(W/m²)",
# )
title_stubs = Dict(
    "lhf" => "Latent heat flux",
    "shf" => "Sensible heat flux",
    "lwu" => "Upward longwave radiation",
    "swu" => "Upward shortwave radiation",
)

include("data_paper_plots.jl")

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
    # fig = CairoMakie.Figure(size = (1800, 400)) # use for single row plotting
    fig = CairoMakie.Figure(size = (1800, 400 * length(short_names)))

    for (idx, short_name) in enumerate(short_names)
        # Plot figures in odd rows, colorbars in even rows
        fig_row = (idx - 1) * 2 + 1

        title_stub = title_stubs[short_name]

        # Access simulation data in the time we want
        sim_var = sim_var_dict[short_name]()
        sim_var_times = [ClimaAnalysis.times(sim_var)[1]]
        kwarg_z = ClimaAnalysis.has_altitude(sim_var) ? Dict(:z => 1) : Dict() # if has altitude, take first layer
        sim_var_sliced = ClimaAnalysis.slice(sim_var; kwarg_z...)
        i = 2 # use second year of simulation
        sim_var_window = ClimaAnalysis.window(
            sim_var_sliced,
            "time",
            left = (i - 1) * 366 * 86400 + 30 * 86400, # 1 year left of year i, in seconds.
            right = i * 366 * 86400, # 1 year right of year i, in seconds
        )
        units_label = "(" * sim_var.attributes["units"] * ")"

        # Access observation data in the time we want
        obs_var = obs_var_dict[short_name](sim_var.attributes["start_date"])
        # ClimaAnalysis.center_longitude!(obs_var, 0) # TODO center lon at 0 to prevent https://github.com/MakieOrg/GeoMakie.jl/issues/290

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

        # obs_var_annual_average below contains annually-averaged data in the provided window
        obs_var_annual_average = ClimaAnalysis.average_time(obs_var_window)

        t = sim_var_times[1]

        # Get colorbar limits to share between heatmaps
        sim_var_no_nans = deepcopy(sim_var_annual_average.data)
        sim_var_no_nans[isnan.(sim_var_no_nans)] .= 0 # TODO this is a hacky way to remove NaNs in data
        sim_extrema = extrema(sim_var_no_nans)
        obs_extrema = extrema(obs_var_annual_average.data)
        # we MUST pass the same clims to the sim plot, obs plot, AND colorbars
        clims = extrema(vcat(sim_extrema..., obs_extrema...))
        clims = (clims[1], floor(clims[2], digits = -1)) # round to nearest 10

        # Plot simulation data heatmap
        sim_title =
            fig_row == 1 ?
            CairoMakie.rich("ClimaLand, annually averaged", fontsize = 18) : "" # title of the figure

        # fig_global_sim = CairoMakie.Figure(size = (600, 400))
        # round_step(x, step) = round(x / step) * step
        viz.heatmap2D_on_globe!(
            fig,
            sim_var_annual_average,
            p_loc = (fig_row, 1), # plot in the first column
            plot_colorbar = false,
            # colorbar_label = "$(sim_var.attributes["long_name"]) $units_label",
            mask = viz.oceanmask(),
            more_kwargs = Dict(
                :mask => ClimaAnalysis.Utils.kwargs(color = :white),
                :plot => ClimaAnalysis.Utils.kwargs(
                    rasterize = true,
                    colorrange = clims,
                ),
                :axis => ClimaAnalysis.Utils.kwargs(
                    title = sim_title,
                    xticklabelsvisible = false, # don't show lat labels
                    yticklabelsvisible = false, # don't show lon labels
                    xgridvisible = false, # don't show lat grid
                    ygridvisible = false, # don't show lon grid
                    height = 235,
                ),
                # :cb => ClimaAnalysis.Utils.kwargs(
                #     colorrange = clims,
                #     ticks = 0:50:round_step(clims[2], 50),
                # ),
            ),
        )
        Makie.Colorbar(
            fig[fig_row + 1, 1],
            # plot,
            label = "$(sim_var.attributes["long_name"]) $units_label",
            vertical = false, # horizontal colorbar
            flipaxis = false, # label underneath colorbar
            height = 15,
            width = 400, # a little smaller
            tellwidth = false, # make colorbar width indep of plot width
            colorrange = clims, # use global min/max across sim and obs plots
            # cb_kwargs...,
        )
        # CairoMakie.save(joinpath(root_path, "$(short_name)_$(t)-annual_avg.pdf"), fig_global_sim)

        # Plot observational data heatmap
        obs_title =
            fig_row == 1 ?
            CairoMakie.rich("ERA5, annually averaged", fontsize = 18) : "" # title of the figure

        # fig_global_obs = CairoMakie.Figure(size = (600, 400))
        viz.heatmap2D_on_globe!(
            fig,
            obs_var_annual_average,
            p_loc = (fig_row, 2), # plot in the second column
            plot_colorbar = false,
            # colorbar_label = "$(sim_var.attributes["long_name"]) $units_label",
            mask = viz.oceanmask(),
            more_kwargs = Dict(
                :mask => ClimaAnalysis.Utils.kwargs(color = :white),
                :plot => ClimaAnalysis.Utils.kwargs(
                    rasterize = true,
                    colorrange = clims,
                ),
                :axis => ClimaAnalysis.Utils.kwargs(
                    title = obs_title,
                    xticklabelsvisible = false, # don't show lat labels
                    yticklabelsvisible = false, # don't show lon labels
                    xgridvisible = false, # don't show lat grid
                    ygridvisible = false, # don't show lon grid
                    height = 235,
                ),
                # :cb => ClimaAnalysis.Utils.kwargs(
                #     colorrange = clims,
                #     ticks = 0:50:round_step(clims[2], 50),
                # ),
            ),
        )
        Makie.Colorbar(
            fig[fig_row + 1, 2],
            # plot,
            label = "$(sim_var.attributes["long_name"]) $units_label",
            vertical = false, # horizontal colorbar
            flipaxis = false, # label underneath colorbar
            height = 15,
            width = 400, # a little smaller
            tellwidth = false, # make colorbar width indep of plot width
            colorrange = clims, # use global min/max across sim and obs plots
            # cb_kwargs...,
        )
        # Makie.colgap!(fig.layout, 0, fig_row) # reduce gap between plot/colorbar pairs
        # CairoMakie.save(joinpath(root_path, "$(short_name)_$(t)-annual_avg-obs.pdf"), fig_global_obs)

        @assert !(plot_bias && plot_seasonal) # only one of these can be true
        ## BIAS PLOT
        if plot_bias
            bias_title =
                fig_row == 1 ?
                CairoMakie.rich("ClimaLand vs ERA5 bias", fontsize = 18) : "" # title of the figure
            bias_colorbar_label = "$(sim_var.attributes["long_name"]): ClimaLand - ERA5 bias $units_label"
            viz.plot_bias_on_globe!(
                fig,
                sim_var_annual_average,
                obs_var_annual_average,
                p_loc = (fig_row, 3), # plot in the third column
                # cmap_extrema = (-50, 50),#compare_vars_biases_plot_extrema[short_name],
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
                    ),
                    :cb => ClimaAnalysis.Utils.kwargs(
                        vertical = false, # horizontal colorbar
                        label = bias_colorbar_label,
                        flipaxis = false, # label underneath colorbar
                        height = 15,
                        width = 400, # a little smaller
                        tellwidth = false, # make colorbar width indep of plot width
                        # colorrange = clims,
                        # ticks = 0:50:round_step(clims[2], 50),
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
                    fontsize = 18,
                ) : "" # title of the figure

            ax = Axis(
                fig[fig_row, 3],
                title = seasonal_title,
                # ylabel = "$units_label",
                ylabel = title_stub * " $units_label",
                # CairoMakie.rich(
                #     title_stub * " $units_label",
                #     fontsize = 18,
                # ),
                # title = CairoMakie.rich(title_stub, fontsize = 18),
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
            # TODO apply shift_to_start_of_previous_month
            CairoMakie.lines!(
                ax,
                sim_var_global_average,
                color = :blue,#RGBf(0.5, 0.5, 0.5),
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

            # CairoMakie.save(
            #     joinpath(root_path, "$(short_name)_global_monthly.pdf"),
            #     fig_seasonal_cycle,
            # )
        end


    end
    CairoMakie.save(joinpath(root_path, "combined_figures_bias.pdf"), fig)
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
