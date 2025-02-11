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
# outdir = "/scratch/clima/slurm-buildkite/climaland-long-runs/3057/climaland-long-runs/snowy_land_longrun_gpu/global_diagnostics/output_active/" # on clima
outdir = "snowy_land_longrun_gpu/global_diagnostics/output_active" # on local
# simdir = ClimaAnalysis.SimDir(outdir)

short_names = ["lhf"]
units_labels = Dict("lhf" => "(J/m²s)")
sim_var_units_labels = Dict("lhf" => "(J/m²s)")
title_stubs = Dict("lhf" => "Latent heat flux")

function get_sim_var_dict(simdir)
    # Dict for loading in simulation data
    sim_var_dict = Dict{String, Any}()

    # convert ET to LHF
    earth_param_set = LP.LandParameters(Float64)
    _LH_v0 = LP.LH_v0(earth_param_set) # J/kg
    sim_var_dict["lhf"] =
        () -> begin
            sim_var_et = get(simdir, short_name = "et") # units (kg/m²s)

            attribs = Dict(
                "short_name" => "lhf",
                "long_name" => "Latent heat flux",
                "units" => "J/m²s",
                "start_date" => sim_var_et.attributes["start_date"],
            )
            sim_var_lhf = ClimaAnalysis.OutputVar(
                attribs,
                sim_var_et.dims,
                sim_var_et.dim_attributes,
                sim_var_et.data * _LH_v0, # J/m²s
            )
            # Shift to the first day and subtract one month as preprocessing
            # sim_var_lhf =
            #     ClimaAnalysis.shift_to_start_of_previous_month(sim_var_lhf)
            return sim_var_lhf
        end

    return sim_var_dict
end

function get_obs_var_dict()
    era5_data_path = joinpath(
        ClimaLand.Artifacts.era5_monthly_averages_2008_folder_path(),
        "era5_monthly_surface_fluxes_200801-200812.nc",
    )
    # era5_monthly_surface_fluxes_200801-200812.nc
    # era5_data_path = joinpath(ClimaLand.Artifacts.era5_monthly_averages_2008_folder_path(lowres = true), "era5_2008_1.0x1.0_lowres.nc")
    # Dict for loading in observational data
    obs_var_dict = Dict{String, Any}()
    obs_var_dict["lhf"] =
        (start_date) -> begin

            obs_var = ClimaAnalysis.OutputVar(
                era5_data_path,
                "mslhf",
                new_start_date = start_date,
                shift_by = Dates.firstdayofmonth,
            ) # units (W/m² = J/m²s)
            obs_var = ClimaAnalysis.replace(obs_var, missing => NaN)

            # Add attributes and make values positive to account for convention difference
            attribs = Dict(
                "short_name" => "lhf",
                "long_name" => "Latent heat flux",
                "units" => "J/m²s",
                "start_date" => start_date,
            )
            obs_var_pos = ClimaAnalysis.OutputVar(
                attribs,
                obs_var.dims,
                obs_var.dim_attributes,
                -1 * obs_var.data,
            )
            # obs_var_pos = ClimaAnalysis.convert_units(
            #     obs_var,
            #     "J/m²s",
            #     conversion_function = (x) -> -x,
            # )
            return obs_var_pos
        end
    return obs_var_dict
end

function get_mask_dict()
    # Dict for loading in masks
    mask_dict = Dict{String, Any}()

    mask_dict["lhf"] =
        (sim_var, obs_var) -> begin
            return ClimaAnalysis.make_lonlat_mask(
                ClimaAnalysis.slice(
                    obs_var,
                    time = ClimaAnalysis.times(obs_var) |> first,
                );
                set_to_val = isnan,
            )
        end
    return mask_dict
end



function make_paper_figures(
    root_path,
    outdir,
    short_names,
    units_labels,
    sim_var_units_labels,
    title_stubs,
)
    # Set up for comparison to data (copied from leaderboard.jl)
    # use sim_var and obs_var together for the seasonal plot because they already have the same units :)
    sim_var_dict = get_sim_var_dict(ClimaAnalysis.SimDir(outdir))
    obs_var_dict = get_obs_var_dict()
    # Set up dict for storing simulation and observational data after processing
    sim_obs_comparsion_dict = Dict()
    mask_dict = get_mask_dict()

    for short_name in [short_names[1]]
        title_stub = title_stubs[short_name]
        units_label = units_labels[short_name]
        sim_var_units_label = sim_var_units_labels[short_name]

        # Access simulation data in the time we want
        sim_var = sim_var_dict[short_name]()
        sim_var_times = [ClimaAnalysis.times(sim_var)[1]]
        kwarg_z = ClimaAnalysis.has_altitude(sim_var) ? Dict(:z => 1) : Dict() # if has altitude, take first layer
        sim_var_sliced = ClimaAnalysis.slice(sim_var; kwarg_z...)
        i = 1 # use first year of simulation
        sim_var_window = ClimaAnalysis.window(
            sim_var_sliced,
            "time",
            left = (i - 1) * 366 * 86400 + 30 * 86400, # 1 year left of year i, in seconds.
            right = i * 366 * 86400, # 1 year right of year i, in seconds
        )

        # Access observation data in the time we want
        obs_var = obs_var_dict[short_name](sim_var.attributes["start_date"])
        obs_var_sliced = ClimaAnalysis.slice(obs_var; kwarg_z...)
        obs_var_window = ClimaAnalysis.window(
            obs_var_sliced,
            "time",
            left = 0, # observation data starts at 2008
            right = 366 * 86400, # 1 year of observation data, in seconds
        )

        # create figure for all 3 plots
        fig = CairoMakie.Figure(size = (1800, 400))

        ## GLOBAL HEATMAP
        # sim_var_annual_average below contains annually-averaged data in the provided window
        sim_var_annual_average = ClimaAnalysis.average_time(sim_var_window)

        # obs_var_annual_average below contains annually-averaged data in the provided window
        obs_var_annual_average = ClimaAnalysis.average_time(obs_var_window)

        t = sim_var_times[1]

        # Get colorbar limits to share between heatmaps
        sim_extrema = extrema(sim_var_annual_average.data)
        obs_extrema = extrema(obs_var_annual_average.data)
        clims = extrema(vcat(sim_extrema..., obs_extrema...))
        # clims = (0, 250) # hardcode for now
        clims = (clims[1], floor(clims[2], digits = -1)) # round to nearest 10
        # @show count(sim_var_annual_average.data .> 250)
        # @show count(obs_var_annual_average.data .> 250)
        # @show count(sim_var_annual_average.data .> 200)
        # @show count(obs_var_annual_average.data .> 200)
        @show clims
        # @assert false

        # Plot simulation data heatmap
        sim_title =
            CairoMakie.rich("ClimaLand, annually averaged", fontsize = 18) # title of the figure
        # fig_global_sim = CairoMakie.Figure(size = (600, 400))
        round_step(x, step) = round(x / step) * step
        viz.heatmap2D_on_globe!(
            fig,
            sim_var_annual_average,
            p_loc = (1, 1), # plot in the first column
            plot_colorbar = true,
            colorbar_label = "Latent heat flux $(units_labels[short_name])",
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
                ),
                # :cb => ClimaAnalysis.Utils.kwargs(
                #     colorrange = clims,
                #     ticks = 0:50:round_step(clims[2], 50),
                # ),
            ),
        )
        Makie.Colorbar(
            fig[2, 1],
            # plot,
            label = "Latent heat flux $(units_labels[short_name])",
            vertical = false, # horizontal colorbar
            flipaxis = false, # label underneath colorbar
            width = 400, # a little smaller
            tellwidth = false, # make colorbar width indep of plot width
            limits = clims, # use global min/max across sim and obs plots
            # cb_kwargs...,
        )
        # CairoMakie.save(joinpath(root_path, "$(short_name)_$(t)-annual_avg.pdf"), fig_global_sim)

        # Plot observational data heatmap
        obs_title = CairoMakie.rich("ERA5, annually averaged", fontsize = 18) # title of the figure
        # fig_global_obs = CairoMakie.Figure(size = (600, 400))
        viz.heatmap2D_on_globe!(
            fig,
            obs_var_annual_average,
            p_loc = (1, 2), # plot in the second column
            plot_colorbar = true,
            colorbar_label = "Latent heat flux $(units_labels[short_name])",
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
                ),
                # :cb => ClimaAnalysis.Utils.kwargs(
                #     colorrange = clims,
                #     ticks = 0:50:round_step(clims[2], 50),
                # ),
            ),
        )
        Makie.Colorbar(
            fig[2, 2],
            # plot,
            label = "Latent heat flux $(units_labels[short_name])",
            vertical = false, # horizontal colorbar
            flipaxis = false, # label underneath colorbar
            width = 400, # a little smaller
            tellwidth = false, # make colorbar width indep of plot width
            limits = clims, # use global min/max across sim and obs plots
            # cb_kwargs...,
        )
        # CairoMakie.save(joinpath(root_path, "$(short_name)_$(t)-annual_avg-obs.pdf"), fig_global_obs)


        ## SEASONAL CYCLE
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
        ax = Axis(
            fig[1, 3],
            title = CairoMakie.rich(
                "ClimaLand vs ERA5 seasonal cycle",
                fontsize = 18,
            ), # title of the figure
            # ylabel = "$sim_var_units_label",
            ylabel = title_stub * " $sim_var_units_label",
            # CairoMakie.rich(
            #     title_stub * " $sim_var_units_label",
            #     fontsize = 18,
            # ),
            # title = CairoMakie.rich(title_stub, fontsize = 18),
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


        CairoMakie.save(joinpath(root_path, "$(short_name)_figures.pdf"), fig)
    end

    # figures = readdir(root_path, join = true)
    # pdfunite() do unite
    #     run(Cmd([unite, figures..., joinpath(root_path, "figures.pdf")]))
    # end
    return nothing
end

make_paper_figures(
    root_path,
    outdir,
    short_names,
    units_labels,
    sim_var_units_labels,
    title_stubs,
)
