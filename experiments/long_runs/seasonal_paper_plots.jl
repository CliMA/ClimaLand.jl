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
root_path = outdir

short_names = ["lhf", "shf", "lwu", "swu"]
title_stubs = Dict(
    "lhf" => "LE",#"Latent heat flux",
    "shf" => "H",#"Sensible heat flux",
    "lwu" => "LW_u",#"Upward longwave radiation",
    "swu" => "SW_u",#"Upward shortwave radiation",
)

include("data_paper_plots.jl")


function make_seasonal_cycle_figure(root_path, outdir, short_names, title_stubs)
    # Set up for comparison to data (copied from leaderboard.jl)
    # use sim_var and obs_var together for the seasonal plot because they already have the same units :)
    sim_var_dict = get_sim_var_dict(ClimaAnalysis.SimDir(outdir))
    obs_var_dict = get_obs_var_dict()

    # create figure for all plots
    # fig = CairoMakie.Figure(size = (1800, 400)) # use for single column plotting
    num_cols = 2
    num_rows = 2
    fig = CairoMakie.Figure(size = (550num_cols, 350num_rows))

    n_vars = length(short_names)
    for (idx, short_name) in enumerate(short_names)
        # Plot figures in odd rows, colorbars in even rows for 1x4 grid plotting
        # fig_row = (idx - 1) * 2 + 1

        # Determine row and column index for 2x2 grid plotting
        row_idx = mod(idx, 2) + 1 # odd var idx in even rows, offset by 1 for title
        col_idx = idx <= fld(n_vars, 2) ? 1 : 2 # first half of vars in first column (floor division)

        title_stub = title_stubs[short_name]

        # Access simulation data in the time we want
        sim_var = sim_var_dict[short_name]()
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

        obs_var_sliced = ClimaAnalysis.slice(obs_var; kwarg_z...)
        obs_var_window = ClimaAnalysis.window(
            obs_var_sliced,
            "time",
            left = 0, # observation data starts at 2008
            right = 366 * 86400, # 1 year of observation data, in seconds
        )

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

        ax = Axis(
            fig[row_idx, col_idx],
            ylabel = title_stub * " $units_label",
            ylabelsize = 17,
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

        # plot model output
        CairoMakie.lines!(
            ax,
            sim_var_global_average,
            color = :blue,
            linewidth = 3,
            label = "ClimaLand",
        )

        # Add comparison to observational data (copied from leaderboard.jl)
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

    save_name = joinpath(root_path, "seasonal_cycle.pdf")
    CairoMakie.save(save_name, fig)
    @show save_name
    return nothing
end

make_seasonal_cycle_figure(root_path, outdir, short_names, title_stubs)
