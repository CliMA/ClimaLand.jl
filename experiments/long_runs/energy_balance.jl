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
# 3057 (prev) -> 3310 (newer)
# outdir = "/scratch/clima/slurm-buildkite/climaland-long-runs/3057/climaland-long-runs/snowy_land_longrun_gpu/global_diagnostics/output_active/" # on clima
# outdir = "/scratch/clima/slurm-buildkite/climaland-long-runs/3310/climaland-long-runs/snowy_land_longrun_gpu/global_diagnostics/output_active/" # on clima
outdir = "snowy_land_longrun_gpu/output_active" # on local
# simdir = ClimaAnalysis.SimDir(outdir)

include("data_paper_plots.jl")

sim_var_dict = get_sim_var_dict(ClimaAnalysis.SimDir(outdir))

var_swd = sim_var_dict["swd"]()
var_swu = sim_var_dict["swu"]()
var_lwd = sim_var_dict["lwd"]()
var_lwu = sim_var_dict["lwu"]()
var_lhf = sim_var_dict["lhf"]()
var_shf = sim_var_dict["shf"]()

# Compute the energy balance
var_net_rad = var_swd - var_swu + var_lwd - var_lwu - var_lhf - var_shf

# Get the year 2 window
i = 2
sim_var_window = ClimaAnalysis.window(
    var_net_rad,
    "time",
    left = (i - 1) * 366 * 86400 + 30 * 86400, # 1 year left of year i, in seconds.
    right = i * 366 * 86400, # 1 year right of year i, in seconds
)

# Plot the energy balance
fig = CairoMakie.Figure(size = (800, 400 ))
fig_row = 1

sim_var_global_average =
    ClimaAnalysis.average_lon(
        ClimaAnalysis.weighted_average_lat(
            ClimaAnalysis.apply_oceanmask(sim_var_window),
        ),
    )

sim_var_time = ClimaAnalysis.average_time(sim_var_global_average)
@show sim_var_time.data

# # fig_seasonal_cycle = CairoMakie.Figure(size = (600, 400))
# seasonal_title =
#     fig_row == 1 ?
#     CairoMakie.rich(
#         "ClimaLand energy balance seasonal cycle",
#         fontsize = 18,
#     ) : "" # title of the figure

# ax = Axis(
#     fig[fig_row, 1],
#     title = seasonal_title,
#     # ylabel = "$units_label",
#     ylabel = "(W/m²)",
#     # CairoMakie.rich(
#     #     title_stub * " $units_label",
#     #     fontsize = 18,
#     # ),
#     # title = CairoMakie.rich(title_stub, fontsize = 18),
#     height = 250,
#     xgridvisible = false,
#     ygridvisible = false,
#     xticks = (
#         1:1:12,
#         [
#             "Jan",
#             "Feb",
#             "Mar",
#             "Apr",
#             "May",
#             "Jun",
#             "Jul",
#             "Aug",
#             "Sep",
#             "Oct",
#             "Nov",
#             "Dev",
#         ],
#     ),
# )
# # [
# # plot model output
# # TODO apply shift_to_start_of_previous_month
# CairoMakie.lines!(
#     ax,
#     sim_var_global_average,
#     color = :blue,#RGBf(0.5, 0.5, 0.5),
#     linewidth = 3,
#     label = "ClimaLand",
# )
# CairoMakie.save(joinpath(root_path, "energy_balance_cl.pdf"), fig)
