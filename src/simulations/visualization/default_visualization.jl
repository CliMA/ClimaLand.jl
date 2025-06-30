import ClimaDiagnostics
import ClimaAnalysis
import ClimaAnalysis.Visualize as viz
using CairoMakie
import GeoMakie
using Dates
import NCDatasets

using Poppler_jll: pdfunite
include("plotting_utils.jl")
include("leaderboard/data_sources.jl")
include("leaderboard/leaderboard.jl")

function make_plots(sim::LandSimulation; outdir = ".", kwargs)
    make_plots(
        model,
        ClimaLand.get_domain(model),
        simulation.diagnostics;
        outdir,
        kwargs...,
    )
end

# How to handle if diagnostics -> held in memory (NonInterpSavingCallback)
# Split up figures function into time series and heatmaps
# Make a plot_var function that takes in time and level, use that
# Pass in short names to plot?
function make_plots(
    model::LandModel,
    domain::SphericalShell,
    diagnostics;
    outdir,
    leaderboard_data_sources = ["ERA5", "ILAMB"],
    make_figures = true,
)
    diagdir = first(diagnostics).output_writer.output_dir
    simdir = ClimaAnalysis.SimDir(diagdir)
    short_names = [d.variable.short_name for d in diagnostics]
    make_figures(outdir, diagdir, short_names)
    diagnostics_folder_path = diagdir
    leaderboard_base_path = outdir
    for data_source in leaderboard_data_sources
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
end

function make_plots(
    model::EnergyHydrology,
    domain::SphericalShelll,
    diagnostics;
    outdir,
    check_conservation = false,
)
    diagdir = first(diagnostics).output_writer.output_dir
    short_names = [d.variable.short_name for d in diagnostics]
    make_figures(outdir, diagdir, short_names)
    ## Conservation
    if check_conservation
        check_conservation(outdir, diagdir)
    end

end
