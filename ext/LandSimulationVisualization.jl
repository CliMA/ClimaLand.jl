module LandSimulationVisualization
import ClimaDiagnostics
import ClimaAnalysis
import ClimaAnalysis.Visualize as viz
using CairoMakie
import GeoMakie
using Dates
import NCDatasets
using ClimaLand
using StatsBase
using Printf
using Poppler_jll: pdfunite

include("land_sim_vis/plotting_utils.jl")
include("land_sim_vis/leaderboard/data_sources.jl")
include("land_sim_vis/leaderboard/leaderboard.jl")
# TODO: How to handle if diagnostics -> held in memory (NonInterpSavingCallback)

"""
    make_leaderboard_plots(sim::ClimaLand.Simulations.LandSimulation; outdir = ".",  leaderboard_data_sources = ["ERA5", "ILAMB"])


Uses the diagnostic output of the `sim` to create leaderboard plots, comparing the output of the simulation to the ``observations"
from ERA5 or ILAMB. In the future, other observations can be included - see the leadboard directory included here for details.
"""
function make_leaderboard_plots(
    sim::ClimaLand.Simulations.LandSimulation;
    outdir = ".",
    leaderboard_data_sources = ["ERA5", "ILAMB"],
)
    make_leaderboard_plots(
        model,
        ClimaLand.get_domain(model),
        simulation.diagnostics;
        outdir,
        leaderboard_data_sources,
    )
end

"""
    make_leaderboard_plots(m, d, diag; nothing, nothing)

A default method for `make_leaderboard_plots` notifying the user that these plots have not been enable for their domain/model combination.
"""
function make_leaderboard_plots(m, d, diag; nothing, nothing)
    @info "Leaderboard plots not configured for your model type and/or your model domain. Model type must be LandModel, and the domain must be global."
end

"""
    function make_leaderboard_plots(
        model::ClimaLand.Simulations.LandModel,
        domain::ClimaLand.Domains.SphericalShell,
        diagnostics,
        outdir;
        leaderboard_data_sources = ["ERA5", "ILAMB"],
)

Creates the leaderboard plots a global LandModel simulations with diagnostics output
as defined by `diagnostics` using `leaderboard_data_sources` as the source of truth; 
saves the output to files in `outdir`.

The first two arguments are used for dispatch only.
"""
function make_leaderboard_plots(
    model::ClimaLand.Simulations.LandModel,
    domain::ClimaLand.Domains.SphericalShell,
    diagnostics;
    outdir,
    leaderboard_data_sources = ["ERA5", "ILAMB"],
)
    # assert that data spans multiple years and is monthly output?
    # check that the short_names include the appropriate variables for the data source?
    diagdir = first(diagnostics).output_writer.output_dir
    simdir = ClimaAnalysis.SimDir(diagdir)
    short_names = [d.variable.short_name for d in diagnostics]
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

"""
     make_heatmaps(
        sim::ClimaLand.Simulations.LandSimulation;
        outdir = ".",
        short_names = nothing,
        date = nothing,
        levels = nothing,
)

Makes heatmaps using the diagnostics output of the `sim` simulation,
specifically for the list of variables `short_names`, at the date
given by `date`, and at the layers defined by `levels`; the output
plots are saved in `outdir`.

Please note that 
- `date` must be a DateTime, and the output closest to this date will be used 
for plotting
- `short_names` can be a string (single variable), a list, or `nothing`, in which
case all possible variables will be plotted
- 2d variables are unaffected by the choice of `levels`, while 3D variables will
be plotted at each level specified. Note that level = 1 corresponds to the top level,
and that the function will error if you provide a level not included in the output.
Passing levels = nothing defaults to plotting surface values of 3D fields.
"""
function make_heatmaps(
    sim::ClimaLand.Simulations.LandSimulation;
    outdir = ".",
    short_names = nothing,
    date = nothing,
    levels = nothing,
)
    make_heatmaps(
        model,
        ClimaLand.get_domain(model),
        simulation.diagnostics;
        outdir,
        short_names,
        date,
        levels,
    )
end

function make_heatmaps(
    model,
    domain::Union{
        ClimaLand.Domains.SphericalShell,
        ClimaLand.Domains.HybridBox,
        ClimaLand.Domains.Plane,
        ClimaLand.Domains.SphericalSurface,
    },
    diagnostics;
    outdir,
    short_names,
    date,
    levels,
) end


function make_heatmaps(
    model,
    domain::Union{ClimaLand.Domains.Point, ClimaLand.Domains.Column},
    diagnostics;
    outdir,
    short_names,
    date,
    levels,
)
    @info "Heatmaps are not available for simualtions with 1d (column) and 0d (point) domains."
end

"""
     make_annual_timeseries(
        sim::ClimaLand.Simulations.LandSimulation;
        outdir = ".",
        short_names = nothing,
)

Makes timeseries of the domain-averaged variable, 
 using the diagnostics output of the `sim` simulation,
specifically for the list of variables `short_names; the output
plots are saved in `outdir`.

Please note that 
- `short_names` can be a string (single variable), a list, or `nothing`, in which
case all possible variables will be plotted
- The top layer of 3D variables is used for plotting.
"""
function make_annual_timeseries(
    sim::ClimaLand.Simulations.LandSimulation;
    outdir = ".",
    short_names = nothing,
)
    make_annual_timeseries(
        model,
        ClimaLand.get_domain(model),
        simulation.diagnostics;
        outdir,
        short_names,
    )
end

function make_annual_timeseries(model, domain, diagnostics; outdir, short_names)
    diagdir = first(diagnostics).output_writer.output_dir
    simdir = ClimaAnalysis.SimDir(diagdir)
    avail_short_names = [d.variable.short_name for d in diagnostics]
    if short_names isa String
        @assert short_names in avail_short_names
        short_names = [short_names]
    elseif short_names isa Nothing
        short_names = avail_short_names
    end
    make_annual_timeseries(outdir, diagdir, short_names)
end

"""
    check_conservation(sim::ClimaLand.Simulations.LandSimulation; outdir = ".")

Creates a plot which assess conservation of energy and water by the simulation;
the outut is saved in `outdir`.
"""
function check_conservation(
    sim::ClimaLand.Simulations.LandSimulation;
    outdir = ".",
)
    check_conservation(
        model,
        ClimaLand.get_domain(model),
        simulation.diagnostics,
        outdir,
    )
end

function check_conservation(m, d, diags, o)
    @info "Conservation checks not configured for your model type and/or your model domain yet. Model type must be EnergyHydrology, and the domain must be global."
end

function check_conservation(
    model::ClimaLand.Soil.EnergyHydrology,
    domain::ClimaLand.Domains.SphericalShell,
    diagnostics,
    outdir,
)
    diagdir = first(diagnostics).output_writer.output_dir
    check_conservation(outdir, diagdir)

end

end
