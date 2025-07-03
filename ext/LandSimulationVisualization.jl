module LandSimulationVisualization
import ClimaLand
import ClimaDiagnostics
import ClimaAnalysis
import ClimaAnalysis.Visualize as viz
using CairoMakie
import GeoMakie
using Dates
import NCDatasets
using Printf
import StatsBase: mean

using Poppler_jll: pdfunite

include("land_sim_vis/plotting_utils.jl")
include("land_sim_vis/leaderboard/data_sources.jl")
include("land_sim_vis/leaderboard/leaderboard.jl")
# How to handle if diagnostics -> held in memory (NonInterpSavingCallback)

function make_leaderboard_plots(sim::LandSimulation; outdir = ".", kwargs)
    make_leaderboard_plots(
        model,
        ClimaLand.get_domain(model),
        simulation.diagnostics,
        outdir;
        kwargs...,
    )
end

function make_leaderboard_plots(m, d, diag, outdir; _...)
    @info "Leaderboard plots not configured for your model type and/or your model domain. Model type must be LandModel, and the domain must be global."
end


function make_leaderboard_plots(
    model::LandModel,
    domain::SphericalShell,
    diagnostics,
    outdir;
    leaderboard_data_sources = ["ERA5", "ILAMB"],
)
    # assert that data spans multiple years and is monthly output?
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

function make_heatmaps(
    sim::LandSimulation;
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
    domain::Union{SphericalShell, HybridBox},
    diagnostics;
    outdir,
    short_names,
    date,
    levels,
)
    diagdir = first(diagnostics).output_writer.output_dir
    simdir = ClimaAnalysis.SimDir(diagdir)
    avail_short_names = [d.variable.short_name for d in diagnostics]
    if short_names isa String
        @assert short_names in avail_short_names
        short_names = [short_name]
    elseif short_names isa Nothing
        short_names = avail_short_names
    end
    make_heatmaps(outdir, diagdir, short_names, date; levels)
end

function make_annual_timeseries(
    sim::LandSimulation;
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


function check_conservation(sim::LandSimulation; outdir = ".")
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
    model::EnergyHydrology,
    domain::SphericalShell,
    diagnostics,
    outdir,
)
    diagdir = first(diagnostics).output_writer.output_dir
    check_conservation(outdir, diagdir)

end



end
