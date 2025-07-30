module LandSimulationVisualizationExt
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

"""
    make_leaderboard_plots(sim::ClimaLand.Simulations.LandSimulation; savedir = ".",  leaderboard_data_sources = ["ERA5", "ILAMB"])


Uses the diagnostic output of the `sim` to create leaderboard plots, comparing the output of the simulation to the ``observations"
from ERA5 or ILAMB. 
"""
function make_leaderboard_plots(
    sim::ClimaLand.Simulations.LandSimulation;
    savedir = ".",
    leaderboard_data_sources = ["ERA5", "ILAMB"],
)
    model = sim.model
    make_leaderboard_plots(
        model,
        ClimaLand.get_domain(model),
        sim.diagnostics;
        savedir,
        leaderboard_data_sources,
    )
end

"""
    make_leaderboard_plots(m, d, diag; nothing, nothing)

A default method for `make_leaderboard_plots` notifying the user that these plots have not been enable for their domain/model combination.
"""
function make_leaderboard_plots(m, d, diag)
    @info "Leaderboard plots not configured for your model type and/or your model domain. Model type must be LandModel, and the domain must be global."
end

"""
    function make_leaderboard_plots(
        model::ClimaLand.LandModel,
        domain::ClimaLand.Domains.SphericalShell,
        diagnostics;
        savedir,
        leaderboard_data_sources = ["ERA5", "ILAMB"],
)

Creates the leaderboard plots a global LandModel simulations with diagnostics output
as defined by `diagnostics` using `leaderboard_data_sources` as the source of truth;
saves the output to files in `savedir`.

The first two arguments are used for dispatch only.
"""
function make_leaderboard_plots(
    model::ClimaLand.LandModel,
    domain::ClimaLand.Domains.SphericalShell,
    diagnostics;
    savedir = ".",
    leaderboard_data_sources = ["ERA5", "ILAMB"],
)
    # assert that data spans multiple years and is monthly output?
    # check that the short_names include the appropriate variables for the data source?
    diagdir = first(diagnostics).output_writer.output_dir
    short_names = [d.variable.short_name for d in diagnostics]
    diagnostics_folder_path = diagdir
    leaderboard_base_path = savedir
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
        savedir = ".",
        short_names = nothing,
        date = nothing,
        levels = nothing,
	plot_name = "figures.pdf"
)

Makes heatmaps using the diagnostics output of the `sim` simulation,
specifically for the list of variables `short_names`, at the date
given by `date`, and at the layers defined by `levels`; the output
plots are saved in `savedir`.

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
    savedir = ".",
    short_names = nothing,
    date = nothing,
    levels = nothing,
    plot_name = "figures.pdf",
)
    model = sim.model
    make_heatmaps(
        model,
        ClimaLand.get_domain(model),
        sim.diagnostics;
        plot_name,
        savedir,
        short_names,
        date,
        levels,
    )
end

function make_heatmaps(
    model::ClimaLand.AbstractModel,
    domain::Union{
        ClimaLand.Domains.SphericalShell,
        ClimaLand.Domains.SphericalSurface,
    },
    diagnostics;
    plot_name = "figures.pdf",
    savedir,
    short_names,
    date,
    levels,
)
    diagdir = first(diagnostics).output_writer.output_dir
    avail_short_names = [d.variable.short_name for d in diagnostics]
    if short_names isa String
        @assert short_names in avail_short_names
        short_names = [short_names]
    elseif short_names isa Nothing
        short_names = avail_short_names
    end
    make_heatmaps(savedir, diagdir, short_names, date; plot_name, levels)
end

function make_heatmaps(
    model::ClimaLand.AbstractModel,
    domain::Union{ClimaLand.Domains.HybridBox, ClimaLand.Domains.Plane},
    diagnostics;
    plot_name = "figures.pdf",
    savedir,
    short_names,
    date,
    levels,
)
    diagdir = first(diagnostics).output_writer.output_dir
    avail_short_names = [d.variable.short_name for d in diagnostics]
    if short_names isa String
        @assert short_names in avail_short_names
        short_names = [short_names]
    elseif short_names isa Nothing
        short_names = avail_short_names
    end
    make_heatmaps(
        savedir,
        diagdir,
        short_names,
        date;
        plot_name,
        levels,
        plot! = viz.heatmap2D!,
        mask = nothing,
        plot_kwargs = Dict(
            :plot => ClimaAnalysis.Utils.kwargs(rasterize = true),
        ),
    )
end


function make_heatmaps(
    model::ClimaLand.AbstractModel,
    domain::Union{ClimaLand.Domains.Point, ClimaLand.Domains.Column},
    diagnostics;
    plot_name,
    savedir,
    short_names,
    date,
    levels,
)
    @info "Heatmaps are not available for simulations with 1d (column) and 0d (point) domains."
end

"""
     make_annual_timeseries(
        sim::ClimaLand.Simulations.LandSimulation;
        savedir = ".",
        short_names = nothing,
	plot_name = "annual_timeseries.pdf",
)

Makes timeseries of the domain-averaged variable,
 using the diagnostics output of the `sim` simulation,
specifically for the list of variables `short_names; the output
plots are saved in `savedir`.

Please note that
- `short_names` can be a string (single variable), a list, or `nothing`, in which
case all possible variables will be plotted
- The top layer of 3D variables is used for plotting.
"""
function make_annual_timeseries(
    sim::ClimaLand.Simulations.LandSimulation;
    savedir = ".",
    short_names = nothing,
    plot_name = "annual_timeseries.pdf",
)
    model = sim.model
    make_annual_timeseries(
        ClimaLand.get_domain(model),
        sim.diagnostics;
        plot_name,
        savedir,
        short_names,
    )
end

function make_annual_timeseries(
    domain::ClimaLand.Domains.AbstractDomain,
    diagnostics;
    plot_name = "annual_timeseries.pdf",
    savedir = ".",
    short_names = nothing,
)
    @info "No method matching make_annual_timeseries for $domain."
end

function make_annual_timeseries(
    domain::Union{
        ClimaLand.Domains.SphericalShell,
        ClimaLand.Domains.SphericalSurface,
    },
    diagnostics;
    plot_name = "annual_timeseries.pdf",
    savedir = ".",
    short_names = nothing,
)
    diagdir = first(diagnostics).output_writer.output_dir
    avail_short_names = [d.variable.short_name for d in diagnostics]
    if short_names isa String
        @assert short_names in avail_short_names
        short_names = [short_names]
    elseif short_names isa Nothing
        short_names = avail_short_names
    end
    make_ocean_masked_annual_timeseries(
        savedir,
        diagdir,
        short_names;
        plot_name,
    )
end

"""
     make_diurnal_timeseries(
        sim::ClimaLand.Simulations.LandSimulation;
        savedir = ".",
        short_names = nothing,
	plot_name = "diurnal_timeseries.pdf",
        comparison_data = nothing,
        spinup_date = sim.start_date
)

Computes and plots the average diurnal cycle of the diagnostic
output of the simulation `sim`,
specifically for the list of variables `short_names; the output
plots are saved in `savedir` under the name `plot_name`. 
Only data after the spinup_date is considered.

The comparison data can optionally be provided as a NamedTuple; the data must 
be labeled with a key equal to the same name as in `short_names`, 
or no comparison will be plotted. The value at each key is the timeseries of
the variable. The timestamp of each is passed as another key `UTC_datetime`,
with values equal to the time at which the observations where made.

Please note that
- `short_names` can be a string (single variable), a list, or `nothing`, in which
case all possible variables will be plotted
- The top layer of 3D variables is used for plotting.
"""
function make_diurnal_timeseries(
    sim::ClimaLand.Simulations.LandSimulation;
    savedir = ".",
    short_names = nothing,
    plot_name = "diurnal_timeseries.pdf",
    comparison_data = nothing,
    spinup_date = sim.start_date,
)
    model = sim.model
    make_diurnal_timeseries(
        ClimaLand.get_domain(model),
        sim.diagnostics,
        sim.start_date;
        plot_name,
        savedir,
        short_names,
        comparison_data,
        spinup_date,
    )
end

function make_diurnal_timeseries(
    domain::ClimaLand.Domains.AbstractDomain,
    diagnostics,
    start_date;
    plot_name = "diurnal_timeseries.pdf",
    savedir = ".",
    short_names = nothing,
    comparison_data = nothing,
    spinup_date = start_date,
)
    @info "No method matching make_diurnal_timeseries for $domain."
end

function make_diurnal_timeseries(
    domain::Union{ClimaLand.Domains.Column, ClimaLand.Domains.Point},
    diagnostics,
    start_date;
    plot_name = "diurnal_timeseries.pdf",
    savedir = ".",
    short_names = nothing,
    comparison_data = nothing,
    spinup_date = start_date,
)
    avail_short_names = [d.variable.short_name for d in diagnostics]
    if short_names isa String
        @assert short_names in avail_short_names
        short_names = [short_names]
    elseif short_names isa Nothing
        short_names = avail_short_names
    end
    diag_ids = [findfirst(sn .== avail_short_names) for sn in short_names]
    make_diurnal_timeseries(
        savedir,
        diagnostics[diag_ids], # To be consistent with global/regional runs, this should be a directory with the saved diagnostics.
        start_date;
        plot_name,
        comparison_data,
        spinup_date,
    )
end

"""
     make_timeseries(
        sim::ClimaLand.Simulations.LandSimulation;
        savedir = ".",
        short_names = nothing,
	plot_name = "variable_timeseries.pdf",
        comparison_data = nothing,
        spinup_date = sim.start_date
)

Computes and plots the timeseries of the diagnostic
output of the simulation `sim`,
specifically for the list of variables `short_names; the output
plots are saved in `savedir` under the name `plot_name`. 
Only data after the spinup_date is considered.

The comparison data can optionally be provided as a NamedTuple; the data must 
be labeled with a key equal to the same name as in `short_names`, 
or no comparison will be plotted. The value at each key is the timeseries of
the variable. The timestamp of each is passed as another key `UTC_datetime`,
with values equal to the time at which the observations where made.

Please note that
- `short_names` can be a string (single variable), a list, or `nothing`, in which
case all possible variables will be plotted
- The top layer of 3D variables is used for plotting.
"""
function make_timeseries(
    sim::ClimaLand.Simulations.LandSimulation;
    savedir = ".",
    short_names = nothing,
    plot_name = "variable_timeseries.pdf",
    comparison_data = nothing,
    spinup_date = sim.start_date,
)
    model = sim.model
    make_timeseries(
        ClimaLand.get_domain(model),
        sim.diagnostics,
        sim.start_date;
        plot_name,
        savedir,
        short_names,
        comparison_data,
        spinup_date,
    )
end

function make_timeseries(
    domain::ClimaLand.Domains.AbstractDomain,
    diagnostics,
    start_date;
    plot_name = "variable_timeseries.pdf",
    savedir = ".",
    short_names = nothing,
    comparison_data = nothing,
    spinup_date = start_date,
)
    @info "No method matching make_timeseries for $domain."
end

function make_timeseries(
    domain::Union{ClimaLand.Domains.Column, ClimaLand.Domains.Point},
    diagnostics,
    start_date;
    plot_name = "variable_timeseries.pdf",
    savedir = ".",
    short_names = nothing,
    comparison_data = nothing,
    spinup_date = start_date,
)
    avail_short_names = [d.variable.short_name for d in diagnostics]
    if short_names isa String
        @assert short_names in avail_short_names
        short_names = [short_names]
    elseif short_names isa Nothing
        short_names = avail_short_names
    end
    diag_ids = [findfirst(sn .== avail_short_names) for sn in short_names]
    make_timeseries(
        savedir,
        diagnostics[diag_ids],# To be consistent with global/regional runs, this should be a directory with the saved diagnostics.
        start_date;
        plot_name,
        comparison_data,
        spinup_date,
    )
end

"""
    check_conservation(sim::ClimaLand.Simulations.LandSimulation; savedir = ".", plot_name = "conservation_figures.pdf")

Creates a plot which assess conservation of energy and water by the simulation;
the outut is saved in `savedir`.
"""
function check_conservation(
    sim::ClimaLand.Simulations.LandSimulation;
    savedir = ".",
    plot_name = "conservation_figures.pdf",
)
    model = sim.model
    check_conservation(
        model,
        ClimaLand.get_domain(model),
        sim.diagnostics,
        savedir,
        plot_name,
    )
end

function check_conservation(m, d, diags, o, pn)
    @info "Conservation checks not configured for your model type and/or your model domain yet. Model type must be EnergyHydrology, and the domain must be global."
end

function check_conservation(
    model::ClimaLand.Soil.EnergyHydrology,
    domain::ClimaLand.Domains.SphericalShell,
    diagnostics,
    savedir,
    plot_name,
)
    diagdir = first(diagnostics).output_writer.output_dir
    check_conservation(savedir, diagdir; plot_name)

end

end
