# Run a single-site ClimaLand simulation using FluxnetSimulations.generic_site_simulation.
#
# Usage:
#   julia --project=experiments experiments/integrated/generic_site/run_generic_site.jl US-MOz
#   julia --project=experiments experiments/integrated/generic_site/run_generic_site.jl US-MOz --local
#   SIM_DURATION_DAYS=14 julia --project=experiments .../run_generic_site.jl DE-Tha
#
# Without `--local`, coordinates are auto-resolved from the FLUXNET2015 metadata
# CSV (requires the `fluxnet2015` artifact, HPC-only). With `--local`, we pull
# coordinates from the per-site Val{} dispatchers — works for the four bundled
# sites (US-MOz, US-NR1, US-Ha1, US-Var) without the FLUXNET2015 artifact.

import ClimaComms
ClimaComms.@import_required_backends
using Dates
using Printf
using Statistics
using ClimaDiagnostics
using ClimaUtilities
using ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.FluxnetSimulations as FluxnetSimulations
using ClimaLand.Simulations: solve!
import ClimaLand.LandSimVis as LandSimVis
using CairoMakie
using ClimaAnalysis
using GeoMakie

site_ID = length(ARGS) >= 1 ? ARGS[1] : "US-MOz"
local_mode = "--local" in ARGS
duration_days = parse(Int, get(ENV, "SIM_DURATION_DAYS", "7"))

const FT = Float64

if local_mode
    site_val = FluxnetSimulations.replace_hyphen(site_ID)
    (; time_offset, lat, long) =
        FluxnetSimulations.get_location(FT, Val(site_val))
    (; atmos_h) = FluxnetSimulations.get_fluxtower_height(FT, Val(site_val))
    coord_kwargs = (; lat, long, time_offset, atmos_h)
else
    coord_kwargs = (;)
end

@info "Running generic_site_simulation" site_ID duration_days local_mode

result = FluxnetSimulations.generic_site_simulation(
    site_ID;
    coord_kwargs...,
    duration = Day(duration_days),
)
@time solve!(result.simulation)

# `time_offset` for comparison-data fetch: reuse what we already resolved if
# possible, otherwise look up via metadata.
time_offset_for_comp = if local_mode
    coord_kwargs.time_offset
else
    FluxnetSimulations.get_site_info(site_ID).time_offset
end

savedir = joinpath(
    pkgdir(ClimaLand),
    "experiments",
    "integrated",
    "generic_site",
    "out",
    site_ID,
)
mkpath(savedir)

comparison_data =
    FluxnetSimulations.get_comparison_data(site_ID, time_offset_for_comp)

LandSimVis.make_diurnal_timeseries(
    result.land_domain,
    result.diags,
    result.start_date;
    savedir,
    short_names = ["gpp", "lhf", "shf", "swu", "lwu"],
    spinup_date = result.start_date + Day(1),
    comparison_data,
)
LandSimVis.make_timeseries(
    result.land_domain,
    result.diags,
    result.start_date;
    savedir,
    short_names = ["swc", "tsoil", "swe"],
    spinup_date = result.start_date + Day(1),
    comparison_data,
)

@info "Wrote plots to $savedir"
