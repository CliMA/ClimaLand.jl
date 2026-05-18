# Copied from experiments/integrated/generic_site/run_generic_site.jl in
# tn/ar/api_singlesite_calibrate

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


function setup_simulation(; site_ID = "US-MOz", duration_days = 7)
    FT = Float64
    local_mode = true

    # local_mode always true
    site_val = FluxnetSimulations.replace_hyphen(site_ID)
    (; time_offset, lat, long) =
        FluxnetSimulations.get_location(FT, Val(site_val))
    (; atmos_h) = FluxnetSimulations.get_fluxtower_height(FT, Val(site_val))
    coord_kwargs = (; lat, long, time_offset, atmos_h)

    @info "Running generic_site_simulation" site_ID duration_days local_mode

    result = FluxnetSimulations.generic_site_simulation(
        site_ID;
        coord_kwargs...,
        duration = Day(duration_days),
    )
    return result.simulation
end
