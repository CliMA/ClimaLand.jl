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


function setup_simulation(; site_ID = "US-MOz", duration_days = 7, space_type = "single")
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
        space_type
    )
    return result.simulation
end


single_simulation = setup_simulation(; site_ID = "US-MOz", duration_days = 7, space_type = "single")
# multi_simulation = setup_simulation(; site_ID = "US-MOz", duration_days = 7, space_type = "multi")

# import ClimaCore
# import ClimaTimeSteppers

# function compare_contents(sim1, sim2)
#     u1 = sim1._integrator.u
#     u2 = sim2._integrator.u
#     p1 = sim1._integrator.p
#     p2 = sim2._integrator.p
#     _compare_contents(u1, u2, [])
#     _compare_contents(p1, p2, [])
#     return nothing
# end

# function _compare_contents(u1, u2, curr_field_names)
#     field_names = propertynames(u1)
#     for field_name in field_names
#         push!(curr_field_names, field_name)
#         _compare_contents(getproperty(u1, field_name), getproperty(u2, field_name), curr_field_names)
#         pop!(curr_field_names)
#     end
# end

# function _compare_contents(u1::Union{ClimaCore.Fields.Field}, u2::Union{ClimaCore.Fields.Field}, field_names)
#     try
#         l2_norm = sum((ClimaCore.Fields.field2array(u1) .- ClimaCore.Fields.field2array(u2)).^2)
#         if !iszero(l2_norm)
#             @info join(field_names, ".") l2_norm
#         end
#     catch
#     end
# end

# compare_contents(single_simulation, multi_simulation)

# for step in 1:100
#     @info step
#     ClimaTimeSteppers.step!(single_simulation._integrator)
#     ClimaTimeSteppers.step!(multi_simulation._integrator)
#     compare_contents(single_simulation, multi_simulation)
# end

solve!(single_simulation)
# solve!(multi_simulation)

# import ClimaCore: Fields
# import JLD2


# for (i, simulation) in enumerate((single_simulation, multi_simulation))
#     varnames = keys(simulation.diagnostics[1].output_writer)

#     var_to_vals = Dict()
#     for varname in varnames
#         diag_vals = values(simulation.diagnostics[1].output_writer[varname])
#         diag_vals = [Fields.field2array(diag_val) for diag_val in diag_vals]
#         diag_vals = vcat(diag_vals...)
#         var_to_vals[varname] = diag_vals
#     end

#     JLD2.save_object("$i.jld2", var_to_vals)
# end
