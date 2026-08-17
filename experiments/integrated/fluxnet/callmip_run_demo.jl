# Full forward run of the default `LandModel` on a 2-column `ColumnEnsemble`,
# with each column driven by a genuinely different FLUXNET site's CalLMIP
# Phase 1b met forcing (see `callmip_forcing.jl`).
#
# Only the forcing, LAI, tower height, and coordinates vary per column; all
# model parameters and initial conditions are the defaults. Diagnostics are
# disabled (postprocessing is not yet correct for N > 1 columns); the final
# state is inspected per column instead.

import ClimaComms
ClimaComms.@import_required_backends

include(joinpath(@__DIR__, "callmip_forcing.jl"))

FT = Float64
context = ClimaComms.context()
ClimaComms.init(context)
toml_dict = LP.create_toml_dict(FT)

site_IDs = ["DK-Sor", "US-NR1"]
hour_offsets_from_UTC = [1, -7] # UTC = local - offset
met_paths = [ClimaLand.Artifacts.callmip_phase1_forcing_path(ID) for ID in site_IDs]
sites = [
    read_callmip_met(path, offset) for
    (path, offset) in zip(met_paths, hour_offsets_from_UTC)
]

start_date = DateTime(2010, 7, 1)
stop_date = start_date + Day(5)
Δt = 450.0

simulation =
    build_callmip_simulation(sites, start_date, stop_date, Δt, toml_dict, FT)
@info "Running LandModel on a ColumnEnsemble" site_IDs start_date stop_date Δt
@time solve!(simulation)

# Inspect the final state per column; the column index is the last dimension
# of `parent(field)`
function column_values(field, col)
    a = parent(field)
    return vec(Array(selectdim(a, ndims(a), col)))
end

Y = simulation._integrator.u
p = simulation._integrator.p
for (i, ID) in enumerate(site_IDs)
    @info "Final state, column $i ($ID)" T_air = column_values(
        p.drivers.T,
        i,
    )[1] SW_d = column_values(p.drivers.SW_d, i)[1] T_canopy = column_values(
        Y.canopy.energy.T,
        i,
    )[1] ϑ_l_top = column_values(Y.soil.ϑ_l, i)[end] T_soil_top = column_values(
        p.soil.T,
        i,
    )[end]
end
