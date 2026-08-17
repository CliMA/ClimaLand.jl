# Demo: load CalLMIP Phase 1b met forcing for two genuinely different FLUXNET
# sites onto a 2-column `ColumnEnsemble` and evaluate the drivers, one value
# per column. See `callmip_forcing.jl` for the loading helpers and
# `callmip_run_demo.jl` for a full forward run.

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

utc_dates, _ = align_sites(sites)
@info "Shared UTC axis" first(utc_dates) last(utc_dates) length(utc_dates)
start_date = DateTime(2010, 1, 1)
@assert first(utc_dates) <= start_date <= last(utc_dates)

# One column per site, built from the files' own coordinates, in file order
land_domain = ColumnEnsemble(;
    zlim = FT.((-2, 0)),
    nelements = 10,
    longlat = [FT.((s.long, s.lat)) for s in sites],
)
surface_space = land_domain.space.surface

forcing = prescribed_forcing_callmip(
    sites,
    surface_space,
    start_date,
    toml_dict,
    FT,
)
LAI = prescribed_lai_callmip(sites, surface_space, start_date)

# Evaluate a few drivers; each destination Field holds one value per column
dest = ClimaCore.Fields.zeros(surface_space)
for (name, tvi) in (
    ("Tair (K)", forcing.atmos.T),
    ("SWdown (W/m^2)", forcing.radiation.SW_d),
    ("rain (m/s)", forcing.atmos.liquid_precip),
    ("LAI (m^2/m^2)", LAI),
)
    evaluate!(dest, tvi, 0.0)
    at_start = Array(ClimaCore.Fields.field2array(dest))
    evaluate!(dest, tvi, 12 * 3600.0)
    at_noon = Array(ClimaCore.Fields.field2array(dest))
    @info name at_start at_noon
end
@info "Per-column tower height (m)" forcing.atmos.h
