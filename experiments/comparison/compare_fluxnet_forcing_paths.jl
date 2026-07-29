# Compare the two forcing-loading paths through a REAL fluxnet simulation:
#
#   old : FLUXNET CSV -> 0-D `TimeVaryingInput(times, values)`
#   new : FLUXNET-style NetCDF -> `DataSource` -> `MultiColumnDataHandler` ->
#         `TimeVaryingInput`  (requires a `PointCloudSpace`, so a `ColumnEnsemble`)
#
# Unlike `compare_forcing_paths.jl`, which drives a bare `EnergyHydrology` column for half a
# day, this is the `run_fluxnet.jl` model: an integrated `LandModel` of soil (with root
# extraction and surface runoff), soil CO2 biogeochemistry, a canopy (two-stream radiative
# transfer, Medlyn conductance, Farquhar photosynthesis, plant hydraulics, big-leaf energy,
# prescribed biomass with MODIS LAI), and snow, started from the site's observed initial
# conditions. That matters for two reasons: the soil's top boundary condition is mediated by
# the canopy rather than bare, and `c_co2` is actually consumed (by photosynthesis) instead
# of merely sitting in the driver cache.
#
# The question this answers is not whether the two paths return the same values -- that is
# settled, and they agree to within 1 ulp of linear-interpolation tie-breaking. It is
# whether a 1-ulp forcing difference STAYS small once it is fed through canopy hydraulics,
# stomatal conductance, and snow accumulation and melt over many steps. So the run is long
# enough to reach the site's first precipitation event, which is also the only way the
# composed `Precip`/`VPD` forcing and the snow model are exercised at all.
#
# The two builds are identical apart from the domain and the forcing object: the same
# parameters, vertical mesh, LAI, and initial conditions. A single-column `ColumnEnsemble`
# reproduces the equivalent `Column` to 1 ulp (see `compare_fluxnet.jl`), so a difference
# here is attributable to the forcing path.
#
# Expect this to take a few minutes: it is ~2000 implicit steps of the full model, twice.
#
# Run: julia --project=.buildkite experiments/comparison/compare_fluxnet_forcing_paths.jl [SITE_ID]
#      (SITE_ID defaults to US-MOz; must be a site with a CSV and a ClimaLand site file)

import ClimaComms
ClimaComms.@import_required_backends
using ClimaCore
using Dates
using DelimitedFiles # loads the FluxnetSimulations extension

using ClimaLand
using ClimaLand.Domains: Column, ColumnEnsemble
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.Simulations: step!
import ClimaLand.FluxnetSimulations as FluxnetSimulations

include(joinpath(@__DIR__, "field_rmse.jl"))
include(joinpath(@__DIR__, "fluxnet_model.jl"))
include(joinpath(@__DIR__, "fluxnet_netcdf.jl"))

const FT = Float64
const SITE_ID = isempty(ARGS) ? "US-MOz" : ARGS[1]

toml_dict = LP.create_toml_dict(FT)
earth_param_set = LP.LandParameters(toml_dict)
thermo_params = LP.thermodynamic_parameters(earth_param_set)

# Site defaults: vertical mesh, location, tower height, timestep, LAI, and all soil/plant
# parameters, shared by both builds.
site = fluxnet_site_setup(FT, SITE_ID)
(; dt, dz_tuple, nelements, zmin, zmax, time_offset, lat, long, atmos_h) = site
(; start_date) = site

# ---------------------------------------------------------------------------
# The two forcings, from the same site data.
# ---------------------------------------------------------------------------
csv_forcing = FluxnetSimulations.prescribed_forcing_fluxnet(
    SITE_ID,
    lat,
    long,
    time_offset,
    atmos_h,
    start_date,
    toml_dict,
    FT,
)

column_domain = Column(;
    zlim = (zmin, zmax),
    nelements = nelements,
    dz_tuple = dz_tuple,
    longlat = (long, lat),
)
ensemble_domain = ColumnEnsemble(;
    zlim = (zmin, zmax),
    nelements = nelements,
    dz_tuple = dz_tuple,
    longlat = (long, lat),
)

met_nc_path = joinpath(mkpath(joinpath(@__DIR__, "output")), "$(SITE_ID)_Met.nc")
isfile(met_nc_path) && rm(met_nc_path)
(; first_precip_date) = write_fluxnet_netcdf(
    met_nc_path,
    SITE_ID,
    lat,
    long,
    time_offset,
    thermo_params,
    FT,
)
@info "Generated NetCDF from CSV" met_nc_path SITE_ID time_offset lat long start_date first_precip_date

nc_forcing = FluxnetSimulations.prescribed_forcing_netcdf(
    met_nc_path,
    ensemble_domain.space.surface,
    atmos_h,
    start_date,
    toml_dict,
    FT;
    hour_offset_from_UTC = time_offset,
)

# ---------------------------------------------------------------------------
# The two builds. `build_fluxnet_sim` (fluxnet_model.jl) is the run_fluxnet.jl model
# parameterized by domain and forcing, so these differ in nothing else.
# ---------------------------------------------------------------------------
sim_csv = build_fluxnet_sim(column_domain, csv_forcing, site, toml_dict)
# sim_nc = build_fluxnet_sim(ensemble_domain, nc_forcing, site, toml_dict)

# ---------------------------------------------------------------------------
# Step both to a series of checkpoints, so any growth of the difference over time is
# visible rather than only its final value.
# ---------------------------------------------------------------------------
# report_state_diffs(sim_csv, sim_nc; col = 1, label = "t = 0 (no steps)")

# precip_seconds =
#     isnothing(first_precip_date) ? nothing :
#     FT(Dates.value(Second(first_precip_date - start_date)))
# isnothing(precip_seconds) && @warn "Site records no precipitation; the composed \
#                                     Precip/VPD forcing and the snow model stay inert."

# checkpoints = if isnothing(precip_seconds)
#     [("1 day", FT(86400)), ("2 days", FT(2 * 86400))]
# else
#     [
#         ("1 day", FT(86400)),
#         ("just before first precip", precip_seconds - 3600),
#         ("at first precip", precip_seconds),
#         ("6 h after first precip", precip_seconds + 6 * 3600),
#     ]
# end

# # A run in which precipitation never became nonzero would report agreement on the composed
# # `Precip`/`VPD` forcing without ever having evaluated it. The drivers hold only the current
# # step's value, which is zero again within hours of a precipitation event, so the maximum is
# # accumulated over every step rather than read off the cache at the end.
# driver_max(sim, name) =
#     maximum(abs, parent(getproperty(sim._integrator.p.drivers, name)))
# precip_maxima = Dict{Symbol, FT}(
#     :P_liq_csv => 0,
#     :P_snow_csv => 0,
#     :P_liq_nc => 0,
#     :P_snow_nc => 0,
# )
# function update_precip_maxima!(maxima, sim_csv, sim_nc)
#     maxima[:P_liq_csv] = max(maxima[:P_liq_csv], driver_max(sim_csv, :P_liq))
#     maxima[:P_snow_csv] = max(maxima[:P_snow_csv], driver_max(sim_csv, :P_snow))
#     maxima[:P_liq_nc] = max(maxima[:P_liq_nc], driver_max(sim_nc, :P_liq))
#     maxima[:P_snow_nc] = max(maxima[:P_snow_nc], driver_max(sim_nc, :P_snow))
#     return nothing
# end
# update_precip_maxima!(precip_maxima, sim_csv, sim_nc)

# let elapsed = FT(0)
#     for (label, target) in checkpoints
#         nsteps = round(Int, (target - elapsed) / dt)
#         for _ in 1:nsteps
#             step!(sim_csv)
#             step!(sim_nc)
#             update_precip_maxima!(precip_maxima, sim_csv, sim_nc)
#         end
#         elapsed += nsteps * dt
#         @assert sum(isnan.(sim_csv._integrator.u)) == 0
#         @assert sum(isnan.(sim_nc._integrator.u)) == 0
#         report_state_diffs(
#             sim_csv,
#             sim_nc;
#             col = 1,
#             label = "$label (t = $elapsed s, $(round(Int, elapsed / dt)) steps)",
#         )
#     end
# end

# @info "Maximum precipitation drivers over every step (m/s)" P_liq_csv =
#     precip_maxima[:P_liq_csv] P_snow_csv = precip_maxima[:P_snow_csv] P_liq_nc =
#     precip_maxima[:P_liq_nc] P_snow_nc = precip_maxima[:P_snow_nc]
