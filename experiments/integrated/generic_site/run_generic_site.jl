# Run the default ClimaLand model at any FLUXNET2015 site and plot it against
# the site observations. Parameters are the spatially varying defaults at the
# site coordinates and LAI is the MODIS climatology, so no site-specific
# configuration is needed.
#
# Usage:
#   julia --project=.buildkite experiments/integrated/generic_site/run_generic_site.jl US-MOz
#
# US-MOz, US-NR1, US-Ha1 and US-Var use the downloadable `fluxnet_sites`
# artifact. Any other site needs the `fluxnet2015` artifact, which is too large
# to be downloadable; see
# https://github.com/CliMA/ClimaArtifacts/tree/main/fluxnet2015 for how to
# download it from fluxnet.org and set it up.

import ClimaComms
ClimaComms.@import_required_backends
using Dates
using ClimaDiagnostics
using ClimaLand
using ClimaLand.Domains: Column
using ClimaLand.Simulations: LandSimulation, solve!
import ClimaLand.Parameters as LP
import ClimaLand.FluxnetSimulations as FluxnetSimulations
import ClimaLand.LandSimVis as LandSimVis
using CairoMakie, ClimaAnalysis, GeoMakie, Printf, Statistics

const FT = Float64
toml_dict = LP.create_toml_dict(FT)

site_ID = length(ARGS) >= 1 ? ARGS[1] : "US-MOz"
site_ID_val = FluxnetSimulations.replace_hyphen(site_ID)
duration = Day(7)

(; time_offset, lat, long) =
    FluxnetSimulations.get_location(FT, Val(site_ID_val))
(; atmos_h) = FluxnetSimulations.get_fluxtower_height(FT, Val(site_ID_val))

# Start at the first row where all forcing variables are observed
(start_date, stop_date) = FluxnetSimulations.get_data_dates(
    site_ID,
    time_offset;
    duration,
    required_columns = FluxnetSimulations.FLUXNET_FORCING_COLUMNS,
)
data_dt = FluxnetSimulations.get_data_dt(site_ID)

# Same timestep and soil layers as the global runs
Δt = 900.0
domain = Column(;
    zlim = (FT(-15), FT(0)),
    nelements = 15,
    dz_tuple = (FT(3), FT(0.05)),
    longlat = (long, lat),
)

forcing = FluxnetSimulations.prescribed_forcing_fluxnet(
    site_ID,
    lat,
    long,
    time_offset,
    atmos_h,
    start_date,
    toml_dict,
    FT,
)
# FLUXNET2015 records start before MODIS, so use the MODIS climatology
LAI = ClimaLand.Canopy.prescribed_climatological_lai_modis(domain.space.surface)

prognostic_land_components = (:canopy, :snow, :soil, :soilco2)
land = LandModel{FT}(
    forcing,
    LAI,
    toml_dict,
    domain,
    Δt;
    prognostic_land_components,
)

set_ic! = FluxnetSimulations.make_set_fluxnet_initial_conditions(
    site_ID,
    start_date,
    time_offset,
    land,
)
diags = ClimaLand.default_diagnostics(
    land,
    start_date;
    output_writer = ClimaDiagnostics.Writers.DictWriter(),
    output_vars = [
        "gpp",
        "er",
        "lhf",
        "shf",
        "swu",
        "lwu",
        "swc",
        "tsoil",
        "swe",
    ],
    reduction_period = data_dt == 3600 ? :hourly : :halfhourly,
)
simulation = LandSimulation(
    start_date,
    stop_date,
    Δt,
    land;
    set_ic!,
    updateat = Second(data_dt),
    diagnostics = diags,
)
@time solve!(simulation)

savedir = joinpath(@__DIR__, "out", site_ID)
mkpath(savedir)
comparison_data = FluxnetSimulations.get_comparison_data(site_ID, time_offset)
LandSimVis.make_diurnal_timeseries(
    domain,
    diags,
    start_date;
    savedir,
    short_names = ["gpp", "er", "lhf", "shf", "swu", "lwu"],
    spinup_date = start_date + Day(1),
    comparison_data,
)
LandSimVis.make_timeseries(
    domain,
    diags,
    start_date;
    savedir,
    short_names = ["swc", "tsoil", "swe"],
    spinup_date = start_date + Day(1),
    comparison_data,
)
@info "Wrote plots to $savedir"
