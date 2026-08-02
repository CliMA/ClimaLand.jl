## Global driver run for the prognostic-carbon equilibrium map (stage 5).
##
## The project goal names *spatial pattern*, which twenty columns cannot test.
## This produces the global driver climatology the offline pool integrator needs,
## on the same 1x1 degree lat-lon grid the ILAMB biomass products use.
##
## The carbon model is deliberately **not** run here. Phase-1 coupling is
## one-way - the pools consume GPP and LAI and feed nothing back - so the drivers
## are identical with and without them, and a CARBON=0 run costs less while
## making rule 1 true by construction rather than by inspection.
##
## Monthly means are sufficient: `check_monthly_drivers.jl` compares
## daily-driven against monthly-driven equilibria at all twenty battery sites and
## finds a median difference of 0.4% and a maximum of 3.3%, against a 3.4x spread
## among the observational products themselves.
##
## Usage (see run_global_driver.pbs):
##   julia --project=.buildkite global_driver.jl

import ClimaComms
ClimaComms.@import_required_backends
using ClimaUtilities.ClimaArtifacts
import ClimaUtilities
import ClimaUtilities.TimeManager: ITime, date
import ClimaParams as CP
using ClimaCore
using ClimaLand
using ClimaLand.Snow
using ClimaLand.Soil
using ClimaLand.Canopy
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.Simulations: LandSimulation, solve!
using Dates

const FT = Float64

# Julia block-buffers stderr when it is redirected to a file, so `@info` from a
# batch job can sit unwritten for the whole run: an empty log then looks
# identical to a hung job. Every stage marker is flushed explicitly, which is
# the only thing that makes this script diagnosable from outside.
stage(msg) = (@info "STAGE $msg"; flush(stderr); flush(stdout))

stage("julia up, packages loaded")

context = ClimaComms.context()
ClimaComms.init(context)
device = ClimaComms.device()

outroot = get(
    ENV,
    "GLOBAL_OUTDIR",
    "/glade/derecho/scratch/arenchon/claude/prognostic_carbon/global_drivers",
)
outdir = ClimaUtilities.OutputPathGenerator.generate_output_path(outroot)

start_date = DateTime(get(ENV, "START", "2008-03-01"))
stop_date = DateTime(get(ENV, "STOP", "2010-03-01"))
Δt = parse(FT, get(ENV, "DT", "900.0"))

# `HybridBox` takes nelements as (nx, ny, nz) where x spans 360 degrees of
# longitude and y spans 180 of latitude, so ClimaLand's shipped default
# (180, 360, 15) is 2 degrees in longitude by 0.5 in latitude - NOT 1x1. The
# default is kept because the other loop's global runs use it, and changing the
# model grid would make these drivers inconsistent with that calibration.
nelem_lon = parse(Int, get(ENV, "NELEM_LON", "180"))
nelem_lat = parse(Int, get(ENV, "NELEM_LAT", "360"))
stage("building domain (nelem lon=$(nelem_lon) lat=$(nelem_lat))")
domain = ClimaLand.Domains.global_box_domain(
    FT;
    context,
    mask_threshold = FT(0.99),
    nelements = (nelem_lon, nelem_lat, 15),
)
stage("domain built; reading parameters")
toml_dict = LP.create_toml_dict(FT)
surface_space = domain.space.surface

stage("setting up ERA5 forcing (regridding can be slow on a cold cache)")

atmos, radiation = ClimaLand.prescribed_forcing_era5(
    start_date,
    stop_date,
    surface_space,
    toml_dict,
    FT;
    max_wind_speed = 25.0,
    context,
)

stage("forcing ready; constructing LandModel")

# No prescribed LAI argument, so the LandModel constructor selects the prognostic
# optimal-LAI model - the same configuration the battery runs under.
model = LandModel{FT}(
    (; atmos, radiation),
    toml_dict,
    domain,
    Δt;
    prognostic_land_components = (:canopy, :lake, :snow, :soil, :soilco2),
)

# Exactly the columns `read_driver_record` consumes, and no more: every extra
# global field is gigabytes of NetCDF for nothing. `tair` is carried because the
# `ct` diagnostic is NaN-masked wherever there is no canopy and needs a fallback.
stage("LandModel constructed; setting up diagnostics")
driver_vars = ["gpp", "lai", "rd", "ct", "tair", "fc3", "pra"]
# Output on a regular 1x1 degree lat-lon grid, which is what the ILAMB biomass
# products use and what ERA5 aligns to. Constructing the writer without
# `horizontal_pts` samples the model's own anisotropic element grid instead,
# which silently produced 2-degree longitude by 0.5-degree latitude output.
out_longs = collect(range(-180.0; length = 360, step = 1.0))
out_lats = collect(range(-90.0; length = 180, step = 1.0))
diagnostics = ClimaLand.default_diagnostics(
    model,
    start_date;
    output_writer = ClimaLand.Diagnostics.NetCDFWriter(
        surface_space,
        outdir;
        start_date,
        horizontal_pts = (out_longs, out_lats),
    ),
    reduction_period = :monthly,
    output_vars = driver_vars,
)

simulation = LandSimulation(start_date, stop_date, Δt, model; diagnostics)

@info "global driver run" start_date stop_date Δt device outdir
@info "resolution" nelements = domain.nelements
@info "driver vars" driver_vars

stage("starting solve")
solve!(simulation)

@info "GLOBAL_DRIVER_DONE" outdir
