# # Global run of soil model

# The code sets up and runs the soil model for on a spherical domain,
# using ERA5 data. In this simulation, we have
# turned lateral flow off because horizontal boundary conditions and the
# land/sea mask are not yet supported by ClimaCore.

# Simulation Setup
# Number of spatial elements: 101 in horizontal, 15 in vertical
# Soil depth: 50 m
# Simulation duration: 2-20 years, based on LONGER_RUN setting
# Timestep: 450 s
# Timestepper: ARS111
# Fixed number of iterations: 3
# Jacobian update: every new Newton iteration
# Atmos forcing update: every 3 hours
import ClimaComms
ClimaComms.@import_required_backends
import ClimaTimeSteppers as CTS
using ClimaUtilities.ClimaArtifacts

using ClimaDiagnostics
using ClimaUtilities

import ClimaUtilities.TimeVaryingInputs: LinearInterpolation, PeriodicCalendar
import ClimaUtilities.ClimaArtifacts: @clima_artifact
import ClimaUtilities.TimeManager: ITime
import ClimaParams as CP
using ClimaCore
using ClimaLand
using ClimaLand.Soil
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.Simulations: LandSimulation, solve!
using Dates

using CairoMakie, GeoMakie, ClimaAnalysis
import ClimaLand.LandSimVis as LandSimVis

const FT = Float64;
# If you want to do a very long run locally, you can enter `export
# LONGER_RUN=""` in the terminal and run this script. If you want to do a very
# long run on Buildkite manually, then make a new build and pass `LONGER_RUN=""`
# as an environment variable. In both cases, the value of `LONGER_RUN` does not
# matter.
const LONGER_RUN = haskey(ENV, "LONGER_RUN") ? true : false
context = ClimaComms.context()
ClimaComms.init(context)
device = ClimaComms.device()
device_suffix = device isa ClimaComms.CPUSingleThreaded ? "cpu" : "gpu"
root_path = "soil_longrun_$(device_suffix)"
diagnostics_outdir = joinpath(root_path, "global_diagnostics")
outdir =
    ClimaUtilities.OutputPathGenerator.generate_output_path(diagnostics_outdir)
time_interpolation_method =
    LONGER_RUN ? LinearInterpolation() : LinearInterpolation(PeriodicCalendar())

# If not LONGER_RUN, run for 2 years; note that the forcing from 2008 is repeated.
# If LONGER run, run for 20 years, with the correct forcing each year.
start_date = LONGER_RUN ? DateTime(2000) : DateTime(2008)
stop_date = LONGER_RUN ? DateTime(2020) : DateTime(2010)
Δt = 450.0
nelements = (101, 15)
domain = ClimaLand.Domains.global_domain(FT; context, nelements)
params = LP.LandParameters(FT)
# Forcing data
if LONGER_RUN
    era5_ncdata_path =
        ClimaLand.Artifacts.find_era5_year_paths(start_date, stop_date; context)
else
    era5_ncdata_path =
        ClimaLand.Artifacts.era5_land_forcing_data2008_path(; context)
end
forcing = ClimaLand.prescribed_forcing_era5(
    era5_ncdata_path,
    domain.space.surface,
    start_date,
    params,
    FT;
    max_wind_speed = 25.0,
    time_interpolation_method,
)
model = ClimaLand.Soil.EnergyHydrology{FT}(domain, forcing, params)
diagnostics = ClimaLand.default_diagnostics(
    model,
    start_date;
    output_writer = ClimaDiagnostics.Writers.NetCDFWriter(
        domain.space.subsurface,
        outdir;
        start_date,
    ),
    conservation = true,
    conservation_period = Day(10),
)
simulation =
    LandSimulation(start_date, stop_date, Δt, model; outdir, diagnostics)

@info "Run: Global Soil Model"
@info "Resolution: $nelements"
@info "Timestep: $Δt s"
@info "Start Date: $start_date"
@info "Stop Date: $stop_date"
ClimaLand.Simulations.solve!(simulation)

short_names = ["swc", "sie", "si", "et"]
LandSimVis.make_annual_timeseries(simulation; savedir = root_path, short_names)
LandSimVis.make_heatmaps(
    simulation;
    savedir = root_path,
    date = stop_date,
    short_names,
)
LandSimVis.check_conservation(simulation; savedir = root_path)
