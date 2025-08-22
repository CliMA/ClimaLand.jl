# # Global run of land model

# The code sets up and runs ClimaLand v1, which
# includes soil, canopy, and snow, on a spherical domain,
# using ERA5 data as forcing. In this simulation, we have
# turned lateral flow off because horizontal boundary conditions and the
# land/sea mask are not yet supported by ClimaCore.

# Simulation Setup
# Number of spatial elements: 101 in horizontal, 15 in vertical
# Soil depth: 50 m
# Simulation duration: 730 d
# Timestep: 450 s
# Timestepper: ARS111
# Fixed number of iterations: 3
# Jacobian update: every new Newton iteration
# Atmos forcing update: every 3 hours

import ClimaComms
ClimaComms.@import_required_backends
using ClimaUtilities.ClimaArtifacts
import ClimaUtilities.TimeManager: ITime, date

import ClimaDiagnostics
import ClimaUtilities

import ClimaUtilities.TimeVaryingInputs:
    TimeVaryingInput, LinearInterpolation, PeriodicCalendar
import ClimaUtilities.ClimaArtifacts: @clima_artifact
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
root_path = "snowy_land_pmodel_longrun_$(device_suffix)"
diagnostics_outdir = joinpath(root_path, "global_diagnostics")
outdir =
    ClimaUtilities.OutputPathGenerator.generate_output_path(diagnostics_outdir)

function setup_model(FT, start_date, stop_date, Δt, domain, earth_param_set)
    earth_param_set = LP.LandParameters(toml_dict)
    era5_time_interpolation_method =
        LONGER_RUN ? LinearInterpolation() :
        LinearInterpolation(PeriodicCalendar())
    surface_domain = ClimaLand.Domains.obtain_surface_domain(domain)
    surface_space = domain.space.surface
    # Forcing data
    if LONGER_RUN
        era5_ncdata_path = ClimaLand.Artifacts.find_era5_year_paths(
            start_date,
            stop_date;
            context,
        )
    else
        era5_ncdata_path =
            ClimaLand.Artifacts.era5_land_forcing_data2008_path(; context)
    end
    atmos, radiation = ClimaLand.prescribed_forcing_era5(
        era5_ncdata_path,
        surface_space,
        start_date,
        earth_param_set,
        FT;
        max_wind_speed = 25.0,
        time_interpolation_method = era5_time_interpolation_method,
    )
    forcing = (; atmos, radiation)

    # Read in LAI from MODIS data
    LAI = ClimaLand.prescribed_lai_modis(surface_space, start_date, stop_date)

    # Overwrite some defaults for the canopy model
    # Energy model
    energy_args = Canopy.BigLeafEnergyParameters(toml_dict)
    energy = Canopy.BigLeafEnergyModel{FT}(energy_args)

    # Plant hydraulics
    conductivity_model = Canopy.PlantHydraulics.Weibull(toml_dict)
    retention_model = Canopy.PlantHydraulics.LinearRetentionCurve(toml_dict)
    hydraulics = Canopy.PlantHydraulicsModel{FT}(
        surface_domain,
        LAI,
        toml_dict;
        conductivity_model,
        retention_model,
    )

    # Roughness lengths
    h_canopy = hydraulics.compartment_surfaces[end]
    z_0m = FT(0.13) * h_canopy
    z_0b = FT(0.1) * z_0m

    ground = ClimaLand.PrognosticGroundConditions{FT}()
    canopy_forcing = (; atmos, radiation, ground)

    # Construct the P model manually since it is not a default
    photosynthesis = PModel{FT}()
    conductance = PModelConductance{FT}()

    canopy = ClimaLand.Canopy.CanopyModel{FT}(
        surface_domain,
        canopy_forcing,
        LAI,
        toml_dict;
        prognostic_land_components = (:canopy, :snow, :soil, :soilco2),
        energy,
        hydraulics,
        photosynthesis,
        conductance,
        z_0m,
        z_0b,
    )

    # Snow model setup
    # Set β = 0 in order to regain model without density dependence
    α_snow = Snow.ZenithAngleAlbedoModel(toml_dict)
    horz_degree_res =
        sum(ClimaLand.Domains.average_horizontal_resolution_degrees(domain)) / 2 # mean of resolution in latitude and longitude, in degrees
    scf = Snow.WuWuSnowCoverFractionModel(toml_dict, horz_degree_res)
    snow = Snow.SnowModel(
        FT,
        surface_domain,
        forcing,
        toml_dict,
        Δt;
        prognostic_land_components = (:canopy, :snow, :soil, :soilco2),
        α_snow,
        scf,
    )

    # Construct the land model with all default components except for snow
    land = LandModel{FT}(forcing, LAI, toml_dict, domain, Δt; snow, canopy)
    return land
end

# If not LONGER_RUN, run for 2 years; note that the forcing from 2008 is repeated.
# If LONGER run, run for 19 years, with the correct forcing each year.
# Note that since the Northern hemisphere's winter season is defined as DJF,
# we simulate from and until the beginning of
# March so that a full season is included in seasonal metrics.
start_date = LONGER_RUN ? DateTime("2000-03-01") : DateTime("2008-03-01")
stop_date = LONGER_RUN ? DateTime("2019-03-01") : DateTime("2010-03-01")
Δt = 450.0
nelements = (101, 15)
domain = ClimaLand.Domains.global_domain(
    FT;
    context,
    nelements,
    mask_threshold = FT(0.99),
)
default_params_filepath =
    joinpath(pkgdir(ClimaLand), "toml", "default_parameters.toml")
toml_dict = LP.create_toml_dict(FT, default_params_filepath)
model = setup_model(FT, start_date, stop_date, Δt, domain, toml_dict)

# Construct the PModel callback
pmodel_cb = ClimaLand.make_PModel_callback(FT, start_date, Δt, model.canopy)
nancheck_cb = ClimaLand.NaNCheckCallback(
    Dates.Month(6),
    start_date,
    ITime(Δt, epoch = start_date),
    mask = ClimaLand.Domains.landsea_mask(ClimaLand.get_domain(model)),
)
report_cb = ClimaLand.ReportCallback(1000)
user_callbacks = (pmodel_cb, nancheck_cb, report_cb)

simulation =
    LandSimulation(start_date, stop_date, Δt, model; outdir, user_callbacks)
@info "Run: Global Soil-Canopy-Snow Model"
@info "Resolution: $nelements"
@info "Timestep: $Δt s"
@info "Start Date: $start_date"
@info "Stop Date: $stop_date"
ClimaLand.Simulations.solve!(simulation)

LandSimVis.make_annual_timeseries(simulation; savedir = root_path)
LandSimVis.make_heatmaps(simulation; savedir = root_path, date = stop_date)
LandSimVis.make_leaderboard_plots(simulation; savedir = root_path)
