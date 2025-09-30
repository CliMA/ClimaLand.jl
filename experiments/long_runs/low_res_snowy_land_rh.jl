# # Global run of land model at low resolution

# The code sets up and runs ClimaLand v1, which
# includes soil, canopy, and snow, on a spherical domain,
# using low resolution ERA5 data as forcing.

# Simulation Setup
# Number of spatial elements: 30 in horizontal, 15 in vertical
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
context = ClimaComms.context()
ClimaComms.init(context)
device = ClimaComms.device()
device_suffix = device isa ClimaComms.CPUSingleThreaded ? "cpu" : "gpu"
root_path = "lowres_snowy_land_longrun_rh_$(device_suffix)"
diagnostics_outdir = joinpath(root_path, "global_diagnostics")
outdir =
    ClimaUtilities.OutputPathGenerator.generate_output_path(diagnostics_outdir)

function setup_model(FT, start_date, stop_date, Δt, domain, toml_dict)
    earth_param_set = LP.LandParameters(toml_dict)
    surface_domain = ClimaLand.Domains.obtain_surface_domain(domain)
    surface_space = domain.space.surface
    # Forcing data
    era5_ncdata_path =
        ClimaLand.Artifacts.era5_land_forcing_data2008_lowres_path(; context)
    atmos, radiation = ClimaLand.prescribed_perturbed_rh_era5(
        era5_ncdata_path,
        surface_space,
        start_date,
        earth_param_set,
        -0.2,
        FT;
        max_wind_speed = 25.0,
    )
    forcing = (; atmos, radiation)

    # Read in LAI from MODIS data
    LAI = ClimaLand.Canopy.prescribed_lai_modis(
        surface_space,
        start_date,
        stop_date,
    )

    # Overwrite some defaults for the canopy model
    # Plant hydraulics
    retention_model = Canopy.PlantHydraulics.LinearRetentionCurve(toml_dict)
    hydraulics = Canopy.PlantHydraulicsModel{FT}(
        surface_domain,
        toml_dict;
        retention_model,
    )

    # Roughness lengths
    h_canopy = hydraulics.compartment_surfaces[end]
    z_0m = FT(0.13) * h_canopy
    z_0b = FT(0.1) * z_0m

    ground = ClimaLand.PrognosticGroundConditions{FT}()
    canopy_forcing = (; atmos, radiation, ground)

    canopy = ClimaLand.Canopy.CanopyModel{FT}(
        surface_domain,
        canopy_forcing,
        LAI,
        toml_dict;
        prognostic_land_components = (:canopy, :snow, :soil, :soilco2),
        hydraulics,
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
# Note that since the Northern hemisphere's winter season is defined as DJF,
# we simulate from and until the beginning of
# March so that a full season is included in seasonal metrics.
start_date = DateTime("2008-03-01")
stop_date = DateTime("2010-03-02")
Δt = 450.0
nelements = (30, 15)
domain = ClimaLand.Domains.global_domain(
    FT;
    context,
    nelements,
    mask_threshold = FT(0.99),
)
toml_dict = LP.create_toml_dict(FT)
model = setup_model(FT, start_date, stop_date, Δt, domain, toml_dict)
user_callbacks = (
    ClimaLand.NaNCheckCallback(
        Dates.Month(6),
        ITime(0, epoch = start_date);
        dt = ITime(Δt, epoch = start_date),
        mask = ClimaLand.Domains.landsea_mask(ClimaLand.get_domain(model)),
    ),
    ClimaLand.ReportCallback(Month(1), ITime(0, epoch = start_date)),
)
simulation =
    LandSimulation(start_date, stop_date, Δt, model; user_callbacks, outdir)
@info "Run: Global Soil-Canopy-Snow Model"
@info "Resolution: $nelements"
@info "Timestep: $Δt s"
@info "Start Date: $start_date"
@info "Stop Date: $stop_date"
CP.log_parameter_information(toml_dict, joinpath(root_path, "parameters.toml"))
ClimaLand.Simulations.solve!(simulation)

LandSimVis.make_annual_timeseries(simulation; savedir = root_path)
LandSimVis.make_heatmaps(simulation; savedir = root_path, date = stop_date)
