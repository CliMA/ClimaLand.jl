# # Regional run of full land model

# The code sets up and runs the soil/canopy/snow model for 6 hours on a small
# region of the globe in Southern California,
# using ERA5 data. In this simulation, we have
# turned lateral flow off because horizontal boundary conditions and the
# land/sea mask are not yet supported by ClimaCore.

# Simulation Setup
# Number of spatial elements: 10x10 in horizontal, 15 in vertical
# Soil depth: 50 m
# Simulation duration: 2 years
# Timestep: 450 s
# Timestepper: ARS111
# Fixed number of iterations: 3
# Jacobian update: every new timestep
# Atmos forcing update: every 3 hours
import ClimaComms
ClimaComms.@import_required_backends
import ClimaTimeSteppers as CTS
using ClimaUtilities.ClimaArtifacts
import ClimaUtilities

import ClimaUtilities.TimeVaryingInputs:
    TimeVaryingInput, LinearInterpolation, PeriodicCalendar
import ClimaUtilities.TimeManager: ITime, date
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
import Interpolations
using CairoMakie, GeoMakie, ClimaAnalysis
import ClimaLand.LandSimVis as LandSimVis

const FT = Float64;
context = ClimaComms.context()
ClimaComms.init(context)
device = ClimaComms.device()
device_suffix = device isa ClimaComms.CPUSingleThreaded ? "cpu" : "gpu"
root_path = "california_longrun_$(device_suffix)"
diagnostics_outdir = joinpath(root_path, "regional_diagnostics")
outdir =
    ClimaUtilities.OutputPathGenerator.generate_output_path(diagnostics_outdir)

function setup_model(FT, context, start_date, Δt, domain, earth_param_set)
    surface_domain = ClimaLand.Domains.obtain_surface_domain(domain)
    surface_space = domain.space.surface

    # Forcing data
    era5_ncdata_path =
        ClimaLand.Artifacts.era5_land_forcing_data2008_path(; context)
    atmos, radiation = ClimaLand.prescribed_forcing_era5(
        era5_ncdata_path,
        surface_space,
        start_date,
        earth_param_set,
        FT;
        max_wind_speed = 25.0,
        time_interpolation_method = LinearInterpolation(PeriodicCalendar()),
        interpolation_method = Interpolations.Linear(),
    )
    forcing = (; atmos, radiation)

    # Read in LAI from MODIS data
    modis_lai_ncdata_path = ClimaLand.Artifacts.modis_lai_multiyear_paths(;
        context = context,
        start_date,
        end_date = stop_date,
    )
    LAI = ClimaLand.prescribed_lai_modis(
        modis_lai_ncdata_path,
        surface_space,
        start_date;
        time_interpolation_method = LinearInterpolation(),
        interpolation_method = Interpolations.Linear(),
    )

    # Overwrite some defaults for the canopy model
    # Energy model
    ac_canopy = FT(2.5e3)
    energy = Canopy.BigLeafEnergyModel{FT}(; ac_canopy)

    # Plant hydraulics
    a = FT(0.2 * 0.0098) # 1/m
    retention_model = Canopy.PlantHydraulics.LinearRetentionCurve{FT}(a)
    hydraulics =
        Canopy.PlantHydraulicsModel{FT}(surface_domain, LAI; retention_model)

    # Roughness lengths
    h_canopy = hydraulics.compartment_surfaces[end]
    z0_m = FT(0.13) * h_canopy
    z0_b = FT(0.1) * z0_m

    ground = ClimaLand.PrognosticGroundConditions{FT}()
    canopy_forcing = (; atmos, radiation, ground)

    canopy = ClimaLand.Canopy.CanopyModel{FT}(
        surface_domain,
        canopy_forcing,
        LAI,
        earth_param_set;
        prognostic_land_components = (:canopy, :snow, :soil, :soilco2),
        energy,
        hydraulics,
        z_0m = z0_m,
        z_0b = z0_b,
    )

    # Snow model setup
    # Set β = 0 in order to regain model without density dependence
    α_snow = Snow.ZenithAngleAlbedoModel(
        FT(0.64),
        FT(0.06),
        FT(2);
        β = FT(0.4),
        x0 = FT(0.2),
    )
    snow = Snow.SnowModel(
        FT,
        surface_domain,
        forcing,
        earth_param_set,
        Δt;
        prognostic_land_components = (:canopy, :snow, :soil, :soilco2),
        α_snow,
    )

    # Construct the land model with all default components except for snow
    land =
        LandModel{FT}(forcing, LAI, earth_param_set, domain, Δt; snow, canopy)
    return land
end

start_date = DateTime(2008)
stop_date = DateTime(2010)
Δt = 450.0
nelements = (10, 10, 15)
@info "Run: Regional Soil-Canopy-Snow Model"
@info "Resolution: $nelements"
@info "Timestep: $Δt s"
@info "Start Date: $start_date"
@info "Stop Date: $stop_date"

radius = FT(6378.1e3)
depth = FT(50)
center_long, center_lat = FT(-117.59736), FT(34.23375)
delta_m = FT(200_000) # in meters
domain = ClimaLand.Domains.HybridBox(;
    xlim = (delta_m, delta_m),
    ylim = (delta_m, delta_m),
    zlim = (-depth, FT(0)),
    nelements = nelements,
    longlat = (center_long, center_lat),
    dz_tuple = FT.((10.0, 0.05)),
)
params = LP.LandParameters(FT)
model = setup_model(FT, context, start_date, Δt, domain, params)
simulation = LandSimulation(start_date, stop_date, Δt, model; outdir)
ClimaLand.Simulations.solve!(simulation)

LandSimVis.make_heatmaps(simulation; savedir = root_path, date = stop_date)
