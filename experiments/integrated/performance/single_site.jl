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

const FT = Float32;

device = ClimaComms.device()
device_suffix = device isa ClimaComms.CPUSingleThreaded ? "cpu" : "gpu"
root_path = "snowy_land_pmodel_longrun_$(device_suffix)"
diagnostics_outdir = joinpath(root_path, "global_diagnostics")
outdir =
    ClimaUtilities.OutputPathGenerator.generate_output_path(diagnostics_outdir)

function setup_model(
    ::Type{FT},
    start_date,
    stop_date,
    Δt,
    domain,
    toml_dict,
) where {FT}
    surface_domain = ClimaLand.Domains.obtain_surface_domain(domain)
    surface_space = domain.space.surface
    # Forcing data - high resolution
    atmos, radiation = ClimaLand.prescribed_forcing_era5(
        start_date,
        stop_date,
        surface_space,
        toml_dict,
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

    ground = ClimaLand.PrognosticGroundConditions{FT}()
    canopy_forcing = (; atmos, radiation, ground)
    prognostic_land_components = (:canopy, :snow, :soil, :soilco2)

    # Construct the P model manually since it is not a default
    photosynthesis = PModel{FT}(domain, toml_dict)
    conductance = PModelConductance{FT}(toml_dict)
    # Use the soil moisture stress function based on soil moisture only
    soil_moisture_stress =
        ClimaLand.Canopy.PiecewiseMoistureStressModel{FT}(domain, toml_dict)
    biomass = ClimaLand.Canopy.PrescribedBiomassModel{FT}(
        domain,
        LAI,
        toml_dict;
        height = ClimaLand.Canopy.clm_canopy_height(
            surface_space;
            max_height = atmos.h * FT(0.9),
        ),
    )
    canopy = ClimaLand.Canopy.CanopyModel{FT}(
        surface_domain,
        canopy_forcing,
        LAI,
        toml_dict;
        prognostic_land_components,
        photosynthesis,
        conductance,
        soil_moisture_stress,
        biomass,
    )

    # Snow model setup
    # Set β = 0 in order to regain model without density dependence
    α_snow = Snow.ZenithAngleAlbedoModel(toml_dict)
    horz_degree_res = FT(1)
    scf = Snow.WuWuSnowCoverFractionModel(toml_dict, horz_degree_res)
    snow = Snow.SnowModel(
        FT,
        surface_domain,
        forcing,
        toml_dict,
        Δt;
        prognostic_land_components,
        α_snow,
        scf,
    )

    # Construct the land model with all default components except for snow
    land = LandModel{FT}(
        forcing,
        LAI,
        toml_dict,
        domain,
        Δt;
        prognostic_land_components,
        snow,
        canopy,
    )
    return land
end
start_date = DateTime("2008-03-01")
stop_date = DateTime("2009-03-01")
Δt = 450.0
longlat = FT.((85,-80));
zlim = FT.((-15,0))
nelements = 15
dz_tuple = FT.((3, 0.05))
domain =
    ClimaLand.Domains.Column(; zlim, longlat, nelements, dz_tuple)
toml_dict = LP.create_toml_dict(FT)

land = setup_model(FT, start_date, stop_date, Δt, domain, toml_dict)
first_save_date = start_date#stop_date - Year(1)
saving_cb = ClimaLand.NonInterpSavingCallback(start_date, stop_date, Second(Δt); first_save_date)
sv = saving_cb.affect!.saved_values;
t0 = ITime(0, Second(1), start_date)
tf = ITime(Second(stop_date - start_date).value, Second(1), start_date)
nancheck_cb = ClimaLand.NaNCheckCallback(
                isnothing(t0.epoch) ? div((tf - t0), 10) : Dates.Month(1),
                t0;
                dt = ITime(Δt))
user_callbacks = (nancheck_cb, saving_cb, ClimaLand.ReportCallback(div((tf - t0), 10), t0))
simulation = LandSimulation(start_date, stop_date, Δt, land; outdir, user_callbacks)
ClimaLand.Simulations.solve!(simulation)

