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

ClimaComms.@import_required_backends
const FT = Float64;
context = ClimaComms.context()
ClimaComms.init(context)
device = ClimaComms.device()
device_suffix = device isa ClimaComms.CPUSingleThreaded ? "cpu" : "gpu"
root_path = "snowy_land_pmodel_desert_site"
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
        context,
        use_lowres_forcing = true,
    )
    forcing = (; atmos, radiation)

    # Read in LAI from MODIS data
    LAI = ClimaLand.Canopy.prescribed_lai_modis(
        surface_space,
        start_date,
        stop_date,
    )

    prognostic_land_components = (:canopy, :snow, :soil, :soilco2)
    # Construct the land model with all default components except for snow
    land = LandModel{FT}(
        forcing,
        LAI,
        toml_dict,
        domain,
        Δt;
        prognostic_land_components,
    )
    return land
end

start_date = DateTime("2000-09-01")
stop_date = DateTime("2001-09-01")
Δt = 450.0
longlat = FT.((-70.0, 53.0))
zlim = FT.((-15, 0))
nelements = 15
dz_tuple = FT.((3, 0.05))
domain = ClimaLand.Domains.Column(; zlim, longlat, nelements, dz_tuple);
toml_dict = LP.create_toml_dict(FT)

model = setup_model(FT, start_date, stop_date, Δt, domain, toml_dict);
diagnostics = ClimaLand.default_diagnostics(
    model,
    start_date,
    outdir;
    reduction_period = :daily,
    output_vars = ["tsoil", "swc", "hr", "sco2", "so2", "scd", "sod", "scms"],
);
simulation =
    LandSimulation(start_date, stop_date, Δt, model; outdir, diagnostics);

@info "Run: Global Soil-Canopy-Snow-SoilCO2 Model"
@info "Resolution: $(domain.nelements)"
@info "Timestep: $Δt s"
@info "Start Date: $start_date"
@info "Stop Date: $stop_date"
CP.log_parameter_information(toml_dict, joinpath(root_path, "parameters.toml"))
ClimaLand.Simulations.solve!(simulation);

LandSimVis.make_timeseries(simulation; savedir = root_path)
