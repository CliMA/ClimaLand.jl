using Dates

const FT = Float32;
const setup = Dict(
    "output_tag" => "first_try",
    "use_neural_albedo" => true,
    "use_neural_depth" => true,
    "use_sfc_temp" => true,
    "max_wind_speed" => FT(25.0),
    "snow_min_density_param" => FT(300),
    "start_date" => DateTime("2008-03-01"),
    "stop_date" => DateTime("2008-06-01"),
    "dt" => FT(450),
    "output_vars" =>
        ["snd", "swe", "snowc", "salb", "swa", "snalb", "galb"], #salb = soil alb, swa = p.snow.α_sfc
)

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

using CairoMakie, GeoMakie, ClimaAnalysis
import ClimaLand.LandSimVis as LandSimVis

using Flux, JLD2, StaticArrays, InteractiveUtils, Adapt
NS = Base.get_extension(ClimaLand, :ConstrainedNeuralModelExt).NeuralSnow;

context = ClimaComms.context()
ClimaComms.init(context)
device = ClimaComms.device()
device_suffix = device isa ClimaComms.CPUSingleThreaded ? "cpu" : "gpu"

output_dir = joinpath("/home/acharbon/thesis_outputs", setup["output_tag"])

toml_dict = LP.create_toml_dict(FT)
Δt = setup["dt"]

domain =
    ClimaLand.Domains.global_box_domain(FT; context, mask_threshold = FT(0.99))

surface_domain = ClimaLand.Domains.obtain_surface_domain(domain)
surface_space = domain.space.surface

atmos, radiation = ClimaLand.prescribed_forcing_era5(
    setup["start_date"],
    setup["stop_date"],
    surface_space,
    toml_dict,
    FT;
    max_wind_speed = setup["max_wind_speed"],
    context,
    #use_lowres_forcing = true,
)
forcing = (; atmos, radiation)

LAI = ClimaLand.Canopy.prescribed_lai_modis(
    surface_space,
    setup["start_date"],
    setup["stop_date"],
)

ground = ClimaLand.PrognosticGroundConditions{FT}()
canopy_forcing = (; atmos, radiation, ground)
prognostic_land_components = (:canopy, :snow, :soil, :soilco2)

photosynthesis = PModel{FT}(domain, toml_dict)
conductance = PModelConductance{FT}(toml_dict)

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

#what about HTESSEL?
if setup["use_neural_albedo"]
    α_snow = NS.NeuralAlbedoModel(toml_dict, domain.space.surface, Δt = Δt)
else
    α_snow = Snow.ZenithAngleAlbedoModel(toml_dict)
end

#what about Anderson76?
if setup["use_neural_depth"]
    density = NS.NeuralDepthModel(toml_dict, Δt = Δt)
else
    density = Snow.MinimumDensityModel(setup["snow_min_density_param"])
end

horz_degree_res =
    sum(ClimaLand.Domains.average_horizontal_resolution_degrees(domain)) / 2 # mean of resolution in latitude and longitude, in degrees
scf = Snow.WuWuSnowCoverFractionModel(toml_dict, horz_degree_res)

if setup["use_sfc_temp"]
    surf_temp = Snow.EquilibriumGradientTemperatureModel{FT}()
else
    surf_temp = Snow.BulkSurfaceTemperatureModel{FT}()
end

snow = Snow.SnowModel(
    FT,
    surface_domain,
    forcing,
    toml_dict,
    Δt;
    prognostic_land_components,
    α_snow,
    scf,
    density,
    surf_temp,
)

model = LandModel{FT}(
    forcing,
    LAI,
    toml_dict,
    domain,
    Δt;
    prognostic_land_components,
    snow,
    canopy,
)

diagnostics = ClimaLand.default_diagnostics(
    model,
    setup["start_date"];
    output_writer = ClimaLand.Diagnostics.default_output_writer(
        domain,
        setup["start_date"],
        output_dir,
    ),
    output_vars = setup["output_vars"],
    reduction_period = :daily,
    reduction_type = :average,
)

simulation = LandSimulation(
    setup["start_date"],
    setup["stop_date"],
    Δt,
    model;
    outdir = output_dir,
    diagnostics,
)

parameter_log_path = joinpath(output_dir, "parameters.toml")

isdir(output_dir) || mkdir(output_dir)
CP.log_parameter_information(toml_dict, parameter_log_path)
ClimaLand.Simulations.solve!(simulation)
close_output_writers(diagnostics)

#can you save/output the SYPD and wall time for each simulation to an ouptut file?
#^^ pipe the run with a "*> path/to/output_file.txt" and move that file into the appropriate directory