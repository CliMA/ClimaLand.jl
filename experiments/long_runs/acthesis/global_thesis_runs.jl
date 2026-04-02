using Dates

setup_only = false
on_local = false
use_col = false
debug_mode = false
dump_err_state = false
stop_if_err = false
gather_diagnostics = true
col_lon_lat = (-17, 64) #only used if use_col is true
const FT = Float64;
const setup = Dict(
    "output_tag" => "all_global",
    "use_neural_albedo" => true,
    "use_neural_depth" => true,
    "use_sfc_temp" => true,
    "max_wind_speed" => FT(25.0),
    "snow_min_density_param" => FT(200), #only used if not using neural models
    "start_date" => DateTime("2000-03-01"), #earliest we have is Jan 1st 1979
    "stop_date" => DateTime("2020-06-01"), #latest we have is 11/7/2024
    "dt" => FT(450),
    "output_vars" => [
        "snd",
        "swe",
        "usnow",
        "snowc",
        "salb",
        #"swa",
        "snalb",
        "galb",
        #"rn",
        "shf",
        "lhf",
        "rnir",
        "rpar",
        #"pcflux",
        "apeflux",
        "esflux",
        #"tsoil",
        "swn",
        "lwn",
        #"swu",
        #"lwu",
        #"ct",
        #"snow",
        #"tair",
        #"lai",
        #"sai",
    ], #salb = soil alb, swa = p.snow.α_sfc
)

@info "\n SETTINGS LOCKED!: $(setup["output_tag"]), ($(setup["use_neural_depth"]))\n"
print("\n SETTINGS LOCKED!: $(setup["output_tag"]), ($(setup["use_neural_depth"]))\n")

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
import ClimaLand.Simulations: LandSimulation, solve!, step!
using ClimaTimeSteppers

using CairoMakie, GeoMakie, ClimaAnalysis
import ClimaLand.LandSimVis as LandSimVis

using Flux, JLD2, StaticArrays, InteractiveUtils, Adapt
NS = Base.get_extension(ClimaLand, :ConstrainedNeuralModelExt).NeuralSnow;

context = ClimaComms.context()
ClimaComms.init(context)
device = ClimaComms.device()
device_suffix = device isa ClimaComms.CPUSingleThreaded ? "cpu" : "gpu"
include("debug_tools.jl")

output_dir = joinpath("/home/acharbon/outputs", setup["output_tag"])

toml_dict = LP.create_toml_dict(FT)
Δt = setup["dt"]

@info "Setting up Simulation!"
print("ALL CODE LOADED; SETUP BEGINS\n")

if use_col
    domain = ClimaLand.Domains.Column(;
        zlim = (-FT(15), FT(0.0)),
        nelements = 15,
        longlat = FT.(col_lon_lat),
        dz_tuple = (FT(3.0), FT(0.05)),
    )
else
    domain = ClimaLand.Domains.global_box_domain(
        FT;
        context,
        mask_threshold = FT(0.99),
    )
end

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
    use_lowres_forcing = on_local,
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

if use_col
    horz_degree_res = FT(1)
else
    horz_degree_res =
        sum(ClimaLand.Domains.average_horizontal_resolution_degrees(domain)) / 2 # mean of resolution in latitude and longitude, in degrees
end
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

output_writer = ClimaLand.Diagnostics.default_output_writer(
    domain,
    setup["start_date"],
    output_dir,
)

if gather_diagnostics
    diagnostics = ClimaLand.default_diagnostics(
        model,
        setup["start_date"];
        output_writer,
        output_vars = setup["output_vars"],
        reduction_period = :daily,
        reduction_type = :average,
    )
else
    diagnostics = nothing
end

simulation = LandSimulation(
    setup["start_date"],
    setup["stop_date"],
    Δt,
    model;
    outdir = output_dir,
    diagnostics,
)

parameter_log_path = joinpath(output_dir, "parameters.toml")

if !on_local
    isdir(output_dir) || mkdir(output_dir)
    CP.log_parameter_information(toml_dict, parameter_log_path)
end

@info "Beginning Simulation!"
print("SIMULATION BEGINS\n")
if !debug_mode
    ClimaLand.Simulations.solve!(simulation)
    @info "Simulation Complete!"
    print("SIMULATION COMPLETE!")
else
    const sim_mask = ClimaLand.Domains.landsea_mask(ClimaLand.get_domain(model))
    const exclude_list =
        (:soilco2, :bidiag_matrix_scratch, :full_bidiag_matrix_scratch)
    old_Y = deepcopy(simulation._integrator.u) #my_to_gpu() makes broadcast-setting the two error?
    old_p = deepcopy(simulation._integrator.p) #my_to_gpu() makes broadcast-setting the two error?
    stop_t = Second(setup["stop_date"] .- setup["start_date"]).value
    state = simulation._integrator
    if setup_only
        error("Stopping here; setup_only = true.\n")
    end
    err_found = [false]
    while simulation._integrator.t.counter < stop_t
        if !err_found[1] &&
           nans_exist(simulation, sim_mask, excl = exclude_list)
            @info "******NaNs DISCOVERED! Writing previous and NaN'd state; at time $(simulation._integrator.t.counter)******"
            (!on_local && dump_err_state) &&
                write_present_state(simulation, old_Y, old_p, output_dir)
            err_found[1] = true
            stop_if_err && break
        #=elseif any_negatives(simulation._integrator.p.canopy.radiative_transfer.nir.refl, sim_mask) || any_negatives(simulation._integrator.p.canopy.radiative_transfer.nir.refl, sim_mask)
            @info "NEGATIVE REFL CONSTANTS!"
            !on_local &&
                write_present_state(simulation, old_Y, old_p, output_dir)
            break
        elseif any_negatives(simulation._integrator.p.α_sfc, sim_mask)
            @info "NEGATIVE SWA!"
            !on_local &&
                write_present_state(simulation, old_Y, old_p, output_dir)
            break =#
        else
            @. old_Y = simulation._integrator.u
            set_fields!(old_p, simulation._integrator.p)
            ClimaLand.Simulations.step!(simulation)
        end
    end
end

#=
if use_col && !on_local
    jldsave("./saved_column_data.jld2", ow = output_writer.dict)
    colstate = Dict(Symbol(split(q, "_")[1]) => [my_to_cpu(parent(v))[1] for (k, v) in output_writer.dict[q]] for q in keys(output_writer.dict))
    t = [Second(q.counter) + q.epoch for q in collect(keys(output_writer.dict["snd_1d_average"]))]
    colstate = NamedTuple(colstate)
    (t = t, colstate...)
    jldsave("./saved_column.jld2", colseries = colstate)
end
=#