using Dates

#= issues: (time index, lon, lat, z, SWE, scf, total_albedo):
 (DateTime("2019-11-03T00:00:00"), 128.0, -73.0, 7.470959, 3.4211667, 1.0, 0.7661099)
 (DateTime("2014-03-16T00:00:00"), 147.0, -6.0, 0.0001744235, 0.00017585493, 0.0029724056, 0.07051761) #why is z<SWE??
 (DateTime("2014-05-31T00:00:00"), -79.0, -2.0, 0.0003941978, 0.00039459218, 0.0067059896, 0.083248325) #why is z<SWE??
 (DateTime("2010-10-27T00:00:00"), -77.0, 3.0, 0.0, 0.0, 0.0, 0.06903501)
 (DateTime("2002-09-07T00:00:00"), -17.0, 64.0, 0.0005595295, 0.0005580454, 0.009502905, 0.04070939)
=#
use_col = false
debug_mode = true
col_lon_lat = (-17, 64) #only used if use_col is true
const FT = Float32;
const setup = Dict(
    "output_tag" => "global_debug_1",
    "use_neural_albedo" => true,
    "use_neural_depth" => true,
    "use_sfc_temp" => true,
    "max_wind_speed" => FT(25.0),
    "snow_min_density_param" => FT(300), #only used if not using neural models
    "start_date" => DateTime("2000-06-01"), #earliest we have is Jan 1st 1979
    "stop_date" => DateTime("2003-12-01"), #latest we have is 7th 2024
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

output_dir = joinpath("/home/acharbon/thesis_outputs", setup["output_tag"])

toml_dict = LP.create_toml_dict(FT)
Δt = setup["dt"]

@info "Setting up Simulation!"

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

diagnostics = ClimaLand.default_diagnostics(
    model,
    setup["start_date"];
    output_writer,
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

const sim_mask = ClimaLand.Domains.landsea_mask(ClimaLand.get_domain(model))

function check_nans(x::ClimaCore.Fields.FieldVector, mask; excl = [])
    return any(
        check_nans(getproperty(x, p), mask) for
        p in propertynames(x) if !(p in excl)
    )
end

function check_nans(x::NamedTuple, mask; excl = [])
    return any(check_nans(x[p], mask) for p in propertynames(x) if !(p in excl))
end

function check_nans(x::ClimaCore.Fields.Field, mask)
    return check_nans(x, axes(x), mask)
end

function check_nans(
    x::ClimaCore.Fields.Field,
    space::Union{
        ClimaCore.Spaces.AbstractSpectralElementSpace,
        ClimaCore.Spaces.AbstractPointSpace,
    },
    mask,
)
    return count_nans(x, mask)
end

function check_nans(
    x::ClimaCore.Fields.Field,
    space::Union{
        ClimaCore.Spaces.ExtrudedFiniteDifferenceSpace,
        ClimaCore.Spaces.FiniteDifferenceSpace,
    },
    mask,
)
    return count_nans(ClimaLand.Domains.top_center_to_surface(x), mask)
end

function count_nans(x::ClimaCore.Fields.Field, mask)
    if isnothing(mask)
        return count(isnan, parent(x)) > 0
    else
        return mapreduce(
            (s, m) -> m != 0 && isnan(s),
            |,
            parent(x),
            parent(mask),
        )
    end
end

@inline get_tendency_object(::Val{ClimaTimeSteppers.ARS111()}, simulation) =
    simulation._integrator.cache.T_exp[1]

function nans_exist(simulation, mask; excl = [])
    nans_in_Y = check_nans(simulation._integrator.u, mask, excl = excl)
    nans_in_p = check_nans(simulation._integrator.p, mask, excl = excl)
    nans_in_dY = check_nans(
        get_tendency_object(Val(simulation.timestepper.name), simulation),
        mask,
        excl = excl,
    )
    return nans_in_Y || nans_in_p || nans_in_dY
end

function write_present_state(simulation, old_Y, outdir)
    curr_p = simulation._integrator.p
    curr_Y = simulation._integrator.u
    curr_dY = get_tendency_object(Val(simulation.timestepper.name), simulation)
    JLD2.jldsave(
        joinpath(outdir, "integrator_state.jld2"),
        p = curr_p,
        Y = curr_Y,
        old_Y = old_Y,
        dY = curr_dY,
    )
end

parameter_log_path = joinpath(output_dir, "parameters.toml")

isdir(output_dir) || mkdir(output_dir)
CP.log_parameter_information(toml_dict, parameter_log_path)

const exclude_list = [:soilco2]
old_Y = deepcopy(simulation._integrator.u)
@info "Beginning Simulation!"
if !debug_mode
    ClimaLand.Simulations.solve!(simulation)
    @info "Simulation Complete!"
else
    stop_t = Second(setup["stop_date"] .- setup["start_date"]).value
    while simulation._integrator.t.counter < stop_t
        if nans_exist(simulation, sim_mask, excl = exclude_list)
            @info "NaNs DISCOVERED!"
            write_present_state(simulation, old_Y, outputp_dir)
            break
        else
            @. old_Y = simulation._integrator.u
            ClimaLand.Simulations.step!(simulation)
        end
    end
    state = simulation._integrator
end
