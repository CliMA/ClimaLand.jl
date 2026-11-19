# # Global run of land model

# The code sets up and runs ClimaLand's default land model
# globally, using ERA5 data as forcing. In this simulation, we have
# turned lateral flow off because horizontal boundary conditions
# are not yet supported by ClimaCore.

# Simulation Setup
# Number of spatial elements: 360 latitude, 180 longitude, 15 vertically
# Soil depth: 15 m
# Simulation duration: set by stop and start date passed as CLA
# LAI: Either prescribed (MODIS climatological (modis_c) or read from MODIS data directly (modis);
# note that climatological required for years outside of 2000-2019) or prognostic (prog)
# Timestep: 900 s
# Atmos forcing update: every 3 hours with ERA5.

import ClimaComms
ClimaComms.@import_required_backends
import ClimaParams as CP
using ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.Simulations: LandSimulation, solve!
using Dates
using CairoMakie, GeoMakie, ClimaAnalysis
import ClimaLand.LandSimVis as LandSimVis
include("parse_cla_snowyland.jl")
const FT = Float64;
context = ClimaComms.context()
ClimaComms.init(context)
device = ClimaComms.device()
device_suffix = device isa ClimaComms.CPUSingleThreaded ? "cpu" : "gpu"
parsed_args = parse_commandline()
start_date = parsed_args["start_date"]
start_date = parsed_args["stop_date"]
LAI_MODE = parsed_args["lai"]
ILAMB_RUN = parsed_args["ilamb"] == true
spinup_months = 12
lai_suffix = ""
@assert lai_mode == "modis" || lai_mode == "modis_c" || lai_mode == "prog"
if lai_mode == "modis" # Years we have data available
    @assert start_date >= DateTime(2000)
    @assert stop_date <= DateTime(2020)
end
if LAI_MODE == "prog"
    spinup_months = 24
    lai_suffix = "_opt_lai"
end
@assert stop_date >= start_date + Month(spinup_months)

root_path = "snowy_land_pmodel$(lai_suffix)_longrun_$(device_suffix)"
diagnostics_outdir = joinpath(root_path, "global_diagnostics")
outdir =
    ClimaUtilities.OutputPathGenerator.generate_output_path(diagnostics_outdir)

Δt = 900.0
domain =
    ClimaLand.Domains.global_box_domain(FT; context, mask_threshold = FT(0.99))
toml_dict = LP.create_toml_dict(FT)
atmos, radiation = ClimaLand.prescribed_forcing_era5(
    start_date,
    stop_date,
    domain.space.surface,
    toml_dict,
    FT;
    max_wind_speed = 25.0,
    context,
)
forcing = (; atmos, radiation)
prognostic_land_components = (:canopy, :lake, :snow, :soil, :soilco2)
if lai_mode != "prog"
    if LAI_MODE == "modis"
        LAI = ClimaLand.Canopy.prescribed_lai_modis(
            domain.space.surface,
            start_date,
            stop_date,
        )
    elseif LAI_MODE = "modis_c"
        LAI = ClimaLand.Canopy.prescribed_lai_modis(
            domain.space.surface,
            start_date,
            stop_date,
        )
    end
    model = LandModel{FT}(
        forcing,
        LAI,
        toml_dict,
        domain,
        Δt;
        prognostic_land_components,
    )
else
    # The LandModel constructor uses the prognostic LAI model if no
    # prescribed LAI is passed; the carbon pools wrap it.
    biomass = ClimaLand.Canopy.PrognosticCarbonModel{FT}(
        ClimaLand.Canopy.ZhouOptimalLAIModel{FT}(
            ClimaLand.Domains.obtain_surface_domain(domain),
            toml_dict,
        ),
            toml_dict,
    )
    model = LandModel{FT}(
        forcing,
        toml_dict,
        domain,
        Δt;
        prognostic_land_components,
        biomass,
    )
end

simulation = LandSimulation(start_date, stop_date, Δt, model; outdir)
@info "Run: Global Soil-Canopy-Snow-Lake-SoilCO2 Model"
@info "LAI: $(lai_mode)"
@info "Resolution: $(domain.nelements)"
@info "Timestep: $Δt s"
@info "Start Date: $start_date"
@info "Stop Date: $stop_date"
CP.log_parameter_information(toml_dict, joinpath(root_path, "parameters.toml"))
# The pools do not feed back on GPP, so their equilibrium under the mean drivers of a
# year is analytic (`Canopy.equilibrium_carbon_pools`) and two simulated years replace
# decades of wood growth. The pools are set to it under the first year, then under the
# second, which follows the transient of the initial state; the leaderboards leave out
# both years
if LAI_MODE == "prog"
    for spinup_years in 1:2
        ClimaLand.Simulations.spin_up_carbon_pools!(
            simulation,
            start_date + Year(spinup_years),
        )
    end
end
ClimaLand.Simulations.solve!(simulation)

LandSimVis.make_annual_timeseries(simulation; savedir = root_path)
LandSimVis.make_heatmaps(simulation; savedir = root_path, date = stop_date)
leaderboard_data_sources = ["ERA5", "FlagshipCarbonMetrics"]
LAI_MODE == "prog" && push!(leaderboard_data_sources, "FlagshipVegetationMetrics")
LandSimVis.make_leaderboard_plots(
    simulation;
    savedir = root_path,
    leaderboard_data_sources,
    spin_up_months,
)

if ILAMB_RUN
    include("../ilamb/ilamb_conversion.jl")
    make_compatible_with_ILAMB(
        joinpath(root_path, "global_diagnostics", "output_active"),
        joinpath(root_path, "global_diagnostics", "ILAMB_diagnostics"),
    )
end
