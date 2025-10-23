"""
This runs a standalone canopy model driven by prescribed atmosphere and soil measured
at a Fluxnet site (Vaira Ranch; US-Var). Furthermore, this demonstrates the use of
a simple piecewise-linear soil moisture stress function to downweight photosynthesis and
stomatal conductance.
"""

import SciMLBase
import ClimaComms
ClimaComms.@import_required_backends
using ClimaCore
import ClimaParams as CP
using Dates
using Statistics
using ClimaDiagnostics
using ClimaUtilities

using ClimaLand
using ClimaLand.Domains: Point
using ClimaLand.Canopy
using ClimaLand.Canopy.PlantHydraulics
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.Simulations: LandSimulation, solve!

import ClimaLand.FluxnetSimulations as FluxnetSimulations
using CairoMakie, ClimaAnalysis, GeoMakie, Printf, StatsBase
import ClimaLand.LandSimVis as LandSimVis

const FT = Float64
toml_dict = LP.create_toml_dict(FT)
climaland_dir = pkgdir(ClimaLand)
prognostic_land_components = (:canopy)
site_ID = "US-Var"
site_ID_val = FluxnetSimulations.replace_hyphen(site_ID)

# Get the default values for this site's domain, location, and parameters
(; dz_tuple, nelements, zmin, zmax) =
    FluxnetSimulations.get_domain_info(FT, Val(site_ID_val))
(; time_offset, lat, long) =
    FluxnetSimulations.get_location(FT, Val(site_ID_val))
(; atmos_h) = FluxnetSimulations.get_fluxtower_height(FT, Val(site_ID_val))
(;
    soil_ν,
    soil_K_sat,
    soil_S_s,
    soil_vg_n,
    soil_vg_α,
    θ_r,
    ν_ss_quartz,
    ν_ss_om,
    ν_ss_gravel,
    z_0m_soil,
    z_0b_soil,
    soil_ϵ,
    soil_α_PAR,
    soil_α_NIR,
    Ω,
    χl,
    α_PAR_leaf,
    λ_γ_PAR,
    τ_PAR_leaf,
    α_NIR_leaf,
    τ_NIR_leaf,
    ϵ_canopy,
    ac_canopy,
    g1,
    Drel,
    g0,
    Vcmax25,
    SAI,
    f_root_to_shoot,
    K_sat_plant,
    ψ63,
    Weibull_param,
    a,
    conductivity_model,
    retention_model,
    plant_ν,
    plant_S_s,
    rooting_depth,
    n_stem,
    n_leaf,
    h_leaf,
    h_stem,
    h_canopy,
    z_0m,
    z_0b,
) = FluxnetSimulations.get_parameters(FT, Val(site_ID_val))

# Construct the ClimaLand domain to run the simulation on
domain = Point(; z_sfc = FT(0.0), longlat = (long, lat))

# Set up the timestepping information for the simulation
dt = Float64(450) # 7.5 minutes

# This reads in the data from the flux tower site and creates
# the atmospheric and radiative driver structs for the model
(start_date, stop_date) = FluxnetSimulations.get_data_dates(
    site_ID,
    time_offset;
    construct_prescribed_ground = true,
    duration = Year(1),
)
soil_driver_args = (
    α_PAR = soil_α_PAR,
    α_NIR = soil_α_NIR,
    ϵ = soil_ϵ,
    ν = soil_ν,
    θ_r = θ_r,
    hydrology_cm = ClimaLand.Soil.vanGenuchten{FT}(;
        α = soil_vg_α,
        n = soil_vg_n,
    ),
)
forcing = FluxnetSimulations.prescribed_forcing_fluxnet(
    site_ID,
    lat,
    long,
    time_offset,
    atmos_h,
    start_date,
    toml_dict,
    FT;
    construct_prescribed_ground = true,
    soil_driver_args = soil_driver_args,
)

# Now we set up the standalone canopy model
# Set up radiative transfer
radiation_parameters = (;
    Ω,
    G_Function = CLMGFunction(χl),
    α_PAR_leaf,
    τ_PAR_leaf,
    α_NIR_leaf,
    τ_NIR_leaf,
)
radiative_transfer =
    Canopy.TwoStreamModel{FT}(domain, toml_dict; radiation_parameters, ϵ_canopy)

# Set up conductance
conductance = PModelConductance{FT}(toml_dict)

# Set up photosynthesis
photosynthesis = PModel{FT}(domain, toml_dict)

# Set up plant hydraulics
hydraulics = Canopy.PlantHydraulicsModel{FT}(
    domain,
    toml_dict;
    n_stem,
    n_leaf,
    h_stem,
    h_leaf,
    ν = plant_ν,
    S_s = plant_S_s,
    conductivity_model,
    retention_model,
)

# Set up biomass model
# Read in LAI from MODIS data
LAI =
    ClimaLand.prescribed_lai_modis(domain.space.surface, start_date, stop_date);
# Get the maximum LAI at this site over the first year of the simulation
maxLAI = FluxnetSimulations.get_maxLAI_at_site(
    start_date,
    lat,
    long;
    ncd_path = modis_lai_ncdata_path[1],
);
RAI = maxLAI * f_root_to_shoot
height = h_stem + h_leaf
biomass =
    Canopy.PrescribedBiomassModel{FT}(; LAI, SAI, RAI, rooting_depth, height)

# Set up the energy model
energy = Canopy.BigLeafEnergyModel{FT}(toml_dict; ac_canopy)

# Combine the components into a CanopyModel
canopy = CanopyModel{FT}(
    domain,
    forcing,
    LAI,
    toml_dict;
    z_0m,
    z_0b,
    radiative_transfer,
    photosynthesis,
    conductance,
    hydraulics,
)

# This makes a set_ic! function that we can pass into the LandSimulation
set_ic! = FluxnetSimulations.make_set_fluxnet_initial_conditions(
    site_ID,
    start_date,
    time_offset,
    canopy,
);

# Now we set up the diagnostics. We'll write them to memory, then at the end of the
# simulation we can write them to a netCDF file
dict_writer = ClimaDiagnostics.Writers.DictWriter()

output_vars = [
    "sif",      # solar induced fluorescence
    "ra",       # autotrophic respiration
    "gs",       # stomatal conductance
    "gpp",      # gross primary productivity
    "msf",      # moisture stress factor
    "vcmax25",   # maximum carboxylation capacity at 25°C
]

diagnostics = ClimaLand.default_diagnostics(
    canopy,
    start_date;
    output_writer = dict_writer,
    output_vars = output_vars,
    reduction_period = :halfhourly,
)

# This gets the dt period of the forcing data, in seconds
# Then, setting updateat specifies that the forcing should be updated at this period
data_dt = Second(FluxnetSimulations.get_data_dt(site_ID, time_offset))
updateat = Array(start_date:data_dt:stop_date);

# Now we can construct the simulation object and solve it.
simulation = LandSimulation(
    start_date,
    stop_date,
    Δt,
    canopy;
    set_ic!,
    updateat,
    user_callbacks = (),
    diagnostics,
)

@time solve!(simulation)

ClimaLand.Diagnostics.close_output_writers(diagnostics)
comparison_data = FluxnetSimulations.get_comparison_data(site_ID, time_offset)
savedir = joinpath(
    pkgdir(ClimaLand),
    "experiments/standalone/Vegetation/$(site_ID)/out",
)
mkpath(savedir)

# get the modeled GPP
time_vector, gpp_model =
    ClimaLand.Diagnostics.diagnostic_as_vectors(dict_writer, "gpp_30m_average")
time_vector = DateTime.(time_vector)

_, msf_model =
    ClimaLand.Diagnostics.diagnostic_as_vectors(dict_writer, "msf_30m_average")

overlapping_times = in.(comparison_data.UTC_datetime, Ref(Set(time_vector)))
gpp_obs = comparison_data.gpp[overlapping_times]

# take the daily means
gpp_model_daily = mean(reshape(gpp_model, 48, 365); dims = 1)[:]
gpp_obs_daily = mean(reshape(gpp_obs, 48, 365); dims = 1)[:]
msf_daily = mean(reshape(msf_model, 48, 365); dims = 1)[:]
day_of_year = range(1, 365)

# plot GPP comparison
p1 = plot(
    day_of_year,
    gpp_model_daily;
    title = "GPP",
    label = "Modeled GPP",
    xlabel = "Day of year",
    ylabel = "GPP (mol/m²/s)",
    framestyle = :box,
    grid = :y,
)
plot!(p1, day_of_year, gpp_obs_daily, label = "Observed GPP")

p2 = plot(
    day_of_year,
    msf_daily;
    title = "Soil moisture stress",
    label = "βm",
    xlabel = "Day of year",
    ylabel = "βm (unitless)",
    legend = :topright,
    framestyle = :box,
    grid = :y,
)

plt = plot(p1, p2; layout = (2, 1), link = :x, size = (900, 550), dpi = 300)
savefig(plt, joinpath(savedir, "gpp_comparison.png"))
