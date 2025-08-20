"""
This runs a standalone canopy model driven by prescribed atmosphere and soil measured
at a Fluxnet site (Vaira Ranch; US-Var). Furthermore, this demonstrates the use of 
a simple piecewise-linear soil moisture stress function to downweight photosynthesis and 
stomatal conductance. 
"""

import SciMLBase
import ClimaTimeSteppers as CTS
import ClimaComms
ClimaComms.@import_required_backends
using ClimaCore
import ClimaParams as CP
using Dates
using Plots
using Statistics
using ClimaDiagnostics
using ClimaUtilities

using ClimaLand
using ClimaLand.Domains: Point
using ClimaLand.Snow
using ClimaLand.Soil
using ClimaLand.Soil.Biogeochemistry
using ClimaLand.Canopy
using ClimaLand.Canopy.PlantHydraulics
import ClimaLand
import ClimaLand.Parameters as LP

import ClimaLand.FluxnetSimulations as FluxnetSimulations
import ClimaLand.Simulations: LandSimulation, solve!

const FT = Float64
earth_param_set = LP.LandParameters(FT)
climaland_dir = pkgdir(ClimaLand)
prognostic_land_components = (:canopy, :snow, :soil, :soilco2)

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
    z0_m,
    z0_b,
) = FluxnetSimulations.get_parameters(FT, Val(site_ID_val))

# Construct the ClimaLand domain to run the simulation on
canopy_domain = Point(; z_sfc = FT(0.0), longlat = (long, lat))

# Set up the timestepping information for the simulation
t0 = Float64(0)
Δt = Float64(450) # 7.5 minutes
N_spinup_days = 15
N_days = N_spinup_days + 340
tf = Float64(t0 + N_days * 3600 * 24)
t_spinup = Float64(t0 + N_spinup_days * 3600 * 24)

timestepper = CTS.ARS111();
ode_algo = CTS.IMEXAlgorithm(
    timestepper,
    CTS.NewtonsMethod(
        max_iters = 3,
        update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
    ),
);

# This reads in the data from the flux tower site and creates
# the atmospheric and radiative driver structs for the model
(start_date, stop_date) = FluxnetSimulations.get_data_dates(
    site_ID,
    time_offset;
    construct_prescribed_soil = true,
    duration = Year(1),
)
soil_driver_args = (
    α_PAR = FT(0.2),
    α_NIR = FT(0.4),
    ϵ = FT(0.99),
    ν = soil_ν,
    θ_r = θ_r,
    hydrology_cm = vanGenuchten{FT}(; α = soil_vg_α, n = soil_vg_n),
)
(; atmos, radiation, soil) = FluxnetSimulations.prescribed_forcing_fluxnet(
    site_ID,
    lat,
    long,
    time_offset,
    atmos_h,
    start_date,
    earth_param_set,
    FT;
    construct_prescribed_soil = true,
    soil_driver_args = soil_driver_args,
)
forcing = (; atmos, radiation, ground = soil)


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
    Canopy.TwoStreamModel{FT}(canopy_domain; radiation_parameters, ϵ_canopy)

# Set up conductance
conductance = PModelConductance{FT}()

# Set up photosynthesis
photosynthesis = PModel{FT}()

# Set up the soil moisture stress function. We'll use a simple linear piecewise function with
# two thresholds θ_c and θ_w
soil_moisture_stress_params = PiecewiseMoistureStressParameters(
    FT;
    θ_c = FT(0.25),
    θ_w = FT(0.13),
    c = FT(1.0),
    β0 = FT(1.0),
)
soil_moisture_stress =
    PiecewiseMoistureStressModel{FT}(soil_moisture_stress_params)

# Set up plant hydraulics
# Read in LAI from MODIS data
surface_space = canopy_domain.space.surface
modis_lai_ncdata_path = ClimaLand.Artifacts.modis_lai_multiyear_paths(;
    context = ClimaComms.context(surface_space),
    start_date,
    stop_date,
)
LAI = ClimaLand.prescribed_lai_modis(
    modis_lai_ncdata_path,
    surface_space,
    start_date,
);
# Get the maximum LAI at this site over the first year of the simulation
maxLAI =
    FluxnetSimulations.get_maxLAI_at_site(modis_lai_ncdata_path[1], lat, long);
RAI = maxLAI * f_root_to_shoot
hydraulics = Canopy.PlantHydraulicsModel{FT}(
    canopy_domain,
    LAI;
    n_stem,
    n_leaf,
    h_stem,
    h_leaf,
    SAI,
    RAI,
    ν = plant_ν,
    S_s = plant_S_s,
    conductivity_model,
    retention_model,
    rooting_depth,
)

# Set up the energy model
energy = Canopy.BigLeafEnergyModel{FT}(; ac_canopy)

# Combine the components into a CanopyModel
canopy = CanopyModel{FT}(
    canopy_domain,
    forcing,
    LAI,
    earth_param_set;
    z_0m = z0_m,
    z_0b = z0_b,
    radiative_transfer,
    photosynthesis,
    conductance,
    soil_moisture_stress,
    hydraulics,
)

# Get the start date in local time, which is useful for plotting diagnostics
# We'll use this time to index the saved variables
start_date_local =
    start_date - FluxnetSimulations.hour_offset_to_period(time_offset)

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
    average_period = :halfhourly,
)

# This gets the dt period of the forcing data, in seconds
# Then, setting updateat specifies that the forcing should be updated at this period
data_dt = Second(FluxnetSimulations.get_data_dt(site_ID))
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
