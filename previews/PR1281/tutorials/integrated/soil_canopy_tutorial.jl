# # Coupling the CliMA Canopy and Soil Hydraulics Models

# In the
# [`previous tutorial`](@ref https://clima.github.io/ClimaLand.jl/dev/generated/soil_plant_hydrology_tutorial/),
# we demonstrated how to run the canopy model in
# standalone mode using prescribed values for the inputs of soil hydraulics
# into the canopy hydraulics model. However, ClimaLand has the built-in capacity
# to couple the canopy model with a soil physics model and timestep the two
# simulations together to model a canopy-soil system. This tutorial
# demonstrates how to setup and run a coupled simulation, again using
# initial conditions, atmospheric and radiative flux conditions, and canopy
# properties observed at the US-MOz flux tower, a flux tower located within an
# oak-hickory forest in Ozark, Missouri, USA. See [Wang et al. 2021]
# (https://doi.org/10.5194/gmd-14-6741-2021) for details on the site and
# parameters.


# In ClimaLand, the coupling of the canopy and soil models is done by
# pairing the inputs and outputs which between the two models so that they match.
# For example, the root extraction of the canopy hydraulics model, which acts as a
# boundary flux for the plant system, is paired with a source term for root extraction in
# the soil model, so that the flux of water from the soil into the roots is
# equal and factored into both models. This pairing is done automatically in the
# constructor for a
# [`SoilCanopyModel`](https://clima.github.io/ClimaLand.jl/dev/APIs/ClimaLand/#LSM-Model-Types-and-methods)
# so that a
# user needs only specify the necessary arguments for each of the component
# models, and the two models will automatically be paired into a coupled
# simulation.

# # Preliminary Setup

# Load External Packages:

import SciMLBase
using Plots
using Statistics
using Dates
using Insolation


# Load CliMA Packages and ClimaLand Modules:

using ClimaCore
import ClimaParams as CP
import ClimaTimeSteppers as CTS
using ClimaLand
using ClimaLand.Domains: Column, obtain_surface_domain
using ClimaLand.Soil
using ClimaLand.Soil.Biogeochemistry
using ClimaLand.Canopy
using ClimaLand.Canopy.PlantHydraulics
import ClimaLand.Simulations: LandSimulation, solve!
import ClimaLand
import ClimaLand.Parameters as LP
using DelimitedFiles
import ClimaLand.FluxnetSimulations as FluxnetSimulations

# Define the floating point precision desired (64 or 32 bit), and get the
# parameter set holding constants used across CliMA Models:

const FT = Float32;
earth_param_set = LP.LandParameters(FT);

# First provide some information about the site
site_ID = "US-MOz"
# Timezone (offset from UTC in hrs)
time_offset = 7
start_date = DateTime(2010) + Hour(time_offset) + Day(150) # start mid year
# Specify the time range and dt value over which to perform the simulation.
N_days = 100
end_date = start_date + Day(N_days)
dt = Float64(30)

# Site latitude and longitude
lat = FT(38.7441) # degree
long = FT(-92.2000) # degree

# Height of the sensor at the site
atmos_h = FT(32) # m


# - We will be using prescribed atmospheric and radiative drivers from the
#   US-MOz tower, which we read in here. We are using prescribed
#   atmospheric and radiative flux conditions, but it is also possible to couple
#   the simulation with atmospheric and radiative flux models. We also
# read in the observed LAI and let that vary in time in a prescribed manner.

# Forcing data
(; atmos, radiation) = FluxnetSimulations.prescribed_forcing_fluxnet(
    site_ID,
    lat,
    long,
    time_offset,
    atmos_h,
    start_date,
    earth_param_set,
    FT,
);

# Setup the domain for the model:

nelements = 10
zmin = FT(-2)
zmax = FT(0)
domain =
    Column(; zlim = (zmin, zmax), nelements = nelements, longlat = (long, lat));

# # Setup the Coupled Canopy and Soil Physics Model

# We want to simulate the canopy-soil system together, so the model type
# [`SoilCanopyModel`](https://clima.github.io/ClimaLand.jl/dev/APIs/ClimaLand/#LSM-Model-Types-and-methods)
# is chosen.
# From the linked documentation, we see that we need to provide the soil model
# type and arguments as well as the canopy model component types, component
# arguments, and the canopy model arguments, so we first need to initialize
# all of these.
prognostic_land_components = (:canopy, :soil, :soilco2);

# For our soil model, we will choose the
# [`EnergyHydrology`](https://clima.github.io/ClimaLand.jl/dev/APIs/Soil/#Soil-Models-2)
# and set up all the necessary arguments. See the
# [tutorial](https://clima.github.io/ClimaLand.jl/dev/generated/Soil/soil_energy_hydrology/)
# on the model for a more detailed explanation of the soil model.

# Define the parameters for the soil model and provide them to the model constructor.

# Soil parameters
soil_ν = FT(0.5) # m3/m3
soil_K_sat = FT(4e-7) # m/s
soil_S_s = FT(1e-3) # 1/m
soil_vg_n = FT(2.05) # unitless
soil_vg_α = FT(0.04) # inverse meters
θ_r = FT(0.067); # m3/m3
ν_ss_quartz = FT(0.1)
ν_ss_om = FT(0.1)
ν_ss_gravel = FT(0.0);

z_0m_soil = FT(0.1)
z_0b_soil = FT(0.1)
soil_ϵ = FT(0.98)
soil_α_PAR = FT(0.2)
soil_α_NIR = FT(0.4)

soil_albedo = ClimaLand.Soil.ConstantTwoBandSoilAlbedo{FT}(;
    PAR_albedo = soil_α_PAR,
    NIR_albedo = soil_α_NIR,
)
runoff_model = ClimaLand.Soil.Runoff.NoRunoff()

soil = EnergyHydrology{FT}(
    domain,
    (; atmos, radiation),
    earth_param_set;
    prognostic_land_components,
    albedo = soil_albedo,
    runoff = runoff_model,
    additional_sources = (ClimaLand.RootExtraction{FT}(),),
    retention_parameters = (;
        ν = soil_ν,
        θ_r,
        K_sat = soil_K_sat,
        hydrology_cm = vanGenuchten{FT}(; α = soil_vg_α, n = soil_vg_n),
    ),
    composition_parameters = (; ν_ss_om, ν_ss_quartz, ν_ss_gravel),
    S_s = soil_S_s,
    z_0m = z_0m_soil,
    z_0b = z_0b_soil,
    emissivity = soil_ϵ,
);

# For the heterotrophic respiration model, see the
# [documentation](https://clima.github.io/ClimaLand.jl/previews/PR214/dynamicdocs/pages/soil_biogeochemistry/microbial_respiration/)
# to understand the parameterisation.
# The domain is defined similarly to the soil domain described above.
soil_organic_carbon =
    ClimaLand.PrescribedSoilOrganicCarbon{FT}(TimeVaryingInput((t) -> 5));
co2_prognostic_soil = Soil.Biogeochemistry.PrognosticMet(soil.parameters);
drivers = Soil.Biogeochemistry.SoilDrivers(
    co2_prognostic_soil,
    soil_organic_carbon,
    atmos,
);
soilco2 = Soil.Biogeochemistry.SoilCO2Model{FT}(domain, drivers);

# Read in prescribed LAI at the site from global MODIS data
surface_space = domain.space.surface;
modis_lai_ncdata_path = ClimaLand.Artifacts.modis_lai_multiyear_paths(;
    start_date = start_date,
    end_date = end_date,
);
LAI = ClimaLand.prescribed_lai_modis(
    modis_lai_ncdata_path,
    surface_space,
    start_date,
);
# Get the maximum LAI at this site over the first year of the simulation
maxLAI =
    FluxnetSimulations.get_maxLAI_at_site(modis_lai_ncdata_path[1], lat, long);

# For a coupled SoilCanopyModel, we provide a flag to the canopy that indicates
#  the ground forcing is prognostic (i.e. the soil model) rather than prescribed.
ground = ClimaLand.PrognosticSoilConditions{FT}();

# Next we need to set up the [`CanopyModel`](https://clima.github.io/ClimaLand.jl/dev/APIs/canopy/Canopy/#Canopy-Model-Structs).
# For more details on the specifics of this model see the previous tutorial.
surface_domain = ClimaLand.Domains.obtain_surface_domain(domain);

# Let's overwrite some default parameters for the canopy model components.
# This involves constructing the components individually and then
# passing them to the canopy model constructor.

# Canopy conductance
conductance = Canopy.MedlynConductanceModel{FT}(surface_domain; g1 = 141);

# Canopy radiative transfer
radiation_parameters = (;
    G_Function = ConstantGFunction(FT(0.5)),
    α_PAR_leaf = 0.1,
    α_NIR_leaf = 0.45,
    τ_PAR_leaf = 0.05,
    τ_NIR_leaf = 0.25,
    Ω = 0.69,
);
radiative_transfer = TwoStreamModel{FT}(surface_domain; radiation_parameters);

# Canopy photosynthesis
photosynthesis_parameters = (; is_c3 = FT(1), Vcmax25 = FT(5e-5));
photosynthesis = FarquharModel{FT}(surface_domain; photosynthesis_parameters);

# Canopy hydraulics
n_stem = Int64(1)
n_leaf = Int64(1)
h_stem = FT(9)
h_leaf = FT(9.5)
f_root_to_shoot = FT(3.5)
plant_ν = FT(0.7)
plant_S_s = FT(1e-2 * 0.0098)
SAI = FT(0.00242)
RAI = (SAI + maxLAI) * f_root_to_shoot;
K_sat_plant = FT(1.8e-8)
ψ63 = FT(-4 / 0.0098)
Weibull_param = FT(4)
conductivity_model =
    PlantHydraulics.Weibull{FT}(K_sat_plant, ψ63, Weibull_param)
hydraulics = Canopy.PlantHydraulicsModel{FT}(
    surface_domain,
    LAI;
    SAI,
    RAI,
    n_stem,
    n_leaf,
    h_stem,
    h_leaf,
    ν = plant_ν,
    S_s = plant_S_s,
    conductivity_model,
    rooting_depth = FT(1),
)

canopy = ClimaLand.Canopy.CanopyModel{FT}(
    surface_domain,
    (; atmos, radiation, ground),
    LAI,
    earth_param_set;
    prognostic_land_components = (:canopy, :soil, :soilco2),
    conductance,
    radiative_transfer,
    photosynthesis,
    hydraulics,
);

# Now we can combine the soil and canopy models into a single combined model.
land = SoilCanopyModel{FT}(soilco2, soil, canopy);

# We need to provide initial conditions for the model:
function set_ic!(Y, p, t0, model)
    Y.soil.ϑ_l = FT(0.4)
    Y.soil.θ_i = FT(0.0)
    T_0 = FT(288.7)
    ρc_s =
        volumetric_heat_capacity.(
            Y.soil.ϑ_l,
            Y.soil.θ_i,
            land.soil.parameters.ρc_ds,
            earth_param_set,
        )
    Y.soil.ρe_int =
        volumetric_internal_energy.(Y.soil.θ_i, ρc_s, T_0, earth_param_set)

    Y.soilco2.C .= FT(0.000412) # set to atmospheric co2, mol co2 per mol air

    canopy_retention_model = canopy.hydraulics.parameters.retention_model
    ψ_stem_0 = FT(-1e5 / 9800)
    ψ_leaf_0 = FT(-2e5 / 9800)

    S_l_ini =
        inverse_water_retention_curve.(
            canopy_retention_model,
            [ψ_stem_0, ψ_leaf_0],
            plant_ν,
            plant_S_s,
        )
    for i in 1:2
        Y.canopy.hydraulics.ϑ_l.:($i) .=
            augmented_liquid_fraction.(plant_ν, S_l_ini[i])
    end
    evaluate!(Y.canopy.energy.T, atmos.T, t0)
end

# Choose the timestepper and solver needed for the problem.
timestepper = CTS.ARS343()
ode_algo = CTS.IMEXAlgorithm(
    timestepper,
    CTS.NewtonsMethod(
        max_iters = 1,
        update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
    ),
);


# Set the callbacks update times, which govern
# how often we save output, and how often we update
# the forcing data ("drivers")

n = 120
saveat = Array(start_date:Second(n * dt):end_date);
sv = (;
    t = Array{DateTime}(undef, length(saveat)),
    saveval = Array{NamedTuple}(undef, length(saveat)),
)
saving_cb = ClimaLand.NonInterpSavingCallback(sv, saveat)
updateat = Array(start_date:Second(1800):end_date);

# Create the LandSimulation object, which will also create and initialize the state vectors,
# the cache, the driver callbacks, and set the initial conditions.
simulation = LandSimulation(
    start_date,
    end_date,
    dt,
    land;
    set_ic! = set_ic!,
    updateat = updateat,
    solver_kwargs = (; saveat = deepcopy(saveat)),
    timestepper = ode_algo,
    user_callbacks = (saving_cb,),
    diagnostics = (),
);

sol = solve!(simulation);

# # Plotting

# Now that we have both a soil and canopy model incorporated together, we will
# show how to plot some model data demonstrating the time series produced from
# each of these models. As before, we may plot the GPP of the system as well
# as transpiration showing fluxes in the canopy.

daily = FT.(sol.t) ./ 3600 ./ 24
model_GPP = [
    parent(sv.saveval[k].canopy.photosynthesis.GPP)[1] for
    k in 1:length(sv.saveval)
]

plt1 = Plots.plot(size = (600, 700));
Plots.plot!(
    plt1,
    daily,
    model_GPP .* 1e6,
    label = "Model",
    xlim = [minimum(daily), maximum(daily)],
    xlabel = "days",
    ylabel = "GPP [μmol/mol]",
);

# Transpiration plot:

T = [
    parent(sv.saveval[k].canopy.turbulent_fluxes.transpiration)[1] for
    k in 1:length(sv.saveval)
]
T = T .* (1e3 * 24 * 3600)

plt2 = Plots.plot(size = (500, 700));
Plots.plot!(
    plt2,
    daily,
    T,
    label = "Model",
    xlim = [minimum(daily), maximum(daily)],
    xlabel = "days",
    ylabel = "Vapor Flux [mm/day]",
);

# Show the two plots together:

Plots.plot(plt1, plt2, layout = (2, 1));

# Save the output:
savefig("ozark_canopy_flux_test.png");
# ![](ozark_canopy_flux_test.png)

# Now, we will plot the augmented volumetric liquid water fraction at different
# depths in the soil over the course of the simulation.

plt1 = Plots.plot(size = (500, 700));
ϑ_l_10 = [parent(sol.u[k].soil.ϑ_l)[end] for k in 1:1:length(sol.t)]
plt1 = Plots.plot(
    daily,
    ϑ_l_10,
    label = "10 cm",
    xlabel = "Days",
    ylabel = "SWC [m/m]",
    xlim = [minimum(daily), maximum(daily)],
    size = (500, 700),
    margins = 10Plots.mm,
    color = "blue",
);

plot!(
    plt1,
    daily,
    [parent(sol.u[k].soil.ϑ_l)[end - 1] for k in 1:1:length(sol.t)],
    label = "20cm",
    color = "red",
);

plot!(
    plt1,
    daily,
    [parent(sol.u[k].soil.ϑ_l)[end - 2] for k in 1:1:length(sol.t)],
    label = "30cm",
    color = "purple",
);

# Save the output:
savefig("ozark_soil_test.png");
# ![](ozark_soil_test.png)

# And now to demonstrate the coupling of the soil and canopy models we will plot
# the water fluxes from the soil up into the plant hydraulic system:

root_stem_flux = [
    sum(sv.saveval[k].root_extraction) .* (1e3 * 3600 * 24) for
    k in 1:length(sol.t)
]
plt1 = Plots.plot(
    daily,
    root_stem_flux,
    label = "soil-root-stem water flux",
    ylabel = "Water flux[mm/day]",
    xlim = [minimum(daily), maximum(daily)],
    size = (500, 700),
    margins = 10Plots.mm,
);

# And save the output
savefig("ozark_soil_plant_flux.png");
# ![](ozark_soil_plant_flux.png)
