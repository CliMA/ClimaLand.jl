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
import ClimaLand
import ClimaLand.Parameters as LP
using DelimitedFiles
FluxnetSimulationsExt =
    Base.get_extension(ClimaLand, :FluxnetSimulationsExt).FluxnetSimulationsExt;

# Define the floating point precision desired (64 or 32 bit), and get the
# parameter set holding constants used across CliMA Models:

const FT = Float32;
earth_param_set = LP.LandParameters(FT);

# - We will be using prescribed atmospheric and radiative drivers from the
#   US-MOz tower, which we read in here. We are using prescribed
#   atmospheric and radiative flux conditions, but it is also possible to couple
#   the simulation with atmospheric and radiative flux models. We also
# read in the observed LAI and let that vary in time in a prescribed manner.


# First provide some information about the site
site_ID = "US-MOz"
# Timezone (offset from UTC in hrs)
time_offset = 7
start_date = DateTime(2010) + Hour(time_offset)

# Site latitude and longitude
lat = FT(38.7441) # degree
long = FT(-92.2000) # degree

# Height of the sensor at the site
atmos_h = FT(32) # m

# Forcing data
(; atmos, radiation) = FluxnetSimulationsExt.prescribed_forcing_fluxnet(
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

# Define the parameters for the soil model and provide them to the model
# parameters struct:

# Soil parameters
ν = FT(0.5) # m3/m3
K_sat = FT(4e-7) # m/s
n = FT(2.05) # unitless
α = FT(0.04) # inverse meters
θ_r = FT(0.067); # m3/m3
soil = EnergyHydrology{FT}(
    domain,
    (; atmos, radiation),
    earth_param_set;
    prognostic_land_components,
    additional_sources = (RootExtraction{FT}(),),
    retention_parameters = (;
        ν,
        θ_r,
        K_sat,
        hydrology_cm = vanGenuchten{FT}(; α, n),
    ),
)

# For the heterotrophic respiration model, see the
# [documentation](https://clima.github.io/ClimaLand.jl/previews/PR214/dynamicdocs/pages/soil_biogeochemistry/microbial_respiration/)
# to understand the parameterisation.
# The domain is defined similarly to the soil domain described above.
soil_organic_carbon =
    ClimaLand.PrescribedSoilOrganicCarbon{FT}(TimeVaryingInput((t) -> 5))
co2_prognostic_soil = Soil.Biogeochemistry.PrognosticMet(soil.parameters)
drivers = Soil.Biogeochemistry.SoilDrivers(
    co2_prognostic_soil,
    soil_organic_carbon,
    atmos,
)
soilco2 = Soil.Biogeochemistry.SoilCO2Model{FT}(; domain, drivers)

# Next we need to set up the [`CanopyModel`](https://clima.github.io/ClimaLand.jl/dev/APIs/canopy/Canopy/#Canopy-Model-Structs).
# For more details on the specifics of this model see the previous tutorial.
(; LAI, maxLAI) =
    FluxnetSimulationsExt.prescribed_LAI_fluxnet(site_ID, start_date)
surface_domain = ClimaLand.Domains.obtain_surface_domain(domain)
ground = PrognosticSoilConditions{FT}();
canopy = ClimaLand.Canopy.CanopyModel{FT}(
    surface_domain,
    (; atmos, radiation, ground),
    LAI,
    earth_param_set,
)
land = SoilCanopyModel{FT}(; soil, soilco2, canopy);

# Now we can initialize the state vectors and model coordinates, and initialize
# the explicit/implicit tendencies as usual. The Richard's equation time
# stepping is done implicitly, while the canopy model may be explicitly stepped,
# so we use an IMEX (implicit-explicit) scheme for the combined model.

Y, p, coords = initialize(land);
exp_tendency! = make_exp_tendency(land);
imp_tendency! = make_imp_tendency(land);
jacobian! = make_jacobian(land);
jac_kwargs =
    (; jac_prototype = ClimaLand.FieldMatrixWithSolver(Y), Wfact = jacobian!);

# Select the timestepper and solvers needed for the specific problem. Specify the time range and dt
# value over which to perform the simulation.

t0 = Float64(150 * 3600 * 24)# start mid year
N_days = 100
tf = t0 + Float64(3600 * 24 * N_days)
dt = Float64(30)

# We need to provide initial conditions for the soil and canopy hydraulics
# models:
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
canopy_ν = canopy.hydraulics.parameters.ν
canopy_S_s = canopy.hydraulics.parameters.S_s
ψ_leaf_0 = FT(-2e5 / 9800)

S_l_ini = inverse_water_retention_curve(
    canopy_retention_model,
    ψ_leaf_0,
    canopy_ν,
    canopy_S_s,
)
Y.canopy.hydraulics.ϑ_l.:1 .= augmented_liquid_fraction(canopy_ν, S_l_ini);
evaluate!(Y.canopy.energy.T, atmos.T, t0)

n = 120
saveat = Array(t0:(n * dt):tf)

timestepper = CTS.ARS343()
ode_algo = CTS.IMEXAlgorithm(
    timestepper,
    CTS.NewtonsMethod(
        max_iters = 1,
        update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
    ),
);

# Now set the initial values for the cache variables for the combined soil and plant model.

set_initial_cache! = make_set_initial_cache(land)
set_initial_cache!(p, Y, t0);

# Set the callbacks, which govern
# how often we save output, and how often we update
# the forcing data ("drivers")

sv = (;
    t = Array{Float64}(undef, length(saveat)),
    saveval = Array{NamedTuple}(undef, length(saveat)),
)
saving_cb = ClimaLand.NonInterpSavingCallback(sv, saveat)
model_drivers = ClimaLand.get_drivers(land)
updatefunc = ClimaLand.make_update_drivers(model_drivers)
updateat = Array(t0:1800:tf)
driver_cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc)
cb = SciMLBase.CallbackSet(driver_cb, saving_cb);

# Carry out the simulation
prob = SciMLBase.ODEProblem(
    CTS.ClimaODEFunction(
        T_exp! = exp_tendency!,
        T_imp! = SciMLBase.ODEFunction(imp_tendency!; jac_kwargs...),
        dss! = ClimaLand.dss!,
    ),
    Y,
    (t0, tf),
    p,
);
sol = SciMLBase.solve(
    prob,
    ode_algo;
    dt = dt,
    callback = cb,
    adaptive = false,
    saveat = saveat,
);

# # Plotting

# Now that we have both a soil and canopy model incorporated together, we will
# show how to plot some model data demonstrating the time series produced from
# each of these models. As before, we may plot the GPP of the system as well
# as transpiration showing fluxes in the canopy.

daily = sol.t ./ 3600 ./ 24
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
