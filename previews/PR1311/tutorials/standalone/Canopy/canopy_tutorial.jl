# # Introduction to the Canopy Model

# This tutorial shows how to instantiate and run a simulation of the
# canopy biophysics model in ClimaLand. A
# [`CanopyModel`](https://clima.github.io/ClimaLand.jl/dev/APIs/canopy/Canopy/#Canopy-Model-Structs)
# including all component
# models is initialized, then an example simulation is run. The initial conditions,
# atmospheric and radiative flux conditions, and canopy properties are set up
# to match those observed at the US-MOz flux tower, a flux tower located within
# an oak-hickory forest in Ozark, Missouri, USA. See [Wang et al. 2021](https://doi.org/10.5194/gmd-14-6741-2021)
# for details on the site and canopy parameters.

# The canopy biophysics model in ClimaLand combines a photosynthesis model with a
# canopy radiative transfer scheme, plant hydraulics model, and stomatal
# conductance model, placing them under either prescribed or simulated (as in a
# full Earth System Model) atmospheric and radiative flux conditions.

# ClimaLand supports either Beer-Lambert law or a Two-Stream model for radiative
# transfer. For this tutorial, we will use the Beer-Lambert law,
# in which the intensity of light absorbed is a negative exponential function of
# depth in the canopy and an exinction coefficient determined by optical depth.

# The model of photosynthesis in Clima Land is the Farquhar Model in which GPP is
# calculated based on C3 and C4 photosynthesis, which determines potential
# leaf-level photosynthesis.

# The plant hydraulics model in ClimaLand solves for the water content within
# bulk root-stem-canopy system using Richards equation discretized into an
# arbitrary number of layers. The water content is related to the water
# potential using a retention curve relationship, and the water potential is
# used to simulate the effect moisture stress has on transpiration and GPP.

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
using StaticArrays
using ClimaLand
using ClimaLand.Domains: Point
using ClimaLand.Canopy
using ClimaLand.Canopy.PlantHydraulics
import ClimaLand
import ClimaLand.Parameters as LP
using DelimitedFiles
import ClimaLand.FluxnetSimulations as FluxnetSimulations

# Define the floating point precision desired (64 or 32 bit), and get the
# parameter set holding constants used across CliMA Models:

const FT = Float32;
earth_param_set = LP.LandParameters(FT);

# First provide some information about the site:
# Timezone (offset from UTC in hrs)
time_offset = 7
start_date = DateTime(2010) + Hour(time_offset)

# Site latitude and longitude
lat = FT(38.7441) # degree
long = FT(-92.2000) # degree

# Height of the sensor at the site
atmos_h = FT(32)
# Site ID
site_ID = "US-MOz";

# # Setup the Canopy Model

# We want to simulate a vegetative canopy in standalone mode, without coupling
# the canopy to atmospheric or soil physics models, so we choose a
# [`CanopyModel`](https://clima.github.io/ClimaLand.jl/dev/APIs/canopy/Canopy/#Canopy-Model-Structs).
# Here we will use the default parameterizations and parameters for ease of setting up
# the model, but these can be overridden by constructing and passing canopy
# components to the `CanopyModel` constructor. This will be explored in a later tutorial.

# Now we define the parameters of the model domain. These values are needed
# by some of the component models. Here we are performing a 1-dimensional
# simulation in a `Point` domain and will use
# single stem and leaf compartments, but for 2D simulations, the parameters of
# the [`domain`](https://clima.github.io/ClimaLand.jl/dev/APIs/shared_utilities/#Domains)
# would change.
domain = Point(; z_sfc = FT(0.0), longlat = (long, lat));

# Select a time range to perform time stepping over, and a dt. As usual,
# the timestep depends on the problem you are solving, the accuracy of the
# solution required, and the timestepping algorithm you are using.
t0 = 0.0
N_days = 364
tf = t0 + 3600 * 24 * N_days
dt = 225.0;

# We will be using prescribed atmospheric and radiative drivers from the
# US-MOz tower, which we read in here. We are using prescribed
# atmospheric and radiative flux conditions, but it is also possible to couple
# the simulation with atmospheric and radiative flux models.
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


# For this canopy, we are running in standalone mode, which means we need to
# use a prescribed soil driver, defined as follows:
ψ_soil = FT(0.0)
T_soil = FT(298.0)
ground = PrescribedGroundConditions{FT}(;
    α_PAR = FT(0.2),
    α_NIR = FT(0.4),
    T = TimeVaryingInput(t -> T_soil),
    ψ = TimeVaryingInput(t -> ψ_soil),
    ϵ = FT(0.99),
);
forcing = (; atmos, radiation, ground);

# Now we read in time-varying LAI from a global MODIS dataset.
surface_space = domain.space.surface;
modis_lai_ncdata_path = ClimaLand.Artifacts.modis_lai_multiyear_paths(
    start_date = start_date + Second(t0),
    end_date = start_date + Second(t0) + Second(tf),
)
LAI = ClimaLand.prescribed_lai_modis(
    modis_lai_ncdata_path,
    surface_space,
    start_date,
);
# Get the maximum LAI at this site over the first year of the simulation
maxLAI =
    FluxnetSimulations.get_maxLAI_at_site(modis_lai_ncdata_path[1], lat, long);

# Construct radiative transfer model, overwriting some default parameters.
radiation_parameters = (;
    Ω = FT(0.69),
    G_Function = ConstantGFunction(FT(0.5)),
    α_PAR_leaf = FT(0.1),
    α_NIR_leaf = FT(0.45),
    τ_PAR_leaf = FT(0.05),
    τ_NIR_leaf = FT(0.25),
)
radiative_transfer =
    Canopy.TwoStreamModel{FT}(domain; radiation_parameters, ϵ_canopy = FT(0.99));

# Construct canopy hydraulics with 1 stem and 1 leaf compartment.
# By default, the model is constructed with a single leaf compartment and no stem.
n_stem = Int64(1)
n_leaf = Int64(1)
h_stem = FT(9)
h_leaf = FT(9.5)
SAI = FT(0.00242)
f_root_to_shoot = FT(3.5)
RAI = FT((SAI + maxLAI) * f_root_to_shoot)
plant_ν = FT(0.7)
plant_S_s = FT(1e-2 * 0.0098)
K_sat_plant = FT(1.8e-8)
ψ63 = FT(-4 / 0.0098)
Weibull_param = FT(4)
conductivity_model =
    PlantHydraulics.Weibull{FT}(K_sat_plant, ψ63, Weibull_param)
hydraulics = Canopy.PlantHydraulicsModel{FT}(
    domain,
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
    rooting_depth = FT(1),
);

# Construct the conductance model.
conductance = Canopy.MedlynConductanceModel{FT}(domain; g1 = FT(141));

# Construct the photosynthesis model.
photosynthesis_parameters = (; is_c3 = FT(1), Vcmax25 = FT(5e-5))
photosynthesis = Canopy.FarquharModel{FT}(domain; photosynthesis_parameters)


# Set up the canopy model using defaults for all parameterizations and parameters,
# except for the hydraulics model defined above.
canopy = ClimaLand.Canopy.CanopyModel{FT}(
    domain,
    forcing,
    LAI,
    earth_param_set;
    hydraulics,
    radiative_transfer,
    conductance,
    photosynthesis,
);

# Initialize the state vectors and obtain the model coordinates, then get the
# explicit time stepping tendency that updates auxiliary and prognostic
# variables that are stepped explicitly.

Y, p, coords = ClimaLand.initialize(canopy);
exp_tendency! = make_exp_tendency(canopy);
imp_tendency! = make_imp_tendency(canopy);
jacobian! = make_jacobian(canopy);
jac_kwargs =
    (; jac_prototype = ClimaLand.FieldMatrixWithSolver(Y), Wfact = jacobian!);

# Provide initial conditions for the canopy hydraulics model
(; retention_model, ν, S_s) = canopy.hydraulics.parameters
ψ_stem_0 = FT(-1e5 / 9800)
ψ_leaf_0 = FT(-2e5 / 9800)

S_l_ini =
    inverse_water_retention_curve.(
        retention_model,
        [ψ_stem_0, ψ_leaf_0],
        ν,
        S_s,
    )
for i in 1:2
    Y.canopy.hydraulics.ϑ_l.:($i) .= augmented_liquid_fraction.(ν, S_l_ini[i])
end
evaluate!(Y.canopy.energy.T, atmos.T, t0)

# Initialize the cache variables for the canopy using the initial
# conditions and initial time.

set_initial_cache! = make_set_initial_cache(canopy)
set_initial_cache!(p, Y, t0);

# Allocate the struct which stores the saved auxiliary state
# and create the callback which saves it at each element in saveat.

n = 16
saveat = Array(t0:(n * dt):tf)
sv = (;
    t = Array{Float64}(undef, length(saveat)),
    saveval = Array{NamedTuple}(undef, length(saveat)),
)
saving_cb = ClimaLand.NonInterpSavingCallback(sv, saveat);

# Create the callback function which updates the forcing variables,
# or drivers.
updateat = Array(t0:1800:tf)
model_drivers = ClimaLand.get_drivers(canopy)
updatefunc = ClimaLand.make_update_drivers(model_drivers)
driver_cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc)
cb = SciMLBase.CallbackSet(driver_cb, saving_cb);


# Select a timestepping algorithm and setup the ODE problem.
timestepper = CTS.ARS111();
ode_algo = CTS.IMEXAlgorithm(
    timestepper,
    CTS.NewtonsMethod(
        max_iters = 6,
        update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
    ),
);

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

# Now, we can solve the problem and store the model data in the saveat array,
# using [`SciMLBase.jl`](https://github.com/SciML/SciMLBase.jl) and
# [`ClimaTimeSteppers.jl`](https://github.com/CliMA/ClimaTimeSteppers.jl).

sol = SciMLBase.solve(prob, ode_algo; dt = dt, callback = cb, saveat = saveat);

# # Create some plots

# We can now plot the data produced in the simulation. For example, GPP:

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
savefig("ozark_standalone_canopy_test.png");
# ![](ozark_standalone_canopy_test.png)
