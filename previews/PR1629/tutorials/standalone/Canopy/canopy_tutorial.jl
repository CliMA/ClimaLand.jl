# # Introduction to the Canopy Model

# This tutorial shows how to instantiate and run a simulation of the
# canopy biophysics model in ClimaLand. A
# [`CanopyModel`](@ref "Canopy Model and Parameters")
# including all component
# models is initialized, then an example simulation is run. The initial conditions,
# atmospheric and radiative flux conditions, and canopy properties are set up
# to match those observed at the US-MOz flux tower, a flux tower located within
# an oak-hickory forest in Ozark, Missouri, USA. See [Wang2021](@citet)
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
import ClimaLand.Domains
using ClimaLand.Canopy
using ClimaLand.Canopy.PlantHydraulics
import ClimaLand.Simulations: LandSimulation, solve!
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaDiagnostics
using DelimitedFiles
import ClimaLand.FluxnetSimulations as FluxnetSimulations
using CairoMakie, ClimaAnalysis, GeoMakie, Printf, StatsBase
import ClimaLand.LandSimVis as LandSimVis

# Define the floating point precision desired (64 or 32 bit), and get the
# parameter set holding constants used across CliMA Models:

const FT = Float32;
toml_dict = LP.create_toml_dict(FT);

# We will use prescribed atmospheric and radiative forcing from the
# US-MOz tower.
site_ID = "US-MOz";
site_ID_val = FluxnetSimulations.replace_hyphen(site_ID);
# Get the latitude and longitude in degrees, as well as the
# time offset in hours of local time from UTC
(; time_offset, lat, long) =
    FluxnetSimulations.get_location(FT, Val(site_ID_val));
# Get the height of the sensors in m
(; atmos_h) = FluxnetSimulations.get_fluxtower_height(FT, Val(site_ID_val));
# Set a start and stop date of the simulation in UTC, as well as
# a timestep in seconds
start_date = DateTime("2010-05-01", "yyyy-mm-dd")
stop_date = DateTime("2010-09-01", "yyyy-mm-dd")
dt = 450.0

# Site latitude and longitude
lat = FT(38.7441) # degree
long = FT(-92.2000) # degree

# # Setup the Canopy Model

# We want to simulate a vegetative canopy in standalone mode, without coupling
# the canopy to atmospheric or soil physics models, so we choose a
# [`CanopyModel`](@ref "Canopy Model and Parameters").
# Here we will use the default parameterizations and parameters for ease of setting up
# the model, but these can be overridden by constructing and passing canopy
# components to the `CanopyModel` constructor. This will be explored in a later tutorial.

# Now we define the parameters of the model domain. These values are needed
# by some of the component models. Here we are performing a 1-dimensional
# simulation in a `Point` domain and will use
# single stem and leaf compartments, but for 2D simulations, the parameters of
# the [`domain`](@ref "Domain Tutorial")
# would change.
domain = ClimaLand.Domains.Point(; z_sfc = FT(0.0), longlat = (long, lat));

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
    toml_dict,
    FT,
);


# For this canopy, we are running in standalone mode, which means we need to
# use a prescribed soil driver, defined as follows:
θ_soil = FT(0.47)
T_soil = FT(298.0)
ground = PrescribedGroundConditions{FT}(;
    α_PAR = FT(0.2),
    α_NIR = FT(0.4),
    T = TimeVaryingInput(t -> T_soil),
    θ = TimeVaryingInput(t -> θ_soil),
    ϵ = FT(0.99),
);
forcing = (; atmos, radiation, ground);

# Now we read in time-varying LAI from a global MODIS dataset.
surface_space = domain.space.surface;
LAI =
    ClimaLand.Canopy.prescribed_lai_modis(surface_space, start_date, stop_date);
# Get the maximum LAI at this site over the first year of the simulation
maxLAI = FluxnetSimulations.get_maxLAI_at_site(start_date, lat, long);

# Construct radiative transfer model, overwriting some default parameters.
radiation_parameters = (;
    Ω = FT(0.69),
    G_Function = ConstantGFunction(FT(0.5)),
    α_PAR_leaf = FT(0.1),
    α_NIR_leaf = FT(0.45),
    τ_PAR_leaf = FT(0.05),
    τ_NIR_leaf = FT(0.25),
)
radiative_transfer = Canopy.TwoStreamModel{FT}(
    domain,
    toml_dict;
    radiation_parameters,
    ϵ_canopy = FT(0.99),
);

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
    toml_dict;
    n_stem,
    n_leaf,
    h_stem,
    h_leaf,
    ν = plant_ν,
    S_s = plant_S_s,
    conductivity_model,
);
rooting_depth = FT(1)
height = FT(18)
biomass =
    Canopy.PrescribedBiomassModel{FT}(; LAI, SAI, RAI, rooting_depth, height);

# Construct the conductance model.
conductance = Canopy.MedlynConductanceModel{FT}(domain, toml_dict; g1 = FT(141));

# Construct the photosynthesis model.
photosynthesis_parameters = (; is_c3 = FT(1), Vcmax25 = FT(5e-5))
photosynthesis =
    Canopy.FarquharModel{FT}(domain, toml_dict; photosynthesis_parameters);


# Set up the canopy model using defaults for all parameterizations and parameters,
# except for the hydraulics model defined above.
canopy = ClimaLand.Canopy.CanopyModel{FT}(
    domain,
    forcing,
    LAI,
    toml_dict;
    hydraulics,
    radiative_transfer,
    conductance,
    photosynthesis,
    biomass,
);

# Provide initial conditions for the canopy hydraulics model
function set_ic!(Y, p, t0, model)
    (; retention_model, ν, S_s) = model.hydraulics.parameters
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
        Y.canopy.hydraulics.ϑ_l.:($i) .=
            augmented_liquid_fraction.(ν, S_l_ini[i])
    end
    evaluate!(Y.canopy.energy.T, atmos.T, t0)
end

# Set up the diagnostics writer, which will save model variables
# throughout the course of the simulation.
diag_writer = ClimaDiagnostics.Writers.DictWriter();
diagnostics = ClimaLand.Diagnostics.default_diagnostics(
    canopy,
    start_date;
    output_vars = ["gpp", "trans"],
    output_writer = diag_writer,
    reduction_period = :hourly,
);



# Select a timestepping algorithm and setup the ODE problem.
timestepper = CTS.ARS111();
ode_algo = CTS.IMEXAlgorithm(
    timestepper,
    CTS.NewtonsMethod(
        max_iters = 6,
        update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
    ),
);

# Create the LandSimulation object, which will also create and initialize the state vectors,
# the cache, the driver callbacks, and set the initial conditions.
simulation = LandSimulation(
    start_date,
    stop_date,
    dt,
    canopy;
    set_ic!,
    updateat = Second(1800),
    timestepper = ode_algo,
    user_callbacks = (),
    diagnostics,
);

# Now we can solve the simulation, which will run the model forward in time.
sol = solve!(simulation);

# # Create some plots

# We can now plot the data produced in the simulation. For example, GPP and transpiration:
LandSimVis.make_diurnal_timeseries(
    simulation;
    short_names = ["gpp", "trans"],
    plot_stem_name = "canopy",
)
# ![](gpp_canopy.png)
# ![](trans_canopy.png)
