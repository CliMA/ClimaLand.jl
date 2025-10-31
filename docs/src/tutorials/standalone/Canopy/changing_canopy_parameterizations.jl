# # Changing Canopy Parameterizations
# In [Default Canopy](@ref "Default Canopy Model Tutorial"), we ran a simple canopy model simulation
# using all of the default parameterizations and parameters.
# ClimaLand provides multiple options for many parameterizations;
# in this tutorial, we will demonstrate how to change a canopy model parameterization.

# Specifically, we'll use non-default parameterizations for two canopy components:
# - energy model: from the default `BigLeafEnergyModel` to `PrescribedCanopyTempModel`
# - radiative transfer: from the default `TwoStreamModel` to `BeerLambertModel`

# First import the Julia packages we'll need.
import ClimaParams as CP
using ClimaUtilities.TimeVaryingInputs: TimeVaryingInput
using ClimaLand
using ClimaLand.Domains
using ClimaLand.Canopy
import ClimaLand.Simulations: LandSimulation, solve!
import ClimaLand.Parameters as LP
using Dates
import ClimaDiagnostics
using CairoMakie, ClimaAnalysis, GeoMakie, Printf, StatsBase
import ClimaLand.LandSimVis as LandSimVis

# Choose a floating point precision, and get the parameter set,
# which holds constants used across CliMA models.
FT = Float32
toml_dict = LP.create_toml_dict(FT);

# We will run this simulation on a point domain at a lat/lon location
# near Yellowstone National Park.
# This is different from the soil example, which ran on a column, because the
# soil model has depth and the canopy model does not.
longlat = FT.((-110.6, 44.6))
domain = Domains.Point(; z_sfc = FT(0.0), longlat);
surface_space = domain.space.surface;

# We choose the initial and final simulation times as DatesTimes, and a timestep in seconds.
start_date = DateTime(2008);
stop_date = start_date + Second(60 * 60 * 72);
dt = 900.0;

# Whereas the soil model takes in 2 forcing objects (atmosphere and radiation),
# the canopy takes in 3 (atmosphere, radiation, and ground). Here we read in the
# first two from ERA5 data, and specify that the following ground conditions will
# be prescribed: emissivity, albedo, temperature, and soil moisture.
# We also set up a constant leaf area index (LAI); for an example reading LAI from
# MODIS data, please see the [canopy tutorial](@ref "Default Canopy Model Tutorial").
# This differs from the soil example because we have the extra inputs of the ground conditions and LAI.
atmos, radiation = ClimaLand.prescribed_forcing_era5(
    start_date,
    stop_date,
    surface_space,
    toml_dict,
    FT;
    use_lowres_forcing = true,
);
ground = PrescribedGroundConditions{FT}();
LAI = TimeVaryingInput((t) -> FT(1.0));

# Now, we can create the canopy model.

# First, let's set up the `PrescribedCanopyTempModel` parameterization.
# This parameterization acts as a flag to signal that the atmosphere temperature
# from the prescribed forcing should be used as the canopy temperature.
# Since this only involves retrieving the temperature from the atmosphere forcing,
# we do not need to provide any additional parameters, and this setup is straightforward.
energy = Canopy.PrescribedCanopyTempModel{FT}();

# Now let's set up the `BeerLambertModel` radiative transfer parameterization.
# We'll explore three different ways to construct this parameterization.
# The first way is to use the `BeerLambertModel` constructor with the default parameters.
# This method requires the simulation domain as it reads in radiation
# parameters from a map of CLM parameters by default.
radiative_transfer = Canopy.BeerLambertModel{FT}(domain, toml_dict);

# Alternatively, we could use the same constructor but provide custom values for
# each of the parameters.
# For example, we might want to use the global maps for albedo, but
# explore how the results change when we assume a spherical
# distribution of leaves (G = 0.5) and no clumping (Ω = 1).
G_Function = Canopy.ConstantGFunction(FT(0.5)); # leaf angle distribution value 0.5
Ω = 1; # clumping index
α_PAR_leaf = 0.1;
α_NIR_leaf = 0.4;
radiation_parameters = (; G_Function, Ω, α_PAR_leaf, α_NIR_leaf);
radiative_transfer =
    Canopy.BeerLambertModel{FT}(domain, toml_dict; radiation_parameters);

# If you want to overwrite all of the parameters, you do not need to use
# global maps, and therefore do not need the domain.
# In a case like this, we may want to use a different `BeerLambertModel` constructor,
# which takes the parameters object directly:
radiative_transfer_parameters = Canopy.BeerLambertParameters(
    toml_dict;
    G_Function,
    Ω,
    α_PAR_leaf,
    α_NIR_leaf,
);
radiative_transfer = Canopy.BeerLambertModel(radiative_transfer_parameters);
# Note these parameters must all be scalars (or a single instance of a G\_Function) or fields of floats and a field of a G\_Function.

# All three of these methods will construct a `BeerLambertModel`;
# the first will use the default parameters, while the second and third use the custom parameters.
# The method you choose will depend on your use case and whether you want to
# use the default parameters or provide custom ones.
# The set of available constructors for all ClimaLand models can be found in the "APIs" section of the documentation.

# Set up optimal LAI model
lai_model =
    Canopy.OptimalLAIModel{FT}(Canopy.OptimalLAIParameters{FT}(toml_dict));

# Now we can create the `CanopyModel` model with the specified energy and radiative transfer
# parameterizations passed as [keyword arguments](@extref Julia Keyword-Arguments).
model = Canopy.CanopyModel{FT}(
    domain,
    (; atmos, radiation, ground),
    LAI,
    toml_dict;
    energy,
    radiative_transfer,
    lai_model,
);

# Define a function to set initial conditions for the prognostic variables.
# Since these are specific to the model physics, the contents here differ from
# the soil example, but the function structure remains the same.
# The variables initialized here are described in the Model Equations section
# of the documentation.
# Note that here we previously set the initial condition for the canopy energy model's
# prognostic temperature. Since we're using the `PrescribedCanopyTempModel`,
# we do not need to set it here.
function set_ic!(Y, p, t0, model)
    ψ_leaf_0 = FT(-2e5 / 9800)
    (; retention_model, ν, S_s) = model.hydraulics.parameters
    S_l_ini = Canopy.PlantHydraulics.inverse_water_retention_curve(
        retention_model,
        ψ_leaf_0,
        ν,
        S_s,
    )
    Y.canopy.hydraulics.ϑ_l.:1 .=
        Canopy.PlantHydraulics.augmented_liquid_fraction.(ν, S_l_ini)
end

# Since we'll want to make some plots, let's set up an object to save the
# model output periodically, as we did for the soil tutorial.
diag_writer = ClimaDiagnostics.Writers.DictWriter();
diagnostics = ClimaLand.Diagnostics.default_diagnostics(
    model,
    start_date;
    output_vars = ["ct", "trans"],
    output_writer = diag_writer,
    reduction_period = :hourly,
);

# Now construct the `LandSimulation` object, which contains the model
# and additional timestepping information. This is identical to the soil example.
simulation = LandSimulation(
    start_date,
    stop_date,
    dt,
    model;
    set_ic!,
    updateat = Second(dt),
    user_callbacks = (),
    diagnostics,
);

# Now we can run the simulation!
solve!(simulation);

# Let's plot some results, for example diurnally averaged canopy temperature and transpiration over time:
LandSimVis.make_diurnal_timeseries(
    simulation;
    short_names = ["ct", "trans"],
    plot_stem_name = "canopy_parameterizations",
);
# ![](ct_canopy_parameterizations.png)
# ![](trans_canopy_parameterizations.png)

# Now you can compare these plots to those generated in the default canopy tutorial.
# How are the results different? How are they the same?

# Atmospheric forcing data citation:
# [Hersbach2020](@citet)
