# # Changing Soil Parameterizations
# In [Getting Started](@ref getting_started.md), we ran a simple soil model simulation
# using all of the default parameterizations and parameters.
# ClimaLand provides multiple options for many parameterizations;
# in this tutorial, we will demonstrate how to change a soil model parameterization.

# Specifically, we'll switch from using the default `CLMTwoBandSoilAlbedo` soil albedo parameterization
# to the `ConstantTwoBandSoilAlbedo` parameterization.
# In both cases, the soil albedo is defined in two bands (PAR and NIR), and can spatially vary or be set to scalar.
# In the more complex `CLMTwoBandSoilAlbedo` scheme, the soil albedo varies temporally due to a dependence
# on soil water content at the surface, via the effective saturation S(θ_sfc): α = α_wet*S + α_dry*(1-S)
# In contrast, the `ConstantTwoBandSoilAlbedo` parameterization does not vary with soil water content (or with time).
# The `ConstantTwoBandSoilAlbedo` parameterization is useful for cases where the soil albedo is
# known to be constant over time, such as in some idealized simulations or when using a fixed albedo value.

# CLM reference: Lawrence, P.J., and Chase, T.N. 2007. Representing a MODIS consistent land surface in the Community Land Model
# (CLM 3.0). J. Geophys. Res. 112:G01023. DOI:10.1029/2006JG000168.

# First we import the necessary Julia packages:
import ClimaParams as CP
using ClimaLand
using ClimaLand.Domains
using ClimaLand.Soil
import ClimaLand.Simulations: LandSimulation, solve!
import ClimaLand.Parameters as LP
import ClimaDiagnostics
using Dates
using CairoMakie, ClimaAnalysis, GeoMakie, Printf, StatsBase
import ClimaLand.LandSimVis as LandSimVis

# Choose a floating point precision, and get the parameter set, which holds constants used across CliMA models.
FT = Float32
earth_param_set = LP.LandParameters(FT);

# We will run this simulation on a column domain with 1 meter depth, at a lat/lon location near Pasadena, California.
zmax = FT(0)
zmin = FT(-1.0)
longlat = FT.((34.1, -118.1));
domain = Domains.Column(; zlim = (zmin, zmax), nelements = 10, longlat);
surface_space = domain.space.surface;

# We choose the initial and final simulation times as DateTimes, and a timestep in seconds.
start_date = DateTime(2008);
end_date = start_date + Second(60 * 60 * 72);
dt = 1000.0;

# The soil model takes in 2 forcing objects, atmosphere and radiation,
# which we read in from ERA5 data.
era5_ncdata_path =
    ClimaLand.Artifacts.era5_land_forcing_data2008_path(; lowres = true);
atmos, radiation = ClimaLand.prescribed_forcing_era5(
    era5_ncdata_path,
    surface_space,
    start_date,
    earth_param_set,
    FT,
);

# Now, we can create the EnergyHydrology model.

# First, let's set up the `ConstantTwoBandSoilAlbedo` parameterization.
# This parameterization requires the PAR and NIR albedo values, which can be scalars
# or fields that vary spatially. Here, we set them to constant values.
# Note that this constructor call is also overriding the default values for
# `PAR_albedo` and `NIR_albedo`, which are 0.2 and 0.4, respectively.
# To use the default values, we would simply call:
# `albedo = Soil.ConstantTwoBandSoilAlbedo{FT}()`
albedo = Soil.ConstantTwoBandSoilAlbedo{FT}(
    PAR_albedo = FT(0.2),
    NIR_albedo = FT(0.4),
);

# Now we can create the `EnergyHydrology` model with the specified albedo parameterization
# passed as a [keyword argument](https://docs.julialang.org/en/v1/manual/functions/#Keyword-Arguments).
model = Soil.EnergyHydrology{FT}(
    domain,
    (; atmos, radiation),
    earth_param_set;
    albedo,
);

# Define a function to set initial conditions for the prognostic variables.
function set_ic!(Y, p, t0, model)
    Y.soil.ϑ_l .= FT(0.24);
    Y.soil.θ_i .= FT(0.0);
    T = FT(290.15);
    ρc_s =
        Soil.volumetric_heat_capacity.(
            Y.soil.ϑ_l,
            Y.soil.θ_i,
            model.parameters.ρc_ds,
            model.parameters.earth_param_set,
        );
    Y.soil.ρe_int .=
        Soil.volumetric_internal_energy.(
            Y.soil.θ_i,
            ρc_s,
            T,
            model.parameters.earth_param_set,
        );
end

# Since we'll want to make some plots, let's set up an object to save the model output periodically.
diag_writer = ClimaDiagnostics.Writers.DictWriter();
diagnostics = ClimaLand.Diagnostics.default_diagnostics(
    model,
    start_date;
    output_vars = ["swc", "tsoil"],
    output_writer = diag_writer,
    average_period = :hourly,
);

# Now construct the LandSimulation object, which contains the model
# and additional timestepping information.
simulation = LandSimulation(
    start_date,
    end_date,
    dt,
    model;
    set_ic!,
    user_callbacks = (),
    diagnostics,
);

# Now we can run the simulation!
solve!(simulation);

# Let's plot some results, for example soil water content and soil temperature over time:
LandSimVis.make_timeseries(
    simulation;
    short_names = ["swc", "tsoil"],
    plot_stem_name = "soil_parameterizations",
);
# ![](swc_soil_parameterizations.png)
# ![](tsoil_soil_parameterizations.png)

# The same plots can be generated for the default case by re-running the simulation
# without providing a custom albedo parameterization.
