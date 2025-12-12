# # Changing Soil Parameterizations
# In [Getting Started](@ref docs/src/getting_started.md), we ran a simple soil model simulation
# using all of the default parameterizations and parameters.
# ClimaLand provides multiple options for many parameterizations;
# in this tutorial, we will demonstrate how to create a soil model with a
# non-default soil model parameterization.

# Specifically, we'll switch from using the default `CLMTwoBandSoilAlbedo` soil albedo parameterization
# to the `ConstantTwoBandSoilAlbedo` parameterization.
# In both cases, the soil albedo is defined in two bands (PAR and NIR), and can spatially vary or be set to scalar.
# In the more complex `CLMTwoBandSoilAlbedo` scheme, the soil albedo varies temporally due to a dependence
# on soil water content at the surface, via the effective saturation S(θ\_sfc): α = α\_wet*S + α\_dry*(1-S)
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
import ClimaLand.Parameters as LP
using Dates

# Choose a floating point precision, and get the parameter set, which holds constants used across CliMA models.
FT = Float32
toml_dict = LP.create_toml_dict(FT);

# We will run this simulation on a column domain with 1 meter depth, at a lat/lon location near Pasadena, California.
zmax = FT(0)
zmin = FT(-1.0)
longlat = FT.((-118.1, 34.1));
domain = Domains.Column(; zlim = (zmin, zmax), nelements = 10, longlat);
surface_space = domain.space.surface;

# We choose the start_date, which is required to setup the forcing,
# which in turn is required by the model.
start_date = DateTime(2008);

# The soil model takes in 2 forcing objects, atmosphere and radiation,
# which we read in from ERA5 data.
era5_ncdata_path =
    ClimaLand.Artifacts.era5_land_forcing_data2008_path(; lowres = true);
atmos, radiation = ClimaLand.prescribed_forcing_era5(
    era5_ncdata_path,
    surface_space,
    start_date,
    toml_dict,
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
model =
    Soil.EnergyHydrology{FT}(domain, (; atmos, radiation), toml_dict; albedo);
# That's it! Now that you have the model, you can create
# the simulation, solve it, and make your plots like you would for
# any other ClimaLand simulation.
