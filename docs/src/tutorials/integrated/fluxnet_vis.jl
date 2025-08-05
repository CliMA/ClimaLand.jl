# # Visualizing the output of a Fluxnet simulation

# We current support some default plotting functions for Fluxnet simulations.
# Please let us know if you would like to add other functionality; we
# are also implementing a method to write the simulation data to
# file so that you can use your plotting methods of choice.

# First, take a look at the fluxnet simulation tutorials for [SoilCanopyModel](docs/src/tutorials/integrated/soil_canopy_fluxnet_tutorial.md) or [LandModel](docs/src/tutorials/integrated/snowy_land_fluxnet_tutorial.md).
# To use the default plotting tools (`LandSimVis`),
# you need to start by importing the packages required:

# ``` julia
# using ClimaLand
# using CairoMakie, ClimaAnalysis, GeoMakie, Poppler_jll, Printf, StatsBase
# import ClimaLand.LandSimVis as LandSimVis;
# ```

# We assume you have setup your `simulation` and the solve step has finished.
# Two key functions are supported, `make_diurnal_timeseries` and
# `make_timeseries`. Both require the `simulation`:

# ``` julia
# LandSimVis.make_diurnal_timeseries(simulation);
# LandSimVis.make_timeseries(simulation);
# ```

# These make plots of each variable saved as a diagnostic during the
# simulation. Only the top layer of vertically resolved variables is
# plotted.
# To plot a subset of the variables, pass in a list of ClimaLand
# Diagnostic short_names. You can also optionally pass in a spinup
# date: only data after this date will be included in plotting. You
# can also pass in a directory and plot name to overwrite the defaults:

# ```julia
# LandSimVis.make_timeseries(
#    simulation;
#    savedir = path_to_dir,
#    plot_name = "myplot.pdf",
#    short_names = ["swc", "si", "swe"],
#    spinup_date = start_date + Day(20))
# ```

# To plot comparison data, you first need to get the comparison data
# ```julia
# comparison_data = FluxnetSimulations.get_comparison_data(site_ID, time_offset)
# LandSimVis.make_timeseries(
#    simulation;
#    short_names = ["swc", "si", "swe"],
#    comparison_data)
# ```

# The diurnal timeseries plotted is an average diurnal timeseries. We
# interpolate the data or simulation output to an hourly interval, and
# then compute the average of each hour of the day over the entire data.
