# # Site-level simulations with the bundled ERA5 artifact

# ClimaLand can run as a single column anywhere on Earth, and it ships a global ERA5
# forcing artifact that downloads automatically. So you can simulate a point of your
# choosing with no data preparation at all -- no Copernicus account, no preprocessing.
#
# The catch is resolution. The bundled artifact is **8° × 8° and covers only 2008**, and
# `prescribed_forcing_era5` regrids it by *nearest neighbour*. Your column therefore
# inherits the climate of whichever cell centre is closest, which can be hundreds of
# kilometres away. This tutorial shows how to choose a site where that is acceptable, how
# to check, and what the model does at two deliberately contrasting locations.
#
# For a site where the 8° cell will not do, the
# site-level ERA5 tutorial covers supplying your own
# high-resolution data. For forcing with no data dependency whatsoever, see the
# [idealized tutorial](@ref "An idealized year-long land model simulation").

# # Preliminary setup

using Dates
import ClimaComms
ClimaComms.@import_required_backends
import ClimaParams as CP
using ClimaLand
using ClimaLand.Domains: Column
import ClimaLand.Parameters as LP
import ClimaLand.Simulations: LandSimulation, solve!
import ClimaUtilities.TimeManager: date
using Printf, Statistics
using CairoMakie, ClimaAnalysis, GeoMakie
import ClimaLand.LandSimVis as LandSimVis;

const FT = Float64;
context = ClimaComms.context();
toml_dict = LP.create_toml_dict(FT);

# # Choosing a site the 8° grid can represent
#
# The artifact's cell centres sit at longitudes 0, 8, 16, … 352 and latitudes −90, −82,
# … −2, 6, 14, 22, 30, … Two conditions make a site usable:
#
# 1. **The site coincides with a cell centre**, so nearest-neighbour regridding introduces
#    no displacement at all.
# 2. **The surrounding 8° region is climatically homogeneous**, so the cell average is a
#    meaningful description of the site rather than a blend of unlike climates.
#
# Both matter. Here is what the artifact actually contains at four candidate sites,
# read directly out of the file:
#
# | site | snaps to | displacement | 2 m temp | pressure | rain | snow |
# |:---|:---|---:|---:|---:|---:|---:|
# | Central Amazon (−64, −2) | (−64, −2) | 0.0° | 23.3 °C | 95.0 kPa | 2791 mm | 0 |
# | Boreal Manitoba (−96, 54) | (−96, 54) | 0.0° | −6.3 °C | 98.7 kPa | 508 mm | 266 mm |
# | Congo basin (24, −2) | (24, −2) | 0.0° | 24.4 °C | 95.4 kPa | 1369 mm | 0 |
# | Sahel, Niger (8, 14) | (8, 14) | 0.0° | 28.5 °C | 94.4 kPa | **41 mm** | 0 |
#
# The Sahel row shows why condition 2 is not optional. The site sits exactly on a cell
# centre, yet the cell spans roughly 10°N–18°N, which at 8°E is mostly Sahara. The result
# is 41 mm/yr -- a desert, not the Sahel. Failing condition 2 is just as damaging as
# failing condition 1.
#
# !!! warning "What happens when a site fails both conditions"
#     IIT Kanpur (80.23°E, 26.51°N) snaps **3.5° away** to the cell centred at (80°E,
#     30°N) -- the Tibetan Plateau. That cell has an annual mean temperature of
#     **−9.7 °C**, a surface pressure of **52.4 kPa** (about 5000 m elevation), and
#     **292 of its 451 mm/yr as snow**. The correct values for the Indo-Gangetic plain are
#     roughly +26 °C, 99 kPa and no snow. A run at Kanpur's true coordinates with this
#     artifact simulates a frozen alpine column, produces perfectly plausible-looking
#     plots, and never warns you. **Always check temperature and pressure against what you
#     expect for the place you think you are simulating.**

# # A reusable site setup
#
# Everything below is independent of location, so we wrap it in a function and call it
# twice. Note the 15 m soil depth with `nelements = 15`: `LandSimulation` defaults to
# initial conditions from a spun-up global file with that geometry.

function run_site(name, longlat; start_date, stop_date, Δt = 450.0)
    domain = Column(;
        zlim = FT.((-15, 0)),
        nelements = 15,
        dz_tuple = FT.((3, 0.05)),
        longlat,
    )
    surface_space = domain.space.surface

    ## use_lowres_forcing selects the bundled 8-degree 2008 artifact
    forcing = ClimaLand.prescribed_forcing_era5(
        start_date,
        stop_date,
        surface_space,
        toml_dict,
        FT;
        max_wind_speed = 25.0,
        context,
        use_lowres_forcing = true,
    )
    LAI = ClimaLand.Canopy.prescribed_lai_modis(
        surface_space,
        start_date,
        stop_date,
    )
    model = ClimaLand.LandModel{FT}(
        forcing,
        LAI,
        toml_dict,
        domain,
        Δt;
        prognostic_land_components = (:canopy, :snow, :soil, :soilco2),
    )
    output_vars = ["swc", "tsoil", "lhf", "shf", "gpp", "lai", "swe", "et"]
    diagnostics = ClimaLand.default_diagnostics(
        model,
        start_date,
        ".";
        reduction_period = :daily,
        output_vars,
    )
    simulation = LandSimulation(start_date, stop_date, Δt, model; diagnostics)
    solve!(simulation)
    return simulation
end;

# The artifact only contains 2008, so that is the year we simulate. Outside 2008 the same
# year is reused periodically, which is fine for a spin-up but means multi-year runs carry
# no interannual variability.

start_date = DateTime(2008, 1, 1)
stop_date = DateTime(2008, 12, 31);

# `max_wind_speed` clips a known ERA5 artefact: occasional spurious 10 m wind spikes that
# would otherwise produce enormous surface fluxes.

# # Two contrasting sites
#
# The central Amazon is energy-limited, wet year-round and evergreen. Boreal Manitoba is
# temperature-limited, strongly seasonal and snow-covered for much of the year. The same
# model code, the same forcing dataset, and the only difference is two numbers.

amazon = run_site("amazon", FT.((-64.0, -2.0)); start_date, stop_date);

# Each site-year takes on the order of a minute of timestepping.

boreal = run_site("boreal", FT.((-96.0, 54.0)); start_date, stop_date);

# # Results
#
# Seasonal means from an actual run of this tutorial:
#
# | variable | Amazon DJF | Amazon JJA | Boreal DJF | Boreal JJA |
# |:---|---:|---:|---:|---:|
# | soil temperature (K) | 296.0 | 295.4 | 266.3 | 287.5 |
# | soil moisture (m³/m³) | 0.53 | 0.57 | 0.11 | 0.26 |
# | latent heat flux (W/m²) | 91.1 | 85.8 | 0.15 | 74.7 |
# | sensible heat flux (W/m²) | 41.7 | 33.3 | −9.1 | 42.6 |
# | GPP (mol CO₂ m⁻² s⁻¹) | 6.5e-6 | 9.0e-6 | 0.0 | 4.2e-6 |
# | LAI (m²/m²) | 5.23 | 6.05 | 0.23 | 1.92 |
# | snow water equivalent (m) | 0.0 | 0.0 | 0.091 | 3.7e-5 |
#
# Four things worth drawing out:
#
# - **The Amazon has essentially no seasonal temperature cycle** -- 295.4 K against
#   296.0 K -- while the boreal site swings 21 K. Latent heat at the Amazon site barely
#   changes between seasons, because water and energy are both always available: it is
#   energy-limited, and the energy supply is nearly constant.
# - **The boreal site exercises the snow model**, which the tropical site never touches.
#   Snow water equivalent reaches 91 mm in winter and peaks at 150 mm. If you want to
#   exercise snow physics, you need a site like this one.
# - **Boreal sensible heat flux goes negative in winter** (−9.1 W/m²), meaning heat flows
#   *from* the atmosphere *into* the surface. That is the correct behaviour over a cold
#   snow surface under a stable boundary layer, and it is a good sign the surface energy
#   balance is behaving.
# - **Boreal GPP is exactly zero in winter.** The canopy is dormant, LAI has dropped to
#   0.23, and there is no carbon uptake at all -- whereas the Amazon assimilates all year.

# `make_timeseries` writes into `savedir` but does not create it, so make the directories
# first.
mkpath("amazon")
mkpath("boreal")
LandSimVis.make_timeseries(amazon; savedir = "amazon");
LandSimVis.make_timeseries(boreal; savedir = "boreal");

# # Caveats
#
# - **You are simulating an 8° cell average, not a point.** Even when a site sits exactly
#   on a cell centre, the forcing describes a region roughly 900 km across. Comparing
#   these runs against point observations such as a flux tower is not meaningful.
# - **Only 2008 exists.** Longer runs repeat it, so there is no interannual variability.
# - **Spatially-varying parameters are *not* coarse.** Soil texture, CLM canopy properties
#   and MODIS LAI are read at the true coordinates from finer maps, so the parameters are
#   site-specific even though the forcing is not. That mismatch is worth remembering: a
#   rainforest's soil and vegetation parameters driven by an 8° regional-average climate.
# - **Check before you trust.** The Kanpur and Sahel cases above both look entirely
#   reasonable in the output and are both wrong. Reading the forcing values at your site is
#   a few seconds of work and is the only thing standing between you and a confidently
#   wrong simulation.
