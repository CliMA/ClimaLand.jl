# # Evaluating the integrated land model globally

# The [global snow/soil/canopy tutorial](@ref "Snow, soil, canopy") shows how to
# *configure* a global integrated `LandModel`. This tutorial runs one end to end and then
# evaluates it, which is the global counterpart of the
# site-level ERA5 tutorial: same fully integrated model, same low-resolution ERA5 forcing,
# but on the sphere, and finishing with spatial diagnostics rather than a point timeseries.
#
# Everything here runs with artifacts ClimaLand downloads automatically -- no
# user-supplied data.
#
# !!! note "Resolution and cost"
#     At `nelements = (20, 7)` this costs roughly 10 s of wall time per simulated day, so a
#     full year is about an hour. This tutorial runs one month to stay follow-along; the
#     observational leaderboard needs at least a year to be meaningful, and the final
#     section explains how to get there.

# # Preliminary setup

using Dates
import ClimaComms
ClimaComms.@import_required_backends
import ClimaParams as CP
using ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.Simulations: LandSimulation, solve!
using CairoMakie, ClimaAnalysis, GeoMakie
import ClimaLand.LandSimVis as LandSimVis;

const FT = Float64;
context = ClimaComms.context();
toml_dict = LP.create_toml_dict(FT);

# # The global domain
#
# `global_domain` builds a spherical shell with a land/sea mask. `nelements` is
# `(horizontal, vertical)`: the first number sets the spectral element count per cube
# panel, the second the number of soil layers. Production runs use `(101, 15)`; we use a
# much coarser grid so this completes quickly.

nelements = (20, 7)
domain = ClimaLand.Domains.global_domain(FT; context, nelements);

# Unlike a column, a global domain has a land/sea mask, which the default NaN-check
# callback uses to ignore ocean points:

@show ClimaLand.Domains.landsea_mask(domain) !== nothing;

# # Forcing
#
# `use_lowres_forcing = true` selects the 8° × 8° ERA5 artifact for 2008. For a *global*
# run this coarseness is a reasonable trade-off and the artifact is small enough to
# download on demand -- which is exactly the opposite of the site-level case, where an 8°
# cell can place your column hundreds of kilometres away in a different climate entirely.
#
# The data exists only for 2008, so multi-year runs reuse it periodically via the default
# `LinearInterpolation(PeriodicCalendar())`.
#
# !!! warning "If you have set an ERA5 artifact override"
#     Had we omitted `use_lowres_forcing`, this would read the high-resolution artifact --
#     which, if you followed the site-level tutorial and pointed
#     `Overrides.toml` at a small regional box, would silently supply
#     that region's data extrapolated across the whole globe. Comment the override out
#     before global runs, or keep `use_lowres_forcing = true` here.

start_date = DateTime(2008, 1, 1)
stop_date = DateTime(2008, 2, 1)
Δt = 450.0;

forcing = ClimaLand.prescribed_forcing_era5(
    start_date,
    stop_date,
    domain.space.surface,
    toml_dict,
    FT;
    max_wind_speed = 25.0,
    context,
    use_lowres_forcing = true,
);

# Leaf area index comes from MODIS, regridded onto the sphere.

LAI = ClimaLand.Canopy.prescribed_lai_modis(
    domain.space.surface,
    start_date,
    stop_date,
);

# # The model
#
# This is the same call as in the site-level tutorial; only the domain differs. Soil
# texture, CLM canopy and photosynthesis parameters, runoff parameters and albedo are all
# read from global maps and regridded onto the domain, so a sphere needs no more
# configuration than a column.

model = ClimaLand.LandModel{FT}(
    forcing,
    LAI,
    toml_dict,
    domain,
    Δt;
    prognostic_land_components = (:canopy, :snow, :soil, :soilco2),
);

# # Diagnostics
#
# Here the column/global distinction matters. `default_diagnostics` dispatches on the
# domain: for `Column` and `Point` it returns an in-memory `DictWriter`, but for a
# `SphericalShell` it returns a `NetCDFWriter` that interpolates onto a regular
# longitude/latitude grid and writes to `outdir`. So unlike the site case, `outdir` is
# meaningful and output lands on disk.

root_path = "global_land_evaluation"
outdir = joinpath(root_path, "diagnostics")
mkpath(outdir)

# Choose the reduction period shorter than the run. Reductions are written at the *end* of
# each window, so asking for `:monthly` over a one-month run schedules its only output at
# `stop_date` and nothing is ever written -- leaving an empty directory and a confusing
# "variable not found" error downstream.
output_vars = ["swc", "tsoil", "lhf", "shf", "et", "gpp", "swe", "swu", "lwu"]
diagnostics = ClimaLand.default_diagnostics(
    model,
    start_date,
    outdir;
    reduction_period = :daily,
    output_vars,
);

# # Running
#
# The default initial conditions are read from a spun-up global file, which is what the
# 15 m / `nelements[2]`-layer soil column is designed to match. For a global run these
# defaults are exactly the intended path -- they were generated for this purpose.

simulation = LandSimulation(
    start_date,
    stop_date,
    Δt,
    model;
    outdir,
    diagnostics,
);
solve!(simulation);

# # Spatial diagnostics
#
# For global runs the natural visualisation is a map rather than a timeseries.
# `make_heatmaps` plots each diagnostic at a chosen date.

LandSimVis.make_heatmaps(
    simulation;
    date = stop_date,
    savedir = root_path,
);

# What to look for as a sanity check, independent of any observational dataset:
#
# - soil moisture should show wet tropics, dry subtropical deserts and a wet high-latitude
#   band, not a uniform field;
# - snow water equivalent should be confined to high latitudes and high elevation, and in
#   January should be northern-hemisphere heavy;
# - latent heat flux should be largest over the warm, wet tropics and near zero over
#   deserts and ice;
# - GPP should be near zero over deserts, ice and the winter high latitudes.
#
# If any field is uniform, or if land and ocean are not clearly distinguished, suspect the
# land/sea mask or the forcing rather than the physics.

# # Comparing against observations
#
# For global output ClimaLand provides a *leaderboard*: it compares the simulation against
# bundled observational datasets, computes bias and RMSE globally, and produces maps and
# summary plots.
#
# ```julia
# LandSimVis.make_leaderboard_plots(simulation; savedir = root_path)
# ```
#
# This needs a simulation at least a year long -- it compares seasonal climatology, so a
# one-month run gives it almost nothing to work with. To produce a real leaderboard, set
#
# ```julia
# stop_date = DateTime(2009, 1, 1)
# ```
#
# and expect roughly an hour at this resolution.
#
# The leaderboard draws on the same observational products the site-level tutorial uses
# pointwise -- MODIS ET, FLUXCOM GPP, CERES radiation -- which means the interpretive
# cautions carry over unchanged, and are worth restating because a global summary statistic
# hides them even more effectively than a point comparison:
#
# - **Upward longwave will look excellent, and that is close to meaningless.** In an
#   offline run it is dominated by σT⁴ with a surface temperature slaved to the prescribed
#   air temperature. At the site we showed a null model using only prescribed ERA5 `t2m`
#   and no land model at all reproduced r = +0.986 against CERES, versus +0.997 for
#   ClimaLand. Globally the same mechanism applies, and averaging over a large domain makes
#   a high correlation even easier to achieve.
# - **Do not compare prescribed fields.** LAI is forcing here, not a prediction.
# - **Prefer variables the land model actually determines**: ET, GPP, the latent/sensible
#   partitioning, soil moisture, snow.
# - **Observational products have their own failure modes**, which spatial averaging can
#   either hide or amplify. At the Kanpur gridcell, MODIS ET reported near-zero values
#   through the pre-monsoon months, almost certainly retrieval failure rather than reality.
#
# A global leaderboard is good at finding where a model is wrong. It is poor at telling you
# whether a variable was ever a fair test. Pairing a global run with one or two carefully
# examined single-column runs -- as in the
# site-level tutorial -- is a far more reliable way to tell those apart.

# # Summary
#
# The integrated `LandModel` constructor is domain-agnostic: the same call that builds a
# single column builds a global simulation, with every spatially-varying parameter
# regridded to whichever domain you hand it. What changes between the two is mostly
# infrastructure -- a land/sea mask appears, diagnostics switch from an in-memory
# `DictWriter` to a `NetCDFWriter`, and visualisation moves from timeseries to maps -- and
# the appropriate forcing resolution, which is the one choice that is far more forgiving
# globally than it is at a point.
