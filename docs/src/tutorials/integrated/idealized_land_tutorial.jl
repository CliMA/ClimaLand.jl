# # An idealized year-long land model simulation

# The other integrated tutorials drive ClimaLand with real data: tower measurements in the
# [Fluxnet tutorials](@ref "Fluxnet simulations with an integrated soil and canopy model"),
# or reanalysis in the site-level ERA5 tutorial. Both tie you to a download. This tutorial
# takes the opposite approach and drives the full soil-canopy-snow-soilCO₂ model with
# forcing written as **plain functions of time**.
#
# That buys three things:
#
# 1. **No forcing data at all.** Nothing to fetch, nothing to regrid, nothing that can be
#    silently wrong. (Soil and canopy *parameters* are still read from global maps, but
#    those artifacts are small and download automatically.)
# 2. **A full year runs on a laptop**, so you see a complete seasonal cycle rather than a
#    few days.
# 3. **You control the experiment.** Because every driver is a function you wrote, you can
#    halve the rainfall or flatten the seasonal cycle and see exactly what the land model
#    does about it. That is much harder to interpret when the forcing is reanalysis.
#
# The cost is that this is not a real place. Treat it as a controlled numerical experiment,
# not as a simulation of anywhere in particular.
#
# ClimaLand supports this directly through `prescribed_analytic_forcing`, which is used
# throughout the test suite.

# # Preliminary setup

using Dates
import ClimaComms
ClimaComms.@import_required_backends
import ClimaParams as CP
using ClimaCore
using ClimaLand
using ClimaLand.Domains: Column
import ClimaLand.Parameters as LP
import ClimaLand.Simulations: LandSimulation, solve!
import ClimaUtilities.TimeManager: date
import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput
using Printf, Statistics
using CairoMakie, ClimaAnalysis, GeoMakie
import ClimaLand.LandSimVis as LandSimVis;

const FT = Float64;
context = ClimaComms.context();
toml_dict = LP.create_toml_dict(FT);

# # Writing the forcing
#
# We build a mid-latitude northern-hemisphere climate: a warm wet summer and a cool dry
# winter, with a diurnal cycle on top.
#
# !!! warning "Driver functions receive an `ITime`, not a `Float64`"
#     ClimaLand evaluates drivers with an `ITime`, so arithmetic like `mod(t, 86400)`
#     fails. Call `float(t)` first, which returns the time in seconds and also works if a
#     plain number is passed. Most existing uses of `prescribed_analytic_forcing` are
#     *constant* functions that ignore `t` entirely, so they never hit this.

const DAY = 86400.0
const SIG = 5.670374419e-8

secs(t) = float(t)
dayfrac(t) = mod(secs(t), DAY) / DAY
## +1 at midsummer, -1 at midwinter, zero near the equinoxes
seasonal(t) = sin(2π * (secs(t) / DAY - 80) / 365.25)
## a half-sine bump centred on local noon, zero at night
noon(t) = max(0.0, cos(2π * (dayfrac(t) - 0.5)));

# Downward shortwave: a diurnal bump whose amplitude is modulated seasonally.
SW_d(t) = 900.0 * (0.7 + 0.3 * seasonal(t)) * noon(t)

# Air temperature: annual mean 288 K, ±12 K seasonally, ±4 K diurnally, warmest
# mid-afternoon.
T_atmos(t) = 288.0 + 12.0 * seasonal(t) + 4.0 * cos(2π * (dayfrac(t) - 0.6))

# Downward longwave from an effective atmospheric emissivity of 0.8.
LW_d(t) = 0.80 * SIG * T_atmos(t)^4

# Specific humidity, wind and pressure kept deliberately simple.
q_atmos(t) = 0.006 + 0.004 * max(0.0, seasonal(t))
u_atmos(t) = 2.0
P_atmos(t) = 101325.0

# Total precipitation as six-hour events every five days, confined to the warm half of the
# year -- roughly 750 mm/yr.
#
# !!! note "Precipitation sign convention"
#     Precipitation is a volume flux of liquid water in m/s and is **negative downward**,
#     matching what the ERA5 helper produces. A positive value here would remove water.
precip(t) =
    -1.5e-6 *
    max(0.0, seasonal(t)) *
    (mod(secs(t), 5 * DAY) < 6 * 3600 ? 1.0 : 0.0)

# The soil and snow models need liquid and solid precipitation *separately*, so we split
# the total by air temperature. Everything falls as rain above freezing and as snow below.
rain(t) = T_atmos(t) > 273.15 ? precip(t) : 0.0
snow(t) = T_atmos(t) > 273.15 ? 0.0 : precip(t);

# # The domain
#
# Even though the forcing is synthetic, the soil and canopy *parameters* are still looked
# up from global maps, so the column needs a `longlat`. Picking a mid-latitude land point
# gives sensible soil texture and vegetation properties to go with the climate we invented.
#
# !!! warning "Longitude must not be exactly zero"
#     A column at a location is built internally as a degenerate 1×1 box, and a longitude
#     of exactly `0.0` collapses one axis, tripping an `ylim[1] < ylim[2]` assertion. Any
#     non-zero longitude is fine.

longlat = FT.((10.0, 45.0))
domain = Column(;
    zlim = FT.((-15, 0)),
    nelements = 15,
    dz_tuple = FT.((3, 0.05)),
    longlat,
)
surface_space = domain.space.surface;

# # Assembling the drivers
#
# ClimaLand provides `prescribed_analytic_forcing` as a convenience wrapper, and it is
# what the test suite uses. For an integrated model with a canopy and snow, though, it is
# worth building the two driver objects directly, because the helper makes two
# simplifications that are invisible until they corrupt your results:
#
# !!! warning "The helper reuses one precipitation function for both rain and snow"
#     It passes the same `precip` argument as *both* `liquid_precip` and `snow_precip`
#     (`drivers.jl:1606-1607`). With its default `precip = (t) -> 0` that is harmless, but
#     supply a real precipitation function and every rain event also delivers an equal
#     snowfall, regardless of temperature -- doubling the water input. When I first ran
#     this tutorial that way, snow water equivalent spiked to 10 mm during every
#     warm-season storm at 294 K.
#
# !!! warning "The helper omits the solar zenith angle"
#     It constructs `PrescribedRadiativeFluxes` without a `cosθs`. When the zenith angle is
#     absent, ClimaLand sets **both** `cosθs` and the diffuse fraction to `NaN`
#     (`drivers.jl:1297-1300`), on the assumption they will not be used. A canopy *does*
#     use them, so the `NaN` propagates into every flux and the simulation silently returns
#     `NaN` for everything without ever crashing. Fine for the standalone soil and bucket
#     tests that use the helper; not fine here.
#
# Building both explicitly avoids each problem: separate rain and snow inputs, and a real
# zenith angle. Passing `toml_dict` to the radiation lets ClimaLand compute the diffuse
# fraction from its empirical relation.

start_date = DateTime(2008, 1, 1)
stop_date = start_date + Year(1)
Δt = 450.0

atmos = ClimaLand.PrescribedAtmosphere(
    TimeVaryingInput(rain),
    TimeVaryingInput(snow),
    TimeVaryingInput(T_atmos),
    TimeVaryingInput(u_atmos),
    TimeVaryingInput(q_atmos),
    TimeVaryingInput(P_atmos),
    start_date,
    FT(10),                      # measurement height, m
    toml_dict,
)

earth_param_set = LP.LandParameters(toml_dict)
cosθs =
    (t, s) -> ClimaLand.default_cos_zenith_angle(
        t,
        s;
        latitude = ClimaCore.Fields.coordinate_field(surface_space).lat,
        longitude = ClimaCore.Fields.coordinate_field(surface_space).long,
        insol_params = earth_param_set.insol_params,
    )
radiation = ClimaLand.PrescribedRadiativeFluxes(
    FT,
    TimeVaryingInput(SW_d),
    TimeVaryingInput(LW_d),
    start_date;
    cosθs,
    toml_dict,
)
forcing = (; atmos, radiation);

# # Leaf area index
#
# The site-level tutorial reads LAI from MODIS. Here it is another function of time: a bare
# winter canopy that greens up through the warm season.

LAI = TimeVaryingInput(t -> FT(0.5 + 2.0 * max(0.0, seasonal(t))));

# # The model
#
# Identical to the data-driven tutorials -- only the drivers differ.

model = ClimaLand.LandModel{FT}(
    forcing,
    LAI,
    toml_dict,
    domain,
    Δt;
    prognostic_land_components = (:canopy, :snow, :soil, :soilco2),
);

# # Diagnostics and running
#
# On a `Column` the default writer is an in-memory `DictWriter`, so nothing is written to
# disk except the plots.

output_vars =
    ["swc", "tsoil", "lhf", "shf", "et", "gpp", "lai", "swe", "precip", "tair"]
diagnostics = ClimaLand.default_diagnostics(
    model,
    start_date,
    ".";
    reduction_period = :daily,
    output_vars,
)
simulation = LandSimulation(start_date, stop_date, Δt, model; diagnostics);

# The timestepping itself takes one to two minutes for the full year on a laptop, varying
# with ClimaLand version; expect a few minutes in total once compilation and plotting are
# included. Shorten `stop_date` to `start_date + Month(1)` while experimenting.

solve!(simulation);

# # What to look for
#
# Because we wrote the forcing, we know what the answer should qualitatively be, which
# makes this a much sharper test of understanding than comparing against reanalysis. The
# seasonal contrasts below are from an actual run of this tutorial. Expect small
# differences between ClimaLand versions; the ones below are reproducible to the precision
# quoted across v1.10 and v1.12.
#
# | variable | Jun-Aug | Dec-Feb |
# |:---|---:|---:|
# | soil temperature (K) | 294.1 | 279.2 |
# | soil moisture (m³/m³) | 0.36 | 0.31 |
# | latent heat flux (W/m²) | ~133 | 6.6 |
# | GPP (mol CO₂ m⁻² s⁻¹) | ~6.1e-6 | 1.9e-6 |
# | LAI (m²/m²) | 2.11 | 0.50 |
#
# - **Soil temperature** tracks the imposed seasonal air temperature, damped and lagged
#   with depth: the surface layer follows closely, the deep soil barely moves.
# - **Soil moisture** rises with the warm-season rain events and draws down between them.
# - **Latent heat flux** shows a 20-fold seasonal contrast, because radiation, temperature,
#   LAI and water availability were all constructed to peak together.
# - **Snow water equivalent stays negligible.** Because precipitation is partitioned by
#   temperature and rain only falls in the warm half of the year, no snow is ever
#   delivered. What remains is a small residue inherited from the spun-up initial
#   conditions: about 0.35 mm on 1 January, decaying to near zero by April, with a ripple
#   below 0.04 mm thereafter. That is three orders of magnitude smaller than what the
#   single-`precip` helper produces -- if you see storm-driven spikes of several mm in the
#   warm season, the rain/snow split is not being applied.

LandSimVis.make_timeseries(simulation; savedir = ".");

# # Experiments worth running
#
# This is the real value of an idealized setup -- each of these is a one-line change and
# the interpretation is unambiguous:
#
# 1. **Move the rain into winter.** Replace `max(0.0, seasonal(t))` in `precip` with
#    `max(0.0, -seasonal(t))`. Snow should now accumulate and `swe` become non-zero, with
#    a melt pulse into soil moisture in spring.
# 2. **Remove the seasonal cycle.** Set `seasonal(t) = 0.0`. The model should approach a
#    diurnally-periodic steady state; how long that takes tells you the system's memory,
#    which is dominated by the deep soil.
# 3. **Halve the rainfall.** Scale `precip` by `0.5` and watch how much of the reduction
#    shows up in evapotranspiration versus soil storage versus runoff.
# 4. **Vary the location.** Keep the forcing but change `longlat`. Only the soil and
#    vegetation parameters change, isolating the effect of soil texture and plant
#    functional type on the same climate.
#
# # Caveats
#
# The forcing here is crude on purpose: humidity is not consistent with the temperature in
# any thermodynamic sense, there is no weather variability beyond the imposed rain events,
# and the radiation ignores clouds entirely. So the absolute fluxes are not meant to be
# realistic, and this run should not be compared against observational products. What it is
# good for is understanding how the components respond to forcing you fully control, and as
# a reproducible starting point that needs no data and completes in one sitting.
