# # Running ClimaLand at an arbitrary site with ERA5 forcing

# The [Fluxnet tutorials](@ref "Fluxnet simulations with an integrated soil and canopy model")
# run ClimaLand at flux-tower sites, where tower observations supply the forcing. But
# ClimaLand can be run as a single column *anywhere* on the globe, driven by reanalysis
# rather than tower data. This tutorial does that for a site with no tower at all:
# IIT Kanpur (26.51°N, 80.23°E), in the Indo-Gangetic plain of Uttar Pradesh, India --
# an irrigated, double-cropped agricultural landscape with a strong monsoon cycle.
#
# The mechanics of a point run are easy: `Column` takes a `longlat` keyword and every
# spatially-varying default parameter is looked up by regridding global maps to that point.
# The hard parts are getting trustworthy forcing and knowing whether the answer means
# anything. This tutorial spends most of its time there, because that is where the traps
# are:
#
# 1. the low-resolution ERA5 artifact can silently give you a climate from hundreds of
#    kilometres away;
# 2. the CDS download has three separate gotchas that each produce a plausible-looking but
#    wrong file;
# 3. periodic diagnostics are timestamped at the *end* of their averaging window, which
#    quietly shifts any model-observation comparison by one period;
# 4. some "validation" variables are largely inherited from the forcing, so agreement with
#    observations does not mean the land model is right.
#
# !!! note "Data prerequisite"
#     This tutorial needs high-resolution ERA5 covering the site, which you must supply
#     yourself (see "Getting the forcing" below). That is why it is not part of the
#     automated documentation build. The
#     [global tutorial](@ref "Evaluating the integrated land model globally")
#     runs anywhere with no extra data.

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
using NCDatasets, Printf, Statistics
using CairoMakie, ClimaAnalysis, GeoMakie
import ClimaLand.LandSimVis as LandSimVis;

const FT = Float64;
context = ClimaComms.context();
toml_dict = LP.create_toml_dict(FT);

# # The domain
#
# A column at a location is built by passing `longlat`. Note the order: `(longitude,
# latitude)`, with longitude in [-180, 180].
#
# The 15 m depth and `nelements = 15` are not arbitrary. `LandSimulation` defaults to
# initial conditions read from a spun-up global file whose soil column is 15 m deep, so
# matching that depth lets you use those initial conditions directly.

longlat = FT.((80.23, 26.51))
zlim = FT.((-15, 0))
nelements = 15
dz_tuple = FT.((3, 0.05))
domain = Column(; zlim, longlat, nelements, dz_tuple);
surface_space = domain.space.surface;

# The resulting space genuinely carries the site coordinates, which is what allows every
# gridded parameter lookup (soil texture, CLM canopy properties, MODIS LAI) to work:

@show ClimaLand.Domains.get_long(surface_space)
@show ClimaLand.Domains.get_lat(surface_space);

# # Getting the forcing
#
# ClimaLand ships two ERA5 options, and for a specific site the choice matters enormously.
#
# **The low-resolution artifact is 8° × 8°.** It is downloadable and convenient, and
# `prescribed_forcing_era5` regrids with nearest-neighbour interpolation. At 8° the nearest
# cell centre to (80.23°E, 26.51°N) is at **(80°E, 30°N)** -- on the Tibetan Plateau. Read
# straight out of that artifact, that cell has an annual mean 2 m temperature of
# **263.4 K (−9.7 °C)**, a surface pressure of **52.4 kPa** (about 5000 m elevation), and
# **291 of its 450 mm/yr of precipitation falling as snow**. The neighbouring cell at 22°N
# gives 298.8 K, 99.3 kPa and 1550 mm/yr with no snow at all.
#
# So a run at the true coordinates of Kanpur, forced with low-resolution ERA5, simulates a
# frozen alpine column. Nothing errors. The plots look plausible. This is the single most
# important thing to check when moving a column to a new location.
#
# **The high-resolution artifact is 1° × 1° for 1979-2024**, but it has no download URL in
# `Artifacts.toml` -- it exists only on the Caltech cluster. To run this site elsewhere you
# supply the data yourself, and the cleanest way is a standard Julia artifact override, so
# that `prescribed_forcing_era5` picks it up with no source changes.
#
# ### Downloading from the Copernicus CDS
#
# Create a free account, accept the ERA5 licence on the dataset page (requests fail with an
# opaque 403 otherwise), and request `reanalysis-era5-single-levels` at native 0.25°.
# ClimaLand reads exactly ten variables; request the *mean-rate* forms, which is why no
# de-accumulation is needed:
#
# | CDS variable | ClimaLand name |
# |:---|:---|
# | `mean_total_precipitation_rate` | `mtpr` |
# | `mean_snowfall_rate` | `msr` |
# | `10m_u_component_of_wind` | `u10` |
# | `10m_v_component_of_wind` | `v10` |
# | `2m_dewpoint_temperature` | `d2m` |
# | `2m_temperature` | `t2m` |
# | `surface_pressure` | `sp` |
# | `mean_surface_downward_short_wave_radiation_flux` | `msdwswrf` |
# | `mean_surface_direct_short_wave_radiation_flux` | `msdrswrf` |
# | `mean_surface_downward_long_wave_radiation_flux` | `msdwlwrf` |
#
# Because regridding is nearest-neighbour, only the containing cell is ever read, so the
# box only needs to comfortably contain the site. `area = [27.5, 79.0, 25.5, 81.5]`
# (N/W/S/E) is 11 × 9 cells and about 35 MB per year. At 0.25° the nearest cell centre to
# the site is (80.25°E, 26.50°N), within roughly 2 km.
#
# Three gotchas, each of which produces a file that looks fine but is not:
#
# 1. **Cost limits.** A whole year of all ten variables is rejected outright. Empirically
#    5 variables × 6 months is still too large; 5 variables × 4 months is accepted. Split
#    the request into 4-month chunks.
# 2. **Stream splitting.** Instantaneous variables (`u10`, `t2m`, `sp`, ...) and time-mean
#    variables (`mtpr`, `msr`, ...) live on different ERA5 streams. Requesting both at once
#    makes CDS return a *ZIP of two NetCDFs* regardless of `download_format`. Request the
#    two groups separately.
# 3. **Names and axis order.** The time-mean variables come back under GRIB names, not the
#    names ClimaLand reads, and latitude is **descending**:
#
# | CDS returns | rename to |
# |:---|:---|
# | `avg_tprate` | `mtpr` |
# | `avg_tsrwe` | `msr` |
# | `avg_sdswrf` | `msdwswrf` |
# | `avg_sdirswrf` | `msdrswrf` |
# | `avg_sdlwrf` | `msdwlwrf` |
#
# Descending latitude is the one that bites hardest: `InterpolationsRegridder` defaults to
# `dim_increasing = (true, true)` and `prescribed_forcing_era5` hardcodes its
# `regridder_kwargs`, so a descending axis errors inside Interpolations.jl. Reverse
# latitude (and the data along it) during preprocessing. Also drop any singleton `number`
# or `expver` coordinates. Dimension *names* need not be `lon`/`lat` -- only the order
# matters, lon first -- and `valid_time` is accepted as-is.
#
# ### Wiring it in
#
# Write one file per year named exactly `era5_<year>_1.0x1.0.nc` (the resolution in the
# name is only a label) into a directory, then add its content hash to
# `~/.julia/artifacts/Overrides.toml`:
#
# ```toml
# f269a0b057b9f438b4caafdef17da73746310787 = "/path/to/your/era5_india"
# ```
#
# That hash is the `git-tree-sha1` of `forty_yrs_era5_land_forcing_data` in ClimaLand's
# `Artifacts.toml`. Restart Julia afterwards.
#
# !!! warning "The override is depot-wide"
#     `Overrides.toml` rebinds the artifact for your entire Julia depot, not just this
#     project. With a small regional box, any *other* ClimaLand run that asks for
#     high-resolution ERA5 will silently receive your regional data extrapolated outward
#     (longitude is `Periodic`, latitude is `Flat`) rather than erroring. Runs for years you
#     did not download do fail loudly. Comment the line out when working on global runs, or
#     use a dedicated `JULIA_DEPOT_PATH`.

# Simulation period. Ending inside 2010 keeps the forcing to a single year file, because
# `find_era5_year_paths` requires a file for every year in `year(start):year(stop)`.
start_date = DateTime(2010, 1, 1)
stop_date = DateTime(2010, 12, 31)
Δt = 450.0;

# `max_wind_speed` clips a known ERA5 artefact: occasional spurious 10 m wind spikes that
# would otherwise generate enormous surface fluxes and destabilise the run.
atmos, radiation = ClimaLand.prescribed_forcing_era5(
    start_date,
    stop_date,
    surface_space,
    toml_dict,
    FT;
    max_wind_speed = 25.0,
    context,
);
forcing = (; atmos, radiation);

# # Check the forcing before running anything
#
# This is the guard against the plateau problem, and it costs seconds. If you see roughly
# 263 K and 52 kPa, your override is not in effect and you are reading the 8° cell.

era5_dir = ClimaLand.Artifacts.era5_land_forcing_data_forty_years_folder_path()
@show era5_dir
@show readdir(era5_dir);

# Expected for the Indo-Gangetic plain, and what this data actually gives: annual mean air
# temperature **299.3 K (26.2 °C)**, surface pressure **99.2 kPa**, total precipitation
# **1018 mm/yr** with **90% of it falling in July-September**, and **zero snowfall**.

# # Leaf area index
#
# LAI is prescribed from MODIS, which covers 2000-2020 at 1°. At this site it resolves the
# real double-cropping cycle: a *rabi* (winter wheat) peak around February and a *kharif*
# (monsoon rice) peak around September.

LAI = ClimaLand.Canopy.prescribed_lai_modis(
    surface_space,
    start_date,
    stop_date,
);

# # Building the model
#
# The `LandModel` convenience constructor assembles soil, canopy, snow and soil CO₂ with
# defaults appropriate to the domain, and every spatially-varying parameter is regridded to
# the site automatically.
#
# !!! warning "Do not build the canopy separately without matching soil parameters"
#     It is tempting to construct the canopy yourself to choose a moisture-stress model.
#     But `PiecewiseMoistureStressModel{FT}(domain, toml_dict)` defaults its thresholds to
#     the *Gupta* van Genuchten parameters, while `EnergyHydrology` defaults its retention
#     curve to *ROSETTA*. Those disagree (at this site by 0.017 in porosity and 0.034 in
#     residual water content), and `LandModel`'s constructor asserts that they match, so
#     construction fails. Taking the defaults wires the canopy's thresholds to the soil
#     model's own parameters. If you do need to override one, override both consistently.

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
# For `Column` and `Point` domains `default_diagnostics` returns an in-memory `DictWriter`
# rather than a `NetCDFWriter`, which is the recommended choice for column runs. A
# consequence worth knowing: the `outdir` argument is then unused and nothing lands on disk
# except the plots you make yourself.

output_vars = [
    "swc",
    "tsoil",
    "lhf",
    "shf",
    "et",
    "gpp",
    "lai",
    "precip",
    "tair",
    "swe",
    "swu",
    "lwu",
]
diagnostics = ClimaLand.default_diagnostics(
    model,
    start_date,
    "."; # unused for a Column: the writer is a DictWriter
    reduction_period = :daily,
    output_vars,
);

# # Running

simulation =
    LandSimulation(start_date, stop_date, Δt, model; diagnostics);
solve!(simulation);

# A full year at Δt = 450 s takes roughly a minute on a laptop.

# # Does the answer look like the Gangetic plain?
#
# Physically expected behaviour at this site, and what the run produces:
#
# - soil moisture rises sharply at monsoon onset in late June, from 0.123 in June to 0.306
#   in August, plateauing near 0.39 through mid-September;
# - latent heat dominates during the monsoon and the partitioning reverses in the dry
#   season: LHF/SHF is 2.51 over June-September but 0.83 over March-May;
# - snow water equivalent is exactly zero all year;
# - LAI is bimodal, peaking at 1.46 in February (rabi) and 1.33 in September (kharif).

LandSimVis.make_timeseries(simulation; savedir = ".");

# # Comparing against observations
#
# ClimaLand bundles several gridded observational products as artifacts, all of which cover
# this site and period. Three are genuinely independent of our forcing:
#
# - **MODIS ET** (satellite retrieval, 0.5°)
# - **FLUXCOM GPP** (machine-learning upscaling of FLUXNET towers, 0.5°)
# - **CERES EBAF** upward shortwave, upward longwave and albedo (satellite, 1°)
#
# ClimaLand's monthly ERA5 surface dataset also provides `mslhf`/`msshf`, but those come
# from ERA5's own land model driven by the same atmosphere we force with, so comparing to
# them is a model-vs-model benchmark, not validation.
#
# !!! warning "Periodic diagnostics are timestamped at the end of their window"
#     A monthly average stamped 1 February is *January's* mean. Grouping by the raw
#     timestamp shifts every simulated month by one and silently corrupts the comparison.
#     Step back into the window before extracting the month:
#     `month.(date.(t) .- Day(1))`. On this run, fixing that changed the GPP correlation
#     from −0.58 to +0.02 and the upward-longwave correlation from +0.78 to +1.00 -- the
#     difference between reporting an anticorrelation and reporting no correlation.
#
# Results for 2010 at this site, monthly:
#
# | variable | reference | sim | obs | bias | RMSE | r |
# |:---|:---|---:|---:|---:|---:|---:|
# | upward longwave | CERES | 467.7 | 460.9 | +6.8 | 7.9 | +1.00 |
# | upward shortwave | CERES | 43.7 | 31.6 | +12.1 | 15.1 | +0.83 |
# | ET | MODIS | 2.52e-5 | 1.62e-5 | +9.0e-6 | 1.26e-5 | +0.72 |
# | GPP | FLUXCOM | 2.23 | 2.99 | −0.76 | 1.75 | +0.02 |
#
# # Interpreting this honestly
#
# ### Not all agreement is meaningful
#
# Upward longwave correlates with CERES at r = +1.00, which looks like a triumph. It is
# not. ClimaLand computes
#
# ```
# LW_u = (1 - ϵ_canopy) * LW_u_ground + ϵ_canopy * σ * T_canopy^4
# ```
#
# so `LW_u` is dominated by σT⁴ using model temperatures -- but in an offline run the
# surface temperature is closely slaved to the *prescribed* air temperature. Comparing a
# null model that uses nothing but prescribed ERA5 `t2m`, `ϵσT_air⁴`, with no land model at
# all:
#
# | vs CERES upward longwave | r | RMSE | bias |
# |:---|---:|---:|---:|
# | ClimaLand | +0.9965 | 7.87 | +6.77 |
# | null: `ϵσT_air⁴` only | +0.9856 | 13.58 | −8.68 |
#
# The null model reproduces the correlation almost perfectly. Note the mechanism is *not*
# the reflected `(1-ϵ)·LW_d` term inherited from the forcing -- that is only 7.8 W/m², 1.7%
# of the total. It is that `corr(LW_u, prescribed T_air) = +0.985`, while
# `corr(CERES, prescribed T_air) = +0.981`. The land model's real contribution is halving
# the RMSE by getting the surface-to-air temperature offset (+2.5 K) roughly right. That is
# worth something, but upward longwave should not be counted as validation.
#
# The same caution applies to prescribed fields: comparing the model's `lai` diagnostic to
# MODIS LAI is circular, since MODIS LAI *is* the forcing. It is a useful check that
# ingestion worked (r = +0.88 here, limited by monthly averaging), not a test of skill.
#
# ### Two real, diagnosable biases
#
# **Albedo is too bright.** Against CERES the model's surface albedo is 0.1905 versus
# 0.1527, a bias of +0.038, and r is only +0.33, so the seasonal cycle is wrong too. The
# discrepancy peaks in April-May (0.22 versus 0.16), exactly when LAI bottoms out at 0.35
# and bare soil dominates. That accounts for 8.5 of the 12.1 W/m² upward-shortwave bias and
# points at the CLM soil albedo for this gridcell being too bright for Gangetic alluvium.
#
# **The Bowen ratio is too high.** LHF is 21 W/m² low and SHF is 19 W/m² high -- nearly
# compensating, so net radiation stays reasonable, consistent with the good longwave. Both
# errors concentrate in the pre-monsoon dry season when modelled soil moisture bottoms at
# 0.125 and the surface cannot evaporate.
#
# ### GPP phasing is the outstanding problem
#
# Annual mean GPP is within 26% of FLUXCOM, but the seasonal cycle carries no information
# (r = +0.02). The model peaks in late February at 4.0e-6 mol CO₂ m⁻² s⁻¹ and reaches its
# *minimum* in early July, precisely when FLUXCOM climbs to its August peak of 5.5e-6. The
# model does produce the bimodal rabi/kharif structure that MODIS LAI implies, but with the
# amplitudes inverted: it makes rabi dominant where FLUXCOM makes kharif dominant by more
# than a factor of two.
#
# Since the leaf area is prescribed and correct, and soil moisture peaks at 0.39 in August,
# the limitation is elsewhere -- plausibly radiative, with monsoon cloud suppressing PAR, or
# in the CLM photosynthesis parameters assigned to this cell, which may not represent
# irrigated cropland.
#
# ### A caveat that runs through everything
#
# ERA5 has no irrigation, which is the dominant dry-season water flux in this landscape.
# Prescribed MODIS LAI, however, *does* see the irrigated rabi crop. So the canopy is
# productive on soil the model believes is dry. That inconsistency is inherent to
# prescribed-LAI runs in irrigated regions, and it means dry-season fluxes here should not
# be trusted quantitatively.
#
# Also note MODIS ET reports 3.9e-8, 3.7e-7 and 1.0e-6 kg m⁻² s⁻¹ for April, May and June:
# implausibly near zero for irrigated pre-monsoon cropland, and most likely retrieval
# failure. Part of the +56% ET bias is therefore an artefact of the reference, not model
# error. Observational products have their own failure modes.

# # Summary
#
# Running ClimaLand at an arbitrary point is a two-line change to the domain. Making the
# result trustworthy is the work:
#
# 1. verify the forcing at the site before running -- check that temperature and pressure
#    are plausible for the place you think you are simulating;
# 2. take component defaults unless you are prepared to keep cross-component parameters
#    consistent yourself;
# 3. correct for end-of-window diagnostic timestamps before comparing to observations;
# 4. ask, for every variable that agrees well, whether it could have agreed without the
#    land model doing anything -- and prefer variables that could not.
