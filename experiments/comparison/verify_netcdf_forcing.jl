# Verify that site-level FLUXNET-style NetCDF forcing read through the new ClimaUtilities
# multi-column machinery (`DataSource` -> `MultiColumnDataHandler` -> `TimeVaryingInput`)
# reproduces a direct `NCDatasets` read of the same file, on a SINGLE-column
# `ColumnEnsemble`.
#
# The handler is the only supported way to read per-column NetCDF onto a ClimaLand domain:
# it requires a `PointCloudSpace` surface, so a plain `Column` cannot use it and a
# one-point `ColumnEnsemble` is the single-site case. Nothing here is specific to N = 1
# except the domain, so the same check scales to N sites.
#
# Ground truth is a hand-rolled linear interpolation of the raw file, so this isolates the
# read/regrid/time-shift path from any model: agreement means the handler returns the
# file's values, in the file's units after the documented conversions, at the right UTC
# instant. Two independent things are checked - the VALUE (per variable) and the TIME SHIFT
# (local standard time -> UTC), the latter being the part most likely to be silently wrong,
# since a whole-hour error still yields plausible-looking forcing.
#
# Run: julia --project=.buildkite experiments/comparison/verify_netcdf_forcing.jl [MET_NC_PATH] [HOUR_OFFSET_FROM_UTC]

import ClimaComms
ClimaComms.@import_required_backends
using ClimaCore
using Dates
using DelimitedFiles # loads the FluxnetSimulations extension
using NCDatasets
using Printf

using ClimaLand
using ClimaLand.Domains: ColumnEnsemble
import ClimaLand.FluxnetSimulations as FluxnetSimulations
import ClimaLand.Parameters as LP
import ClimaParams as CP

# `snow_precip_fraction` and `hour_offset_to_period` are internal to the extension, and the
# ground truth has to apply exactly the same functions the forcing path does.
const FSExt = Base.get_extension(ClimaLand, :FluxnetSimulationsExt)

const FT = Float64

# Any gap-filled FLUXNET/PLUMBER2-style `*_Met.nc` works. The default is the CH-Dav file in
# the ClimaUtilities worktree; the `callmip_phase1_forcing` artifact used by the multi-site
# driver has no download URL, so it cannot be resolved here.
const MET_NC_PATH =
    length(ARGS) >= 1 ? ARGS[1] :
    "/home/kphan2/worktree/ClimaUtilities.jl/multi-col-intp/CH-Dav_1997-2014_FLUXNET2015_Met.nc"
# Local standard time offset of the site (CH-Dav is CET, +1).
const HOUR_OFFSET_FROM_UTC = length(ARGS) >= 2 ? parse(FT, ARGS[2]) : FT(1)

isfile(MET_NC_PATH) || error("No met file at $MET_NC_PATH")
@info "Verifying NetCDF forcing" MET_NC_PATH HOUR_OFFSET_FROM_UTC

toml_dict = LP.create_toml_dict(FT)
earth_param_set = LP.LandParameters(toml_dict)
thermo_params = LP.thermodynamic_parameters(earth_param_set)

# ---------------------------------------------------------------------------
# Raw read: the ground truth. Dates are shifted local -> UTC exactly as the
# `DataSource` `time_transform` in `prescribed_forcing_netcdf` does.
# ---------------------------------------------------------------------------
raw = NCDataset(MET_NC_PATH, "r") do ds
    local_dates = ds["time"][:]
    utc_dates =
        local_dates .- FSExt.hour_offset_to_period(HOUR_OFFSET_FROM_UTC)
    read_var(name) = FT.(coalesce.(ds[name][1, 1, :], NaN))
    (;
        utc_dates,
        Tair = read_var("Tair"),
        Wind = read_var("Wind"),
        Qair = read_var("Qair"),
        Psurf = read_var("Psurf"),
        SWdown = read_var("SWdown"),
        LWdown = read_var("LWdown"),
        CO2air = read_var("CO2air"),
        Precip = read_var("Precip"),
        VPD = read_var("VPD"),
    )
end
@info "File record" n_times = length(raw.utc_dates) first_utc =
    first(raw.utc_dates) last_utc = last(raw.utc_dates) n_nan_Tair =
    count(isnan, raw.Tair)

# Linear interpolation of a raw series at a UTC instant, matching
# `LinearInterpolation()` in the TimeVaryingInput.
function raw_at(values, date)
    dates = raw.utc_dates
    date <= first(dates) && return values[begin]
    date >= last(dates) && return values[end]
    i = searchsortedlast(dates, date)
    dates[i] == date && return values[i]
    w =
        Dates.value(Dates.Millisecond(date - dates[i])) /
        Dates.value(Dates.Millisecond(dates[i + 1] - dates[i]))
    return values[i] + w * (values[i + 1] - values[i])
end

# ---------------------------------------------------------------------------
# Single-column domain at the file's own coordinates. `MultiColumnDataHandler`
# matches files to columns within 1e-3 degrees, so reading the location from the
# file is what guarantees the match.
# ---------------------------------------------------------------------------
(; long, lat) = FluxnetSimulations.get_location_netcdf(FT, MET_NC_PATH)
@info "Site location from file" long lat

land_domain = ColumnEnsemble(;
    zlim = (FT(-2), FT(0)),
    nelements = 10,
    longlat = (long, lat),
)
surface_space = land_domain.space.surface
@info "Surface space" space_type = typeof(surface_space) is_point_cloud =
    surface_space isa ClimaCore.Spaces.PointCloudSpace

# t = 0 sits a few days into the record so the interpolation stencil and the date
# bookkeeping are exercised away from the endpoints.
start_date = first(raw.utc_dates) + Day(3)
atmos_h = FT(30)

(; atmos, radiation) = FluxnetSimulations.prescribed_forcing_netcdf(
    MET_NC_PATH,
    surface_space,
    atmos_h,
    start_date,
    toml_dict,
    FT;
    hour_offset_from_UTC = HOUR_OFFSET_FROM_UTC,
)

# ---------------------------------------------------------------------------
# Compare handler output against the raw read at a spread of times, including
# half-step offsets so time interpolation is hit and not just snapshot lookup.
# ---------------------------------------------------------------------------
dt_file_seconds = Dates.value(Dates.Second(raw.utc_dates[2] - raw.utc_dates[1]))
query_seconds =
    FT[0, 0.5 * dt_file_seconds, 3600, 86400, 86400 + 1800, 30 * 86400]

scratch = ClimaCore.Fields.zeros(surface_space)
function handler_at(tvi, t)
    ClimaLand.evaluate!(scratch, tvi, t)
    return parent(scratch)[1]
end

# The snow split is a composed variable, so the ground truth has to recompute it from the
# raw Tair/VPD at the same instant.
function snow_frac(date)
    return FSExt.snow_precip_fraction(
        raw_at(raw.Tair, date) - 273.15,
        raw_at(raw.VPD, date);
        thermo_params,
    )
end

# Each entry maps a driver to the raw quantity it should equal, applying the unit
# conversions documented on `prescribed_forcing_netcdf`.
checks = (
    ("T (K)", atmos.T, d -> raw_at(raw.Tair, d)),
    ("u (m/s)", atmos.u, d -> raw_at(raw.Wind, d)),
    ("q (kg/kg)", atmos.q, d -> raw_at(raw.Qair, d)),
    ("P (Pa)", atmos.P, d -> raw_at(raw.Psurf, d)),
    ("SW_d (W/m^2)", radiation.SW_d, d -> raw_at(raw.SWdown, d)),
    ("LW_d (W/m^2)", radiation.LW_d, d -> raw_at(raw.LWdown, d)),
    ("c_co2 (mol/mol)", atmos.c_co2, d -> raw_at(raw.CO2air, d) * 1e-6),
    (
        "liquid_precip (m/s)",
        atmos.liquid_precip,
        d -> -(raw_at(raw.Precip, d) / 1000) * (1 - snow_frac(d)),
    ),
    (
        "snow_precip (m/s)",
        atmos.snow_precip,
        d -> -(raw_at(raw.Precip, d) / 1000) * snow_frac(d),
    ),
)

@info "Handler vs raw read" start_date file_dt_seconds = dt_file_seconds
for (name, tvi, truth) in checks
    for t in query_seconds
        date = start_date + Millisecond(round(Int, 1000 * t))
        got = handler_at(tvi, t)
        want = truth(date)
        scale = max(abs(want), eps(FT))
        @printf(
            "%-22s t=%9.1f s  handler=% .8e  raw=% .8e  abs=%.3e  rel=%.3e\n",
            name,
            t,
            got,
            want,
            abs(got - want),
            abs(got - want) / scale
        )
    end
end

# ---------------------------------------------------------------------------
# Time-shift check on its own. Comparing the handler at t = 0 against the raw
# series interpolated at UTC and at UTC +/- 1h shows where the local -> UTC
# conversion actually landed: the middle row is the one that should match, and
# the neighbours should differ (unless the field is flat there).
# ---------------------------------------------------------------------------
@info "Time-shift discrimination (T at t = 0)"
got_T = handler_at(atmos.T, FT(0))
@printf("  handler: % .8e\n", got_T)
for shift in (-Hour(1), Hour(0), Hour(1))
    want = raw_at(raw.Tair, start_date + shift)
    @printf(
        "  raw at start_date %+dh: % .8e   |handler - raw| = %.3e\n",
        Dates.value(Dates.Hour(shift)),
        want,
        abs(got_T - want)
    )
end
