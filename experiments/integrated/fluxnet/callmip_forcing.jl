# Helpers to load CalLMIP Phase 1b met forcing (PLUMBER2 `*_Met.nc`) for one or
# more FLUXNET sites onto a `ColumnEnsemble`, using the in-memory multi-point
# `TimeVaryingInput(times, matrix, space)` from ClimaUtilities.
#
# Each row of a `(n_sites, n_times)` matrix drives one ensemble column, in the
# same order as the `longlat` vector used to build the domain. All columns
# share one UTC time axis, so the sites' records must overlap in calendar time.

import ClimaCore
using Dates
using NCDatasets
import DelimitedFiles # triggers FluxnetSimulationsExt (snow_precip_fraction)
using ClimaLand
using ClimaLand.Domains: Column, ColumnEnsemble
import ClimaLand.Parameters as LP
using ClimaLand.Simulations: LandSimulation, solve!
import ClimaUtilities.TimeVaryingInputs:
    TimeVaryingInput, LinearInterpolation, evaluate!

const FSExt = Base.get_extension(ClimaLand, :FluxnetSimulationsExt)

"""
    read_callmip_met(met_nc_path, hour_offset_from_UTC)

Read one CalLMIP/PLUMBER2 `*_Met.nc` file fully into memory and return a
NamedTuple with the UTC timestamps (`UTC = local - hour_offset_from_UTC`), the
site coordinates and tower height, and the forcing time series in the file's
native units. The data are assumed gap-filled; missing values are an error.
"""
function read_callmip_met(met_nc_path, hour_offset_from_UTC)
    NCDataset(met_nc_path) do ds
        function read_series(name)
            v = ds[name][1, 1, :]
            any(ismissing, v) &&
                error("$name in $met_nc_path contains missing values")
            return Float64.(v)
        end
        local_dates = ds["time"][:]
        utc_dates =
            local_dates .- FSExt.hour_offset_to_period(hour_offset_from_UTC)
        return (;
            utc_dates,
            long = Float64(only(ds["longitude"])),
            lat = Float64(only(ds["latitude"])),
            atmos_h = Float64(only(ds["reference_height"])),
            Tair = read_series("Tair"),     # K
            Wind = read_series("Wind"),     # m/s
            Qair = read_series("Qair"),     # kg/kg
            Psurf = read_series("Psurf"),   # Pa
            SWdown = read_series("SWdown"), # W/m^2
            LWdown = read_series("LWdown"), # W/m^2
            CO2air = read_series("CO2air"), # ppm
            Precip = read_series("Precip"), # kg/m^2/s
            VPD = read_series("VPD"),       # hPa
            LAI = haskey(ds, "LAI_alternative") ?
                  read_series("LAI_alternative") : nothing,
        )
    end
end

"""
    align_sites(sites)

Trim the sites returned by [`read_callmip_met`](@ref) to their shared UTC time
axis and return `(utc_dates, aligned_sites)`. All sites must be sampled at the
same interval and on the same time grid.
"""
function align_sites(sites)
    dt = sites[1].utc_dates[2] - sites[1].utc_dates[1]
    for s in sites
        @assert all(==(dt), diff(s.utc_dates)) "met data must be equispaced in time"
    end
    t0 = maximum(s -> first(s.utc_dates), sites)
    t1 = minimum(s -> last(s.utc_dates), sites)
    @assert t0 <= t1 "the sites' records do not overlap in time"
    utc_dates = collect(t0:dt:t1)
    n = length(utc_dates)
    slice(v::Vector, rng) = v[rng]
    slice(v, _) = v # scalars and `nothing` pass through
    aligned = map(sites) do s
        i0 = findfirst(==(t0), s.utc_dates)
        isnothing(i0) &&
            error("site time axes are offset by a non-multiple of $dt")
        return map(Base.Fix2(slice, i0:(i0 + n - 1)), s)
    end
    return utc_dates, aligned
end

"""
    prescribed_forcing_callmip(sites, surface_space, start_date, toml_dict, FT;
                               split_precip = true,
                               gustiness = 1,
                               time_interpolation_method = LinearInterpolation())

Construct the `PrescribedAtmosphere` and `PrescribedRadiativeFluxes` for the
`ColumnEnsemble` surface space from the sites' met data, with one column per
site. `sites` must be in the column (`longlat`) order of `surface_space`;
`start_date` is the UTC datetime corresponding to t = 0.

The tower height is read per site from the met files; for a single site it is
kept a scalar, otherwise it becomes a per-column surface Field.
"""
function prescribed_forcing_callmip(
    sites,
    surface_space,
    start_date,
    toml_dict,
    FT;
    split_precip = true,
    gustiness = 1,
    time_interpolation_method = LinearInterpolation(),
)
    utc_dates, aligned = align_sites(sites)
    seconds = Float64.(Dates.value.(utc_dates .- start_date)) ./ 1000

    # Rows of each matrix must follow the column order of the space
    space_long = Array(
        ClimaCore.Fields.field2array(
            ClimaLand.Domains.get_long(surface_space),
        ),
    )
    space_lat = Array(
        ClimaCore.Fields.field2array(ClimaLand.Domains.get_lat(surface_space)),
    )
    @assert length(space_long) == length(aligned) "one site per column is required"
    @assert all(
        isapprox.(space_long, [s.long for s in aligned]; atol = 1e-3),
    ) && all(isapprox.(space_lat, [s.lat for s in aligned]; atol = 1e-3)) "`sites` must be in the column (longlat) order of `surface_space`"

    matrix_tvi(rows) = TimeVaryingInput(
        seconds,
        stack(rows; dims = 1),
        surface_space;
        method = time_interpolation_method,
    )

    atmos_T = matrix_tvi([s.Tair for s in aligned])
    atmos_u = matrix_tvi([s.Wind for s in aligned])
    atmos_q = matrix_tvi([s.Qair for s in aligned])
    atmos_P = matrix_tvi([s.Psurf for s in aligned])
    SW_d = matrix_tvi([s.SWdown for s in aligned])
    LW_d = matrix_tvi([s.LWdown for s in aligned])
    c_co2 = matrix_tvi([s.CO2air .* 1e-6 for s in aligned]) # ppm -> mol/mol

    earth_param_set = LP.LandParameters(toml_dict)
    thermo_params = LP.thermodynamic_parameters(earth_param_set)
    # Precip is a mass flux (kg/m^2/s); convert to a negative (downward) volume
    # flux (m/s), optionally split into rain and snow
    if split_precip
        snow_frac = [
            FSExt.snow_precip_fraction.(s.Tair .- 273.15, s.VPD; thermo_params) for s in aligned
        ]
        atmos_P_liq = matrix_tvi([
            -(s.Precip ./ 1000) .* (1 .- f) for
            (s, f) in zip(aligned, snow_frac)
        ])
        atmos_P_snow = matrix_tvi([
            -(s.Precip ./ 1000) .* f for (s, f) in zip(aligned, snow_frac)
        ])
    else
        atmos_P_liq = matrix_tvi([-(s.Precip ./ 1000) for s in aligned])
        atmos_P_snow = matrix_tvi([zero(s.Precip) for s in aligned])
    end

    heights = FT[s.atmos_h for s in aligned]
    atmos_h =
        length(heights) == 1 ? heights[1] :
        ClimaCore.Fields.array2field(heights, surface_space)

    atmos = ClimaLand.PrescribedAtmosphere(
        atmos_P_liq,
        atmos_P_snow,
        atmos_T,
        atmos_u,
        atmos_q,
        atmos_P,
        start_date,
        atmos_h,
        toml_dict;
        gustiness = FT(gustiness),
        c_co2,
    )

    cos_zenith_angle =
        (t, s) -> ClimaLand.default_cos_zenith_angle(
            t,
            s;
            insol_params = earth_param_set.insol_params,
            longitude = ClimaLand.Domains.get_long(surface_space),
            latitude = ClimaLand.Domains.get_lat(surface_space),
        )
    radiation = ClimaLand.PrescribedRadiativeFluxes(
        FT,
        SW_d,
        LW_d,
        start_date;
        cosθs = cos_zenith_angle,
        toml_dict = toml_dict,
    )
    return (; atmos, radiation)
end

"""
    prescribed_lai_callmip(sites, surface_space, start_date;
                           time_interpolation_method = LinearInterpolation())

Construct a per-column `TimeVaryingInput` of LAI (m^2/m^2) from the sites'
`LAI_alternative` (MODIS) time series.
"""
function prescribed_lai_callmip(
    sites,
    surface_space,
    start_date;
    time_interpolation_method = LinearInterpolation(),
)
    utc_dates, aligned = align_sites(sites)
    any(s -> isnothing(s.LAI), aligned) &&
        error("some met files do not provide LAI_alternative")
    seconds = Float64.(Dates.value.(utc_dates .- start_date)) ./ 1000
    return TimeVaryingInput(
        seconds,
        stack([s.LAI for s in aligned]; dims = 1),
        surface_space;
        method = time_interpolation_method,
    )
end

"""
    build_callmip_simulation(sites, start_date, stop_date, Δt, toml_dict, FT;
                             zlim = FT.((-2, 0)),
                             nelements = 10,
                             domain_type = :ensemble)

Build the default `LandModel` on a `ColumnEnsemble` with one column per site,
driven by the sites' met forcing and LAI, and return the `LandSimulation`.
All model parameters and initial conditions are the defaults; diagnostics are
disabled (postprocessing is not yet correct for N > 1 columns).

Pass `domain_type = :column` (single site only) to use a plain `Column`
instead; its geometry differs from the ensemble's, so results agree with an
ensemble column only up to roundoff.
"""
function build_callmip_simulation(
    sites,
    start_date,
    stop_date,
    Δt,
    toml_dict,
    FT;
    zlim = FT.((-2, 0)),
    nelements = 10,
    domain_type = :ensemble,
)
    land_domain = if domain_type == :ensemble
        ColumnEnsemble(;
            zlim,
            nelements,
            longlat = [FT.((s.long, s.lat)) for s in sites],
        )
    elseif domain_type == :column
        length(sites) == 1 ||
            error("domain_type = :column supports a single site")
        Column(;
            zlim,
            nelements,
            longlat = FT.((sites[1].long, sites[1].lat)),
        )
    else
        error("unknown domain_type $domain_type")
    end
    surface_space = land_domain.space.surface
    forcing = prescribed_forcing_callmip(
        sites,
        surface_space,
        start_date,
        toml_dict,
        FT,
    )
    LAI = prescribed_lai_callmip(sites, surface_space, start_date)
    # Dynamic calls avoid a Julia 1.12 codegen crash ("Unreachable reached")
    # when the large default-kwarg constructors are inferred inside a function
    land = Base.invokelatest(
        LandModel{FT},
        forcing,
        LAI,
        toml_dict,
        land_domain,
        Δt,
    )
    return Base.invokelatest(
        LandSimulation,
        start_date,
        stop_date,
        Δt,
        land;
        diagnostics = (),
    )
end
