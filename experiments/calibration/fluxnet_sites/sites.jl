# Reading the FLUXNET sites: their forcing on the calendar of the simulations,
# and the monthly means of their observed fluxes.
#
# The sites run together as a `ColumnEnsemble`, with one clock, but their
# records cover different years. Each site contributes up to
# `length(CALIBRATION_YEARS)` consecutive calendar years of complete forcing,
# which are simulated as `CALIBRATION_YEARS` (see `source_years`); a site with
# fewer years repeats them, and only their last repetition is compared to
# observations. The LAI is the MODIS climatology, so it does not depend on the
# year either.

import CSV
import DelimitedFiles
import Dates
import JLD2
import LinearAlgebra
import Statistics
import ClimaLand
import ClimaLand.FluxnetSimulations as FluxnetSimulations
import ClimaLand.Parameters as LP
import EnsembleKalmanProcesses as EKP

include(
    joinpath(
        @__DIR__,
        "..",
        "..",
        "integrated",
        "generic_site",
        "list_fluxnet_sites.jl",
    ),
)

const FluxnetSimulationsExt =
    Base.get_extension(ClimaLand, :FluxnetSimulationsExt)

"""FLUXNET column of each observed variable, and the column of its quality flag
(`nothing` for a variable that is not gap-filled)."""
const OBSERVED_COLUMNS = Dict(
    "lhf" => ("LE_CORR", "LE_F_MDS_QC"),
    "shf" => ("H_CORR", "H_F_MDS_QC"),
    "swu" => ("SW_OUT", nothing),
    "lwu" => ("LW_OUT", nothing),
    "gpp" => ("GPP_DT_VUT_REF", "NEE_VUT_REF_QC"),
)

"""Converts an observed variable to the observation units: GPP from
μmol CO2 m^-2 s^-1 to g C m^-2 day^-1. The energy fluxes are in W m^-2."""
const FLUXNET_TO_OBS_UNITS = Dict("gpp" => x -> x * 1e-6 * 12.011 * 86400)

"""Forcing variables saved for the simulations, in model units."""
const FORCING_NAMES = (:T, :P, :u, :q, :SW_d, :LW_d, :c_co2, :P_liq, :P_snow)

"""Time step of the forcing saved for the simulations. Hourly records are
interpolated to it."""
const FORCING_STEP = Dates.Minute(30)

"""Annual mean CO2 (ppm) at Mauna Loa (NOAA GML) at the start of each decade,
interpolated for the sites that do not measure CO2."""
const MAUNA_LOA_CO2 =
    (years = [1990, 2000, 2010, 2020], ppm = [354.4, 369.7, 389.9, 414.2])

const EPOCH = Dates.DateTime(1900)

seconds_since_epoch(t) = Dates.value(t - EPOCH) / 1000

"""
    parse_timestamp(x)

Return the `DateTime` of a FLUXNET timestamp (`yyyymmddHHMM`) read as a number.
"""
function parse_timestamp(x)
    n = round(Int, x)
    return Dates.DateTime(
        div(n, 10^8),
        rem(div(n, 10^6), 100),
        rem(div(n, 10^4), 100),
        rem(div(n, 100), 100),
        rem(n, 100),
    )
end

"""
    site_paths()

Return the IDs of the sites and the paths to their half-hourly or hourly CSV
files: those listed, comma separated, in the environment variable
`FLUXNET_SITES`, as found by `ClimaLand.Artifacts.experiment_fluxnet_data_path`
(US-MOz, US-NR1, US-Ha1 and US-Var come from the downloadable `fluxnet_sites`
artifact), or else every site of the `fluxnet2015` artifact.
"""
function site_paths()
    if haskey(ENV, "FLUXNET_SITES")
        site_IDs = String.(strip.(split(ENV["FLUXNET_SITES"], ',')))
        return [
            site_ID =>
                ClimaLand.Artifacts.experiment_fluxnet_data_path(site_ID) for
            site_ID in site_IDs
        ]
    end
    return [
        site_ID => fluxnet2015_csv_path(site_ID) for
        site_ID in list_fluxnet_sites()
    ]
end

"""
    fluxnet2015_csv_path(site_ID)

Return the path to the half-hourly, or else hourly, FULLSET CSV file of
`site_ID` in the `fluxnet2015` artifact.
"""
function fluxnet2015_csv_path(site_ID)
    root = ClimaLand.Artifacts.fluxnet2015_data_path()
    site_dir = only(
        filter(
            d ->
                startswith(d, "FLX_$(site_ID)_") && isdir(joinpath(root, d)),
            readdir(root),
        ),
    )
    for resolution in ("HH", "HR")
        file = replace(site_dir, "FULLSET" => "FULLSET_$resolution"; count = 1)
        path = joinpath(root, site_dir, "$file.csv")
        isfile(path) && return path
    end
    error("No half-hourly or hourly CSV file in $(joinpath(root, site_dir))")
end

"""
    read_site_csv(path)

Return a dictionary from the FLUXNET columns used here to their values in the
CSV file at `path`, with `NaN` for missing values. Columns absent from the file
are absent from the dictionary.
"""
function read_site_csv(path)
    observed_columns = [
        column for columns in values(OBSERVED_COLUMNS) for
        column in columns if !isnothing(column)
    ]
    wanted = unique([
        "TIMESTAMP_START",
        "TIMESTAMP_END",
        FluxnetSimulations.FLUXNET_FORCING_COLUMNS...,
        observed_columns...,
    ])
    present = intersect(wanted, split(readline(path), ','))
    file = CSV.File(path; select = present, types = Float64)
    to_value(x) = ismissing(x) || x == -9999 ? NaN : Float64(x)
    return Dict(
        name => to_value.(getproperty(file, Symbol(name))) for name in present
    )
end

"""
    interpolate_finite(times, values, query_times)

Linearly interpolate the finite `values` at the increasing `times` to
`query_times`, holding the first and last finite values beyond them.
"""
function interpolate_finite(times, values, query_times)
    valid = isfinite.(values)
    t, v = times[valid], values[valid]
    isempty(t) && error("No finite value to interpolate")
    return map(query_times) do query
        i = searchsortedlast(t, query)
        i == 0 && return v[1]
        i == length(t) && return v[end]
        w = (query - t[i]) / (t[i + 1] - t[i])
        return v[i] + w * (v[i + 1] - v[i])
    end
end

"""
    monthly_observations(data, utc_start, steps_per_day, years)

Return a dictionary mapping each of `SHORT_NAMES` to the `12 × length(years)`
matrix of its monthly means (UTC months) in `years`, in observation units, with
`NaN` for a month without enough data (see `MIN_GOOD_FRACTION` and
`MIN_COVERAGE`). `utc_start` are the starts of the averaging periods of the
rows of `data`.
"""
function monthly_observations(data, utc_start, steps_per_day, years)
    rows_of_month = Dict{Dates.Date, Vector{Int}}()
    for (i, t) in enumerate(utc_start)
        push!(get!(rows_of_month, Dates.Date(Dates.yearmonth(t)...), Int[]), i)
    end
    return Dict(
        map(SHORT_NAMES) do name
            monthly = fill(NaN, 12, length(years))
            value_column, qc_column = OBSERVED_COLUMNS[name]
            has_columns =
                haskey(data, value_column) &&
                (isnothing(qc_column) || haskey(data, qc_column))
            has_columns || return name => monthly
            values = data[value_column]
            good =
                isnothing(qc_column) ? isfinite.(values) :
                isfinite.(values) .& (data[qc_column] .<= 1)
            convert_units = get(FLUXNET_TO_OBS_UNITS, name, identity)
            for (j, year) in enumerate(years), month in 1:12
                rows = get(rows_of_month, Dates.Date(year, month), Int[])
                expected = Dates.daysinmonth(year, month) * steps_per_day
                observed = filter(i -> isfinite(values[i]), rows)
                enough =
                    length(observed) >= MIN_COVERAGE * expected &&
                    count(good[rows]) >= MIN_GOOD_FRACTION * expected
                enough || continue
                monthly[month, j] =
                    convert_units(Statistics.mean(values[observed]))
            end
            name => monthly
        end,
    )
end

"""
    best_years(complete_years, n_observed)

Return the consecutive years among `complete_years`, at most
`length(CALIBRATION_YEARS)`, with the most observed monthly means
(`n_observed[year]`), the latest in a tie, or no year if none is observed.
"""
function best_years(complete_years, n_observed)
    best, best_count = Int[], 0
    for first_year in sort(complete_years)
        years = Int[]
        for year in first_year .+ (0:(length(CALIBRATION_YEARS) - 1))
            year in complete_years || break
            push!(years, year)
        end
        n = sum(n_observed[year] for year in years)
        if n > 0 && n >= best_count
            best, best_count = years, n
        end
    end
    return best
end

"""
    source_years(years)

Return the year of the site's record simulated in each year from
`SPINUP_START`, and whether each of `CALIBRATION_YEARS` is compared to
observations. The site's consecutive `years` are simulated as the last
calibration years and compared to observations; the earlier years, starting
with the spinup year, cycle through them.
"""
function source_years(years)
    n_years, n = length(CALIBRATION_YEARS), length(years)
    sources = [years[mod1(k - (n_years - n), n)] for k in 0:n_years]
    compared = [k > n_years - n for k in 1:n_years]
    return sources, compared
end

"""
    source_time(t, sources)

Return the time of the site's record simulated at time `t`, where `sources` are
the years of the record simulated in each year from `SPINUP_START` (see
`source_years`). February 29 maps to February 28 of a common year.
"""
function source_time(t, sources)
    years_since_start = Dates.year(t) - Dates.year(SPINUP_START)
    k = clamp(years_since_start, 0, length(sources) - 1)
    year = sources[k + 1] + years_since_start - k
    day = Dates.day(t)
    Dates.month(t) == 2 && day == 29 && !Dates.isleapyear(year) && (day = 28)
    return Dates.DateTime(
        year,
        Dates.month(t),
        day,
        Dates.hour(t),
        Dates.minute(t),
    )
end

"""Times of the forcing saved for the simulations."""
forcing_times() = collect(SPINUP_START:FORCING_STEP:STOP_DATE)

"""
    mauna_loa_co2(t)

Return the annual mean CO2 (ppm) at Mauna Loa interpolated at the time `t`.
"""
mauna_loa_co2(t) = only(
    interpolate_finite(
        Float64.(MAUNA_LOA_CO2.years),
        MAUNA_LOA_CO2.ppm,
        [Dates.year(t) + (Dates.dayofyear(t) - 1) / Dates.daysinyear(t)],
    ),
)

"""
    read_site(site_ID, path, thermo_params)

Return the site `site_ID`, from its FLUXNET CSV file at `path`: its location,
tower height, simulated years, forcing at `forcing_times()` in model units, and
observed monthly means over `CALIBRATION_YEARS` (`NaN` where not compared). Return
`nothing` if no year of the site has complete forcing and an observed monthly
mean.
"""
function read_site(site_ID, path, thermo_params)
    site_val = FluxnetSimulations.replace_hyphen(site_ID)
    (; time_offset, lat, long) =
        FluxnetSimulations.get_location(FT, Val(site_val))
    (; atmos_h) = FluxnetSimulations.get_fluxtower_height(FT, Val(site_val))
    data = read_site_csv(path)

    # FLUXNET timestamps are in local standard time
    local_start = parse_timestamp.(data["TIMESTAMP_START"])
    local_end = parse_timestamp.(data["TIMESTAMP_END"])
    data_dt = only(unique(local_end .- local_start))
    utc_start = local_start .- Dates.Minute(round(Int, 60 * time_offset))
    utc_mid = seconds_since_epoch.(utc_start .+ data_dt / 2)
    steps_per_day = Dates.Minute(Dates.Day(1)) ÷ Dates.Minute(data_dt)

    meteorology =
        filter(!=("CO2_F_MDS"), FluxnetSimulations.FLUXNET_FORCING_COLUMNS)
    all(haskey(data, name) for name in meteorology) || return nothing
    complete = reduce(.&, (isfinite.(data[name]) for name in meteorology))
    local_year = Dates.year.(local_start)
    complete_years = filter(unique(local_year)) do year
        expected = Dates.daysinyear(year) * steps_per_day
        count(complete[local_year .== year]) >= MIN_FORCING_COVERAGE * expected
    end
    isempty(complete_years) && return nothing
    monthly =
        monthly_observations(data, utc_start, steps_per_day, complete_years)
    n_observed = Dict(
        year => sum(
            count(isfinite, monthly[name][:, j]) for name in SHORT_NAMES
        ) for (j, year) in enumerate(complete_years)
    )
    years = best_years(complete_years, n_observed)
    isempty(years) && return nothing
    sources, compared = source_years(years)

    query =
        [seconds_since_epoch(source_time(t, sources)) for t in forcing_times()]
    at_forcing_times(name) = interpolate_finite(utc_mid, data[name], query)
    TA, VPD, PA, WS, SW, LW = at_forcing_times.([
        "TA_F",
        "VPD_F",
        "PA_F",
        "WS_F",
        "SW_IN_F",
        "LW_IN_F",
    ],)
    # P_F is accumulated over a record, in mm
    precip_rate = interpolate_finite(
        utc_mid,
        data["P_F"] ./ (1000 * Dates.value(Dates.Second(data_dt))),
        query,
    )
    snow_fraction =
        FluxnetSimulationsExt.snow_precip_fraction.(TA, VPD; thermo_params)
    in_years = in(years).(local_year)
    measures_co2 =
        haskey(data, "CO2_F_MDS") &&
        count(isfinite, data["CO2_F_MDS"][in_years]) >= 0.5 * count(in_years)
    co2 =
        measures_co2 ? at_forcing_times("CO2_F_MDS") :
        mauna_loa_co2.(source_time.(forcing_times(), Ref(sources)))
    forcing = (;
        T = TA .+ 273.15,
        P = PA .* 1000,
        u = WS,
        # A VPD above the saturation vapor pressure gives a negative q
        q = max.(
            FluxnetSimulationsExt.compute_q.(PA, VPD, TA; thermo_params),
            0,
        ),
        SW_d = max.(SW, 0),
        LW_d = LW,
        c_co2 = co2 .* 1e-6,
        P_liq = .-precip_rate .* (1 .- snow_fraction),
        P_snow = .-precip_rate .* snow_fraction,
    )

    observed = Dict(
        map(SHORT_NAMES) do name
            values = fill(NaN, 12, length(CALIBRATION_YEARS))
            for k in findall(compared)
                values[:, k] = monthly[name][
                    :,
                    findfirst(==(sources[k + 1]), complete_years),
                ]
            end
            name => vec(values)
        end,
    )
    return (;
        site_ID,
        long,
        lat,
        time_offset,
        atmos_h,
        data_dt_minutes = Dates.value(Dates.Minute(data_dt)),
        years,
        measures_co2,
        forcing,
        observed,
    )
end

"""
    build_site_data(data_path, sites_path)

Read the sites of `site_paths()` and save, for the sites with a year of
complete forcing and an observed monthly mean, the forcing of the simulations
and the observed monthly means to the JLD2 file `data_path`, and the table of
the sites to the CSV file `sites_path`.
"""
function build_site_data(data_path, sites_path)
    thermo_params =
        LP.thermodynamic_parameters(LP.LandParameters(LP.create_toml_dict(FT)))
    sites = []
    for (site_ID, csv_path) in site_paths()
        site = try
            read_site(site_ID, csv_path, thermo_params)
        catch error
            @warn "Skipping $site_ID, which could not be read" exception =
                error
            nothing
        end
        if isnothing(site)
            @info "Skipping $site_ID: no year with complete forcing and observations"
        else
            push!(sites, site)
        end
    end
    isempty(sites) &&
        error("No site has a year with complete forcing and observations")
    @info "Calibrating on $(length(sites)) sites"

    times = Dates.value.(forcing_times() .- SPINUP_START) ./ 1000
    # One row per site, in the order of the columns
    forcing = NamedTuple{FORCING_NAMES}(
        map(FORCING_NAMES) do name
            Float32.(
                stack(
                    [getproperty(site.forcing, name) for site in sites];
                    dims = 1,
                ),
            )
        end,
    )
    observed = Dict(
        name => stack([site.observed[name] for site in sites]) for
        name in SHORT_NAMES
    )
    mkpath(dirname(data_path))
    JLD2.jldsave(data_path; times, forcing, observed)
    write_sites(sites_path, sites)
    return nothing
end

"""
    write_sites(path, sites)

Write the table of the `sites` to the CSV file at `path`.
"""
function write_sites(path, sites)
    mkpath(dirname(path))
    CSV.write(
        path,
        (;
            site = [site.site_ID for site in sites],
            longitude = [site.long for site in sites],
            latitude = [site.lat for site in sites],
            utc_offset = [site.time_offset for site in sites],
            tower_height = [site.atmos_h for site in sites],
            time_step_minutes = [site.data_dt_minutes for site in sites],
            first_year = [first(site.years) for site in sites],
            last_year = [last(site.years) for site in sites],
            measures_co2 = [site.measures_co2 for site in sites],
        ),
    )
    return nothing
end

"""
    read_sites(path)

Return the site IDs, locations (`(long, lat)`), and tower heights in the CSV
file at `path` written by `write_sites`.
"""
function read_sites(path)
    table = CSV.File(path)
    return (;
        site_IDs = String.(table.site),
        longlat = collect(
            zip(Float64.(table.longitude), Float64.(table.latitude)),
        ),
        tower_height = Float64.(table.tower_height),
    )
end

"""
    make_observation(observed, masks)

Return the `EKP.Observation` of the entries in `masks` of the `observed`
monthly means, stacked in the order of `SHORT_NAMES`, with a diagonal
covariance of variance `NOISE_STD[name]^2`.
"""
function make_observation(observed, masks)
    observations = map(filter(name -> any(masks[name]), SHORT_NAMES)) do name
        mask = masks[name]
        EKP.Observation(
            vec(observed[name])[mask],
            LinearAlgebra.Diagonal(fill(NOISE_STD[name]^2, count(mask))),
            name,
        )
    end
    return EKP.combine_observations(observations)
end
