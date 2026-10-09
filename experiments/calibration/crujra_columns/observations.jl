# Reading the observations at the calibration columns and aggregating monthly
# series (observed or simulated) into seasonal means.

import Dates
import DelimitedFiles
import LinearAlgebra
import NCDatasets
import Statistics
import ClimaLand
import EnsembleKalmanProcesses as EKP

"""
    fluxcom_energy_dir()

Return the directory holding the FLUXCOM latent (`hfls`) and sensible (`hfss`)
heat flux files.

Until the `ilamb_fluxcom_energy` artifact is available, this defaults to a
committed 2008-2010 subset of it. Set `FLUXCOM_ENERGY_DIR` to use the full
record, as `select_columns.jl` needs.
"""
fluxcom_energy_dir() = get(
    ENV,
    "FLUXCOM_ENERGY_DIR",
    joinpath(@__DIR__, "tmp_artifacts", "ilamb_fluxcom_energy"),
)

"""
    ObsSource{F <: Function}

A monthly, regular latitude-longitude NetCDF file holding one observed
variable.
"""
struct ObsSource{F <: Function}
    "Path to the NetCDF file"
    path::String
    "Name of the variable in the file"
    varname::String
    "Converts the file's units to the units used for calibration"
    to_obs_units::F
end

"""
    obs_sources()

Return a dictionary mapping each calibrated short name to its `ObsSource`.

LHF, SHF, and GPP come from FLUXCOM (0.5°), SWU and LWU from CERES EBAF ed4.2
(1°). The energy fluxes are in W m^-2 and GPP in g C m^-2 day^-1.
"""
function obs_sources()
    ilamb(filename) = ClimaLand.Artifacts.ilamb_dataset_path(filename)
    fluxcom_energy(filename) = joinpath(fluxcom_energy_dir(), filename)
    MJ_per_day_to_W(x) = x / 0.0864
    return Dict(
        "lhf" => ObsSource(
            fluxcom_energy("hfls_FLUXCOM_hfls.nc"),
            "hfls",
            MJ_per_day_to_W,
        ),
        "shf" => ObsSource(
            fluxcom_energy("hfss_FLUXCOM_hfss.nc"),
            "hfss",
            MJ_per_day_to_W,
        ),
        "swu" => ObsSource(ilamb("rsus_CERESed4.2_rsus.nc"), "rsus", identity),
        "lwu" => ObsSource(ilamb("rlus_CERESed4.2_rlus.nc"), "rlus", identity),
        "gpp" => ObsSource(ilamb("gpp_FLUXCOM_gpp.nc"), "gpp", identity),
    )
end

"""
    lon_index(lon)
    lat_index(lat)

Return the index of the 1° cell containing `lon` (in degrees, either in
[-180, 180] or [0, 360]) or `lat`. Cell centers are at -179.5:179.5 and
-89.5:89.5.
"""
lon_index(lon) = clamp(floor(Int, mod(lon + 180, 360)) + 1, 1, 360)
lat_index(lat) = clamp(floor(Int, lat + 90) + 1, 1, 180)

"""
    cell_center(i_lon, i_lat)

Return the (longitude, latitude) of the center of a 1° cell.
"""
cell_center(i_lon, i_lat) = (i_lon - 180.5, i_lat - 90.5)

"""
    read_monthly_1deg(source::ObsSource, months)

Return a `360 × 180 × length(months)` array of the variable of `source`, in
calibration units, averaged onto 1° cells (see [`lon_index`](@ref)). `months`
are the first days of the months to read.

A 1° cell is `NaN` unless every source cell it contains is finite, so that only
cells fully covered by land data are kept.
"""
function read_monthly_1deg(source::ObsSource, months)
    NCDatasets.NCDataset(source.path) do ds
        file_months =
            [Dates.Date(Dates.year(t), Dates.month(t)) for t in ds["time"][:]]
        time_indices = map(months) do month
            i = findfirst(==(month), file_months)
            isnothing(i) && error("$(source.path) has no data for $month")
            i
        end
        lons = Float64.(ds["lon"][:])
        lats = Float64.(ds["lat"][:])
        cells_per_degree =
            round(Int, 1 / abs(lons[2] - lons[1])) *
            round(Int, 1 / abs(lats[2] - lats[1]))
        i_lons = lon_index.(lons)
        i_lats = lat_index.(lats)

        out = fill(NaN, 360, 180, length(months))
        sums = zeros(360, 180)
        counts = zeros(Int, 360, 180)
        for (k, t) in enumerate(time_indices)
            slice = coalesce.(ds[source.varname][:, :, t], NaN32)
            fill!(sums, 0)
            fill!(counts, 0)
            for (j, i_lat) in enumerate(i_lats), (i, i_lon) in enumerate(i_lons)
                value = slice[i, j]
                isfinite(value) || continue
                sums[i_lon, i_lat] += value
                counts[i_lon, i_lat] += 1
            end
            full = counts .== cells_per_degree
            view(out, :, :, k)[full] .=
                source.to_obs_units.(sums[full] ./ cells_per_degree)
        end
        return out
    end
end

"""
    month_starts(start_date, stop_date)

Return the first days of the months in `[start_date, stop_date)`.
"""
month_starts(start_date, stop_date) = collect(
    Dates.Date(start_date):Dates.Month(1):(Dates.Date(stop_date) - Dates.Day(1)),
)

"""
    read_columns(path = COLUMNS_FILE)

Return the `longlat` tuples and area weights (summing to one) of the columns
listed in the CSV file at `path`.
"""
function read_columns(path = COLUMNS_FILE)
    table, header = DelimitedFiles.readdlm(path, ',', Float64; header = true)
    vec(header) == ["longitude", "latitude", "weight"] ||
        error("Unexpected header $header in $path")
    longlat = [(table[i, 1], table[i, 2]) for i in axes(table, 1)]
    return (; longlat, weight = table[:, 3])
end

"""
    monthly_obs_at_columns(short_name, months, longlat)

Return a `length(months) × length(longlat)` matrix of the observed `short_name`
in the 1° cell containing each column.
"""
function monthly_obs_at_columns(short_name, months, longlat)
    data = read_monthly_1deg(obs_sources()[short_name], months)
    return reduce(
        hcat,
        [data[lon_index(lon), lat_index(lat), :] for (lon, lat) in longlat],
    )
end

"""
    season_starts(start_date, stop_date)

Return the first days of the meteorological seasons (DJF, MAM, JJA, SON) that
lie entirely within `[start_date, stop_date)`. `start_date` must be the first
day of a season.
"""
function season_starts(start_date, stop_date)
    start = Dates.Date(start_date)
    (Dates.day(start) == 1 && Dates.month(start) in (3, 6, 9, 12)) ||
        error("$start_date is not the first day of a season")
    return collect(
        start:Dates.Month(3):(Dates.Date(stop_date) - Dates.Month(3)),
    )
end

"""
    seasonal_means(monthly, months, seasons)

Average the rows of `monthly`, dated by the first days of `months`, over each
season in `seasons` (see [`season_starts`](@ref)). Any missing (`NaN`) month
makes the seasonal mean `NaN`.
"""
function seasonal_means(monthly::AbstractMatrix, months, seasons)
    return reduce(
        vcat,
        map(seasons) do season
            rows = map(0:2) do k
                row = findfirst(==(season + Dates.Month(k)), months)
                isnothing(row) &&
                    error("No data for $(season + Dates.Month(k))")
                row
            end
            Statistics.mean(monthly[rows, :], dims = 1)
        end,
    )
end

"""
    observed_seasonal_means(longlat)

Return a dictionary mapping each of `SHORT_NAMES` to the `seasons × columns`
matrix of its observed seasonal means over the calibration period.
"""
function observed_seasonal_means(longlat)
    months = month_starts(CALIBRATION_START, STOP_DATE)
    seasons = season_starts(CALIBRATION_START, STOP_DATE)
    return Dict(
        short_name => seasonal_means(
            monthly_obs_at_columns(short_name, months, longlat),
            months,
            seasons,
        ) for short_name in SHORT_NAMES
    )
end

"""
    entry_weights(weight, n_seasons)

Return the area weight, normalized to a mean of one over the columns, of each
entry of a vectorized `n_seasons × columns` matrix.
"""
entry_weights(weight, n_seasons) = vec(
    repeat(permutedims(weight .* (length(weight) / sum(weight))), n_seasons),
)

"""
    make_observation(observed, masks, weight)

Return the `EKP.Observation` of the entries in `masks` of the `observed`
seasonal means, stacked in the order of `SHORT_NAMES`. The covariance is
diagonal, with variance `NOISE_STD[name]^2` divided by the relative area
weight of the column.
"""
function make_observation(observed, masks, weight)
    n_seasons = size(observed[first(SHORT_NAMES)], 1)
    relative_weight = entry_weights(weight, n_seasons)
    return EKP.combine_observations(
        map(SHORT_NAMES) do name
            mask = masks[name]
            EKP.Observation(
                vec(observed[name])[mask],
                LinearAlgebra.Diagonal(
                    NOISE_STD[name]^2 ./ relative_weight[mask],
                ),
                name,
            )
        end,
    )
end
