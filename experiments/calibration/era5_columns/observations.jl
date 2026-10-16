# Reading the observations at the points of the ERA5 1° grid and aggregating
# monthly series (observed or simulated) into seasonal means.

import Dates
import LinearAlgebra
import NCDatasets
import Statistics
import ClimaLand
import EnsembleKalmanProcesses as EKP

# The observations are read on the points of the ERA5 1° grid, excluding the
# poles: longitudes -180:179 and latitudes -89:89, indexed by `point_index`.
const POINT_LONS = -180:179
const POINT_LATS = -89:89

"""
    point_index(lon, lat)

Return the indices in `POINT_LONS` and `POINT_LATS` of a point of the ERA5 1°
grid.
"""
function point_index(lon, lat)
    (isinteger(lon) && isinteger(lat) && abs(lat) < 90) ||
        error("($lon, $lat) is not a point of the ERA5 1° grid")
    return (mod(Int(lon) + 180, 360) + 1, Int(lat) + 90)
end

"""ERA5 variable of each energy flux, and the sign that makes it positive
upward."""
const ERA5_VARIABLES = Dict(
    "lhf" => ("mslhf", -1),
    "shf" => ("msshf", -1),
    "swu" => ("msuwswrf", 1),
    "lwu" => ("msuwlwrf", 1),
)

"""
    read_months(ds, varname, months)

Return the monthly time slices of the variable `varname` (dimensions
longitude, latitude, time) of the NetCDF dataset `ds` for the months whose
first days are `months`, with `NaN` for missing values.
"""
function read_months(ds, varname, months)
    file_months =
        [Dates.Date(Dates.year(t), Dates.month(t)) for t in ds["time"][:]]
    indices = map(months) do month
        i = findfirst(==(month), file_months)
        isnothing(i) && error("No data for $month")
        i
    end
    first_index, last_index = extrema(indices)
    data = ds[varname][:, :, first_index:last_index]
    return Float64.(coalesce.(data[:, :, indices .- (first_index - 1)], NaN))
end

"""
    read_era5_monthly(short_name, months)

Return a `360 × 179 × length(months)` array of the ERA5 monthly mean of the
energy flux `short_name` (W m^-2, positive upward) at the points of the ERA5 1°
grid. `months` are the first days of the months to read.
"""
function read_era5_monthly(short_name, months)
    varname, sign = ERA5_VARIABLES[short_name]
    path = joinpath(
        ClimaLand.Artifacts.era5_monthly_averages_single_level_path(),
        "era5_monthly_averages_surface_single_level_197901-202410.nc",
    )
    NCDatasets.NCDataset(path) do ds
        (ds["longitude"][:] == 0:359 && ds["latitude"][:] == -90:90) ||
            error("Unexpected ERA5 grid in $path")
        i_lons = [mod(lon, 360) + 1 for lon in POINT_LONS]
        j_lats = [lat + 91 for lat in POINT_LATS]
        return sign .* read_months(ds, varname, months)[i_lons, j_lats, :]
    end
end

"""
    read_gosif_monthly(months)

Return a `360 × 179 × length(months)` array of the GOSIF GPP (g C m^-2 day^-1)
of the `inversion_nee` artifact at the points of the ERA5 1° grid. `months`
are the first days of the months to read.

The GOSIF 1° cells are centered half a degree from the points, so the value at
a point is the area-weighted mean of the four cells around it: the mean over
the 1° cell centered on the point. It is `NaN` unless all four are observed.
"""
function read_gosif_monthly(months)
    path = ClimaLand.Artifacts.inversion_nee_dataset_path()
    NCDatasets.NCDataset(path) do ds
        lons, lats = ds["lon"][:], ds["lat"][:]
        (lons == -179.5:179.5 && lats == -89.5:89.5) ||
            error("Unexpected GOSIF grid in $path")
        ds["gpp"].attrib["units"] == "g C m-2 month-1" ||
            error("Unexpected GOSIF units in $path")
        cells =
            read_months(ds, "gpp", months) ./
            reshape(Dates.daysinmonth.(months), 1, 1, :)
        out = zeros(length(POINT_LONS), length(POINT_LATS), length(months))
        for (j, lat) in enumerate(POINT_LATS), (i, lon) in enumerate(POINT_LONS)
            # Cells centered at lon ± 0.5 and lat ± 0.5
            i_cells = mod1.((lon + 180, lon + 181), 360)
            j_cells = (lat + 90, lat + 91)
            weights = [cosd(lats[jc]) for ic in i_cells, jc in j_cells]
            out[i, j, :] = sum(
                weights[a, b] .* cells[i_cells[a], j_cells[b], :] for
                a in 1:2, b in 1:2
            )
            out[i, j, :] ./= sum(weights)
        end
        return out
    end
end

"""
    read_monthly(short_name, months)

Return a `360 × 179 × length(months)` array of the monthly mean observations
of `short_name` at the points of the ERA5 1° grid: GOSIF for `gpp`, ERA5 for
the energy fluxes. `months` are the first days of the months to read.
"""
read_monthly(short_name, months) =
    short_name == "gpp" ? read_gosif_monthly(months) :
    read_era5_monthly(short_name, months)

"""
    month_starts(start_date, stop_date)

Return the first days of the months in `[start_date, stop_date)`.
"""
month_starts(start_date, stop_date) = collect(
    Dates.Date(start_date):Dates.Month(1):(Dates.Date(stop_date) - Dates.Day(1)),
)

"""
    monthly_obs_at_columns(short_name, months, longlat)

Return a `length(months) × length(longlat)` matrix of the observed `short_name`
at the columns, which are points of the ERA5 1° grid.
"""
function monthly_obs_at_columns(short_name, months, longlat)
    data = read_monthly(short_name, months)
    return reduce(
        hcat,
        [data[point_index(lon, lat)..., :] for (lon, lat) in longlat],
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
