# Select the columns simulated by the CRUJRA column-ensemble calibration and
# write them to `columns.csv`.
#
# Candidates are the 1° land cells observed throughout `CLIMATOLOGY_YEARS` (see
# `select_columns` for the two strata). Each candidate is described by the
# observed seasonal climatology (DJF, MAM, JJA, SON) of the calibrated variables,
# and the candidates are grouped by area-weighted k-means. Each cluster
# contributes the cell closest to its centroid, weighted by the cluster's share
# of the candidates' area, so the columns span the observed climates of global
# land and their weighted mean approximates a global land mean.
#
# This needs the full FLUXCOM record (not the committed 2008-2010 subset):
#   FLUXCOM_ENERGY_DIR=<ilamb_fluxcom_energy artifact> \
#       julia --project=.buildkite experiments/calibration/crujra_columns/select_columns.jl

import Dates
import Random
import Statistics
import DelimitedFiles

include(joinpath(@__DIR__, "config.jl"))
include(joinpath(@__DIR__, "observations.jl"))

"""Years of the climatology used to describe the candidate cells, the overlap
of CERES (from March 2000) and FLUXCOM GPP (until December 2013)."""
const CLIMATOLOGY_YEARS = 2001:2013

const N_COLUMNS = 200

"""
    seasonal_climatology(short_name, years)

Return a `360 × 180 × 4` array of the DJF, MAM, JJA, and SON means of
`short_name` over `years`, on 1° cells. A cell is `NaN` if any month is
missing.
"""
function seasonal_climatology(short_name, years)
    months = month_starts(
        Dates.DateTime(first(years), 1, 1),
        Dates.DateTime(last(years) + 1, 1, 1),
    )
    data = read_monthly_1deg(obs_sources()[short_name], months)
    season_of_month(m) = mod(div(Dates.month(m), 3), 4) + 1
    return cat(
        [
            Statistics.mean(
                data[:, :, findall(m -> season_of_month(m) == s, months)],
                dims = 3,
            ) for s in 1:4
        ]...;
        dims = 3,
    )
end

"""
    weighted_kmeans(X, w, k; rng, max_iterations = 500)

Cluster the columns of `X` (features × points) into `k` clusters minimizing the
`w`-weighted within-cluster sum of squares, with k-means++ initialization.
Return the cluster assignment of each point.
"""
function weighted_kmeans(X, w, k; rng, max_iterations = 500)
    n = size(X, 2)
    sqdist(c, centers) = vec(sum(abs2, X .- centers[:, c], dims = 1))

    # k-means++: draw each new center with probability ∝ w * (distance to the
    # nearest center)^2
    centers = zeros(size(X, 1), k)
    centers[:, 1] = X[:, findfirst(cumsum(w) .>= rand(rng) * sum(w))]
    nearest = sqdist(1, centers)
    for c in 2:k
        p = cumsum(w .* nearest)
        centers[:, c] = X[:, findfirst(p .>= rand(rng) * last(p))]
        nearest = min.(nearest, sqdist(c, centers))
    end

    assignments = zeros(Int, n)
    for iteration in 1:max_iterations
        distances = reduce(hcat, [sqdist(c, centers) for c in 1:k])
        new_assignments = [argmin(row) for row in eachrow(distances)]
        if new_assignments == assignments
            @info "k-means converged after $iteration iterations"
            return assignments
        end
        assignments = new_assignments
        for c in 1:k
            members = findall(==(c), assignments)
            if isempty(members)
                # Reseed an empty cluster at the worst-represented point
                worst = argmax(w .* minimum(distances, dims = 2)[:])
                centers[:, c] = X[:, worst]
            else
                centers[:, c] = X[:, members] * w[members] / sum(w[members])
            end
        end
    end
    @warn "k-means did not converge after $max_iterations iterations"
    return assignments
end

"""
    interior_land()

Return a `360 × 180` mask of the 1° cells that are land in ClimaLand's 1° land
sea mask, as are their eight neighbors.
"""
function interior_land()
    land = NCDatasets.NCDataset(
        ClimaLand.Artifacts.landseamask_file_path(; resolution = "1deg"),
    ) do ds
        Array(ds["landsea"]) .== 1
    end
    size(land) == (360, 180) || error("Unexpected land sea mask size")
    return [
        all(
            land[mod1(i + di, 360), clamp(j + dj, 1, 180)] for
            di in -1:1, dj in -1:1
        ) for i in 1:360, j in 1:180
    ]
end

"""
    cluster_cells(climatologies, short_names, mask, n_columns, rng)

Cluster the 1° cells in `mask` by the seasonal climatology of `short_names` and
return the center (longitude, latitude) of the cell closest to each cluster's
centroid, with the cluster's area (sum of the cosines of the latitudes).
"""
function cluster_cells(climatologies, short_names, mask, n_columns, rng)
    candidates = findall(vec(mask))
    lonlat = [
        cell_center(Tuple(c)...) for
        c in CartesianIndices((360, 180))[candidates]
    ]
    area = [cosd(lat) for (_, lat) in lonlat]

    # Standardize each variable (over all seasons) so that they all contribute
    # to the distance, while keeping the seasonal cycle within a variable
    features = reduce(
        vcat,
        map(short_names) do name
            values = reshape(climatologies[name], 360 * 180, 4)[candidates, :]'
            μ = sum(values .* area') / (4 * sum(area))
            σ = sqrt(sum(abs2.(values .- μ) .* area') / (4 * sum(area)))
            (values .- μ) ./ σ
        end,
    )

    assignments = weighted_kmeans(features, area, n_columns; rng)
    return map(1:n_columns) do c
        members = findall(==(c), assignments)
        centroid = features[:, members] * area[members] / sum(area[members])
        distances = vec(sum(abs2, features[:, members] .- centroid, dims = 1))
        (lonlat[members[argmin(distances)]]..., sum(area[members]))
    end
end

"""
    select_columns(; n_columns = N_COLUMNS, years = CLIMATOLOGY_YEARS, seed = RNG_SEED)

Return the longitudes, latitudes, and area weights of the selected columns.

The candidates are split in two strata, which get columns in proportion to
their area:
- cells where all of `SHORT_NAMES` are observed, clustered on all of them;
- interior land cells equatorward of 60° that FLUXCOM does not cover (mostly
  deserts), clustered on SWU and LWU, the only variables observed there. This
  keeps bare soil in the calibration of the radiative parameters.
"""
function select_columns(;
    n_columns = N_COLUMNS,
    years = CLIMATOLOGY_YEARS,
    seed = RNG_SEED,
)
    climatologies = Dict(
        short_name => seasonal_climatology(short_name, years) for
        short_name in SHORT_NAMES
    )
    observed(names) = reduce(
        .&,
        [
            all(isfinite, climatologies[name], dims = 3)[:, :, 1] for
            name in names
        ],
    )
    radiation_names = ["swu", "lwu"]
    latitudes = [cell_center(i, j)[2] for i in 1:360, j in 1:180]

    fully_observed = observed(SHORT_NAMES)
    radiation_only =
        interior_land() .& .!fully_observed .& observed(radiation_names) .&
        (abs.(latitudes) .< 60)
    area(mask) = sum(cosd.(latitudes[mask]))
    total_area = area(fully_observed) + area(radiation_only)
    n_radiation_only = round(Int, n_columns * area(radiation_only) / total_area)
    @info "Candidate 1° cells" fully_observed = count(fully_observed) radiation_only =
        count(radiation_only) n_radiation_only

    rng = Random.MersenneTwister(seed)
    columns = vcat(
        cluster_cells(
            climatologies,
            SHORT_NAMES,
            fully_observed,
            n_columns - n_radiation_only,
            rng,
        ),
        cluster_cells(
            climatologies,
            radiation_names,
            radiation_only,
            n_radiation_only,
            rng,
        ),
    )
    columns = [(lon, lat, a / total_area) for (lon, lat, a) in columns]
    return sort(columns, by = c -> (c[2], c[1]))
end

if abspath(PROGRAM_FILE) == @__FILE__
    columns = select_columns()
    open(COLUMNS_FILE, "w") do io
        println(io, "longitude,latitude,weight")
        for (lon, lat, weight) in columns
            println(io, "$lon,$lat,$(round(weight, sigdigits = 6))")
        end
    end
    @info "Wrote $(length(columns)) columns to $COLUMNS_FILE"
end
