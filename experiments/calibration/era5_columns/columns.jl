# Selection of the simulated columns, and their CSV file.
#
# The candidates are the points of the ERA5 1° grid inside land (see
# `interior_land`) where all of `SELECTION_NAMES` are observed throughout
# `CLIMATOLOGY_YEARS`, which excludes the ice sheets that GOSIF does not cover.
# Each candidate is described by the observed seasonal climatology (DJF, MAM,
# JJA, SON) of these variables, and the candidates are grouped by
# area-weighted k-means. Each cluster contributes the point closest to its
# centroid, weighted by the cluster's share of the candidates' area, so the
# columns span the observed climates of global land and their weighted mean
# approximates a global land mean.

import DelimitedFiles
import NCDatasets
import Random
import Statistics
import ClimaLand

"""
    seasonal_climatology(short_name, years)

Return a `360 × 179 × 4` array of the DJF, MAM, JJA, and SON means of
`short_name` over `years` at the points of the ERA5 1° grid. A point is `NaN`
if any month is missing.
"""
function seasonal_climatology(short_name, years)
    months = month_starts(
        Dates.DateTime(first(years), 1, 1),
        Dates.DateTime(last(years) + 1, 1, 1),
    )
    data = read_monthly(short_name, months)
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
    interior_land()

Return a `360 × 179` mask of the points of the ERA5 1° grid that are land in
ClimaLand's 1° land sea mask, as are their eight neighbors.
"""
function interior_land()
    path = ClimaLand.Artifacts.landseamask_file_path(; resolution = "1deg")
    land = NCDatasets.NCDataset(path) do ds
        # The mask samples the 60 arcsecond mask near the integer degrees
        (
            all(abs.(ds["lon"][:] .- (-180:179)) .< 0.01) &&
            all(abs.(ds["lat"][:] .- (-90:89)) .< 0.01)
        ) || error("Unexpected land sea mask grid in $path")
        Array(ds["landsea"]) .== 1
    end
    is_land(lon, lat) = land[mod(lon + 180, 360) + 1, clamp(lat + 91, 1, 180)]
    return [
        all(is_land(lon + dlon, lat + dlat) for dlon in -1:1, dlat in -1:1) for
        lon in POINT_LONS, lat in POINT_LATS
    ]
end

"""
    weighted_kmeans(X, w, k; rng, max_iterations = 500)

Cluster the columns of `X` (features × points) into `k` clusters minimizing the
`w`-weighted within-cluster sum of squares (Lloyd's algorithm). Return the
cluster assignment of each point.

The initial centers are drawn by k-means++ (Arthur and Vassilvitskii, 2007,
"k-means++: The advantages of careful seeding"): each new center is a point
drawn with probability proportional to its weight times its squared distance
to the nearest center already drawn. Unlike centers drawn uniformly, this
spreads them over the data, and the expected sum of squares is within a factor
of O(log k) of the optimum.
"""
function weighted_kmeans(X, w, k; rng, max_iterations = 500)
    n = size(X, 2)
    sqdist(c, centers) = vec(sum(abs2, X .- centers[:, c], dims = 1))

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
    select_columns(; n_columns = N_COLUMNS, years = CLIMATOLOGY_YEARS, seed = RNG_SEED)

Return the longitude, latitude, and area weight (summing to one) of the
selected columns, sorted by latitude and longitude.
"""
function select_columns(;
    n_columns = N_COLUMNS,
    years = CLIMATOLOGY_YEARS,
    seed = RNG_SEED,
)
    climatologies = Dict(
        short_name => seasonal_climatology(short_name, years) for
        short_name in SELECTION_NAMES
    )
    observed = reduce(
        .&,
        [
            all(isfinite, climatologies[name], dims = 3)[:, :, 1] for
            name in SELECTION_NAMES
        ],
    )
    candidates = findall(vec(interior_land() .& observed))
    lonlat = [
        (POINT_LONS[c[1]], POINT_LATS[c[2]]) for c in
        CartesianIndices((length(POINT_LONS), length(POINT_LATS)))[candidates]
    ]
    area = [cosd(lat) for (_, lat) in lonlat]
    @info "Selecting $n_columns columns among $(length(candidates)) points"

    # Standardize each variable (over all seasons) so that they all contribute
    # to the distance, while keeping the seasonal cycle within a variable
    features = reduce(
        vcat,
        map(SELECTION_NAMES) do name
            values = reshape(climatologies[name], :, 4)[candidates, :]'
            μ = sum(values .* area') / (4 * sum(area))
            σ = sqrt(sum(abs2.(values .- μ) .* area') / (4 * sum(area)))
            (values .- μ) ./ σ
        end,
    )

    assignments = weighted_kmeans(
        features,
        area,
        n_columns;
        rng = Random.MersenneTwister(seed),
    )
    columns = map(1:n_columns) do c
        members = findall(==(c), assignments)
        centroid = features[:, members] * area[members] / sum(area[members])
        distances =
            vec(sum(abs2, features[:, members] .- centroid, dims = 1))
        (lonlat[members[argmin(distances)]]..., sum(area[members]) / sum(area))
    end
    return sort(columns, by = c -> (c[2], c[1]))
end

"""
    write_columns(path, columns)

Write the longitude, latitude, and area weight of `columns` to the CSV file at
`path`.
"""
function write_columns(path, columns)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "longitude,latitude,weight")
        for (lon, lat, weight) in columns
            println(io, "$lon,$lat,$(round(weight, sigdigits = 6))")
        end
    end
    return nothing
end

"""
    read_columns(path)

Return the `longlat` tuples and area weights of the columns listed in the CSV
file at `path`.
"""
function read_columns(path)
    table, header = DelimitedFiles.readdlm(path, ',', Float64; header = true)
    vec(header) == ["longitude", "latitude", "weight"] ||
        error("Unexpected header $header in $path")
    longlat = [(table[i, 1], table[i, 2]) for i in axes(table, 1)]
    return (; longlat, weight = table[:, 3])
end
