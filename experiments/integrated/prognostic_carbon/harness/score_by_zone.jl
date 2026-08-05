## Biomass skill broken down by broad climate zone.
##
## A global mean hides compensating errors: the model can be unbiased overall
## while being systematically high in one zone and low in another, and the whole
## story of this model is that its error is regionally structured rather than
## uniform. This reports bias and RMSE per zone so that structure is visible in
## the summary rather than only in a map.
##
## Zones are defined from the model's own drivers, not from a map, so they are
## reproducible and carry no extra boundary condition:
##
##   semi-arid   MAP < 0.6 m/yr, at any latitude - taken first, because a dry
##               cell behaves as a dry cell whatever its latitude
##   tropics     |lat| < 23.5
##   temperate   23.5 <= |lat| < 50
##   boreal      |lat| >= 50
##
## The semi-arid class deliberately cuts across the latitude bands. Reporting
## Sahara inside "tropics" would put the model's largest relative errors in the
## same row as closed-canopy rainforest, which is exactly the averaging this
## breakdown exists to undo.
##
## Usage:
##   julia --project=.buildkite score_by_zone.jl <equilibrium_carbon.nc> <driver_outdir>

include(joinpath(@__DIR__, "global_equilibrium.jl"))

const ARID_MAP = FT(0.6)     # m/yr
const TROPIC_LAT = FT(23.5)
const BOREAL_LAT = FT(50.0)

function zone_of(lat, map_)
    isfinite(map_) && map_ < ARID_MAP && return :semiarid
    abs(lat) < TROPIC_LAT && return :tropics
    abs(lat) < BOREAL_LAT && return :temperate
    return :boreal
end

const ZONES = (:tropics, :temperate, :boreal, :semiarid)
const ZONE_LABEL = Dict(
    :tropics => "tropics",
    :temperate => "temperate",
    :boreal => "boreal",
    :semiarid => "semi-arid",
)

function main(ncpath, driverdir)
    ds = NCDatasets.NCDataset(ncpath)
    lons = Array{FT}(Array(ds["lon"][:]))
    lats = Array{FT}(Array(ds["lat"][:]))
    # C_stem, matching every other biomass comparison in this harness.
    M = Array{Union{Missing, FT}}(Array(ds["C_stem"][:, :]))
    close(ds)

    _, _, pra = read_monthly(driverdir, "pra")
    MAP = fill(FT(NaN), length(lons), length(lats))
    for i in eachindex(lons), j in eachindex(lats)
        all(m -> valid(pra[i, j, m]), 1:12) || continue
        MAP[i, j] = abs(sum(FT(pra[i, j, m]) for m in 1:12) / 12)
    end

    for (name, path) in PRODUCTS
        isfile(path) || continue
        olons, olats, O = obs_grid(name, path)
        acc = Dict(z => (FT[], FT[]) for z in ZONES)
        for i in eachindex(lons), j in eachindex(lats)
            valid(M[i, j]) || continue
            o = nearest(olons, olats, O, lons[i], lats[j])
            valid(o) || continue
            z = zone_of(lats[j], MAP[i, j])
            push!(acc[z][1], FT(M[i, j]))
            push!(acc[z][2], FT(o))
        end
        println("\n=== $name ===")
        println(
            rpad("zone", 12),
            rpad("cells", 8),
            rpad("model", 9),
            rpad("obs", 9),
            rpad("bias", 9),
            rpad("RMSE", 9),
            "r",
        )
        for z in ZONES
            m, o = acc[z]
            n = length(m)
            n < 20 && continue
            mb, ob = sum(m) / n, sum(o) / n
            sm = sqrt(sum((m .- mb) .^ 2) / n)
            so = sqrt(sum((o .- ob) .^ 2) / n)
            r = sum((m .- mb) .* (o .- ob)) / n / max(sm * so, eps(FT))
            rmse = sqrt(sum((m .- o) .^ 2) / n)
            println(
                rpad(ZONE_LABEL[z], 12),
                rpad(n, 8),
                rpad(round(mb, digits = 2), 9),
                rpad(round(ob, digits = 2), 9),
                rpad(round(mb - ob, digits = 2), 9),
                rpad(round(rmse, digits = 2), 9),
                round(r, digits = 3),
            )
        end
    end
end

isempty(ARGS) && error(
    "usage: julia score_by_zone.jl <equilibrium_carbon.nc> <driver_outdir>",
)
main(ARGS[1], ARGS[2])
