## Does the global equilibrium map reproduce the single-column results?
##
## The map and the battery reach the same quantity by different routes: the
## battery integrates daily drivers from a column run at one longitude/latitude,
## the map integrates a monthly climatology sampled from a global run at the same
## point. Agreement is not guaranteed by construction, so it is worth checking -
## it exercises the global driver output, the dimension-order handling, the
## land mask, and the climatology averaging all at once, against numbers that
## were produced before any of that code existed.
##
## Disagreement is informative rather than fatal. Two known sources:
##   - the global run samples a 1 degree cell, the battery a single column, and
##     coastal or mountainous sites differ legitimately;
##   - monthly vs daily drivers, worth ~0.4% median (check_monthly_drivers.jl).
## A site that disagrees by more than a factor of two is a bug, not a sampling
## difference, and that is the threshold this reports on.
##
## Usage:
##   julia --project=.buildkite check_map_vs_sites.jl <equilibrium_carbon.nc>

import NCDatasets

const FT = Float64
const TSV = joinpath(@__DIR__, "equilibrium_pools.tsv")
const SITES = joinpath(@__DIR__, "test_sites.csv")

read_sites(csv) = [
    (name = f[4], lon = parse(FT, f[1]), lat = parse(FT, f[2])) for
    f in (split(strip(l), ',') for l in readlines(csv)) if
    length(f) == 4 && f[1] != "lon" && !startswith(f[1], "#")
]

function read_pools(path)
    out = Dict{String, FT}()
    for l in readlines(path)
        (isempty(strip(l)) || startswith(l, '#') || startswith(l, "site")) &&
            continue
        f = split(strip(l), '\t')
        length(f) >= 5 && (out[f[1]] = parse(FT, f[5]))  # C_stem
    end
    return out
end

"""
    nearest_land(lons, lats, A, lon, lat; radius = 2)

Value at the nearest cell, falling back to the nearest finite cell within
`radius` cells. Battery sites sit at round coordinates that can land on a
coastal cell the global land mask excludes; searching a small neighbourhood
keeps such a site comparable instead of silently dropping it.
"""
function nearest_land(lons, lats, A, lon, lat; radius = 2)
    l = lon
    maximum(lons) <= 180 && l > 180 && (l -= 360)
    i0 = argmin(abs.(lons .- l))
    j0 = argmin(abs.(lats .- lat))
    isfinite(A[i0, j0]) && return A[i0, j0], 0
    best, bestd = FT(NaN), typemax(Int)
    for di in (-radius):radius, dj in (-radius):radius
        i, j = i0 + di, j0 + dj
        (1 <= i <= size(A, 1) && 1 <= j <= size(A, 2)) || continue
        isfinite(A[i, j]) || continue
        d = di * di + dj * dj
        d < bestd && (best, bestd = A[i, j], d)
    end
    return best, isfinite(best) ? round(Int, sqrt(bestd)) : -1
end

function main(ncpath)
    sites = read_sites(SITES)
    site_stem = read_pools(TSV)
    ds = NCDatasets.NCDataset(ncpath)
    lons = Array{FT}(Array(ds["lon"][:]))
    lats = Array{FT}(Array(ds["lat"][:]))
    dims = collect(NCDatasets.dimnames(ds["C_stem"]))
    A = Array(ds["C_stem"][:, :])
    # written as (lon, lat) by global_equilibrium.jl, but check rather than trust
    if lowercase(dims[1]) != "lon"
        A = permutedims(A, (2, 1))
    end
    close(ds)

    println(
        rpad("site", 21),
        rpad("battery", 10),
        rpad("map", 10),
        rpad("ratio", 9),
        "note",
    )
    nbad, ncmp = 0, 0
    for s in sites
        haskey(site_stem, s.name) || continue
        b = site_stem[s.name]
        m, shift = nearest_land(lons, lats, A, s.lon, s.lat)
        if !isfinite(m)
            println(
                rpad(s.name, 21),
                rpad(round(b, digits = 2), 10),
                "no land cell within 2 cells",
            )
            continue
        end
        ncmp += 1
        # Both near zero is agreement, not a 0/0 blow-up.
        ratio = (b < 0.05 && m < 0.05) ? 1.0 : m / max(b, 0.01)
        bad = ratio > 2 || ratio < 0.5
        nbad += bad
        println(
            rpad(s.name, 21),
            rpad(round(b, digits = 2), 10),
            rpad(round(m, digits = 2), 10),
            rpad(round(ratio, digits = 2), 9),
            bad ? "OUTSIDE 2x" :
            (shift > 0 ? "ok (shifted $shift cells)" : "ok"),
        )
    end
    println("\n$(ncmp - nbad)/$ncmp sites agree within a factor of two")
    println(
        nbad == 0 ? "MAP CONSISTENT with the single-column battery" :
        "MAP DISAGREES at $nbad site(s) - investigate before using the map",
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    isempty(ARGS) &&
        error("usage: julia check_map_vs_sites.jl <equilibrium_carbon.nc>")
    main(ARGS[1])
end
