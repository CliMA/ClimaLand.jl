## Per-site equilibrium seeds, sampled from the global equilibrium map.
##
## `check_drift.jl` asks whether the coupled model holds the offline steady
## state, which needs each site started at its own equilibrium. Sampling the
## global map is the way to get that now the map carries the seasonality limit:
## a per-site offline spinup cannot, because the water deficit it needs is built
## from a gridded precipitation climatology rather than a single column.
##
## Column 7 is `D_annual`. Seeding it matters as much as the pools: it fills over
## a two-year memory from an aseasonal start, so a drift run that leaves it at
## zero under-suppresses woody allocation for most of its length and grows wood
## the equilibrium does not have - which looks like offline/coupled disagreement
## and is only spin-up.
##
## Usage:
##   julia --project=.buildkite seeds_from_map.jl <equilibrium_carbon.nc> <driver_outdir>

include(joinpath(@__DIR__, "global_equilibrium.jl"))

const SITES = joinpath(@__DIR__, "test_sites.csv")

read_sites(csv) = [
    (name = f[4], lon = parse(FT, f[1]), lat = parse(FT, f[2])) for
    f in (split(strip(l), ',') for l in readlines(csv)) if
    length(f) == 4 && f[1] != "lon" && !startswith(f[1], "#")
]

function nearest_land(lons, lats, A, lon, lat; radius = 2)
    l = lon
    maximum(lons) <= 180 && l > 180 && (l -= 360)
    i0 = argmin(abs.(lons .- l))
    j0 = argmin(abs.(lats .- lat))
    valid(A[i0, j0]) && return FT(A[i0, j0])
    best, bestd = FT(NaN), typemax(Int)
    for di in (-radius):radius, dj in (-radius):radius
        i, j = i0 + di, j0 + dj
        (1 <= i <= size(A, 1) && 1 <= j <= size(A, 2)) || continue
        valid(A[i, j]) || continue
        d = di * di + dj * dj
        d < bestd && (best, bestd = FT(A[i, j]), d)
    end
    return best
end

function main(ncpath, driverdir)
    sites = read_sites(SITES)
    ds = NCDatasets.NCDataset(ncpath)
    lons = Array{FT}(Array(ds["lon"][:]))
    lats = Array{FT}(Array(ds["lat"][:]))
    fld(n) = Array{Union{Missing, FT}}(Array(ds[n][:, :]))
    leaf, stem, root = fld("C_leaf"), fld("C_stem"), fld("C_root")
    cveg = fld("cVeg")
    close(ds)

    # The map stores no C_sugar; recover it as the residual of cVeg.
    _, _, P = read_monthly(driverdir, "precip")
    _, _, T = read_monthly(driverdir, "tair")
    toml_dict = LP.create_toml_dict(FT)
    pc = Canopy.PrognosticCarbonParameters(toml_dict)
    D = annual_deficit_grid(lons, lats, T, FT(0.6), P) ./ 1000  # mm -> m

    println("site\tbiome\tC_sugar\tC_leaf\tC_stem\tC_root\tD_annual")
    for s in sites
        l = nearest_land(lons, lats, leaf, s.lon, s.lat)
        st = nearest_land(lons, lats, stem, s.lon, s.lat)
        r = nearest_land(lons, lats, root, s.lon, s.lat)
        v = nearest_land(lons, lats, cveg, s.lon, s.lat)
        d = nearest_land(lons, lats, D, s.lon, s.lat)
        all(isfinite, (l, st, r, v)) || continue
        sugar = max(v - l - st - r, zero(FT))
        println(
            s.name,
            "\t-\t",
            round(sugar, digits = 5),
            "\t",
            round(l, digits = 5),
            "\t",
            round(st, digits = 5),
            "\t",
            round(r, digits = 5),
            "\t",
            isfinite(d) ? round(d, digits = 5) : "",
        )
    end
end

isempty(ARGS) && error(
    "usage: julia seeds_from_map.jl <equilibrium_carbon.nc> <driver_outdir>",
)
main(ARGS[1], ARGS[2])
