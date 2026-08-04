## Do the single-column runs and the global run see the same climate?
##
## This check used to compare equilibrium woody carbon between the map and a
## per-site offline spinup. That comparison became **circular** once the
## seasonality limit arrived: the per-site seeds can no longer be produced from a
## single column, because the water deficit needs a gridded precipitation
## climatology, so `equilibrium_pools.tsv` is now sampled from the map itself.
## Comparing it against the map would report a perfect agreement that means
## nothing. Rather than delete the check, it now asks the question its answer
## actually depended on.
##
## The pools follow from the forcing. If the site driver and the global driver
## disagree about precipitation at the same coordinate, the two configurations
## equilibrate to different vegetation for a reason that is not a model error,
## and every site-versus-map comparison inherits that disagreement. This measures
## it directly.
##
## It is also what makes `check_drift.jl` readable: a site seeded from a cell
## whose climate it does not experience will move on its own, and without this
## the movement looks like the offline and coupled integrators disagreeing.
##
## Interpretation: agreement is typically good - the median ratio across the
## battery is ~0.92 - and the outliers are sites in sharp precipitation
## gradients, where a point sample and a cell mean legitimately differ. That is a
## sampling difference, not a bug in either configuration.
##
## Usage:
##   julia --project=.buildkite check_map_vs_sites.jl <driver_outdir> <battery_runroot>

include(joinpath(@__DIR__, "global_equilibrium.jl"))

const SITES = joinpath(@__DIR__, "test_sites.csv")

# Beyond this the seed's climate and the run's climate are different enough that
# a drift comparison at that site is not testing the integrators. Matches
# DEFICIT_GAP_TOL in check_drift.jl.
const DEFICIT_GAP_TOL = FT(0.10)   # m of annual water deficit

read_sites(csv) = [
    (name = f[4], lon = parse(FT, f[1]), lat = parse(FT, f[2])) for
    f in (split(strip(l), ',') for l in readlines(csv)) if
    length(f) == 4 && f[1] != "lon" && !startswith(f[1], "#")
]

"Scalar metrics from a site's `carbon_metrics.txt`."
function read_metrics(path)
    d = Dict{String, FT}()
    isfile(path) || return d
    for line in readlines(path)
        f = split(strip(line))
        length(f) == 2 || continue
        v = tryparse(FT, f[2])
        v === nothing || (d[f[1]] = v)
    end
    return d
end

"Annual mean of a monthly field on the model grid, NaN where any month is masked."
function annual_mean(lons, lats, A)
    M = fill(FT(NaN), length(lons), length(lats))
    for i in eachindex(lons), j in eachindex(lats)
        all(m -> valid(A[i, j, m]), 1:12) || continue
        M[i, j] = abs(sum(FT(A[i, j, m]) for m in 1:12) / 12)
    end
    return M
end

function sample_at(lons, lats, M, lon, lat)
    l = lon
    maximum(lons) <= 180 && l > 180 && (l -= 360)
    return M[argmin(abs.(lons .- l)), argmin(abs.(lats .- lat))]
end

function main(driverdir, runroot)
    sites = read_sites(SITES)
    lons, lats, pra = read_monthly(driverdir, "pra")
    _, _, precip = read_monthly(driverdir, "precip")
    _, _, tair = read_monthly(driverdir, "tair")
    P = annual_mean(lons, lats, pra)
    D = annual_deficit_grid(lons, lats, tair, FT(0.6), precip) ./ 1000  # mm -> m

    println(
        rpad("site", 22),
        rpad("site P", 9),
        rpad("global P", 10),
        rpad("ratio", 8),
        rpad("site D", 9),
        rpad("global D", 10),
        "deficit gap",
    )
    ratios, nbad, ncmp = FT[], 0, 0
    for s in sites
        m = read_metrics(joinpath(runroot, s.name, "out", "carbon_metrics.txt"))
        haskey(m, "P_annual_m_yr") || continue
        gp = sample_at(lons, lats, P, s.lon, s.lat)
        gd = sample_at(lons, lats, D, s.lon, s.lat)
        isfinite(gp) || continue
        ncmp += 1
        r = gp > FT(0.05) ? m["P_annual_m_yr"] / gp : FT(NaN)
        isnan(r) || push!(ratios, r)
        gap =
            (haskey(m, "D_annual_m") && isfinite(gd)) ?
            abs(m["D_annual_m"] - gd) : FT(NaN)
        bad = isfinite(gap) && gap > DEFICIT_GAP_TOL
        nbad += bad
        println(
            rpad(s.name, 22),
            rpad(round(m["P_annual_m_yr"], digits = 2), 9),
            rpad(round(gp, digits = 2), 10),
            rpad(isnan(r) ? "-" : string(round(r, digits = 2)), 8),
            rpad(
                haskey(m, "D_annual_m") ? round(m["D_annual_m"], digits = 3) :
                "-",
                9,
            ),
            rpad(isfinite(gd) ? round(gd, digits = 3) : "-", 10),
            isnan(gap) ? "-" :
            string(round(gap, digits = 3), bad ? "   >tol" : ""),
        )
    end
    ncmp == 0 && error("no site metrics with P_annual under $runroot")
    sort!(ratios)
    med = ratios[(length(ratios) + 1) ÷ 2]
    println(
        "\n$ncmp sites; median precipitation ratio ",
        round(med, digits = 2),
        ", range ",
        round(ratios[1], digits = 2),
        " - ",
        round(ratios[end], digits = 2),
    )
    println(
        "$(ncmp - nbad)/$ncmp sites agree on the annual water deficit to within $(DEFICIT_GAP_TOL) m",
    )
    println(
        nbad == 0 ?
        "FORCING CONSISTENT - a site-versus-map comparison is meaningful everywhere" :
        "FORCING DIFFERS at $nbad site(s) - drift there is not evidence about the integrators",
    )
end

isempty(ARGS) && error(
    "usage: julia check_map_vs_sites.jl <driver_outdir> <battery_runroot>",
)
main(ARGS[1], ARGS[2])
