## Does the coupled model hold the offline steady state? (stage 5 verification)
##
## Every equilibrium in the global map came from the offline integrator. The
## coupled model has never actually been asked to carry those pools, so the
## whole map rests on an assumption: that `step_pools` reproduces the fluxes
## `update_carbon_fluxes!` computes inside ClimaLand. This tests it directly.
##
## A battery run seeded at each site's own 400-year offline equilibrium should
## sit still. Drift means the two integrations disagree, and the map is wrong by
## whatever that difference compounds to over centuries.
##
## One caveat this check cannot see on its own. Since the seeds are sampled from
## the global equilibrium map, the seed carries the *grid cell's* climate while the
## run experiences the *site's*. Those usually agree - median precipitation ratio
## 0.92 over 15 sites - but in a sharp precipitation gradient they do not, and a
## site seeded from a drier cell than it actually experiences will grow wood and
## report large drift that has nothing to do with the integrators. Check
## `P_annual_m_yr` and `D_annual_m` in the site metrics against the seed before
## reading a large drift as offline/coupled disagreement.
##
## What counts as still: the run is 2 years against stem turnover of 30 years and
## more, so even a perfectly consistent pair drifts a little - the seeded state
## is the equilibrium of the *recycled climatology*, while the coupled run sees
## two particular years of weather. A few percent over two years is that; tens of
## percent is a real inconsistency.
##
## Usage:  julia --project=.buildkite check_drift.jl <battery_runroot>

const FT = Float64
const SEEDS = joinpath(@__DIR__, "equilibrium_pools.tsv")

# A site is only a fair test of the integrators if the climate it runs under
# matches the climate its seed came from. The seeds are sampled from the global
# equilibrium map, so they carry the grid cell's water deficit while the run
# experiences the site's; in a sharp precipitation gradient those differ enough
# to move the equilibrium on their own. Sites beyond this gap are reported
# separately rather than counted as integrator disagreement.
const DEFICIT_GAP_TOL = FT(0.10)   # m of annual water deficit

function read_seeds(path)
    out = Dict{String, NTuple{5, FT}}()
    for l in readlines(path)
        (isempty(strip(l)) || startswith(l, '#') || startswith(l, "site")) &&
            continue
        f = split(strip(l), '\t')
        length(f) >= 6 && (
            out[f[1]] = (
                parse(FT, f[3]),
                parse(FT, f[4]),
                parse(FT, f[5]),
                parse(FT, f[6]),
                # Column 7 is D_annual, absent from older TSVs.
                length(f) >= 7 && !isempty(f[7]) ? parse(FT, f[7]) : FT(NaN),
            )
        )
    end
    return out
end

function read_metrics(path)
    d = Dict{String, FT}()
    for l in readlines(path)
        f = split(strip(l))
        length(f) >= 2 || continue
        v = tryparse(FT, f[2])
        v === nothing || (d[f[1]] = v)
    end
    return d
end

function main(runroot)
    seeds = read_seeds(SEEDS)
    println(
        rpad("site", 21),
        rpad("C_stem seed", 13),
        rpad("C_stem final", 14),
        rpad("drift", 10),
        rpad("cVeg drift", 12),
        "deficit gap",
    )
    worst, worst_site = 0.0, ""
    n, nbad = 0, 0
    for (name, s) in sort(collect(seeds); by = first)
        m = joinpath(runroot, name, "out", "carbon_metrics.txt")
        isfile(m) || continue
        d = read_metrics(m)
        haskey(d, "C_stem_kgC_m2") || continue
        stem0, stem1 = s[3], d["C_stem_kgC_m2"]
        veg0 = sum(s[1:4])
        veg1 = get(d, "cVeg_kgC_m2", NaN)
        # Sites with no vegetation at all are exactly zero on both sides; a
        # relative drift there is 0/0, not an error.
        rel = stem0 > 0.05 ? (stem1 - stem0) / stem0 : 0.0
        vrel = veg0 > 0.05 ? (veg1 - veg0) / veg0 : 0.0
        n += 1
        bad = abs(rel) > 0.10
        nbad += bad
        if abs(rel) > abs(worst)
            worst, worst_site = rel, name
        end
        println(
            rpad(name, 21),
            rpad(round(stem0, digits = 3), 13),
            rpad(round(stem1, digits = 3), 14),
            rpad(string(round(100 * rel, digits = 1), "%"), 10),
            rpad(string(round(100 * vrel, digits = 1), "%"), 12),
            rpad(
                haskey(d, "D_annual_m") && length(s) >= 5 ?
                string(round(abs(d["D_annual_m"] - s[5]), digits = 3)) : "-",
                9,
            ),
            bad ? "  >10%" : "",
        )
    end
    n == 0 && error("no carbon_metrics.txt under $runroot")
    println(
        "\n$(n - nbad)/$n sites hold their offline equilibrium to within 10% over 2 years",
    )
    println("worst: $worst_site at $(round(100 * worst, digits = 1))%")
    println(
        nbad == 0 ?
        "DRIFT OK - the offline integrator reproduces the coupled model, so the global map stands" :
        "DRIFT SUSPECT at $nbad site(s) - check the deficit gap column before \
concluding the integrations disagree: a site whose run sees a different climate \
from its seed will move for that reason alone",
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    isempty(ARGS) && error("usage: julia check_drift.jl <battery_runroot>")
    main(ARGS[1])
end
