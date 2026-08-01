## Are monthly-mean drivers good enough to spin the pools up on? (stage 5)
##
## The 20-site battery records drivers daily. A *global* equilibrium map cannot
## afford daily output at every land column, and the global runs that already
## exist write monthly means. Before building the global evaluation on monthly
## drivers, this establishes whether that coarsening changes the answer.
##
## It should not be assumed either way. `Rm` carries a Q10 and a saturating
## substrate ramp, and allocation is nonlinear in sugar, so by Jensen's
## inequality the equilibrium of the averaged drivers is *not* the average of
## the equilibria. The question is only whether the bias is small compared with
## the 3.4x spread among the observational products it will be judged against.
##
## Usage:  julia --project=.buildkite check_monthly_drivers.jl <battery_runroot>

include(joinpath(@__DIR__, "offline_spinup.jl"))

import ClimaLand.Parameters as LP
import ClimaLand.Canopy

read_sites(csv) = [
    (name = f[4], biome = f[3]) for
    f in (split(strip(l), ',') for l in readlines(csv)) if
    length(f) == 4 && f[1] != "lon" && !startswith(f[1], "#")
]

"""
    to_monthly(d)

Collapse a daily driver record to 12 monthly means, then re-expand it to a
365-day record holding each month's mean over its days. Re-expanding rather than
integrating on a monthly timestep keeps the integrator identical between the two
arms, so the only thing under test is the loss of sub-monthly variability.
"""
function to_monthly(d)
    n = length(first(values(d)))
    # Day-of-year to month, on a no-leap calendar. The record is whole years of
    # recycled forcing, so a fixed calendar is exact rather than approximate.
    lens = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    month_of = Int[]
    for (m, L) in enumerate(lens), _ in 1:L
        push!(month_of, m)
    end
    out = Dict{String, Vector{FT}}()
    for (c, v) in d
        sums = zeros(FT, 12)
        cnts = zeros(Int, 12)
        for j in 1:n
            m = month_of[mod1(j, 365)]
            sums[m] += v[j]
            cnts[m] += 1
        end
        means = [cnts[m] > 0 ? sums[m] / cnts[m] : FT(0) for m in 1:12]
        out[c] = [means[month_of[mod1(j, 365)]] for j in 1:n]
    end
    return out
end

function main(runroot)
    sites = read_sites(joinpath(@__DIR__, "test_sites.csv"))
    toml_dict = LP.create_toml_dict(FT)
    p = Canopy.PrognosticCarbonParameters(toml_dict)

    println(
        rpad("site", 21),
        rpad("biome", 27),
        rpad("C_stem daily", 14),
        rpad("C_stem monthly", 16),
        "rel. diff",
    )
    diffs = FT[]
    for s in sites
        f = joinpath(runroot, s.name, "out", "driver_record.csv")
        isfile(f) || continue
        d = read_driver_record(f)
        kw = (; years = 400, T_ref_tau = p.T_ref_τ_stem, q_tau = p.q_τ_stem)
        daily, _ = spinup(d, p; kw...)
        monthly, _ = spinup(to_monthly(d), p; kw...)
        # Relative to the daily answer, which is the reference by construction.
        rel =
            daily[3] > FT(0.01) ? (monthly[3] - daily[3]) / daily[3] :
            monthly[3] - daily[3]
        push!(diffs, abs(rel))
        println(
            rpad(s.name, 21),
            rpad(s.biome, 27),
            rpad(round(daily[3], digits = 3), 14),
            rpad(round(monthly[3], digits = 3), 16),
            string(round(100 * rel, digits = 1), "%"),
        )
    end

    isempty(diffs) && error("no driver records under $runroot")
    sort!(diffs)
    println(
        "\nmedian |rel diff| = ",
        round(100 * diffs[max(1, div(length(diffs), 2))], digits = 1),
        "%,  max = ",
        round(100 * maximum(diffs), digits = 1),
        "%",
    )
    println(
        maximum(diffs) < 0.15 ?
        "MONTHLY OK - the coarsening is small next to the 3.4x product spread" :
        "MONTHLY NOT OK - the global evaluation needs finer driver output",
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    isempty(ARGS) &&
        error("usage: julia check_monthly_drivers.jl <battery_runroot>")
    main(ARGS[1])
end
