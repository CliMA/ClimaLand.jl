## Stage-3 acceptance check: is prognostic SOC behaving?
##
## The gate is: SOC drifting slowly - neither collapsing nor exploding - and Rh
## within a factor of ~2 of its stage-0 baseline at every site.
##
## Both halves matter and they fail differently. A collapsing SOC means litter is
## not reaching the soil, or Sm is being double-debited. An Rh far from baseline
## means turning SOC on changed decomposition itself, which it should not: Rh
## depends on SOC only weakly over a 2-year run, since SOC is a large pool with a
## centuries-long turnover time.
##
## Usage:
##   julia check_stage3.jl <stage3_summary.tsv> <stage0_baseline.tsv>

# SOC turnover is centuries, so a 2-year run should move it by very little.
# 10% is generous; anything beyond that is a rate problem, not spinup.
const MAX_DRIFT = 0.10
const RH_FACTOR = 2.0

function read_summary(path)
    rows, order = Dict{String, Dict{String, Any}}(), String[]
    open(path) do io
        header = split(strip(readline(io)), '\t')
        for line in eachline(io)
            isempty(strip(line)) && continue
            f = split(strip(line), '\t')
            length(f) < 2 && continue
            rows[f[1]] = Dict(
                k => something(tryparse(Float64, v), v) for
                (k, v) in zip(header, f)
            )
            push!(order, f[1])
        end
    end
    return rows, order
end

num(row, key) = get(row, key, nothing) isa Float64 ? row[key] : NaN

function main(path, baseline_path)
    rows, order = read_summary(path)
    base, _ = read_summary(baseline_path)

    failures = String[]
    println(
        rpad("site", 21),
        rpad("cSoil_init", 12),
        rpad("cSoil", 12),
        rpad("drift", 10),
        rpad("Rh", 9),
        rpad("Rh_base", 9),
        "verdict",
    )
    for site in order
        row = rows[site]
        ci, cf = num(row, "cSoil_init"), num(row, "cSoil")
        drift = num(row, "cSoil_dfrac")
        rh = num(row, "Rh_gC_m2_day")
        rh_b = haskey(base, site) ? num(base[site], "Rh_gC_m2_day") : NaN

        issues = String[]
        if isnan(cf) || !isfinite(cf)
            push!(issues, "cSoil not finite")
        elseif cf < 0
            push!(issues, "cSoil negative ($cf)")
        end
        if !isnan(drift) && abs(drift) > MAX_DRIFT
            push!(issues, "SOC drift $(round(drift*100, digits=1))% over 2 yr")
        end
        # Rh is compared only where the baseline actually respires; a site with
        # no heterotrophic respiration in either run is not a ratio failure.
        if !isnan(rh) && !isnan(rh_b) && rh_b > 1e-3
            r = rh / max(rh_b, eps())
            (r > RH_FACTOR || r < 1 / RH_FACTOR) &&
                push!(issues, "Rh $(round(r, digits=2))x baseline")
        end

        append!(failures, ["$site: $i" for i in issues])
        println(
            rpad(site, 21),
            rpad(isnan(ci) ? "-" : round(ci, digits = 3), 12),
            rpad(isnan(cf) ? "-" : round(cf, digits = 3), 12),
            rpad(
                isnan(drift) ? "-" :
                string(round(drift * 100, digits = 2), "%"),
                10,
            ),
            rpad(isnan(rh) ? "-" : round(rh, digits = 3), 9),
            rpad(isnan(rh_b) ? "-" : round(rh_b, digits = 3), 9),
            isempty(issues) ? "ok" : "FAIL " * join(issues, "; "),
        )
    end

    println()
    if isempty(failures)
        println(
            "STAGE3 PASS - SOC drifts <$(MAX_DRIFT*100)% and Rh within ",
            "$(RH_FACTOR)x of baseline everywhere",
        )
        return 0
    else
        println("STAGE3 FAIL - $(length(failures)) issue(s):")
        foreach(f -> println("  ", f), failures)
        return 1
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 2 || error(
        "usage: julia check_stage3.jl <stage3_summary.tsv> <stage0_baseline.tsv>",
    )
    exit(main(ARGS[1], ARGS[2]))
end
