## Rule-1 check: the carbon model must not change GPP or LAI.
##
## Diffs a CARBON=1 battery summary against the committed CARBON=0 baseline,
## site by site. Phase 1 is one-way coupled, so any difference beyond round-off
## is a bug in the carbon model by definition, not a tuning question.
##
## Usage:
##   julia check_rule1.jl <carbon_summary.tsv> <baseline_summary.tsv>
##
## Exits non-zero if any site fails, so a battery script can gate on it.

const REL_TOL = 1e-10  # round-off; phase 1 should reproduce the baseline exactly

function read_summary(path)
    rows = Dict{String, Dict{String, Any}}()
    open(path) do io
        header = split(strip(readline(io)), '\t')
        for line in eachline(io)
            isempty(strip(line)) && continue
            fields = split(strip(line), '\t')
            length(fields) < 2 && continue
            row = Dict{String, Any}()
            for (k, v) in zip(header, fields)
                row[k] = something(tryparse(Float64, v), v)
            end
            rows[fields[1]] = row
        end
    end
    return rows
end

function main(carbon_path, baseline_path)
    carbon = read_summary(carbon_path)
    baseline = read_summary(baseline_path)

    sites = sort(collect(intersect(keys(carbon), keys(baseline))))
    isempty(sites) &&
        error("no sites in common between $carbon_path and $baseline_path")

    missing_sites = setdiff(keys(baseline), keys(carbon))
    if !isempty(missing_sites)
        println(
            "NOTE: $(length(missing_sites)) baseline site(s) absent from the ",
            "carbon run, not checked: ",
            join(sort(collect(missing_sites)), ", "),
        )
    end

    failures = String[]
    println(
        rpad("site", 22),
        rpad("field", 16),
        rpad("baseline", 22),
        rpad("carbon", 22),
        "rel_diff",
    )
    for site in sites
        for field in ("GPP_gC_m2_day", "LAI")
            b = get(baseline[site], field, nothing)
            c = get(carbon[site], field, nothing)
            (b isa Float64 && c isa Float64) || continue
            scale = max(abs(b), abs(c))
            rel = scale == 0 ? abs(c - b) : abs(c - b) / scale
            flag = rel <= REL_TOL ? "ok" : "FAIL"
            rel > REL_TOL && push!(failures, "$site/$field rel=$rel")
            println(
                rpad(site, 22),
                rpad(field, 16),
                rpad(b, 22),
                rpad(c, 22),
                rel,
                "  ",
                flag,
            )
        end
    end

    println()
    if isempty(failures)
        println(
            "RULE1 PASS - GPP and LAI unchanged at all $(length(sites)) sites ",
            "(rel tol $REL_TOL)",
        )
        return 0
    else
        println("RULE1 FAIL - $(length(failures)) violation(s):")
        foreach(f -> println("  ", f), failures)
        return 1
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 2 || error(
        "usage: julia check_rule1.jl <carbon_summary.tsv> <baseline_summary.tsv>",
    )
    exit(main(ARGS[1], ARGS[2]))
end
