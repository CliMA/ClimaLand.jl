## Stage-2 acceptance check: is the pool-based autotrophic respiration sane?
##
## The gate is: annual NPP/GPP in 0.3-0.6 at the forest and grassland sites, and
## Ra neither collapsing to zero nor exceeding GPP in the annual mean anywhere.
##
## Also tests the falsifiable prediction recorded before the run: at a site whose
## pools are all zero (the true deserts), pool-based Ra must be ~0. The JULES
## term cannot do this - it respires prescribed SAI and RAI, which do not vanish
## over bare ground - so this is the single number that shows the replacement
## did what it was for.
##
## Usage:
##   julia check_stage2.jl <stage2_summary.tsv> [<baseline_summary.tsv>]
##
## Exits non-zero if the acceptance criteria are not met.

const BAND = (0.3, 0.6)
# Biomes the gate names explicitly. Deserts, tundra and savanna are reported but
# not gated: the band is a forest/grassland expectation.
const GATED = (
    "tropical_rainforest",
    "temperate_deciduous_forest",
    "temperate_grassland",
    "temperate_grassland_C4",
    "boreal_forest",
)

function read_summary(path)
    rows = Dict{String, Dict{String, Any}}()
    order = String[]
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
    base =
        isnothing(baseline_path) ? nothing : first(read_summary(baseline_path))

    failures = String[]
    println(
        rpad("site", 21),
        rpad("biome", 27),
        rpad("GPP", 9),
        rpad("Ra", 9),
        rpad("NPP/GPP", 9),
        isnothing(base) ? "" : rpad("Ra_base", 9),
        "verdict",
    )
    for site in order
        row = rows[site]
        gpp, ra, biome =
            num(row, "GPP_gC_m2_day"), num(row, "Ra_gC_m2_day"), row["biome"]
        frac = gpp > 0 ? 1 - ra / gpp : NaN
        gated = biome in GATED

        verdict = if gpp == 0
            # No photosynthesis at all: Ra must be ~0, or the model is
            # respiring carbon it never fixed. This is the prediction.
            # 1e-3 g C m^-2 day^-1 is 0.36 g C m^-2 yr^-1: physically nil, and
            # ~800x below the 0.813 the JULES term produced at these sites. A
            # tighter threshold flags decaying seeded pools, which is noise.
            ra <= 1e-3 ? "PREDICTION OK (Ra~0 with GPP=0)" :
            (
                push!(failures, "$site: Ra=$ra with GPP=0"); "FAIL Ra>0 with GPP=0"
            )
        elseif ra > gpp
            push!(failures, "$site: Ra ($ra) exceeds GPP ($gpp)")
            "FAIL Ra>GPP"
        elseif ra <= 0
            push!(failures, "$site: Ra collapsed to $ra")
            "FAIL Ra<=0"
        elseif gated && !(BAND[1] <= frac <= BAND[2])
            push!(
                failures,
                "$site ($biome): NPP/GPP=$(round(frac, digits=3)) outside $BAND",
            )
            "FAIL band"
        elseif gated
            "ok"
        else
            "ok (not gated)"
        end

        ra_base =
            isnothing(base) || !haskey(base, site) ? "" :
            rpad(round(num(base[site], "Ra_gC_m2_day"), digits = 3), 9)
        println(
            rpad(site, 21),
            rpad(biome, 27),
            rpad(round(gpp, digits = 3), 9),
            rpad(round(ra, digits = 3), 9),
            rpad(isnan(frac) ? "n/a" : round(frac, digits = 3), 9),
            ra_base,
            verdict,
        )
    end

    println()
    if isempty(failures)
        println(
            "STAGE2 PASS - Ra sane at every site; gated biomes inside $BAND",
        )
        return 0
    else
        println("STAGE2 FAIL - $(length(failures)) issue(s):")
        foreach(f -> println("  ", f), failures)
        return 1
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    isempty(ARGS) && error(
        "usage: julia check_stage2.jl <stage2_summary.tsv> [<baseline.tsv>]",
    )
    exit(main(ARGS[1], length(ARGS) > 1 ? ARGS[2] : nothing))
end
