## Is there ANY constant stem allocation that works? (stage 5)
##
## Stage 4 found forests under-built and grasslands over-built simultaneously,
## which is the symptom MODEL.md §2.3 names as the trigger for letting `f_stem`
## vary with climate. Before making that change, this establishes the stronger
## claim the fallback actually rests on: that no single `(f_stem, τ_stem)` pair
## satisfies both, so the failure is structural rather than a bad choice of
## constant.
##
## The offline integrator makes this cheap - each site equilibrates in well under
## a second, so a full grid over 20 sites costs seconds rather than the
## node-hours a coupled sweep would.
##
## Usage:
##   julia sweep_allocation.jl <battery_runroot> [--years N]

include(joinpath(@__DIR__, "offline_spinup.jl"))

import ClimaLand.Parameters as LP
import ClimaLand.Canopy

# MODEL.md §7 targets, kg C m^-2. Only the classes the specification actually
# names are gated; the rest are reported but not scored, because inventing a
# target would let the sweep "pass" on a number nobody chose.
const TARGETS = Dict(
    "tropical_rainforest" => (10.0, 20.0),
    "temperate_deciduous_forest" => (5.0, 15.0),
    "boreal_forest" => (5.0, 15.0),
    "temperate_grassland" => (0.3, 1.0),
    "temperate_grassland_C4" => (0.3, 1.0),
    "desert" => (0.0, 1.0),
)

read_sites(csv) = [
    (name = f[4], biome = f[3]) for
    f in (split(strip(l), ',') for l in readlines(csv)) if
    length(f) == 4 && f[1] != "lon" && !startswith(f[1], "#")
]

"""
    equilibrium_cVeg(records, params; years)

Spins up every site offline and returns its equilibrium cVeg.
"""
function equilibrium_cVeg(records, p; years = 400, map_half = 0.0, map_n = 2.0)
    out = Dict{String, Float64}()
    stem = Dict{String, Float64}()
    for (name, d) in records
        pools, _ = spinup(d, p; years, map_half, map_n)
        out[name] = sum(pools)
        stem[name] = pools[3]
    end
    return out, stem
end

"""
    score(cveg, sites)

Fraction of gated sites whose equilibrium cVeg falls inside its biome target,
reported separately for woody and herbaceous classes. Keeping them apart is the
point: a single number would hide the conflict this sweep exists to expose.
"""
function score(cveg, sites)
    woody_ok, woody_n, herb_ok, herb_n = 0, 0, 0, 0
    for s in sites
        haskey(TARGETS, s.biome) || continue
        haskey(cveg, s.name) || continue
        lo, hi = TARGETS[s.biome]
        ok = lo <= cveg[s.name] <= hi
        if occursin("grassland", s.biome) || s.biome == "desert"
            herb_n += 1
            herb_ok += ok
        else
            woody_n += 1
            woody_ok += ok
        end
    end
    return (; woody_ok, woody_n, herb_ok, herb_n)
end

function main(runroot; years = 400)
    sites = read_sites(joinpath(@__DIR__, "test_sites.csv"))
    records = Tuple{String, Any}[]
    for s in sites
        f = joinpath(runroot, s.name, "out", "driver_record.csv")
        isfile(f) && push!(records, (s.name, read_driver_record(f)))
    end
    isempty(records) && error("no driver records under $runroot")
    @info "loaded $(length(records)) driver records"

    toml_dict = LP.create_toml_dict(Float64)
    println(
        rpad("map_half", 11),
        rpad("f_stem_c3", 13),
        rpad("woody in band", 15),
        rpad("herb in band", 14),
        "C_stem at key sites",
    )
    # map_half = 0 is the constant-f_stem control, so the mechanism is always
    # compared against the behaviour it is meant to replace.
    for map_half in (0.0, 0.4, 0.7, 1.0, 1.5)
        for f_stem in (0.4, 0.5)
            p = Canopy.PrognosticCarbonParameters(toml_dict; f_stem_c3 = f_stem)
            cveg, stem = equilibrium_cVeg(records, p; years, map_half)
            sc = score(cveg, sites)
            println(
                rpad(map_half == 0 ? "off" : map_half, 11),
                rpad(f_stem, 13),
                rpad("$(sc.woody_ok)/$(sc.woody_n)", 15),
                rpad("$(sc.herb_ok)/$(sc.herb_n)", 14),
                "stem gp/pampas/siberia ",
                round(get(stem, "us_great_plains", NaN), digits = 2),
                " / ",
                round(get(stem, "pampas_argentina", NaN), digits = 2),
                " / ",
                round(get(stem, "central_siberia", NaN), digits = 2),
            )
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    isempty(ARGS) && error("usage: julia sweep_allocation.jl <battery_runroot>")
    yi = findfirst(==("--years"), ARGS)
    main(ARGS[1]; years = yi === nothing ? 400 : parse(Float64, ARGS[yi + 1]))
end
