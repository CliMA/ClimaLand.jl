## Can rainfall *seasonality* separate wet savanna from wet forest?
##
## The binned diagnostic (`bin_by_obs.jl`) localised the model's overprediction:
## it is not spread across treeless land, it sits in the ~6% of treeless cells
## whose mean annual precipitation is above the woody-fraction half-point. Those
## cells - Cerrado, the Sahel fringe, the Pampas, miombo - receive forest-sized
## rainfall, so no function of MAP can ever remove wood from them. Five climate
## predictors have now failed on exactly them.
##
## But MAP is an annual total, and the classical discriminator of the tropical
## forest/savanna boundary is not the total: it is how the rain is distributed
## through the year. A cell receiving 1500 mm evenly supports closed forest; the
## same 1500 mm delivered in five months, followed by a seven-month drought that
## dries the fuel and carries fire, supports savanna. Dry-season length and
## maximum cumulative water deficit are the standard measures.
##
## This asks whether that separation is actually present in the data, *before*
## proposing any model change. The test is a contingency table restricted to the
## cells the MAP ramp cannot reach:
##
##   among wet cells (MAP >= map_half), do the ones the observations say are
##   treeless have a longer dry season than the ones they say are forest?
##
## If the two distributions overlap, seasonality is a sixth failed predictor and
## the residual is not reachable by any function of precipitation at all. If they
## separate, it justifies making the woody fraction depend on seasonality as well
## as amount - still a global constant, still a climate driver, still no PFT.
##
## Precipitation is taken from GPCC rather than from the model so that the test
## is of the real world, not of the forcing. If it succeeds, the model-side
## version must be recomputed from ERA5.
##
## Usage:
##   julia --project=.buildkite check_seasonality.jl <equilibrium_carbon.nc> <driver_outdir>

include(joinpath(@__DIR__, "global_equilibrium.jl"))

# PR_PATH, DRY_MM_PER_MONTH, pr_climatology and annual_deficit come from
# global_equilibrium.jl, which also applies the predictor this script evaluates.

"""
    mcwd(p)

Maximum cumulative water deficit (mm, positive). Runs the monthly water balance
`P - 100` around the year twice so the accumulation starts from the wet season
whatever month the record begins in, which is what makes it comparable between
hemispheres.
"""
function mcwd(p)
    w, worst = FT(0), FT(0)
    for pass in 1:2, m in 1:12
        w = min(w + p[m] - DRY_MM_PER_MONTH, FT(0))
        pass == 2 && (worst = min(worst, w))
    end
    return -worst
end

quantiles(v, qs) =
    isempty(v) ? fill(NaN, length(qs)) :
    (
        s = sort(v); [
            s[clamp(round(Int, q * length(s)), 1, length(s))] for q in qs
        ]
    )

function seasonality_main(ncpath, driverdir)
    ds = NCDatasets.NCDataset(ncpath)
    lons = Array{FT}(Array(ds["lon"][:]))
    lats = Array{FT}(Array(ds["lat"][:]))
    cVeg = Array{Union{Missing, FT}}(Array(ds["cVeg"][:, :]))
    close(ds)

    _, _, pra = read_monthly(driverdir, "pra")
    MAP = fill(FT(NaN), size(pra, 1), size(pra, 2))
    for i in axes(pra, 1), j in axes(pra, 2)
        s, n = FT(0), 0
        for m in 1:12
            valid(pra[i, j, m]) && (s += FT(pra[i, j, m]); n += 1)
        end
        n == 12 && (MAP[i, j] = abs(s / 12))
    end

    plons, plats, P = pr_climatology(PR_PATH)
    toml_dict = LP.create_toml_dict(FT)
    half = Canopy.PrognosticCarbonParameters(toml_dict).map_half_woody
    @info "wet-cell threshold" map_half = half

    for (name, path) in PRODUCTS
        isfile(path) || continue
        olons, olats, O = obs_grid(name, path)
        # Two classes among the cells the MAP ramp cannot reach.
        tree_d, tree_w, tree_m, tree_a = FT[], FT[], FT[], FT[]
        grass_d, grass_w, grass_m, grass_a = FT[], FT[], FT[], FT[]
        for i in eachindex(lons), j in eachindex(lats)
            valid(cVeg[i, j]) || continue
            valid(MAP[i, j]) && MAP[i, j] >= half || continue
            o = nearest(olons, olats, O, lons[i], lats[j])
            valid(o) || continue
            p = [
                nearest(plons, plats, view(P, :, :, m), lons[i], lats[j])
                for m in 1:12
            ]
            all(valid, p) || continue
            pv = FT.(p)
            d, w, mm = FT(dry_months(pv)), mcwd(pv), FT(cVeg[i, j])
            a = annual_deficit(pv)
            if o < 0.5
                push!(grass_d, d)
                push!(grass_w, w)
                push!(grass_m, mm)
                push!(grass_a, a)
            elseif o > 5.0
                push!(tree_d, d)
                push!(tree_w, w)
                push!(tree_m, mm)
                push!(tree_a, a)
            end
        end
        (isempty(grass_d) || isempty(tree_d)) && continue
        qd_g, qd_t = quantiles(grass_d, [0.25, 0.5, 0.75]),
        quantiles(tree_d, [0.25, 0.5, 0.75])
        qw_g, qw_t = quantiles(grass_w, [0.25, 0.5, 0.75]),
        quantiles(tree_w, [0.25, 0.5, 0.75])
        qa_g, qa_t = quantiles(grass_a, [0.25, 0.5, 0.75]),
        quantiles(tree_a, [0.25, 0.5, 0.75])
        println("\n=== $name : wet cells only (MAP >= $half m/yr) ===")
        println(
            rpad("class", 20),
            rpad("cells", 8),
            rpad("dry months q25/50/75", 24),
            rpad("MCWD mm q25/50/75", 24),
            rpad("annual deficit mm", 22),
            "model cVeg",
        )
        for (lab, n, qd, qw, qa, mv) in (
            ("observed treeless", length(grass_d), qd_g, qw_g, qa_g, grass_m),
            ("observed forest", length(tree_d), qd_t, qw_t, qa_t, tree_m),
        )
            println(
                rpad(lab, 20),
                rpad(n, 8),
                rpad(join(round.(qd, digits = 1), "/"), 24),
                rpad(join(round.(Int, qw), "/"), 24),
                rpad(join(round.(Int, qa), "/"), 22),
                round(sum(mv) / length(mv), digits = 2),
            )
        end
        # A threshold is only useful if one class sits mostly on each side of it.
        # Reported for the median of the treeless class, which is the best case
        # a single global cut can do.
        # Skill of the best single global cut on each candidate, scored the same
        # way for all three so they are comparable: the cut is placed at the
        # treeless median and read as "how much of each class falls above it".
        # A useless predictor puts the same share of both classes above the cut.
        for (lab, g, t, q) in (
            ("dry months", grass_d, tree_d, qd_g[2]),
            ("MCWD mm", grass_w, tree_w, qw_g[2]),
            ("annual deficit mm", grass_a, tree_a, qa_g[2]),
        )
            hit = 100 * count(>=(q), g) / length(g)
            fa = 100 * count(>=(q), t) / length(t)
            println(
                rpad("  cut $lab >= $(round(q, digits = 1))", 34),
                rpad("treeless $(round(Int, hit))%", 16),
                rpad("forest $(round(Int, fa))%", 14),
                "separation $(round(Int, hit - fa)) pts",
            )
        end
    end
end

isempty(ARGS) && error(
    "usage: julia check_seasonality.jl <equilibrium_carbon.nc> <driver_outdir>",
)
seasonality_main(ARGS[1], ARGS[2])
