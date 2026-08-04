## Where does the model go wrong, and is it production or turnover?
##
## The global totals hide the actual error. Globally the model's vegetation
## carbon turnover time is close to the observational estimate, which would
## suggest the pools are built the right way - but the global mean is an average
## over cells that are individually wrong in opposite directions. This bins every
## land cell by how much woody carbon the *observations* report there and
## reports, per bin, the model/observed ratio together with the two factors whose
## product it is:
##
##     cVeg = tau_veg * NPP
##
## so a bin where the ratio is high because NPP is high is a photosynthesis
## problem, and a bin where it is high because tau is long is an allocation
## problem - carbon put into wood, which turns over in decades, where it should
## have gone into leaves and roots, which turn over in months. Those need
## different fixes, and the pool ratio alone cannot tell them apart.
##
## Binning by the observation rather than by the model is deliberate: binning by
## the model would put every cell the model wrongly forested into the "forest"
## bin and hide exactly the error being looked for.
##
## Usage:
##   julia --project=.buildkite bin_by_obs.jl <equilibrium_carbon.nc>

include(joinpath(@__DIR__, "global_equilibrium.jl"))

# Edges in kg C m^-2 of observed woody carbon, descending. The lowest bin is the
# one that matters: it is land the observations say carries essentially no wood.
const EDGES = [(10.0, Inf), (5.0, 10.0), (2.0, 5.0), (0.5, 2.0), (0.0, 0.5)]

bin_label(lo, hi) =
    hi == Inf ? "> $(lo)" : (lo == 0.0 ? "< $(hi)" : "$(lo)-$(hi)")

mean_or_nan(v) = isempty(v) ? NaN : sum(v) / length(v)

function bin_main(ncpath, driverdir = nothing)
    ds = NCDatasets.NCDataset(ncpath)
    lons = Array{FT}(Array(ds["lon"][:]))
    lats = Array{FT}(Array(ds["lat"][:]))
    fld(n) =
        haskey(ds, n) ? Array{Union{Missing, FT}}(Array(ds[n][:, :])) : nothing
    cVeg = fld("cVeg")
    tau = fld("tau_veg")
    npp = fld("npp")
    stem, leaf, root = fld("C_stem"), fld("C_leaf"), fld("C_root")
    close(ds)
    cVeg === nothing && error("$ncpath has no cVeg")

    # Mean annual precipitation, if the driver run is given. MAP is the sole
    # input to the woody-fraction ramp, so it decides whether that ramp can reach
    # a cell at all: no exponent removes wood from a cell whose rainfall is above
    # the half-point, however treeless the observations say it is.
    MAP = nothing
    if driverdir !== nothing
        _, _, pra = read_monthly(driverdir, "pra")
        MAP = fill(FT(NaN), size(pra, 1), size(pra, 2))
        for i in axes(pra, 1), j in axes(pra, 2)
            s, n = FT(0), 0
            for m in 1:12
                valid(pra[i, j, m]) && (s += FT(pra[i, j, m]); n += 1)
            end
            # `pra` is the model's own running annual precipitation in m/yr -
            # already the quantity the ramp consumes, not a rate to integrate.
            # `offline_spinup.jl` feeds it to `woody_fraction` unconverted.
            n == 12 && (MAP[i, j] = abs(s / 12))
        end
    end

    for (name, path) in PRODUCTS
        isfile(path) || continue
        olons, olats, O = obs_grid(name, path)
        println("\n=== $name ===")
        println(
            rpad("obs bin", 10),
            rpad("cells", 8),
            rpad("obs", 8),
            rpad("model", 8),
            rpad("ratio", 8),
            rpad("tau_veg", 9),
            rpad("NPP", 7),
            rpad("stem", 7),
            rpad("leaf", 7),
            rpad("root", 7),
            # cVeg carries a fourth pool. Non-structural carbon has no litter
            # flux at all - it leaves only by respiration and allocation - so it
            # inflates both cVeg and the diagnosed tau_veg without ever becoming
            # litter, and it is invisible unless named.
            rpad("sugar", 7),
            MAP === nothing ? "" : "MAP    wet%  dry|wet  wetshare",
        )
        for (lo, hi) in EDGES
            om, mm, tm, nm, pm = FT[], FT[], FT[], FT[], FT[]
            sm, lm, rm, gm = FT[], FT[], FT[], FT[]
            # Split the bin's own cells by whether the ramp can act on them. A
            # bin mean hides this completely: a handful of wet cells carrying
            # forest-sized biomass outweighs hundreds of dry cells carrying
            # almost none, so "the bin is 10x too high" can be true while the
            # overwhelming majority of its cells are already right.
            drym, wetm = FT[], FT[]
            for i in eachindex(lons), j in eachindex(lats)
                m = cVeg[i, j]
                valid(m) || continue
                o = nearest(olons, olats, O, lons[i], lats[j])
                valid(o) || continue
                (lo <= o < hi) || continue
                push!(om, FT(o))
                push!(mm, FT(m))
                tau === nothing ||
                    (valid(tau[i, j]) && push!(tm, FT(tau[i, j])))
                npp === nothing ||
                    (valid(npp[i, j]) && push!(nm, FT(npp[i, j])))
                MAP === nothing ||
                    (valid(MAP[i, j]) && push!(pm, FT(MAP[i, j])))
                if MAP !== nothing && valid(MAP[i, j])
                    push!(MAP[i, j] >= 0.8 ? wetm : drym, FT(m))
                end
                if stem !== nothing && valid(stem[i, j])
                    push!(sm, FT(stem[i, j]))
                    push!(lm, FT(leaf[i, j]))
                    push!(rm, FT(root[i, j]))
                    push!(
                        gm,
                        FT(m) - FT(stem[i, j]) - FT(leaf[i, j]) -
                        FT(root[i, j]),
                    )
                end
            end
            isempty(om) && continue
            ō, m̄ = mean_or_nan(om), mean_or_nan(mm)
            println(
                rpad(bin_label(lo, hi), 10),
                rpad(length(om), 8),
                rpad(round(ō, digits = 2), 8),
                rpad(round(m̄, digits = 2), 8),
                rpad(round(m̄ / max(ō, 0.01), digits = 1), 8),
                rpad(round(mean_or_nan(tm), digits = 1), 9),
                rpad(round(mean_or_nan(nm), digits = 2), 7),
                rpad(round(mean_or_nan(sm), digits = 2), 7),
                rpad(round(mean_or_nan(lm), digits = 2), 7),
                rpad(round(mean_or_nan(rm), digits = 2), 7),
                rpad(round(mean_or_nan(gm), digits = 2), 7),
                # "wet%" is the share of the bin the MAP ramp cannot touch: cells
                # already above the half-point, where the woody fraction is >= 0.5
                # no matter how sharp the exponent is made.
                isempty(pm) ? "" :
                string(
                    round(mean_or_nan(pm), digits = 2),
                    "   ",
                    round(Int, 100 * count(>=(0.8), pm) / length(pm)),
                    "%   ",
                    round(mean_or_nan(drym), digits = 2),
                    "|",
                    round(mean_or_nan(wetm), digits = 2),
                    "   ",
                    round(
                        Int,
                        100 * sum(wetm; init = 0.0) / max(
                            sum(wetm; init = 0.0) + sum(drym; init = 0.0),
                            1e-9,
                        ),
                    ),
                    "%",
                ),
            )
        end
    end
end

isempty(ARGS) &&
    error("usage: julia bin_by_obs.jl <equilibrium_carbon.nc> [driver_outdir]")
bin_main(ARGS[1], length(ARGS) > 1 ? ARGS[2] : nothing)
