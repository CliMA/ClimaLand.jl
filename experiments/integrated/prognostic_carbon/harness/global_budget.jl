## Global carbon budget from the equilibrium map.
##
## The pools alone cannot say whether the model is right for the right reason: a
## forest of the correct size can be built from twice the production and twice
## the turnover. This reports the *fluxes* that produced the pools, each against
## an independent observational constraint the biomass products do not supply:
##
##   GPP        ~ 120-130 Pg C/yr   (GOSIF, FLUXCOM; the calibration target)
##   NPP        ~ 55-65 Pg C/yr     (roughly half of GPP across syntheses)
##   litterfall = NPP at equilibrium, by construction
##   tau_veg    ~ 10-25 yr globally  (cVeg/litterfall; Carvalhais et al. 2014
##                                    report a much longer *ecosystem* tau that
##                                    includes soil, so the two are not the same
##                                    number and should not be compared)
##
## `tau_veg` is the one diagnostic here that is not a total. It says how long
## carbon stays in living tissue, and it is the quantity a turnover-time error
## shows up in directly: too much biomass from correct production means tau is
## too long, and that is a different fix from too much production.
##
## Land area comes from the grid, not from a land-fraction field: the driver run
## masked ocean at a 0.99 threshold, so a surviving cell is essentially all land.
##
## Usage:
##   julia --project=.buildkite global_budget.jl <equilibrium_carbon.nc>

import NCDatasets

const FT = Float64
const R_EARTH = FT(6.371e6)
const PG_PER_KG = FT(1e-12)

"""
    cell_areas(lons, lats)

Area (m^2) of each cell of a regular lon-lat grid, from the cell edges implied by
the spacing. Latitude weighting is `sin` of the edges rather than `cos` of the
centre, which matters at the poles where the two differ by percent.
"""
function cell_areas(lons, lats)
    dlon = abs(lons[2] - lons[1])
    dlat = abs(lats[2] - lats[1])
    A = Array{FT}(undef, length(lons), length(lats))
    for (j, lat) in enumerate(lats)
        lo = deg2rad(max(lat - dlat / 2, -90))
        hi = deg2rad(min(lat + dlat / 2, 90))
        band = R_EARTH^2 * deg2rad(dlon) * abs(sin(hi) - sin(lo))
        A[:, j] .= band
    end
    return A
end

valid(x) = x !== missing && isfinite(x)

function total(A, area)
    s = zero(FT)
    for i in eachindex(A)
        valid(A[i]) && (s += FT(A[i]) * area[i])
    end
    return s * PG_PER_KG
end

function land_area(mask, area)
    s = zero(FT)
    for i in eachindex(mask)
        valid(mask[i]) && (s += area[i])
    end
    return s
end

function main(path)
    ds = NCDatasets.NCDataset(path)
    lons = Array{FT}(Array(ds["lon"][:]))
    lats = Array{FT}(Array(ds["lat"][:]))
    get2d(n) = haskey(ds, n) ? Array(ds[n][:, :]) : nothing
    area = cell_areas(lons, lats)

    cVeg = get2d("cVeg")
    cVeg === nothing && error("$path has no cVeg")
    aland = land_area(cVeg, area)

    println("grid $(length(lons)) x $(length(lats))")
    println(
        "land area  ",
        round(aland / 1e12, digits = 1),
        " x 10^12 m2 (observed ice-free land ~130)",
    )
    println()

    println(rpad("pool", 14), rpad("Pg C", 10), "share")
    stocks = FT[]
    for n in ("C_stem", "C_root", "C_leaf")
        A = get2d(n)
        A === nothing && continue
        push!(stocks, total(A, area))
    end
    veg = total(cVeg, area)
    for (n, v) in zip(("C_stem", "C_root", "C_leaf"), stocks)
        println(
            rpad(n, 14),
            rpad(round(v, digits = 1), 10),
            string(round(100 * v / veg, digits = 1), "%"),
        )
    end
    println(rpad("cVeg", 14), rpad(round(veg, digits = 1), 10), "100%")
    println()

    have_flux = get2d("litterfall") !== nothing
    if !have_flux
        println("no flux fields in this file - rerun global_equilibrium.jl")
        close(ds)
        return
    end

    println(rpad("flux", 14), rpad("Pg C/yr", 10), "expected")
    for (n, exp) in (
        ("gpp", "120-130 (GOSIF)"),
        ("ra", "~half of GPP"),
        ("npp", "55-65"),
        ("litterfall", "= NPP at equilibrium"),
    )
        A = get2d(n)
        A === nothing && continue
        println(rpad(n, 14), rpad(round(total(A, area), digits = 1), 10), exp)
    end

    lit = total(get2d("litterfall"), area)
    println()
    println(
        "tau_veg (global) ",
        round(veg / lit, digits = 1),
        " yr   = cVeg / litterfall",
    )

    # An area-weighted mean of a per-cell ratio is not the ratio of the totals;
    # both are reported because they answer different questions - the first is
    # the typical column, the second the whole biosphere.
    tv = get2d("tau_veg")
    if tv !== nothing
        num, den = zero(FT), zero(FT)
        for i in eachindex(tv)
            valid(tv[i]) || continue
            num += FT(tv[i]) * area[i]
            den += area[i]
        end
        println(
            "tau_veg (area-weighted mean of cells) ",
            round(num / den, digits = 1),
            " yr",
        )
    end
    close(ds)
end

if abspath(PROGRAM_FILE) == @__FILE__
    isempty(ARGS) &&
        error("usage: julia global_budget.jl <equilibrium_carbon.nc>")
    main(ARGS[1])
end
