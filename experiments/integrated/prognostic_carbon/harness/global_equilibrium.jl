## Global equilibrium carbon map from monthly driver output (stage 5).
##
## Twenty columns test magnitude; only a map tests *spatial pattern*, which is
## half of what the project set out to produce. This integrates the pools to
## steady state independently in every land column of the 1x1 degree driver run
## `global_driver.jl` writes, then scores the result against the gridded ILAMB
## woody-biomass products on their own grids.
##
## Why offline rather than coupled: stem turnover is decades, and with the
## temperature scaling it reaches centuries in the boreal, so a coupled global
## run would need millennia of integration to equilibrate. The pools are driven
## one-way, so integrating them offline against a driver climatology is exact
## given that climatology - not an approximation of the coupled answer, but the
## same answer for less compute.
##
## Monthly drivers are justified by `check_monthly_drivers.jl`: median 0.4% and
## maximum 3.3% difference against daily-driven equilibria at the battery sites.
##
## Usage:
##   julia --project=.buildkite -t auto global_equilibrium.jl <driver_outdir> [--years N]

include(joinpath(@__DIR__, "offline_spinup.jl"))

import NCDatasets
import ClimaLand.Parameters as LP
import ClimaLand.Canopy

# The driver record is a 365-day climatology, so this is the length of the cycle
# the fluxes are averaged over, not a calendar convention.
const SECONDS_PER_YEAR = FT(365 * 86400)

const OBS_DIR = "/glade/campaign/cesm/community/lmwg/diag/ILAMB/DATA/biomass"
const PRODUCTS = [
    ("XuSaatchi", joinpath(OBS_DIR, "XuSaatchi2021", "XuSaatchi.nc")),
    ("Thurner", joinpath(OBS_DIR, "Thurner", "biomass_0.5x0.5.nc")),
    ("ESACCI", joinpath(OBS_DIR, "ESACCI", "biomass.nc")),
    ("GEOCARBON", joinpath(OBS_DIR, "GEOCARBON", "biomass.nc")),
    ("Saatchi2011", joinpath(OBS_DIR, "Saatchi2011", "biomass_0.5x0.5.nc")),
    ("USForest", joinpath(OBS_DIR, "USForest", "biomass_0.5x0.5.nc")),
]

# netCDF classic missing-value sentinel: finite, and not `missing`, so it passes
# both of the obvious checks and enters a mean as ~1e36 if not rejected here.

# Carbon fraction of each product. Four of the six report **dry biomass**, not
# carbon - readable only in `long_name`, not in `units`, which is how an earlier
# version of this file got it wrong and made the model look unbiased against
# ESACCI and too low against GEOCARBON. 0.5 is the standard woody carbon
# fraction, and it reconciles GEOCARBON with the ILAMB benchmark total (774 Pg
# as-read x 0.5 = 387 Pg, against ILAMB's 364 Pg).
#   XuSaatchi   "annual carbon density ... live woody vegetation"  -> carbon
#   Thurner     "Carbon Mass in Vegetation"                        -> carbon
#   ESACCI      "Above-ground biomass"                             -> biomass
#   GEOCARBON   "above_ground_biomass"                             -> biomass
#   Saatchi2011 "above- and below-ground live biomass"             -> biomass
#   USForest    "US 48-States and Alaska Forest Biomass"           -> biomass
const CARBON_FRACTION = Dict(
    "XuSaatchi" => 1.0,
    "Thurner" => 1.0,
    "ESACCI" => 0.5,
    "GEOCARBON" => 0.5,
    "Saatchi2011" => 0.5,
    "USForest" => 0.5,
)

const SENTINEL = 1e30
valid(x) = !ismissing(x) && isfinite(x) && abs(x) < SENTINEL

const DAYS_IN_MONTH = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

"""
    as_lon_lat_time(ds, name)

Load a variable as `A[lon, lat, time]` whatever order it is stored in.

This is not defensive boilerplate. ClimaLand's own `NetCDFWriter` produces
`(time, lon, lat)` in Julia's index order while every ILAMB biomass product is
`(lon, lat, time)`, so a single hard-coded convention silently transposes one of
the two - and a transposed field is not an error, just a wrong map. Dispatching
on the dimension names removes the choice.
"""
function as_lon_lat_time(ds, name)
    v = ds[name]
    dims = collect(NCDatasets.dimnames(v))
    A = Array(v[ntuple(_ -> Colon(), length(dims))...])
    want = ("lon", "lat", "time")
    perm = Int[]
    for w in want
        k = findfirst(d -> lowercase(d) == w, dims)
        k === nothing || push!(perm, k)
    end
    length(perm) == length(dims) || error(
        "$name has dimensions $dims; expected lon/lat and optionally time",
    )
    A = permutedims(A, perm)
    # A product with a single map and no time axis still indexes uniformly.
    return ndims(A) == 2 ? reshape(A, size(A, 1), size(A, 2), 1) : A
end

"""
    month_of_year_index()

Day-of-year to month on a no-leap calendar. The drivers are a whole number of
recycled years, so a fixed calendar is exact rather than approximate.
"""
function month_of_year_index()
    m = Vector{Int}(undef, 365)
    j = 1
    for (mo, L) in enumerate(DAYS_IN_MONTH), _ in 1:L
        m[j] = mo
        j += 1
    end
    return m
end

"""
    read_monthly(dir, var)

Reads `<var>_1M_average.nc` and returns `(lons, lats, A)` with `A[lon, lat,
month]` averaged over whole years into a 12-month climatology. Averaging to a
climatology rather than integrating against the raw series is deliberate: the
driver run is short, and a climatology removes its particular weather while
keeping the seasonal cycle the pools actually respond to.
"""
function read_monthly(dir, var)
    path = joinpath(dir, "$(var)_1M_average.nc")
    isfile(path) || error("missing driver file $path")
    ds = NCDatasets.NCDataset(path)
    try
        lons = Array{FT}(coalesce.(Array(ds["lon"][:]), NaN))
        lats = Array{FT}(coalesce.(Array(ds["lat"][:]), NaN))
        raw = as_lon_lat_time(ds, var)
        nt = size(raw, 3)
        nyr = div(nt, 12)
        nyr >= 1 || error("$var has only $nt monthly records; need >= 12")
        clim = fill(FT(NaN), size(raw, 1), size(raw, 2), 12)
        for m in 1:12
            for i in axes(raw, 1), j in axes(raw, 2)
                s = FT(0)
                n = 0
                for y in 0:(nyr - 1)
                    v = raw[i, j, y * 12 + m]
                    if valid(v)
                        s += FT(v)
                        n += 1
                    end
                end
                clim[i, j, m] = n > 0 ? s / n : FT(NaN)
            end
        end
        return lons, lats, clim
    finally
        close(ds)
    end
end

"""
    obs_grid(path)

Time-mean of an observational product on its own grid, in kg C m^-2.
"""
function obs_grid(name, path)
    ds = NCDatasets.NCDataset(path)
    try
        vn = first([
            k for k in keys(ds) if occursin("biomass", lowercase(k)) ||
            occursin("cveg", lowercase(k))
        ])
        v = ds[vn]
        f =
            (occursin("kg", get(v.attrib, "units", "")) ? FT(1) : FT(0.1)) *
            FT(CARBON_FRACTION[name])
        lons = Array{FT}(coalesce.(Array(ds["lon"][:]), NaN))
        lats = Array{FT}(coalesce.(Array(ds["lat"][:]), NaN))
        A = as_lon_lat_time(ds, vn)
        B = fill(FT(NaN), size(A, 1), size(A, 2))
        for i in axes(A, 1), j in axes(A, 2)
            s = FT(0)
            n = 0
            for t in axes(A, 3)
                x = A[i, j, t]
                if valid(x)
                    s += FT(x)
                    n += 1
                end
            end
            B[i, j] = n > 0 ? f * s / n : FT(NaN)
        end
        return lons, lats, B
    finally
        close(ds)
    end
end

"""
    nearest(lons, lats, A, lon, lat)

Nearest-neighbour sample, handling either a -180..180 or a 0..360 convention.
Nearest-neighbour rather than interpolation on purpose: the products carry hard
land/ocean edges, and interpolating across them would invent biomass in the sea.
"""
function nearest(lons, lats, A, lon, lat)
    l = lon
    if maximum(lons) > 180 && l < 0
        l += 360
    elseif maximum(lons) <= 180 && l > 180
        l -= 360
    end
    i = argmin(abs.(lons .- l))
    j = argmin(abs.(lats .- lat))
    return A[i, j]
end

# `map_half` selects the precipitation-based woody fraction. `nothing` means
# "use the model's own parameter", which is the only correct default: this script
# exists to reproduce what the model does, and a hard-coded 0 would silently
# disable a mechanism the model has enabled.
const PR_PATH = "/glade/campaign/cesm/community/lmwg/diag/ILAMB/DATA/pr/GPCCv2018/pr.nc"

# The conventional threshold for a "dry" month in tropical ecology, and the ET
# proxy in the standard maximum-cumulative-water-deficit definition.
const DRY_MM_PER_MONTH = FT(100)

"""
    pr_climatology(path)

Monthly precipitation climatology in mm/month on the product's own grid.
"""
function pr_climatology(path)
    ds = NCDatasets.NCDataset(path)
    try
        lons = Array{FT}(coalesce.(Array(ds["lon"][:]), NaN))
        lats = Array{FT}(coalesce.(Array(ds["lat"][:]), NaN))
        v = ds["pr"]
        units = get(v.attrib, "units", "")
        A = as_lon_lat_time(ds, "pr")
        nt = size(A, 3)
        nyr = div(nt, 12)
        clim = fill(FT(NaN), size(A, 1), size(A, 2), 12)
        for m in 1:12
            # Days in the month convert a rate to a monthly total; using 30 for
            # every month would bias the dry-season count near the threshold.
            dpm = FT(DAYS_IN_MONTH[m])
            f =
                occursin("mm d-1", units) ? dpm :
                occursin("kg m-2 s-1", units) ? dpm * FT(86400) :
                error("unhandled precipitation units '$units'")
            for i in axes(A, 1), j in axes(A, 2)
                s, n = FT(0), 0
                for y in 0:(nyr - 1)
                    x = A[i, j, y * 12 + m]
                    valid(x) && (s += FT(x); n += 1)
                end
                n > 0 && (clim[i, j, m] = f * s / n)
            end
        end
        return lons, lats, clim
    finally
        close(ds)
    end
end

"""
    annual_deficit(p)

Total annual water deficit (mm): the monthly shortfall below the threshold,
summed. Unlike MCWD this needs no running maximum and no memory of where the wet
season started - it is a pointwise function of precipitation accumulated over a
year, i.e. exactly the machinery `P_annual` already uses in the model. If it
separates as well as MCWD, it is the form to implement.
"""
annual_deficit(p) = sum(max(DRY_MM_PER_MONTH - x, zero(FT)) for x in p)

"Number of months below the dry threshold."
dry_months(p) = count(<(DRY_MM_PER_MONTH), p)


"""
    annual_deficit_grid(lons, lats)

Annual water deficit (mm) sampled onto the model grid. Kept here rather than in
`check_seasonality.jl` so the diagnostic that evaluated this predictor and the
equilibrium run that applies it cannot drift apart.
"""
function annual_deficit_grid(lons, lats)
    plons, plats, P = pr_climatology(PR_PATH)
    D = fill(FT(NaN), length(lons), length(lats))
    for i in eachindex(lons), j in eachindex(lats)
        p = ntuple(
            m -> nearest(plons, plats, view(P, :, :, m), lons[i], lats[j]),
            12,
        )
        all(valid, p) && (D[i, j] = annual_deficit(FT.(p)))
    end
    return D
end

# `deficit_half` (mm/yr) switches on the seasonality limit tested by
# `check_seasonality.jl`: woody allocation is further reduced where the annual
# water deficit is large, which is the one predictor found so far that separates
# wet savanna from wet forest. Precipitation for it comes from GPCC, so this is a
# feasibility test of the predictor, not of the model - a model-side version must
# read ERA5 and be re-verified.
function main(
    dir;
    years = 400,
    map_half = nothing,
    map_n = nothing,
    q_map = 1.0,
    deficit_half = nothing,
    deficit_n = 4.0,
)
    mo = month_of_year_index()
    vars = ("gpp", "rd", "ct", "tair", "fc3", "pra")
    @info "reading monthly driver climatologies" dir vars
    data = Dict{String, Any}()
    lons = lats = nothing
    for v in vars
        lons, lats, data[v] = read_monthly(dir, v)
    end
    nlon, nlat = length(lons), length(lats)
    @info "grid" nlon nlat threads = Threads.nthreads()

    toml_dict = LP.create_toml_dict(FT)
    p = Canopy.PrognosticCarbonParameters(toml_dict)
    mh = map_half === nothing ? p.map_half_woody : FT(map_half)
    # The exponent controls how abruptly wood disappears below the half-point.
    # It matters far more than it looks: stem carbon turns over 15-30x slower
    # than leaf and root, so even a few percent of allocation to wood dominates
    # the equilibrium pool, and a gentle ramp leaves woody carbon standing in
    # cells the observations say carry none.
    mn = map_n === nothing ? p.n_map_woody : FT(map_n)
    @info "woody fraction" map_half = mh n = mn
    Dfield =
        deficit_half === nothing ? nothing : annual_deficit_grid(lons, lats)
    deficit_half === nothing || @info "seasonality limit" deficit_half deficit_n

    C_stem = fill(FT(NaN), nlon, nlat)
    C_leaf = fill(FT(NaN), nlon, nlat)
    C_root = fill(FT(NaN), nlon, nlat)
    cVeg = fill(FT(NaN), nlon, nlat)
    # `spinup` already averages these over the last driver cycle; carrying them
    # out costs nothing and they are what the pools cannot be judged by alone.
    # Litterfall is the flux the soil carbon pool consumes, and cVeg/litterfall
    # is the vegetation carbon turnover time - a quantity with its own
    # observational estimates, independent of the biomass products.
    gpp_eq = fill(FT(NaN), nlon, nlat)
    npp_eq = fill(FT(NaN), nlon, nlat)
    ra_eq = fill(FT(NaN), nlon, nlat)
    litter_eq = fill(FT(NaN), nlon, nlat)
    tau_veg = fill(FT(NaN), nlon, nlat)

    Threads.@threads for j in 1:nlat
        # One scratch record per thread-iteration, reused across longitudes so
        # the sweep does not allocate 365-element vectors 65000 times.
        d = Dict{String, Vector{FT}}(v => Vector{FT}(undef, 365) for v in vars)
        for i in 1:nlon
            # A column is land if GPP is defined there in every month. The
            # driver run masks ocean to NaN, so this is the run's own land mask
            # rather than a second one that could disagree with it.
            ok = true
            for m in 1:12
                isfinite(data["gpp"][i, j, m]) || (ok = false; break)
            end
            ok || continue
            for v in vars, k in 1:365
                x = data[v][i, j, mo[k]]
                # `ct` is NaN-masked wherever leaf+stem area index is zero; with
                # no canopy the canopy temperature relaxes to air temperature.
                if !isfinite(x)
                    x =
                        v == "ct" ? data["tair"][i, j, mo[k]] :
                        v == "pra" ? FT(0) : FT(NaN)
                end
                d[v][k] = x
            end
            any(v -> any(!isfinite, d[v]), vars) && continue
            # Declines from 1 in an aseasonal climate toward 0 as the
            # dry-season deficit grows, so a cell with forest rainfall delivered
            # in five months is no longer treated as a forest cell.
            ws = FT(1)
            if Dfield !== nothing && isfinite(Dfield[i, j])
                ws = 1 / (1 + (Dfield[i, j] / FT(deficit_half))^FT(deficit_n))
            end
            pools, fl = spinup(
                d,
                p;
                years,
                map_half = mh,
                map_n = mn,
                w_scale = ws,
                q_map,
                T_ref_tau = p.T_ref_τ_stem,
                q_tau = p.q_τ_stem,
            )
            C_leaf[i, j], C_stem[i, j], C_root[i, j] =
                pools[2], pools[3], pools[4]
            cVeg[i, j] = sum(pools)
            gpp_eq[i, j] = fl.GPP * SECONDS_PER_YEAR
            ra_eq[i, j] = fl.Ra * SECONDS_PER_YEAR
            npp_eq[i, j] = (fl.GPP - fl.Ra) * SECONDS_PER_YEAR
            litter_eq[i, j] = fl.litter * SECONDS_PER_YEAR
            # Undefined where nothing grows; 0/0 would report a spurious age.
            tau_veg[i, j] =
                fl.litter > 0 ? cVeg[i, j] / (fl.litter * SECONDS_PER_YEAR) :
                FT(NaN)
        end
    end

    nland = count(isfinite, C_stem)
    @info "equilibrated" nland years

    out = joinpath(
        dir,
        (
            map_half === nothing &&
            map_n === nothing &&
            q_map == 1 &&
            deficit_half === nothing
        ) ? "equilibrium_carbon.nc" :
        deficit_half === nothing ?
        "equilibrium_carbon_mh$(mh)_n$(mn)_q$(q_map).nc" :
        "equilibrium_carbon_dh$(deficit_half)_dn$(deficit_n).nc",
    )
    NCDatasets.NCDataset(out, "c") do ds
        NCDatasets.defDim(ds, "lon", nlon)
        NCDatasets.defDim(ds, "lat", nlat)
        NCDatasets.defVar(ds, "lon", lons, ("lon",))
        NCDatasets.defVar(ds, "lat", lats, ("lat",))
        for (n, A, l) in (
            ("C_stem", C_stem, "equilibrium woody carbon"),
            ("C_leaf", C_leaf, "equilibrium leaf carbon"),
            ("C_root", C_root, "equilibrium root carbon"),
            ("cVeg", cVeg, "equilibrium total vegetation carbon"),
        )
            v = NCDatasets.defVar(ds, n, A, ("lon", "lat"))
            v.attrib["units"] = "kg m-2"
            v.attrib["long_name"] = l
        end
        for (n, A, l) in (
            ("gpp", gpp_eq, "equilibrium gross primary production"),
            ("npp", npp_eq, "equilibrium net primary production"),
            ("ra", ra_eq, "equilibrium autotrophic respiration"),
            ("litterfall", litter_eq, "equilibrium litter carbon flux"),
        )
            v = NCDatasets.defVar(ds, n, A, ("lon", "lat"))
            v.attrib["units"] = "kg m-2 yr-1"
            v.attrib["long_name"] = l
        end
        v = NCDatasets.defVar(ds, "tau_veg", tau_veg, ("lon", "lat"))
        v.attrib["units"] = "yr"
        v.attrib["long_name"] = "vegetation carbon turnover time (cVeg/litterfall)"
    end
    @info "wrote $out"

    # Score against every product on the cells where that product has data.
    println(
        "\n",
        rpad("product", 14),
        rpad("n cells", 10),
        rpad("model mean", 12),
        rpad("obs mean", 11),
        rpad("bias", 9),
        "spatial r",
    )
    for (name, path) in PRODUCTS
        isfile(path) || continue
        olons, olats, O = try
            obs_grid(name, path)
        catch e
            @warn "could not read $name" exception = e
            continue
        end
        xs, ys = FT[], FT[]
        for j in 1:nlat, i in 1:nlon
            isfinite(C_stem[i, j]) || continue
            o = nearest(olons, olats, O, lons[i], lats[j])
            isfinite(o) || continue
            push!(xs, C_stem[i, j])
            push!(ys, FT(o))
        end
        n = length(xs)
        n > 10 || continue
        mx, my = sum(xs) / n, sum(ys) / n
        sx = sqrt(sum((xs .- mx) .^ 2) / n)
        sy = sqrt(sum((ys .- my) .^ 2) / n)
        r =
            (sx > 0 && sy > 0) ? sum((xs .- mx) .* (ys .- my)) / (n * sx * sy) :
            NaN
        println(
            rpad(name, 14),
            rpad(n, 10),
            rpad(round(mx, digits = 2), 12),
            rpad(round(my, digits = 2), 11),
            rpad(round(mx - my, digits = 2), 9),
            round(r, digits = 3),
        )
    end
    println(
        "\nSpatial correlation is the quantity of interest here: the products \
         disagree on magnitude by a median factor of ~3.4x, so a bias against \
         any one of them is weaker evidence than whether the pattern lines up.",
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    isempty(ARGS) &&
        error("usage: julia -t auto global_equilibrium.jl <driver_outdir>")
    yi = findfirst(==("--years"), ARGS)
    mi = findfirst(==("--map-half"), ARGS)
    qi = findfirst(==("--q-map"), ARGS)
    ni = findfirst(==("--n-map"), ARGS)
    di = findfirst(==("--deficit-half"), ARGS)
    dni = findfirst(==("--deficit-n"), ARGS)
    main(
        ARGS[1];
        years = yi === nothing ? 400 : parse(FT, ARGS[yi + 1]),
        map_half = mi === nothing ? nothing : parse(FT, ARGS[mi + 1]),
        map_n = ni === nothing ? nothing : parse(FT, ARGS[ni + 1]),
        q_map = qi === nothing ? 1.0 : parse(FT, ARGS[qi + 1]),
        deficit_half = di === nothing ? nothing : parse(FT, ARGS[di + 1]),
        deficit_n = dni === nothing ? 4.0 : parse(FT, ARGS[dni + 1]),
    )
end
