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
function obs_grid(path)
    ds = NCDatasets.NCDataset(path)
    try
        vn = first([
            k for k in keys(ds) if occursin("biomass", lowercase(k)) ||
            occursin("cveg", lowercase(k))
        ])
        v = ds[vn]
        f = occursin("kg", get(v.attrib, "units", "")) ? FT(1) : FT(0.1)
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
function main(dir; years = 400, map_half = nothing)
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
    @info "woody fraction" map_half = mh n = p.n_map_woody

    C_stem = fill(FT(NaN), nlon, nlat)
    C_leaf = fill(FT(NaN), nlon, nlat)
    C_root = fill(FT(NaN), nlon, nlat)
    cVeg = fill(FT(NaN), nlon, nlat)

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
            pools, _ = spinup(
                d,
                p;
                years,
                map_half = mh,
                T_ref_tau = p.T_ref_τ_stem,
                q_tau = p.q_τ_stem,
            )
            C_leaf[i, j], C_stem[i, j], C_root[i, j] =
                pools[2], pools[3], pools[4]
            cVeg[i, j] = sum(pools)
        end
    end

    nland = count(isfinite, C_stem)
    @info "equilibrated" nland years

    out = joinpath(
        dir,
        map_half === nothing ? "equilibrium_carbon.nc" :
        "equilibrium_carbon_maphalf$(mh).nc",
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
            obs_grid(path)
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
    main(
        ARGS[1];
        years = yi === nothing ? 400 : parse(FT, ARGS[yi + 1]),
        map_half = mi === nothing ? nothing : parse(FT, ARGS[mi + 1]),
    )
end
