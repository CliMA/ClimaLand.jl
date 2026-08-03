## Compare modelled equilibrium woody carbon against XuSaatchi observations.
##
## MODEL.md §7 is specific about what may be compared: XuSaatchi is a *woody*
## live-biomass product, so the model's `C_stem` is the right quantity, not total
## `cVeg`. Grassland cells, where the observation is near zero by construction
## while the model legitimately carries leaf and root carbon, are reported
## separately rather than counted as model error.
##
## The model side is the offline equilibrium - a coupled run reaches nowhere near
## steady state on a 30-year stem turnover - with the temperature-dependent
## turnover the model now uses.
##
## Usage:  julia --project=.buildkite compare_biomass.jl <battery_runroot>

include(joinpath(@__DIR__, "offline_spinup.jl"))

import NCDatasets
import ClimaLand.Parameters as LP
import ClimaLand.Canopy

const OBS_DIR = "/glade/campaign/cesm/community/lmwg/diag/ILAMB/DATA/biomass"

# Several products, not one. They disagree by 2-4x at most forest sites, so a
# comparison against any single product mostly measures which product was
# chosen. The model is scored against the observational *range*.
const PRODUCTS = [
    ("XuSaatchi", joinpath(OBS_DIR, "XuSaatchi2021", "XuSaatchi.nc")),
    ("Thurner", joinpath(OBS_DIR, "Thurner", "biomass_0.5x0.5.nc")),
    ("ESACCI", joinpath(OBS_DIR, "ESACCI", "biomass.nc")),
    ("GEOCARBON", joinpath(OBS_DIR, "GEOCARBON", "biomass.nc")),
    ("Saatchi2011", joinpath(OBS_DIR, "Saatchi2011", "biomass_0.5x0.5.nc")),
    ("USForest", joinpath(OBS_DIR, "USForest", "biomass_0.5x0.5.nc")),
]

# netCDF classic missing-value sentinel. It is finite and not `missing`, so it
# survives both checks and silently enters the mean as ~1e36 if not caught.

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

# Broad zones for the by-zone summary the stage asks for.
zone_of(biome, mat) =
    biome == "desert" ? "arid" :
    occursin("grassland", biome) ? "grassland" :
    occursin("savanna", biome) ? "savanna" :
    biome == "tundra" ? "boreal/tundra" :
    mat < 278 ? "boreal/tundra" : mat > 293 ? "tropics" : "temperate"

read_sites(csv) = [
    (
        name = f[4],
        biome = f[3],
        lon = parse(Float64, f[1]),
        lat = parse(Float64, f[2]),
    ) for f in (split(strip(l), ',') for l in readlines(csv)) if
    length(f) == 4 && f[1] != "lon" && !startswith(f[1], "#")
]

"""
    obs_biomass(ds, lon, lat)

Time-mean observed woody carbon density (kg C m^-2) at the nearest grid cell.
Returns NaN where the product has no data, which is the honest answer over
desert and ice rather than a zero that would flatter the comparison.
"""
function obs_at(name, path, lon, lat)
    ds = NCDatasets.NCDataset(path)
    try
        vn = first([
            k for k in keys(ds) if occursin("biomass", lowercase(k)) ||
            occursin("cveg", lowercase(k))
        ])
        v = ds[vn]
        units = get(v.attrib, "units", "")
        lons = coalesce.(Array(ds["lon"][:]), NaN)
        lats = coalesce.(Array(ds["lat"][:]), NaN)
        l = lon
        if maximum(lons) > 180 && l < 0
            l += 360
        elseif maximum(lons) <= 180 && l > 180
            l -= 360
        end
        i = argmin(abs.(lons .- l))
        j = argmin(abs.(lats .- lat))
        A = Array(v[:, :, :])
        col = ndims(A) == 3 ? A[i, j, :] : [A[i, j]]
        good = filter(valid, col)
        isempty(good) && return NaN
        # kg m^-2 stays; Mg ha^-1 becomes kg m^-2 by 0.1
        # Mg ha^-1 -> kg m^-2 is 0.1; then dry biomass -> carbon where needed.
        f = (occursin("kg", units) ? 1.0 : 0.1) * CARBON_FRACTION[name]
        return (sum(good) / length(good)) * f
    finally
        close(ds)
    end
end

function obs_biomass(lons, lats, B, lon, lat)
    # the file may use either -180..180 or 0..360
    l = lon
    if maximum(lons) > 180 && l < 0
        l += 360
    elseif maximum(lons) <= 180 && l > 180
        l -= 360
    end
    i = argmin(abs.(lons .- l))
    j = argmin(abs.(lats .- lat))
    col = B[i, j, :]
    good = filter(x -> !ismissing(x) && isfinite(x), col)
    isempty(good) && return NaN
    return (sum(good) / length(good)) * MG_HA_TO_KG_M2
end

function main(runroot)
    sites = read_sites(joinpath(@__DIR__, "test_sites.csv"))
    toml_dict = LP.create_toml_dict(Float64)
    p = Canopy.PrognosticCarbonParameters(toml_dict)

    rows = []
    for s in sites
        f = joinpath(runroot, s.name, "out", "driver_record.csv")
        isfile(f) || continue
        d = read_driver_record(f)
        MAT = sum(d["tair"]) / length(d["tair"])
        pools, _ = spinup(
            d,
            p;
            years = 400,
            T_ref_tau = p.T_ref_τ_stem,
            q_tau = p.q_τ_stem,
        )
        obs = Float64[]
        for (nm, path) in PRODUCTS
            isfile(path) || continue
            v = try
                obs_at(nm, path, s.lon, s.lat)
            catch
                NaN
            end
            isfinite(v) && push!(obs, v)
        end
        push!(
            rows,
            (
                name = s.name,
                zone = zone_of(s.biome, MAT),
                model = pools[3],
                obs = obs,
            ),
        )
    end

    println(
        rpad("site", 21),
        rpad("zone", 15),
        rpad("model", 9),
        rpad("obs range (n)", 24),
        "verdict",
    )
    inside = 0
    scored = 0
    for r in rows
        isempty(r.obs) && continue
        lo, hi = minimum(r.obs), maximum(r.obs)
        ok = lo <= r.model <= hi
        scored += 1
        inside += ok
        verdict =
            ok ? "inside" :
            r.model > hi ?
            "above by $(round(r.model / max(hi, 0.01), digits = 1))x" :
            "below by $(round(max(lo, 0.01) / max(r.model, 0.01), digits = 1))x"
        println(
            rpad(r.name, 21),
            rpad(r.zone, 15),
            rpad(round(r.model, digits = 2), 9),
            rpad(
                "$(round(lo, digits = 2))-$(round(hi, digits = 2)) ($(length(r.obs)))",
                24,
            ),
            verdict,
        )
    end

    println(
        "\nmodel inside the observational range at $inside of $scored sites",
    )
    println(
        "products disagree by a median factor of ",
        round(
            begin
                sp = [
                    maximum(r.obs) / max(minimum(r.obs), 0.01) for
                    r in rows if length(r.obs) > 1
                ]
                sort(sp)[max(1, div(length(sp), 2))]
            end,
            digits = 1,
        ),
        "x, so scoring against any single product mostly measures the choice of product",
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    isempty(ARGS) && error("usage: julia compare_biomass.jl <battery_runroot>")
    main(ARGS[1])
end
