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

const OBS_PATH = "/glade/campaign/cesm/community/lmwg/diag/ILAMB/DATA/biomass/XuSaatchi2021/XuSaatchi.nc"
const MG_HA_TO_KG_M2 = 0.1   # 1 Mg ha^-1 = 0.1 kg m^-2

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

    ds = NCDatasets.NCDataset(OBS_PATH)
    lons = Array(ds["lon"][:])
    lats = Array(ds["lat"][:])
    # NCDatasets already returns Julia-order (lon, lat, time) for a CF file
    # declared (time, lat, lon); permuting again would transpose it back.
    B = Array(ds["biomass"][:, :, :])
    lons = coalesce.(lons, NaN)
    lats = coalesce.(lats, NaN)
    close(ds)

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
        push!(
            rows,
            (
                name = s.name,
                biome = s.biome,
                zone = zone_of(s.biome, MAT),
                model = pools[3],           # C_stem: the woody pool
                obs = obs_biomass(lons, lats, B, s.lon, s.lat),
            ),
        )
    end

    println(
        rpad("site", 21),
        rpad("zone", 15),
        rpad("model C_stem", 14),
        rpad("XuSaatchi", 12),
        "diff",
    )
    for r in rows
        d = r.model - r.obs
        println(
            rpad(r.name, 21),
            rpad(r.zone, 15),
            rpad(round(r.model, digits = 2), 14),
            rpad(isnan(r.obs) ? "no data" : round(r.obs, digits = 2), 12),
            isnan(r.obs) ? "" : string(d > 0 ? "+" : "", round(d, digits = 2)),
        )
    end

    println(
        "\n",
        rpad("zone", 16),
        rpad("n", 4),
        rpad("bias", 10),
        rpad("RMSE", 10),
        "mean obs",
    )
    for z in (
        "tropics",
        "temperate",
        "boreal/tundra",
        "savanna",
        "grassland",
        "arid",
    )
        sel = [r for r in rows if r.zone == z && isfinite(r.obs)]
        isempty(sel) && continue
        diffs = [r.model - r.obs for r in sel]
        bias = sum(diffs) / length(diffs)
        rmse = sqrt(sum(abs2, diffs) / length(diffs))
        mobs = sum(r.obs for r in sel) / length(sel)
        println(
            rpad(z, 16),
            rpad(length(sel), 4),
            rpad(string(bias > 0 ? "+" : "", round(bias, digits = 2)), 10),
            rpad(round(rmse, digits = 2), 10),
            round(mobs, digits = 2),
        )
    end

    # The woody comparison is only meaningful where the observation is woody.
    woody = [
        r for r in rows if isfinite(r.obs) && !(r.zone in ("grassland", "arid"))
    ]
    if !isempty(woody)
        diffs = [r.model - r.obs for r in woody]
        println(
            "\nwoody zones only (n=",
            length(woody),
            "):  bias ",
            round(sum(diffs) / length(diffs), digits = 2),
            "   RMSE ",
            round(sqrt(sum(abs2, diffs) / length(diffs)), digits = 2),
            "   kg C m^-2",
        )
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    isempty(ARGS) && error("usage: julia compare_biomass.jl <battery_runroot>")
    main(ARGS[1])
end
