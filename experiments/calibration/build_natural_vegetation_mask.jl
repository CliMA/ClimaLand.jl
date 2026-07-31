# Build the natural-vegetation mask used to restrict the stage-3 LAI
# calibration to undisturbed vegetation (see apply_natural_vegetation_mask in
# generate_observations.jl).
#
# Why this exists: the optimal-LAI model simulates natural vegetation responding
# to climate. It has no crops, no harvest, no land management and no fire, so
# grid cells dominated by those processes score the model against an LAI cycle it
# cannot represent. Zhou et al. (2025) exclude cropland-dominated, snow/ice and
# non-vegetated areas from their own evaluation for the same reason.
#
# Two criteria, both from data already on Derecho:
#
#   1. CLM `PCT_NATVEG` >= NATVEG_MIN_PCT. The natural-vegetation *landunit*
#      fraction, so a single threshold drops cropland, urban, lake, wetland and
#      glacier cells at once. Zhou threshold cropland cover at >50% from MODIS
#      MCD12C1; that product is not on this system and cannot be downloaded, and
#      PCT_NATVEG covers a superset of the same intent.
#   2. GFED4.1s mean annual burned fraction <= BURNED_MAX_PCT_PER_YEAR. Removes
#      habitually-burning savanna (N. Australia, Cerrado, the African savanna
#      belt), where MODIS sees post-fire LAI collapse and the model cannot.
#
# Deserts and other bare-but-natural cells are deliberately KEPT: they are
# undisturbed, and the model's over-prediction in arid regions is a documented
# residual that belongs in the score rather than masked out of it.
#
# Two caveats worth repeating in any write-up that uses this mask:
#   - the CLM surface dataset is `simyr2000`, so cropland expansion between 2000
#     and the calibration year is not represented;
#   - GFED4.1s ends in 2016, so the fire criterion is necessarily climatological
#     ("this cell burns often") rather than specific to the calibrated year.
#
# The mask is static for a given grid. Regenerate it whenever `nelements`
# changes. It is not committed (`*.jld2` is gitignored).
#
# Usage:
#   julia --project=.buildkite experiments/calibration/build_natural_vegetation_mask.jl [<output_path>]
#
# Source paths default to the Derecho ClimaArtifacts mirror and can be
# overridden with the CLM_SURFDATA_PATH and GFED_BURNTAREA_PATH env vars.

import ClimaAnalysis
import ClimaComms
import ClimaLand
import JLD2
import NCDatasets
import Statistics

include(joinpath(@__DIR__, "observation_utils.jl"))

const DEFAULT_MASK_PATH = joinpath(@__DIR__, "natural_vegetation_mask.jld2")

const CLM_SURFDATA_PATH = get(
    ENV,
    "CLM_SURFDATA_PATH",
    "/glade/campaign/univ/ucit0011/ClimaArtifacts/artifacts/clm_data/surfdata_0.9x1.25_16pfts__CMIP6_simyr2000_c170616.nc",
)
const GFED_BURNTAREA_PATH = get(
    ENV,
    "GFED_BURNTAREA_PATH",
    "/glade/campaign/univ/ucit0011/ClimaArtifacts/artifacts/ILAMB/DATA/burntArea/GFED4.1S/burntArea.nc",
)

"""Minimum natural-vegetation landunit fraction (%) for a cell to be kept."""
const NATVEG_MIN_PCT = 80.0

"""Maximum mean annual burned fraction (% yr^-1) for a cell to be kept."""
const BURNED_MAX_PCT_PER_YEAR = 5.0

"""Years of GFED4.1s (1997-2016) averaged into the burned-fraction climatology."""
const BURNED_YEARS = 2001:2016

"""
    to_pm180(lon)

Return `lon` in [-180, 180), the convention of the calibration grid.
"""
to_pm180(lon) = lon >= 180 ? lon - 360 : lon

"""
    regrid(target_lats, target_lons, src_lats, src_lons, data)

Regrid `data` (indexed `[lon, lat]`, as NCDatasets returns a `(lat, lon)`
variable) onto the target grid by averaging every source cell whose center falls
inside the target cell, falling back to the nearest source cell where none does.
Averaging keeps fractional fields meaningful when the source is finer than the
target (GFED at 0.5 degrees); the fallback covers the coarser source (CLM at
0.9x1.25 degrees). NaNs are ignored in the average.
"""
function regrid(target_lats, target_lons, src_lats, src_lons, data)
    src_lons = to_pm180.(src_lons)
    dlat = abs(target_lats[2] - target_lats[1])
    dlon = abs(target_lons[2] - target_lons[1])

    # For each source cell, the target cell it falls into (0 = outside).
    assign(src, target, d) =
        map(src) do s
            i = findfirst(t -> t - d / 2 <= s < t + d / 2, target)
            isnothing(i) ? 0 : i
        end
    lat_bin = assign(src_lats, target_lats, dlat)
    lon_bin = assign(src_lons, target_lons, dlon)

    sums = zeros(length(target_lats), length(target_lons))
    counts = zeros(Int, length(target_lats), length(target_lons))
    for (jsrc, j) in enumerate(lat_bin), (isrc, i) in enumerate(lon_bin)
        (j == 0 || i == 0) && continue
        v = data[isrc, jsrc]
        (ismissing(v) || isnan(v)) && continue
        sums[j, i] += v
        counts[j, i] += 1
    end

    out = fill(NaN, length(target_lats), length(target_lons))
    for j in eachindex(target_lats), i in eachindex(target_lons)
        if counts[j, i] > 0
            out[j, i] = sums[j, i] / counts[j, i]
        else
            jn = argmin(abs.(src_lats .- target_lats[j]))
            in_ = argmin(abs.(src_lons .- target_lons[i]))
            v = data[in_, jn]
            out[j, i] = (ismissing(v) || isnan(v)) ? NaN : v
        end
    end
    return out
end

"""
    natural_vegetation_fraction(lats, lons)

Return the CLM natural-vegetation landunit percentage on the target grid.
"""
function natural_vegetation_fraction(lats, lons)
    isfile(CLM_SURFDATA_PATH) || error(
        "CLM surface dataset not found at $CLM_SURFDATA_PATH. Set CLM_SURFDATA_PATH.",
    )
    NCDatasets.NCDataset(CLM_SURFDATA_PATH) do ds
        src_lats = Array(ds["LATIXY"][:, :])[1, :]
        src_lons = Array(ds["LONGXY"][:, :])[:, 1]
        natveg = Array(ds["PCT_NATVEG"][:, :])
        return regrid(lats, lons, src_lats, src_lons, natveg)
    end
end

"""
    burned_fraction(lats, lons)

Return the GFED4.1s mean annual burned fraction (% yr^-1) over `BURNED_YEARS` on
the target grid. Cells with no GFED data are returned as 0 (treated as
unburned): they are outside the product's land footprint and are removed by the
ocean mask anyway.
"""
function burned_fraction(lats, lons)
    isfile(GFED_BURNTAREA_PATH) || error(
        "GFED4.1s burned-area file not found at $GFED_BURNTAREA_PATH. Set GFED_BURNTAREA_PATH.",
    )
    NCDatasets.NCDataset(GFED_BURNTAREA_PATH) do ds
        src_lats = Array(ds["lat"][:])
        src_lons = Array(ds["lon"][:])
        # The file starts in 1997-01 and is monthly; select whole years.
        first_year = 1997
        i1 = (first(BURNED_YEARS) - first_year) * 12 + 1
        i2 = (last(BURNED_YEARS) - first_year + 1) * 12
        ba = Array(ds["burntArea"][:, :, i1:i2])
        # _FillValue is -999; NCDatasets returns it as missing, but guard anyway.
        ba = map(x -> (ismissing(x) || x < -900) ? NaN : Float64(x), ba)
        annual =
            dropdims(sum(replace(ba, NaN => 0.0), dims = 3), dims = 3) ./
            length(BURNED_YEARS)
        # Keep NaN where the product has no data at all, so `regrid` can
        # distinguish "never burns" from "not observed".
        observed = dropdims(any(.!isnan.(ba), dims = 3), dims = 3)
        annual[.!observed] .= NaN
        regridded = regrid(lats, lons, src_lats, src_lons, annual)
        return map(x -> isnan(x) ? 0.0 : x, regridded)
    end
end

"""
    build_natural_vegetation_mask(nelements)

Return an `OutputVar` on the calibration grid that is `1.0` where the cell is
natural vegetation that rarely burns and `0.0` elsewhere.
"""
function build_natural_vegetation_mask(nelements)
    lats, lons = get_lat_lon_from_resolution(nelements)
    natveg = natural_vegetation_fraction(lats, lons)
    burned = burned_fraction(lats, lons)

    keep = (natveg .>= NATVEG_MIN_PCT) .& (burned .<= BURNED_MAX_PCT_PER_YEAR)
    data = map(k -> k ? 1.0 : 0.0, keep)

    attribs = Dict(
        "short_name" => "natural_vegetation",
        "long_name" => "Natural undisturbed vegetation mask (1 where kept)",
        "natveg_min_pct" => NATVEG_MIN_PCT,
        "burned_max_pct_per_year" => BURNED_MAX_PCT_PER_YEAR,
    )
    dim_attribs = Dict(
        "latitude" => Dict("units" => "degrees_north"),
        "longitude" => Dict("units" => "degrees_east"),
    )
    dims = Dict("latitude" => lats, "longitude" => lons)
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)
    return var, natveg, burned
end

if abspath(PROGRAM_FILE) == @__FILE__
    output_path = length(ARGS) >= 1 ? ARGS[1] : DEFAULT_MASK_PATH
    nelements = (180, 360, 15)
    maskvar, natveg, burned = build_natural_vegetation_mask(nelements)

    # Report the removal fractions over land only, so the numbers are meaningful.
    lats, lons = get_lat_lon_from_resolution(nelements)
    ocean_mask = make_ocean_mask(nelements)
    land_probe = ClimaAnalysis.OutputVar(
        Dict("short_name" => "probe"),
        Dict("latitude" => lats, "longitude" => lons),
        Dict(
            "latitude" => Dict("units" => "degrees_north"),
            "longitude" => Dict("units" => "degrees_east"),
        ),
        ones(length(lats), length(lons)),
    )
    land = .!isnan.(ocean_mask(land_probe).data)
    weights = repeat(cosd.(lats), 1, length(lons))
    area_frac(m) = sum(weights[land .& m]) / sum(weights[land])

    @info "Natural-vegetation mask" nelements NATVEG_MIN_PCT BURNED_MAX_PCT_PER_YEAR
    @info "Removed by each criterion (area-weighted fraction of land)" natveg =
        area_frac(natveg .< NATVEG_MIN_PCT) fire =
        area_frac(burned .> BURNED_MAX_PCT_PER_YEAR) combined =
        area_frac(maskvar.data .== 0.0)
    @info "Land cells kept" kept = count(land .& (maskvar.data .== 1.0)) total =
        count(land)
    JLD2.save_object(output_path, maskvar)
    @info "Saved mask" output_path
end
