# Area-weighted global comparison of modeled LAI (from a snowy_land_pmodel run) to the
# MODIS LAI climatology. Reports area-weighted (cos-latitude) global RMSE and bias over
# the intersection of finite land cells — the flux/area-relevant skill metric.
#
# Both fields are on the same 1x1 lon/lat grid; MODIS is a 12-month climatology and the
# model output is monthly, so each is reduced to its time-mean before comparison.
#
# Usage: julia --project=.buildkite experiments/long_runs/compare_global_lai_modis.jl \
#            <run_dir> [modis_climatology.nc]
import NCDatasets as NCD
using Statistics

run_dir = ARGS[1]
modis_path =
    length(ARGS) >= 2 ? ARGS[2] :
    joinpath(
        homedir(),
        ".julia/artifacts/cd070348ab70f5d1b165c814f70e3822c0173eed/modis_lai_climatology.nc",
    )
model_file = joinpath(
    run_dir,
    "global_diagnostics",
    "output_active",
    "lai_1M_average.nc",
)

# Months of model spin-up to skip before scoring (0 for the MODIS climatology). Set via
# ENV for a multi-year run, e.g. SPINUP_MONTHS=24 uses only the final year of a 3-yr run.
const SPINUP_MONTHS = parse(Int, get(ENV, "SPINUP_MONTHS", "0"))

# Load `lai` and reduce over its time axis (identified by dimension name) to (lat, lon).
# `skip` drops that many leading time steps (model spin-up) before averaging.
function load_time_mean_lai(path; skip = 0)
    ds = NCD.Dataset(path)
    v = ds["lai"]
    dims = NCD.dimnames(v)
    data = Array{Float64}(coalesce.(v[:, :, :], NaN))
    tax = findfirst(d -> occursin("time", lowercase(d)), dims)
    if skip > 0 && skip < size(data, tax)
        data = selectdim(data, tax, (skip + 1):size(data, tax)) |> collect
    end
    tm = dropdims(_nanmean(data, tax); dims = tax)  # remaining axes in original order
    remaining = deleteat!(collect(dims), tax)
    latax = findfirst(d -> occursin("lat", lowercase(d)), remaining)
    lat = Array{Float64}(ds["lat"][:])
    lon = Array{Float64}(ds["lon"][:])
    close(ds)
    # Return as (lat, lon) regardless of stored order.
    field = latax == 1 ? tm : permutedims(tm, (2, 1))
    return field, lat, lon
end

# Time-mean ignoring NaN; an all-NaN cell (e.g. ocean, masked at every step) stays NaN
# so it defines the land domain rather than leaking in as a spurious zero.
function _nanmean(a, dim)
    s = sum(x -> isfinite(x) ? x : 0.0, a; dims = dim)
    n = sum(isfinite.(a); dims = dim)
    return ifelse.(n .== 0, NaN, s ./ n)
end

model, lat_m, lon_m = load_time_mean_lai(model_file; skip = SPINUP_MONTHS)
obs, lat_o, lon_o = load_time_mean_lai(modis_path)

# Remap MODIS (lat, lon) onto the model grid by nearest neighbour. The two 1-degree
# grids are offset by half a cell (model on integer nodes, MODIS on half-degree centres),
# a sub-gridscale shift for a global skill metric. Longitude is treated periodically.
nearest(grid, x) = argmin(abs.(grid .- x))
lat_idx = [nearest(lat_o, l) for l in lat_m]
lon_idx = [nearest(lon_o, l) for l in lon_m]
obs_on_model = obs[lat_idx, lon_idx]

d = model .- obs_on_model
valid = isfinite.(d)
w = cosd.(lat_m) .* ones(1, length(lon_m))          # cos-latitude area weight (lat, lon)
w[.!valid] .= 0.0
W = sum(w)
bias = sum(w .* ifelse.(valid, d, 0.0)) / W
rmse = sqrt(sum(w .* ifelse.(valid, d .^ 2, 0.0)) / W)
gm_model = sum(w .* ifelse.(valid, model, 0.0)) / W
gm_obs = sum(w .* ifelse.(valid, obs_on_model, 0.0)) / W

println("GLOBAL_LAI_RMSE $rmse")
println("GLOBAL_LAI_BIAS $bias")
println("GLOBAL_LAI_MODEL_MEAN $gm_model")
println("GLOBAL_LAI_MODIS_MEAN $gm_obs")
println("N_VALID_CELLS $(count(valid))")
