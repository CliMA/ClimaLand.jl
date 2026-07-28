# Build the model-LAI validity mask used to align the MODIS-LAI observation with
# the simulated-LAI footprint (see preprocess_single_obs_var in
# generate_observations.jl).
#
# Why this exists: the ClimaLand `lai` diagnostic is NaN over a broader set of
# grid cells than the observation's `make_ocean_mask` removes (interpolated
# diagnostics are not mask-aware, so coastline / partially-ocean cells the mask
# keeps still come back NaN from the simulation). Those cells stay finite in the
# MODIS obs, so the forward-model G ensemble ends up with NaNs at obs-valid
# cells. With TransformUnscented (UKI) the innovation `y - mean(G)` then goes
# NaN and every parameter update collapses to NaN. Masking the observation by
# the simulated-LAI validity footprint removes the mismatch at the source.
#
# The footprint is structural (identical across ensemble members with valid
# parameters), so it is derived once from a reference forward-model run and
# saved. Regenerate it whenever `nelements` (the model grid) changes.
#
# Usage:
#   julia --project=.buildkite experiments/calibration/build_model_lai_validity_mask.jl \
#       <reference_simdir> [<output_path>]
#
# <reference_simdir> is a global_diagnostics output directory from a valid
# forward-model run, e.g.
#   <cal_output>/iteration_001/member_001/global_diagnostics/output_active

import ClimaAnalysis
import ClimaLand
import JLD2

include(
    joinpath(
        pkgdir(ClimaLand),
        "experiments",
        "calibration",
        "observation_utils.jl",
    ),
)
include(
    joinpath(
        pkgdir(ClimaLand),
        "ext",
        "land_sim_vis",
        "leaderboard",
        "data_sources.jl",
    ),
)

const DEFAULT_MASK_PATH = joinpath(
    pkgdir(ClimaLand),
    "experiments",
    "calibration",
    "model_lai_validity_mask.jld2",
)

"""
    build_model_lai_validity_mask(reference_simdir)

Return an `OutputVar` on the diagnostics grid that is `1.0` where the simulated
`lai` diagnostic is finite and `0.0` where it is NaN, taking the union of NaN
locations across all seasons (the structural footprint). The variable is
preprocessed with the same steps used for the simulation side of the
observation map (`preprocess_sim_var` + `average_season_across_time`) so its
grid matches the resampled observation before windowing.
"""
function build_model_lai_validity_mask(reference_simdir)
    simdir = ClimaAnalysis.SimDir(reference_simdir)
    var = get(simdir, "lai")
    var = preprocess_sim_var(var)
    var = ClimaAnalysis.average_season_across_time(var, ignore_nan = false)

    tname = ClimaAnalysis.time_name(var)
    tidx = var.dim2index[tname]
    nan_union =
        reduce((a, b) -> a .| b, eachslice(isnan.(var.data); dims = tidx))

    lats = ClimaAnalysis.latitudes(var)
    lons = ClimaAnalysis.longitudes(var)
    # nan_union follows `var`'s spatial axis order (lon, lat); OutputVar expects
    # the data laid out as (latitude, longitude), so permute to match the dims.
    lat_idx = var.dim2index[ClimaAnalysis.latitude_name(var)]
    lon_idx = var.dim2index[ClimaAnalysis.longitude_name(var)]
    spatial_nan = lat_idx < lon_idx ? nan_union : permutedims(nan_union, (2, 1))
    validity = map(isnan_cell -> isnan_cell ? 0.0 : 1.0, spatial_nan)

    attribs = Dict(
        "short_name" => "model_lai_validity",
        "long_name" => "Simulated LAI validity mask (1 where sim LAI is finite)",
    )
    dim_attribs = Dict(
        "latitude" => Dict("units" => "degrees_north"),
        "longitude" => Dict("units" => "degrees_east"),
    )
    dims = Dict("latitude" => lats, "longitude" => lons)
    return ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, validity)
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) >= 1 || error(
        "Usage: julia --project=.buildkite build_model_lai_validity_mask.jl <reference_simdir> [<output_path>]",
    )
    reference_simdir = ARGS[1]
    output_path = length(ARGS) >= 2 ? ARGS[2] : DEFAULT_MASK_PATH
    maskvar = build_model_lai_validity_mask(reference_simdir)
    n_valid = count(==(1.0), maskvar.data)
    n_total = length(maskvar.data)
    @info "Built model LAI validity mask" reference_simdir n_valid n_total
    JLD2.save_object(output_path, maskvar)
    @info "Saved mask" output_path
end
