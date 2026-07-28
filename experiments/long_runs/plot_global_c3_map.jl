# Global maps of the C3 fraction: dynamic model (optimal-LAI C3/C4 competition, the `fc3`
# diagnostic) vs the prescribed CLM `c3_proportion` map, and their difference. Assesses
# how the dynamic C3/C4 competition performs against the standard prescribed distribution.
# Three stacked panels; also prints the area-weighted global-mean C3 fraction + bias/RMSE.
#
# Usage: julia --project=.buildkite experiments/long_runs/plot_global_c3_map.jl \
#            <run_dir> [out.png] [veg_map.nc]
import NCDatasets as NCD
import CairoMakie as MK

run_dir = ARGS[1]
out_png = length(ARGS) >= 2 ? ARGS[2] : joinpath(run_dir, "c3_fraction_map.png")
veg_map =
    length(ARGS) >= 3 ? ARGS[3] :
    joinpath(
        homedir(),
        ".julia/artifacts/6284ddefbc7937d9c1fb68fa731ff3f00b68e917/vegetation_properties_map.nc",
    )
model_file = joinpath(
    run_dir,
    "global_diagnostics",
    "output_active",
    "fc3_1M_average.nc",
)
const SPINUP_MONTHS = parse(Int, get(ENV, "SPINUP_MONTHS", "0"))

# Time-mean model fc3 -> (lat, lon); ocean cells (NaN through the run) stay NaN.
function load_fc3(path; skip = 0)
    ds = NCD.Dataset(path)
    v = ds["fc3"]
    dims = NCD.dimnames(v)
    data = Array{Float64}(coalesce.(v[:, :, :], NaN))
    tax = findfirst(d -> occursin("time", lowercase(d)), dims)
    laax = findfirst(d -> occursin("lat", lowercase(d)), dims)
    loax = findfirst(d -> occursin("lon", lowercase(d)), dims)
    if skip > 0 && skip < size(data, tax)
        data = collect(selectdim(data, tax, (skip + 1):size(data, tax)))
    end
    s = sum(x -> isfinite(x) ? x : 0.0, data; dims = tax)
    n = sum(isfinite.(data); dims = tax)
    tm = dropdims(ifelse.(n .== 0, NaN, s ./ n); dims = tax)
    lat = Array{Float64}(ds["lat"][:])
    lon = Array{Float64}(ds["lon"][:])
    close(ds)
    field = laax < loax ? tm : permutedims(tm, (2, 1))  # -> (lat, lon)
    return field, lat, lon
end

model, lat_m, lon_m = load_fc3(model_file; skip = SPINUP_MONTHS)

# CLM prescribed C3 proportion, nearest-remapped to the model grid.
dso = NCD.Dataset(veg_map)
obs_raw = Array{Float64}(coalesce.(dso["c3_proportion"][:, :], NaN))  # (lat, lon)
lat_o = Array{Float64}(dso["lat"][:])
lon_o = Array{Float64}(dso["lon"][:])
close(dso)
# c3_proportion stored (lat, lon)? guard by matching sizes.
obs_ll =
    size(obs_raw, 1) == length(lat_o) ? obs_raw : permutedims(obs_raw, (2, 1))

nearest(grid, x) = argmin(abs.(grid .- x))
lon_m360 = [l < 0 && maximum(lon_o) > 180 ? l + 360 : l for l in lon_m]
obs = obs_ll[
    [nearest(lat_o, l) for l in lat_m],
    [nearest(lon_o, l) for l in lon_m360],
]

d = model .- obs
valid = isfinite.(model) .& isfinite.(obs)          # model land mask ∩ obs
w = cosd.(lat_m) .* ones(1, length(lon_m))
w[.!valid] .= 0.0
W = sum(w)
gm_model = sum(w .* ifelse.(valid, model, 0.0)) / W
gm_obs = sum(w .* ifelse.(valid, obs, 0.0)) / W
bias = gm_model - gm_obs
rmse = sqrt(sum(w .* ifelse.(valid, d .^ 2, 0.0)) / W)
println("GLOBAL_C3_MODEL_MEAN $gm_model")
println("GLOBAL_C3_OBS_MEAN $gm_obs")
println("GLOBAL_C3_BIAS $bias")
println("GLOBAL_C3_RMSE $rmse")
println("N_VALID_CELLS $(count(valid))")

fig = MK.Figure(size = (1050, 1200))
function panel(row, field, title, cmap, crange, cblabel)
    ax = MK.Axis(
        fig[row, 1],
        title = title,
        xlabel = "longitude",
        ylabel = "latitude",
    )
    hm = MK.heatmap!(
        ax,
        lon_m,
        lat_m,
        permutedims(field, (2, 1));
        colormap = cmap,
        colorrange = crange,
        nan_color = (:gray, 0.12),
    )
    MK.Colorbar(fig[row, 2], hm, label = cblabel)
end
panel(
    1,
    model,
    "Dynamic C3 fraction — model (optimal-LAI C3/C4 competition, $(basename(run_dir)))",
    :viridis,
    (0, 1),
    "fraction C3",
)
panel(
    2,
    obs,
    "Prescribed C3 fraction — CLM c3_proportion (data)",
    :viridis,
    (0, 1),
    "fraction C3",
)
panel(3, d, "Difference (model − CLM)", :RdBu, (-1, 1), "Δ fraction C3")
MK.save(out_png, fig)
println("WROTE $out_png")
