# Global map of modeled-minus-MODIS annual-mean LAI for a snowy_land_pmodel run.
# Renders a diverging (blue = model too low, red = model too high) PNG over land.
#
# Usage: julia --project=.buildkite experiments/long_runs/plot_global_lai_bias_map.jl \
#            <run_dir> [out.png] [modis_climatology.nc]
import NCDatasets as NCD
import CairoMakie as MK

run_dir = ARGS[1]
out_png =
    length(ARGS) >= 2 ? ARGS[2] : joinpath(run_dir, "lai_minus_modis_map.png")
modis_path =
    length(ARGS) >= 3 ? ARGS[3] :
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

function load_time_mean_lai(path)
    ds = NCD.Dataset(path)
    v = ds["lai"]
    dims = NCD.dimnames(v)
    data = Array{Float64}(coalesce.(v[:, :, :], NaN))
    tax = findfirst(d -> occursin("time", lowercase(d)), dims)
    s = sum(x -> isfinite(x) ? x : 0.0, data; dims = tax)
    n = sum(isfinite.(data); dims = tax)
    tm = dropdims(ifelse.(n .== 0, NaN, s ./ n); dims = tax)
    remaining = deleteat!(collect(dims), tax)
    latax = findfirst(d -> occursin("lat", lowercase(d)), remaining)
    lat = Array{Float64}(ds["lat"][:])
    lon = Array{Float64}(ds["lon"][:])
    close(ds)
    field = latax == 1 ? tm : permutedims(tm, (2, 1))  # -> (lat, lon)
    return field, lat, lon
end

model, lat_m, lon_m = load_time_mean_lai(model_file)
obs, lat_o, lon_o = load_time_mean_lai(modis_path)

nearest(grid, x) = argmin(abs.(grid .- x))
obs_on_model =
    obs[[nearest(lat_o, l) for l in lat_m], [nearest(lon_o, l) for l in lon_m]]
d = model .- obs_on_model                    # (lat, lon); NaN over ocean (model mask)

# Plot as (lon, lat) with a symmetric diverging scale.
fig = MK.Figure(size = (1100, 560))
ax = MK.Axis(
    fig[1, 1],
    title = "Modeled − MODIS annual-mean LAI  ($(basename(run_dir)))",
    xlabel = "longitude",
    ylabel = "latitude",
)
hm = MK.heatmap!(
    ax,
    lon_m,
    lat_m,
    permutedims(d, (2, 1));
    colormap = :RdBu,
    colorrange = (-3, 3),
    nan_color = (:gray, 0.15),
)
MK.Colorbar(fig[1, 2], hm, label = "ΔLAI (model − MODIS)")
MK.save(out_png, fig)
println("WROTE $out_png")
