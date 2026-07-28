# Global LAI seasonal cycle: area-weighted global-mean modeled vs MODIS LAI for each
# calendar month, plus per-month spatial bias and RMSE over land. Answers "do we capture
# the global seasonality?" Writes a PNG (two panels) and prints the monthly table.
#
# Usage: julia --project=.buildkite experiments/long_runs/plot_global_lai_seasonal_cycle.jl \
#            <run_dir> [out.png] [modis_climatology.nc]
import NCDatasets as NCD
import CairoMakie as MK
import Dates

run_dir = ARGS[1]
out_png =
    length(ARGS) >= 2 ? ARGS[2] : joinpath(run_dir, "lai_seasonal_cycle.png")
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

# Simulation start for the model's time axis (units are raw "s" since sim start with no
# reference date, so calendar months are derived from this base). Matches the
# snowy_land_pmodel default (2008-03-01); override with ARGS[4] as "YYYY-MM-DD".
const BASE_DATE =
    length(ARGS) >= 4 ? Dates.DateTime(ARGS[4]) : Dates.DateTime(2008, 3, 1)
# Model spin-up months to drop before building the climatology (ENV, default 0).
const SPINUP_MONTHS = parse(Int, get(ENV, "SPINUP_MONTHS", "0"))

# Calendar month of each time value: MODIS decodes to real datetimes; the model gives
# raw seconds-since-sim-start (Float64), converted via BASE_DATE.
_months(tv) =
    eltype(tv) <: Real ?
    [Dates.month(BASE_DATE + Dates.Second(round(Int, t))) for t in tv] :
    [Dates.month(t) for t in tv]

# Load `lai` as a monthly climatology (12, lat, lon), binning each time step by its
# calendar month; ocean/masked cells (all-NaN through the year) stay NaN. `positional`
# assigns months 1..12 by index (for a 12-step climatology whose stored timestamps drift,
# e.g. the MODIS product) instead of decoding the time axis.
function load_monthly_climatology(path; positional = false, skip = 0)
    ds = NCD.Dataset(path)
    v = ds["lai"]
    dims = NCD.dimnames(v)
    raw = Array{Float64}(coalesce.(v[:, :, :], NaN))
    tax = findfirst(d -> occursin("time", lowercase(d)), dims)
    laax = findfirst(d -> occursin("lat", lowercase(d)), dims)
    loax = findfirst(d -> occursin("lon", lowercase(d)), dims)
    months0 =
        positional ? [mod1(m, 12) for m in 1:size(raw, tax)] :
        _months(ds["time"][:])
    if skip > 0 && skip < size(raw, tax)                     # drop model spin-up months
        raw = collect(selectdim(raw, tax, (skip + 1):size(raw, tax)))
        months0 = months0[(skip + 1):end]
    end
    months = months0
    lat = Array{Float64}(ds["lat"][:])
    lon = Array{Float64}(ds["lon"][:])
    close(ds)
    nlat, nlon = length(lat), length(lon)
    clim = fill(NaN, 12, nlat, nlon)
    for m in 1:12
        idxs = findall(==(m), months)
        isempty(idxs) && continue
        acc = zeros(nlat, nlon)
        cnt = zeros(Int, nlat, nlon)
        for t in idxs
            sel = Vector{Any}(undef, 3)
            sel[tax] = t
            sel[laax] = Colon()
            sel[loax] = Colon()
            slab = raw[sel...]                       # (lat, lon) or (lon, lat)
            slab = laax < loax ? slab : permutedims(slab, (2, 1))
            fin = isfinite.(slab)
            acc[fin] .+= slab[fin]
            cnt[fin] .+= 1
        end
        clim[m, :, :] = ifelse.(cnt .== 0, NaN, acc ./ max.(cnt, 1))
    end
    return clim, lat, lon
end

model, lat_m, lon_m = load_monthly_climatology(model_file; skip = SPINUP_MONTHS)
obs, lat_o, lon_o = load_monthly_climatology(modis_path; positional = true)

# Nearest-neighbour remap MODIS (0.5° offset) onto the model grid.
nearest(grid, x) = argmin(abs.(grid .- x))
lat_idx = [nearest(lat_o, l) for l in lat_m]
lon_idx = [nearest(lon_o, l) for l in lon_m]
obs_r = obs[:, lat_idx, lon_idx]

# cos-latitude area weights (lat, lon).
w0 = cosd.(lat_m) .* ones(1, length(lon_m))

mon = ["J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"]
model_gm = fill(NaN, 12)
obs_gm = fill(NaN, 12)
bias = fill(NaN, 12)
rmse = fill(NaN, 12)
for m in 1:12
    dm = model[m, :, :]
    om = obs_r[m, :, :]
    valid = isfinite.(dm) .& isfinite.(om)
    w = w0 .* valid
    W = sum(w)
    model_gm[m] = sum(w .* ifelse.(valid, dm, 0.0)) / W
    obs_gm[m] = sum(w .* ifelse.(valid, om, 0.0)) / W
    bias[m] = model_gm[m] - obs_gm[m]
    rmse[m] = sqrt(sum(w .* ifelse.(valid, (dm .- om) .^ 2, 0.0)) / W)
end

println("month  model  MODIS   bias    rmse")
for m in 1:12
    println(
        rpad(mon[m], 6),
        rpad(round(model_gm[m], digits = 3), 7),
        rpad(round(obs_gm[m], digits = 3), 7),
        rpad(round(bias[m], digits = 3), 8),
        round(rmse[m], digits = 3),
    )
end
amp_m = maximum(model_gm) - minimum(model_gm)
amp_o = maximum(obs_gm) - minimum(obs_gm)
println(
    "seasonal amplitude: model ",
    round(amp_m, digits = 3),
    " MODIS ",
    round(amp_o, digits = 3),
)

# First-harmonic (annual) phase: fit x[m] ≈ mean + A·cos(2π(m − φ)/12) and report the
# peak month φ ∈ (0,12]. The model−MODIS phase lag (wrapped to ±6 months) quantifies the
# greenup timing error: negative = model leads (greens up early), positive = model lags.
function harmonic_peak_month(x)
    ω = 2π / 12
    a = sum(x[m] * cos(ω * m) for m in 1:12)
    b = sum(x[m] * sin(ω * m) for m in 1:12)
    φ = atan(b, a) / ω                       # month of the cosine peak, in (−6, 6]
    return mod(φ - 1, 12) + 1                # shift to (0, 12] with m indexing months 1..12
end
peak_m = harmonic_peak_month(model_gm)
peak_o = harmonic_peak_month(obs_gm)
lag = mod(peak_m - peak_o + 6, 12) - 6        # wrap to (−6, 6]; <0 model leads
println(
    "seasonal phase (harmonic peak month): model ",
    round(peak_m, digits = 2),
    " MODIS ",
    round(peak_o, digits = 2),
    "  lag ",
    round(lag, digits = 2),
    " mo (",
    lag < 0 ? "model leads/early" : "model lags/late",
    ")",
)
println(
    "discrete peak month: model ",
    argmax(model_gm),
    " MODIS ",
    argmax(obs_gm),
)

# Two-panel figure: seasonal cycle (lines) + per-month bias & RMSE (bars).
fig = MK.Figure(size = (900, 720))
ax1 = MK.Axis(
    fig[1, 1],
    title = "Global-mean land LAI seasonal cycle  ($(basename(run_dir)))",
    xlabel = "month",
    ylabel = "area-weighted LAI",
    xticks = (1:12, mon),
)
MK.lines!(ax1, 1:12, obs_gm, color = :black, linewidth = 3, label = "MODIS")
MK.scatter!(ax1, 1:12, obs_gm, color = :black)
MK.lines!(
    ax1,
    1:12,
    model_gm,
    color = :seagreen,
    linewidth = 3,
    label = "model",
)
MK.scatter!(ax1, 1:12, model_gm, color = :seagreen)
MK.axislegend(ax1, position = :rb)

ax2 = MK.Axis(
    fig[2, 1],
    title = "per-month bias (model − MODIS) and spatial RMSE over land",
    xlabel = "month",
    ylabel = "LAI",
    xticks = (1:12, mon),
)
MK.barplot!(ax2, 1:12, bias, color = :steelblue, label = "bias")
MK.scatterlines!(
    ax2,
    1:12,
    rmse,
    color = :darkorange,
    linewidth = 2,
    label = "RMSE",
)
MK.hlines!(ax2, [0.0], color = :gray, linestyle = :dash)
MK.axislegend(ax2, position = :rt)
MK.save(out_png, fig)
println("WROTE $out_png")
