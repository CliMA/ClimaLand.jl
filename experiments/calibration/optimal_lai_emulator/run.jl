# Offline emulator of ZhouOptimalLAIModel: reproduces the screening and calibration of
# the moist-season VPD (with snowmelt), the season-scaled leaf cost, the tree leaf
# retention and the climate tree share. See README.md.
#
#     CLM_SURFDATA=/path/to/surfdata_0.9x1.25_16pfts__CMIP6_simyr2000_c170616.nc \
#         julia --project=.buildkite -i experiments/calibration/optimal_lai_emulator/run.jl

import ClimaComms
ClimaComms.@import_required_backends
include(joinpath(@__DIR__, "data.jl"))
include(joinpath(@__DIR__, "model.jl"))
include(joinpath(@__DIR__, "calibrate.jl"))
include(joinpath(@__DIR__, "snow.jl"))

# ---- points, observations and CLM maps ----
cand = [(lo > 180 ? lo - 360 : lo, la) for lo in 0:8:352 for la in -58:8:78]
const PTS = cand
dsS = NCDataset(SURFDATA)
CLM =
    (; lon = Float64.(dsS["LONGXY"][:, 1]), lat = Float64.(dsS["LATIXY"][1, :]))
LANDF = Float64.(coalesce.(dsS["LANDFRAC_PFT"][:, :], 0))
GLAC = Float64.(coalesce.(dsS["PCT_GLACIER"][:, :], 0))
LAKE = Float64.(coalesce.(dsS["PCT_LAKE"][:, :], 0))
close(dsS)
island(lo, la) =
    grid_at(CLM.lon, CLM.lat, LANDF, lo, la) >= 0.9 &&
    grid_at(CLM.lon, CLM.lat, LAKE, lo, la) < 30 &&
    grid_at(CLM.lon, CLM.lat, GLAC, lo, la) < 50
# MODIS averaged over the 1° land cells around each point (ocean is filled with 0)
function modis_monthly_land(points)
    ds = [
        NCDataset(joinpath(MODIS_DIR, "Yuan_et_al_$(y)_1x1.nc")) for
        y in (2008, 2009)
    ]
    lon = ds[1]["lon"][:];
    lat = ds[1]["lat"][:];
    t = vcat((d["time"][:] for d in ds)...)
    lai = cat((coalesce.(d["lai"][:, :, :], NaN) for d in ds)...; dims = 3);
    foreach(close, ds)
    out = fill(NaN, length(points), 12);
    ncell = zeros(Int, length(points))
    for (s, (lon0, lat0)) in enumerate(points)
        cells = [
            (i, j) for i in findall(x -> abs(x - lon0) < 0.6, lon) for
            j in findall(x -> abs(x - lat0) < 0.6, lat) if
            island(lon[i], lat[j])
        ]
        ncell[s] = length(cells);
        isempty(cells) && continue
        series = [mean(lai[i, j, k] for (i, j) in cells) for k in eachindex(t)]
        for m in 1:12
            ms = DateTime(2008, m, 1);
            tm = ms + ((ms + Month(1)) - ms) ÷ 2
            k = findlast(<=(tm), t);
            w = (tm - t[k]) / (t[k + 1] - t[k])
            out[s, m] = (1 - w) * series[k] + w * series[k + 1]
        end
    end
    return out, ncell
end
const OBS, NCELL = modis_monthly_land(PTS)
CLMTREE = [map_at(TREE_FILE, "tree_share", p...) for p in PTS]
CROP = [map_at(CROP_FILE, "crop_fraction", p...) for p in PTS]
glacier = [grid_at(CLM.lon, CLM.lat, GLAC, p...) for p in PTS]
LAND = findall(
    i -> NCELL[i] >= 2 && all(isfinite, OBS[i, :]) && glacier[i] < 50,
    eachindex(PTS),
)
NAT = filter(i -> isnan(CROP[i]) || CROP[i] <= 0.5, LAND)
VEG = filter(i -> isfinite(CLMTREE[i]), NAT)
const WGT = [cosd(p[2]) for p in PTS]
CLIM = point_climates(PTS)
FEAT = [features(c) for c in CLIM]

# ---- moist growing season: warm, and the 30-day water input ≥ PET/2 ----
const C4P = LAIPM.β_c4
function with_vpd(c::PointClimate, f::Features, vpd)
    χ3 = [
        Canopy.optimal_chi(c.T[h], c.P, C_CO2, vpd, LAIPM.β_c3, PC) for
        h in eachindex(c.T)
    ]
    χ4 = [
        Canopy.optimal_chi(c.T[h], c.P, C_CO2, vpd, C4P, PC) for
        h in eachindex(c.T)
    ]
    return PointClimate(
        c.lon,
        c.lat,
        c.month,
        c.T,
        c.precip,
        c.PET,
        c.VPD,
        c.P,
        c.A0c3,
        c.A0c4,
        χ3,
        χ4,
    ),
    Features(
        f.Pa,
        f.PETa,
        f.AI,
        f.gd,
        vpd,
        f.Tgs,
        f.A0c3a,
        f.A0c4a,
        f.dry_all,
        f.dry_warm,
        f.frost,
        f.pcv,
    )
end
CW = Vector{PointClimate}(undef, length(PTS));
FW = Vector{Features}(undef, length(PTS))
# water input of the moist season: precipitation (`:precip`) or rain and the melt of a
# degree-day snow store (`:snow`, as in the model)
function set_input!(kind; ddf = 3.0)
    for i in eachindex(PTS)
        input = kind == :precip ? CLIM[i].precip : water_input(CLIM[i], ddf)
        CW[i], FW[i] = moist_climate(i, input)
    end
end
set_input!(:snow)

# ---- θ = log z_tree, log z_grass, σ, f0_max, humid width of f0, tree retention,
# tree-share logistic (intercept, L_tree, dry warm months, T_growing), maintenance
# share m, log Q10 of the maintenance ----
cfg_of(θ) = Config(
    z_tree = exp(θ[1]),
    z_grass = exp(θ[2]),
    sigma = max(θ[3], 0.3),
    f0_max = clamp(θ[4], 0.3, 1.0),
    w_humid = max(θ[5], 0.0),
    ever_floor = clamp(θ[6], 0.0, 1.0),
    m_maint = length(θ) >= 11 ? clamp(θ[11], 0.0, 1.0) : 0.0,
    q10 = length(θ) >= 12 ? exp(θ[12]) : 2.0,
)
function tree_of(θ, cfg)
    β = θ[7:10]
    return [
        logistic(
            β[1] +
            β[2] * tree_lai_max(CW[i], FW[i], cfg) +
            β[3] * FW[i].dry_warm +
            β[4] * FW[i].Tgs,
        ) for i in eachindex(PTS)
    ]
end
function simulate(θ, idx)
    cfg = cfg_of(θ);
    tree = tree_of(θ, cfg)
    return [emulate(CW[i], FW[i], tree[i], 1.0, 1.0, cfg) for i in idx], tree
end
function losses(θ, idx)
    sims, tree = simulate(θ, idx)
    S = reduce(vcat, (s' for s in sims));
    w = WGT[idx]
    mse_lai = wmean(vec(mean((S .- OBS[idx, :]) .^ 2; dims = 2)), w)
    bias = wmean(vec(mean(S; dims = 2)) .- vec(mean(OBS[idx, :]; dims = 2)), w)
    v = filter(i -> isfinite(CLMTREE[i]), idx)
    return mse_lai, wmean((tree[v] .- CLMTREE[v]) .^ 2, WGT[v]), bias
end
# the defaults of toml/default_parameters.toml
θ_default = [
    log(9.92),
    log(154.0),
    1.09,
    0.65,
    0.604,
    0.5,
    -1.07,
    0.0254,
    -0.372,
    0.0917,
    0.43,
    log(2.0),
]
sims, tree = simulate(θ_default, NAT)
println("default parameters: ", fmt(score(sims, NAT)))

# ---- calibration: weak priors around μ, a penalty on the global bias, and a
# checkerboard of 24° blocks for cross-validation ----
blk(i) = (Int(fld(PTS[i][1] + 180, 24)) + Int(fld(PTS[i][2] + 90, 24))) % 2
FOLD = [filter(i -> blk(i) == b, NAT) for b in (0, 1)]
μ = [
    log(12.227),
    log(100.0),
    1.01,
    0.65,
    0.604,
    0.4,
    -1.08,
    0.109,
    -0.301,
    0.0693,
    0.3,
    log(2.0),
]
sd = [0.5, 0.5, 0.3, 0.1, 0.4, 0.3, 1.0, 0.2, 0.2, 0.05, 0.3, 0.5]
steps = [0.3, 0.3, 0.2, 0.06, 0.3, 0.15, 0.4, 0.08, 0.08, 0.02, 0.15, 0.3]
function calibrate(idx, x0; evals = 800, l0 = losses(x0, idx))
    J(θ) = (
        l = losses(θ, idx);
        l[1] / l0[1] +
        l[2] / l0[2] +
        4 * l[3]^2 / l0[1] +
        0.01 * sum(((θ .- μ) ./ sd) .^ 2)
    )
    θ, _, _ = nelder_mead(J, x0, steps; maxevals = evals)
    θ, Jv, _ = nelder_mead(J, θ, steps ./ 2; maxevals = evals)   # restart
    return θ, Jv
end
# θ, J = calibrate(NAT, θ_default)                         # ≈ 20 min
# θcv = [calibrate(FOLD[b], θ_default; evals = 600)[1] for b in 1:2]
# score(simulate(θcv[1], FOLD[2])[1], FOLD[2])            # held-out block
