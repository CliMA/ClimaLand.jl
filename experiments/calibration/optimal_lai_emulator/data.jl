# Offline emulator of ZhouOptimalLAIModel, data part: hourly low-res ERA5 2008 (recycled)
# at the land points of its 8° grid, MODIS 2008 monthly LAI, and CLM5 maps (tree share,
# evergreen share of trees, crop fraction) used only as fitting/scoring targets.
using ClimaLand, NCDatasets, Statistics, Printf, Dates, LinearAlgebra
import ClimaLand.Canopy
import ClimaLand.Parameters as LP
const FT = Float64
const ART = joinpath(homedir(), ".julia/artifacts")
const ERA5_FILE = joinpath(
    ART,
    "fd826542e6ab1758c8d5a2cc6dbf1f67b47aadee/era5_2008_1.0x1.0_lowres.nc",
)
const MODIS_DIR = joinpath(ART, "81d4bc1b22e94d66eec3d80a2d21b5dd5303bf8f")
const REPO = pkgdir(ClimaLand)
const TREE_FILE = joinpath(REPO, "artifacts/clm_tree_share/clm_tree_share.nc")
const CROP_FILE =
    joinpath(REPO, "artifacts/clm_crop_fraction/clm_crop_fraction.nc")
const SURFDATA = get(ENV, "CLM_SURFDATA", joinpath(@__DIR__, "surfdata.nc"))

const TOML = LP.create_toml_dict(FT)
const EARTH = LP.LandParameters(TOML)
const PM = Canopy.PModel{FT}(
    ClimaLand.Domains.Point(; z_sfc = FT(0), longlat = FT.((-60, -3))),
    TOML,
)
const PC = PM.constants
const LPAR = Canopy.OptimalLAIParameters{FT}(TOML)
const LAIPM = Canopy.optimal_lai_pmodel_parameters(PM.parameters, LPAR)
const ϵC = FT(TOML["canopy_emissivity"])
const λγ = FT(TOML["wavelength_per_PAR_photon"])
const THERMO = LP.thermodynamic_parameters(EARTH)
const TF = LP.T_freeze(EARTH)
const DAY = FT(86400)
const YEAR = 365 * DAY
const C_CO2 = FT(3.9e-4)

# CLM5 evergreen share of trees: natpft 1, 2 (needleleaf evergreen), 4, 5 (broadleaf
# evergreen) over trees 1-8 (0-based; Julia indices are one more)
function clm_evergreen_share()
    ds = NCDataset(SURFDATA)
    pct = Float64.(coalesce.(ds["PCT_NAT_PFT"][:, :, :], 0))
    close(ds)
    trees = dropdims(sum(pct[:, :, 2:9]; dims = 3); dims = 3)
    ever = pct[:, :, 2] .+ pct[:, :, 3] .+ pct[:, :, 5] .+ pct[:, :, 6]
    return ifelse.(trees .> 0, ever ./ trees, NaN)
end

function map_at(file, var, lon0, lat0)
    ds = NCDataset(file)
    lon = Float64.(Array(ds["lon"]));
    lat = Float64.(Array(ds["lat"]))
    x = Float64.(coalesce.(Array(ds[var]), NaN))
    close(ds)
    return grid_at(lon, lat, x, lon0, lat0)
end
grid_at(lon, lat, x, lon0, lat0) = x[
    argmin(abs.(mod.(lon .- lon0 .+ 180, 360) .- 180)),
    argmin(abs.(lat .- lat0)),
]

# MODIS LAI at the 1° ERA5 point (mean of the four cells around it), mid-month 2008
function modis_monthly(points)
    ds = [
        NCDataset(joinpath(MODIS_DIR, "Yuan_et_al_$(y)_1x1.nc")) for
        y in (2008, 2009)
    ]
    lon = ds[1]["lon"][:];
    lat = ds[1]["lat"][:]
    t = vcat((d["time"][:] for d in ds)...)
    lai = cat((coalesce.(d["lai"][:, :, :], NaN) for d in ds)...; dims = 3)
    foreach(close, ds)
    out = fill(NaN, length(points), 12)
    for (s, (lon0, lat0)) in enumerate(points)
        i = findall(x -> abs(x - lon0) < 0.6, lon);
        j = findall(x -> abs(x - lat0) < 0.6, lat)
        series = [
            begin
                v = filter(!isnan, vec(lai[i, j, k]));
                isempty(v) ? NaN : mean(v)
            end for k in eachindex(t)
        ]
        for m in 1:12
            ms = DateTime(2008, m, 1);
            tm = ms + ((ms + Month(1)) - ms) ÷ 2
            k = findlast(<=(tm), t);
            w = (tm - t[k]) / (t[k + 1] - t[k])
            out[s, m] = (1 - w) * series[k] + w * series[k + 1]
        end
    end
    return out
end

# Hourly forcing and the quantities that depend only on it
struct PointClimate
    lon::Float64
    lat::Float64
    month::Vector{Int}
    T::Vector{FT}       # K
    precip::Vector{FT}  # mol m^-2 s^-1
    PET::Vector{FT}
    VPD::Vector{FT}     # Pa
    P::FT               # mean surface pressure, Pa
    A0c3::Vector{FT}    # mol m^-2 s^-1, at fAPAR = 1
    A0c4::Vector{FT}
    χc3::Vector{FT}     # at the growing-season VPD
    χc4::Vector{FT}
end

function point_climates(points)
    ds = NCDataset(ERA5_FILE)
    lon = ds["lon"][:];
    lat = ds["lat"][:]
    month = Dates.month.(ds["time"][:])
    ρ_m = LP.ρ_m_liq(EARTH);
    ρ_l = LP.ρ_cloud_liq(EARTH)
    σ = LP.Stefan(EARTH);
    M_w = LP.molar_mass_water(EARTH)
    TD = ClimaLand.Thermodynamics
    out = PointClimate[]
    for (lon0, lat0) in points
        i = findfirst(==(mod(lon0, 360)), lon);
        j = findfirst(==(lat0), lat)
        g(v) = FT.(ds[v][i, j, :])
        T, Td, P = g("t2m"), g("d2m"), g("sp")
        q = FT.(
            ClimaLand.specific_humidity_from_dewpoint.(Td, T, P, Ref(EARTH)),
        )
        SW, LW = g("msdwswrf"), g("msdwlwrf")
        precip = max.(g("mtpr"), 0) ./ ρ_l .* ρ_m
        PET = [
            Canopy.potential_evaporation(
                SW[h],
                LW[h],
                T[h],
                P[h],
                ϵC,
                σ,
                M_w,
                THERMO,
            ) for h in eachindex(T)
        ]
        VPD = [
            max(TD.vapor_pressure_deficit(THERMO, T[h], P[h], q[h]), zero(FT)) for h in eachindex(T)
        ]
        grow = T .> TF
        vpd_gs = sum(VPD .* grow) / max(sum(grow), 1)
        PPFD = [
            Canopy.compute_PPFD(
                SW[h] / 2,
                λγ,
                PC.lightspeed,
                PC.planck_h,
                PC.N_a,
            ) for h in eachindex(T)
        ]
        p3 = [
            Canopy.compute_A0_and_χ(
                one(FT),
                LAIPM,
                PC,
                EARTH,
                T[h],
                P[h],
                q[h],
                C_CO2,
                PPFD[h],
                one(FT),
                vpd_gs,
            ) for h in eachindex(T)
        ]
        p4 = [
            Canopy.compute_A0_and_χ(
                zero(FT),
                LAIPM,
                PC,
                EARTH,
                T[h],
                P[h],
                q[h],
                C_CO2,
                PPFD[h],
                one(FT),
                vpd_gs,
            ) for h in eachindex(T)
        ]
        push!(
            out,
            PointClimate(
                lon0,
                lat0,
                month,
                T,
                precip,
                PET,
                VPD,
                mean(P),
                [p.A0_c3 for p in p3],
                [p.A0_c4 for p in p4],
                [p.χ for p in p3],
                [p.χ for p in p4],
            ),
        )
    end
    close(ds)
    return out
end
