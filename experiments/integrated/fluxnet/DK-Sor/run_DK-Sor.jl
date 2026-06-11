# DK-Sor FLUXNET site simulation using PModel photosynthesis.
# Site: Sorø Forest, Denmark (DK-Sor)
# Forest type: Temperate deciduous broadleaf (European beech, Fagus sylvatica)
# Lat: 55.49°N  Lon: 11.64°E  UTC offset: +1
# Forcing is read from the local NetCDF file in DK_Sor/.
# Runs from 1997-01-01 to 1998-01-01 (1 year).
# Output is written to experiments/integrated/fluxnet/DK-Sor/out/.
import ClimaComms
ClimaComms.@import_required_backends
using ClimaCore
import ClimaParams as CP
using Dates
using ClimaDiagnostics
using ClimaUtilities

using ClimaLand
using ClimaLand.Domains: Column
using ClimaLand.Snow
using ClimaLand.Soil
using ClimaLand.Soil.Biogeochemistry
using ClimaLand.Canopy
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.Simulations: LandSimulation, solve!

import ClimaLand.FluxnetSimulations as FluxnetSimulations
import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput, evaluate!
import ClimaUtilities.TimeManager: date
using NCDatasets

const FT = Float64
toml_dict = LP.create_toml_dict(FT)

const site_ID     = "DK-Sor"
const site_ID_val = :DK_Sor
const start_date  = DateTime(1997, 1, 1)   # UTC
const stop_date   = DateTime(1998, 1, 1)   # UTC
const dt          = Float64(450)           # 7.5-minute timestep

climaland_dir = pkgdir(ClimaLand)
met_nc_path   = joinpath(climaland_dir, "DK_Sor", "DK-Sor_1997-2014_FLUXNET2015_Met.nc")

prognostic_land_components = (:canopy, :snow, :soil, :soilco2)

# ── Site metadata ──────────────────────────────────────────────────────────────
(; dz_tuple, nelements, zmin, zmax) =
    FluxnetSimulations.get_domain_info(FT, Val(site_ID_val))
(; time_offset, lat, long) =
    FluxnetSimulations.get_location(FT, Val(site_ID_val))
(; atmos_h) = FluxnetSimulations.get_fluxtower_height(FT, Val(site_ID_val))
(;
    soil_ν,
    soil_K_sat,
    soil_S_s,
    soil_vg_n,
    soil_vg_α,
    θ_r,
    ν_ss_quartz,
    ν_ss_om,
    ν_ss_gravel,
    z_0m_soil,
    z_0b_soil,
    soil_ϵ,
    soil_α_PAR,
    soil_α_NIR,
    Ω,
    χl,
    α_PAR_leaf,
    λ_γ_PAR,
    τ_PAR_leaf,
    α_NIR_leaf,
    τ_NIR_leaf,
    ϵ_canopy,
    ac_canopy,
    Drel,
    g0,
    SAI,
    f_root_to_shoot,
    K_sat_plant,
    ψ63,
    Weibull_param,
    a,
    conductivity_model,
    retention_model,
    plant_ν,
    plant_S_s,
    rooting_depth,
    h_canopy,
) = FluxnetSimulations.get_parameters(FT, Val(site_ID_val))

# ── Domain ─────────────────────────────────────────────────────────────────────
land_domain    = Column(;
    zlim     = (zmin, zmax),
    nelements,
    dz_tuple,
    longlat  = (long, lat),
)
surface_domain = ClimaLand.Domains.obtain_surface_domain(land_domain)

# ── Forcing from NetCDF ────────────────────────────────────────────────────────
(; atmos, radiation) = FluxnetSimulations.prescribed_forcing_netcdf(
    met_nc_path, lat, long, time_offset, atmos_h, start_date, toml_dict, FT,
)

# ── LAI from DK-Sor NetCDF ─────────────────────────────────────────────────────
surface_space = land_domain.space.surface
LAI, maxLAI = NCDataset(met_nc_path, "r") do ds
    time_vals  = ds["time"][:]
    lai_data   = Float64.(coalesce.(ds["LAI"][1, 1, :], NaN))
    lai_secs   = Float64[
        Second(t - Hour(time_offset) - start_date).value for t in time_vals
    ]
    valid      = .!isnan.(lai_data)
    tvi        = TimeVaryingInput(lai_secs[valid], lai_data[valid])
    mx         = maximum(lai_data[valid])
    tvi, mx
end
RAI = maxLAI * f_root_to_shoot

# ── Soil ───────────────────────────────────────────────────────────────────────
forcing = (; atmos, radiation)
soil_albedo = Soil.ConstantTwoBandSoilAlbedo{FT}(;
    PAR_albedo = soil_α_PAR,
    NIR_albedo = soil_α_NIR,
)
retention_parameters = (;
    ν        = soil_ν,
    θ_r,
    K_sat    = soil_K_sat,
    hydrology_cm = vanGenuchten{FT}(; α = soil_vg_α, n = soil_vg_n),
)
composition_parameters = (; ν_ss_om, ν_ss_quartz, ν_ss_gravel)

soil = Soil.EnergyHydrology{FT}(
    land_domain,
    forcing,
    toml_dict;
    prognostic_land_components,
    additional_sources = (ClimaLand.RootExtraction{FT}(),),
    albedo             = soil_albedo,
    runoff             = ClimaLand.Soil.Runoff.SurfaceRunoff(),
    retention_parameters,
    composition_parameters,
    S_s   = soil_S_s,
    z_0m  = z_0m_soil,
    z_0b  = z_0b_soil,
    emissivity = soil_ϵ,
)

# ── SoilCO2 ────────────────────────────────────────────────────────────────────
co2_prognostic_soil = Soil.Biogeochemistry.PrognosticMet(soil.parameters)
drivers  = Soil.Biogeochemistry.SoilDrivers(co2_prognostic_soil, atmos)
soilco2  = Soil.Biogeochemistry.SoilCO2Model{FT}(land_domain, drivers, toml_dict)

# ── Canopy (PModel) ────────────────────────────────────────────────────────────
radiation_parameters = (;
    Ω,
    G_Function   = CLMGFunction(χl),
    α_PAR_leaf,
    τ_PAR_leaf,
    α_NIR_leaf,
    τ_NIR_leaf,
)
radiative_transfer = Canopy.TwoStreamModel{FT}(
    surface_domain, toml_dict;
    radiation_parameters,
    ϵ_canopy,
)

photosynthesis  = Canopy.PModel{FT}(surface_domain, toml_dict)
conductance     = Canopy.PModelConductance{FT}(toml_dict; Drel)
soil_moisture_stress = Canopy.PiecewiseMoistureStressModel{FT}(
    land_domain, toml_dict;
    soil_params = (; ν = soil_ν, θ_r),
)

hydraulics = Canopy.PlantHydraulicsModel{FT}(
    surface_domain, toml_dict;
    ν              = plant_ν,
    S_s            = plant_S_s,
    conductivity_model,
    retention_model,
)
biomass = Canopy.PrescribedBiomassModel{FT}(; LAI, SAI, RAI, rooting_depth, height = h_canopy)
energy  = Canopy.BigLeafEnergyModel{FT}(toml_dict; ac_canopy)
ground  = ClimaLand.PrognosticGroundConditions{FT}()

canopy = Canopy.CanopyModel{FT}(
    surface_domain,
    (; atmos, radiation, ground),
    LAI,
    toml_dict;
    prognostic_land_components,
    radiative_transfer,
    photosynthesis,
    conductance,
    soil_moisture_stress,
    hydraulics,
    energy,
    biomass,
)

# ── Snow ───────────────────────────────────────────────────────────────────────
snow = Snow.SnowModel(FT, surface_domain, forcing, toml_dict, dt;
    prognostic_land_components)

# ── Integrated land model ──────────────────────────────────────────────────────
land = LandModel{FT}(canopy, snow, soil, soilco2, nothing)

# ── Initial conditions ─────────────────────────────────────────────────────────
# Read initial air temperature from NetCDF to initialise soil and canopy.
T_init_K = NCDataset(met_nc_path, "r") do ds
    Float64(coalesce(ds["Tair"][1, 1, 1], 283.15))
end

function set_ic!(Y, p, t, model)
    FT = eltype(Y.soil.ρe_int)
    # Soil moisture: 95 % of porosity
    θ_l_0 = soil.parameters.θ_r .+
             (soil.parameters.ν .- soil.parameters.θ_r) .* FT(0.95)
    Y.soil.ϑ_l .= θ_l_0
    Y.soil.θ_i .= FT(0)

    ρc_s = ClimaLand.Soil.volumetric_heat_capacity.(
        Y.soil.ϑ_l, Y.soil.θ_i,
        soil.parameters.ρc_ds,
        soil.parameters.earth_param_set,
    )
    Y.soil.ρe_int .= ClimaLand.Soil.volumetric_internal_energy.(
        Y.soil.θ_i, ρc_s, FT(T_init_K),
        soil.parameters.earth_param_set,
    )

    # Snow: none at start
    Y.snow.S   .= FT(0)
    Y.snow.S_l .= FT(0)
    Y.snow.U   .= FT(0)

    # Canopy energy and hydraulics
    Y.canopy.energy.T .= FT(T_init_K)
    ψ_leaf_0 = FT(-2e5 / 9800)
    hydraulics = model.canopy.hydraulics
    S_l_ini = ClimaLand.Canopy.inverse_water_retention_curve.(
        hydraulics.parameters.retention_model,
        ψ_leaf_0,
        hydraulics.parameters.ν,
        hydraulics.parameters.S_s,
    )
    Y.canopy.hydraulics.ϑ_l .= ClimaLand.Canopy.augmented_liquid_fraction.(
        hydraulics.parameters.ν, S_l_ini)

    # SoilCO2
    Y.soilco2.CO2 .= FT(6e-5)
    Y.soilco2.O2  .= FT(0.08)
    Y.soilco2.SOC .= FT(5.0)
end

# ── Diagnostics ────────────────────────────────────────────────────────────────
output_vars = [
    "sif", "ra", "gs", "gpp", "ct", "swu", "lwu",
    "er",  "hr", "et", "msf", "shf", "lhf", "rn",
    "swe", "swc", "tsoil", "si", "nee",
]
diags = ClimaLand.default_diagnostics(
    land,
    start_date;
    output_writer     = ClimaDiagnostics.Writers.DictWriter(),
    output_vars,
    reduction_period  = :halfhourly,
)

# ── Simulation ─────────────────────────────────────────────────────────────────
simulation = LandSimulation(
    start_date, stop_date, dt, land;
    set_ic!,
    updateat    = dt,
    diagnostics = diags,
)
@time solve!(simulation)

# ── Save output ────────────────────────────────────────────────────────────────
using CairoMakie, ClimaAnalysis, GeoMakie, Printf, Statistics
import ClimaLand.LandSimVis as LandSimVis

savedir = joinpath(pkgdir(ClimaLand),
    "experiments/integrated/fluxnet/$(site_ID)/out")
mkpath(savedir)

# Load observed LE, H, NEE from daily flux NetCDF.
# NEE is gC/m²/d → mol CO2/m²/s (÷ 12 g/mol ÷ 86400 s/d).
# LHF (Qle) and SHF (Qh) are already in W/m².
flux_nc_path = joinpath(climaland_dir,
    "DK_Sor", "DK-Sor_daily_aggregated_1997-2013_FLUXNET2015_Flux.nc")
comparison_data = NCDataset(flux_nc_path, "r") do ds
    flux_times = DateTime.(ds["time"][:])
    lhf_obs = Float64.(coalesce.(ds["Qle_daily"][:], NaN))
    shf_obs = Float64.(coalesce.(ds["Qh_daily"][:],  NaN))
    nee_obs = Float64.(coalesce.(ds["NEE_daily"][:], NaN)) ./ (12.0 * 86400.0)
    mask = (flux_times .>= start_date) .& (flux_times .< stop_date)
    (; UTC_datetime = flux_times[mask],
       lhf = lhf_obs[mask], shf = shf_obs[mask], nee = nee_obs[mask])
end

LandSimVis.make_timeseries(
    land_domain, diags, start_date;
    savedir,
    short_names  = ["nee", "gpp", "shf", "lhf", "swu", "lwu"],
    spinup_date  = start_date + Day(20),
    comparison_data,
)
LandSimVis.make_timeseries(
    land_domain, diags, start_date;
    savedir,
    short_names  = ["swc", "tsoil", "swe"],
    spinup_date  = start_date + Day(20),
)
LandSimVis.make_timeseries(
    land_domain, diags, start_date;
    savedir,
    short_names  = ["hr"],
    spinup_date  = start_date + Day(20),
)

println("Finished. Output in: $savedir")
