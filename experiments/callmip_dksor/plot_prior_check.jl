"""
Plot model prior vs. FLUXNET observations for DK-Sor (3 representative years).

Runs ClimaLand with default parameters for years 2005, 2008, 2010 and plots
monthly NEE, LHF, SHF against FLUXNET observations with RMSE and R².

Output: experiments/callmip_dksor/output_prior_check/prior_vs_obs.png
"""

import ClimaComms
ClimaComms.@import_required_backends

using ClimaLand
using ClimaLand.Domains: Column
using ClimaLand.Snow
using ClimaLand.Soil
using ClimaLand.Soil.Biogeochemistry
using ClimaLand.Canopy
using ClimaCore
using ClimaDiagnostics
using ClimaUtilities
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.Simulations: LandSimulation, solve!
import ClimaLand.FluxnetSimulations as FluxnetSimulations
import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput
import ClimaUtilities.TimeManager: date
import EnsembleKalmanProcesses as EKP     # needed before loading observations.jld2
using NCDatasets
using Dates
using Statistics
using LinearAlgebra
using JLD2
using CairoMakie

# ── Configuration ──────────────────────────────────────────────────────────────
const FT          = Float64
const DT          = Float64(900)
const SITE_ID_VAL = :DK_Sor
const CLIMALAND_DIR = pkgdir(ClimaLand)
const MET_NC_PATH   = joinpath(CLIMALAND_DIR, "DK_Sor",
    "DK-Sor_1997-2014_FLUXNET2015_Met.nc")
const FLUX_NC_PATH  = joinpath(CLIMALAND_DIR, "DK_Sor",
    "DK-Sor_daily_aggregated_1997-2013_FLUXNET2015_Flux.nc")
const OUTDIR        = joinpath(CLIMALAND_DIR,
    "experiments/callmip_dksor/output_prior_check")

const CHECK_YEARS = [2005, 2008, 2010]

const _SITE_LOC    = FluxnetSimulations.get_location(FT, Val(SITE_ID_VAL))
const _SITE_HEIGHT = FluxnetSimulations.get_fluxtower_height(FT, Val(SITE_ID_VAL))
const _SITE_PARAMS = FluxnetSimulations.get_parameters(FT, Val(SITE_ID_VAL))

const lat         = _SITE_LOC.lat
const long        = _SITE_LOC.long
const time_offset = _SITE_LOC.time_offset
const atmos_h     = _SITE_HEIGHT.atmos_h

"""Run one year with default parameters, return (nee_monthly, lhf_monthly, shf_monthly)."""
function run_prior_year(year)
    start_date = DateTime(year, 1, 1)
    stop_date  = DateTime(year + 1, 1, 1)

    local_toml = LP.create_toml_dict(FT)   # no overrides = default parameters

    (; dz_tuple, nelements, zmin, zmax) =
        FluxnetSimulations.get_domain_info(FT, Val(SITE_ID_VAL))
    land_domain = Column(;
        zlim      = (FT(zmin), FT(zmax)),
        nelements, dz_tuple,
        longlat   = (long, lat),
    )
    canopy_domain = ClimaLand.Domains.obtain_surface_domain(land_domain)

    forcing = FluxnetSimulations.prescribed_forcing_netcdf(
        MET_NC_PATH, lat, long, time_offset, atmos_h,
        start_date, local_toml, FT,
    )

    LAI, maxLAI = NCDataset(MET_NC_PATH, "r") do ds
        time_vals = ds["time"][:]
        lai_data  = Float64.(coalesce.(ds["LAI_alternative"][1, 1, :], NaN))
        lai_secs  = Float64[
            Second(t - Hour(time_offset) - start_date).value
            for t in time_vals
        ]
        valid = .!isnan.(lai_data)
        TimeVaryingInput(lai_secs[valid], lai_data[valid]), maximum(lai_data[valid])
    end

    (; soil_ν, θ_r, soil_K_sat, soil_vg_n, soil_vg_α,
       ν_ss_om, ν_ss_quartz, ν_ss_gravel,
       z_0m_soil, z_0b_soil, soil_α_PAR, soil_α_NIR, soil_ϵ,
       Ω, χl, α_PAR_leaf, τ_PAR_leaf, α_NIR_leaf, τ_NIR_leaf,
       ϵ_canopy, ac_canopy, Drel, g0,
       SAI, f_root_to_shoot, rooting_depth, h_canopy,
       conductivity_model, retention_model, plant_ν, plant_S_s) = _SITE_PARAMS

    RAI = maxLAI * f_root_to_shoot

    prognostic_land_components = (:canopy, :snow, :soil, :soilco2)

    soil_albedo = Soil.ConstantTwoBandSoilAlbedo{FT}(;
        PAR_albedo = soil_α_PAR, NIR_albedo = soil_α_NIR)
    retention_parameters = (;
        ν      = soil_ν,
        θ_r,
        K_sat  = soil_K_sat,
        hydrology_cm = vanGenuchten{FT}(; α = soil_vg_α, n = soil_vg_n),
    )
    composition_parameters = (; ν_ss_om, ν_ss_quartz, ν_ss_gravel)

    soil = Soil.EnergyHydrology{FT}(
        land_domain, forcing, local_toml;
        prognostic_land_components,
        additional_sources   = (ClimaLand.RootExtraction{FT}(),),
        albedo               = soil_albedo,
        runoff               = ClimaLand.Soil.Runoff.SurfaceRunoff(),
        retention_parameters, composition_parameters,
        S_s = FT(1e-3), z_0m = z_0m_soil, z_0b = z_0b_soil,
    )

    co2_drivers = Soil.Biogeochemistry.SoilDrivers(
        Soil.Biogeochemistry.PrognosticMet(soil.parameters), forcing.atmos,
    )
    soilco2 = Soil.Biogeochemistry.SoilCO2Model{FT}(
        land_domain, co2_drivers, local_toml,
    )

    radiation_parameters = (;
        Ω, G_Function = CLMGFunction(χl),
        α_PAR_leaf, τ_PAR_leaf, α_NIR_leaf, τ_NIR_leaf,
    )
    radiative_transfer = Canopy.TwoStreamModel{FT}(
        canopy_domain, local_toml; radiation_parameters, ϵ_canopy,
    )

    biomass = Canopy.PrescribedBiomassModel{FT}(
        land_domain, LAI, local_toml;
        SAI = SAI, RAI = RAI, rooting_depth = rooting_depth,
        height = h_canopy,
    )
    canopy = Canopy.CanopyModel{FT}(
        canopy_domain,
        (; atmos = forcing.atmos, radiation = forcing.radiation,
           ground = ClimaLand.PrognosticGroundConditions{FT}()),
        LAI, local_toml;
        prognostic_land_components,
        radiative_transfer,
        photosynthesis       = Canopy.PModel{FT}(canopy_domain, local_toml),
        conductance          = Canopy.PModelConductance{FT}(local_toml; Drel),
        soil_moisture_stress = Canopy.PiecewiseMoistureStressModel{FT}(
            land_domain, local_toml; soil_params = (; ν = soil_ν, θ_r)),
        hydraulics = Canopy.PlantHydraulicsModel{FT}(
            canopy_domain, local_toml;
            ν = plant_ν, S_s = plant_S_s, conductivity_model, retention_model),
        energy  = Canopy.BigLeafEnergyModel{FT}(local_toml; ac_canopy),
        biomass,
    )
    snow = Snow.SnowModel(
        FT, canopy_domain, forcing, local_toml, DT;
        prognostic_land_components,
    )
    land = LandModel{FT}(canopy, snow, soil, soilco2, nothing)

    T_init_K = NCDataset(MET_NC_PATH, "r") do ds
        t_dates = DateTime.(ds["time"][:])
        idx_yr  = findfirst(t -> Dates.year(t) == year, t_dates)
        isnothing(idx_yr) ? 283.15 : Float64(coalesce(ds["Tair"][1, 1, idx_yr], 283.15))
    end

    function set_ic!(Y, p, t, model)
        FT_l = eltype(Y.soil.ρe_int)
        ν_l  = FT_l(soil_ν)
        θr_l = FT_l(θ_r)
        Y.soil.ϑ_l .= θr_l + (ν_l - θr_l) * FT_l(0.95)
        Y.soil.θ_i .= FT_l(0)
        ρc_s = ClimaLand.Soil.volumetric_heat_capacity.(
            Y.soil.ϑ_l, Y.soil.θ_i,
            soil.parameters.ρc_ds, soil.parameters.earth_param_set)
        Y.soil.ρe_int .= ClimaLand.Soil.volumetric_internal_energy.(
            Y.soil.θ_i, ρc_s, FT_l(T_init_K),
            soil.parameters.earth_param_set)
        Y.snow.S .= FT_l(0); Y.snow.S_l .= FT_l(0); Y.snow.U .= FT_l(0)
        Y.canopy.energy.T .= FT_l(T_init_K)
        Y.canopy.hydraulics.ϑ_l .= model.canopy.hydraulics.parameters.ν
        Y.soilco2.CO2 .= FT_l(6e-5)
        Y.soilco2.O2  .= FT_l(0.21)
        SOC_top = FT_l(15.0); SOC_bot = FT_l(0.5)
        τ_soc   = FT_l(1.0 / log(SOC_top / SOC_bot))
        z = ClimaCore.Fields.coordinate_field(axes(Y.soilco2.SOC)).z
        @. Y.soilco2.SOC = SOC_bot + (SOC_top - SOC_bot) * exp(z / τ_soc)
    end

    diags = ClimaLand.default_diagnostics(
        land, start_date, "";
        output_writer    = ClimaDiagnostics.Writers.DictWriter(),
        output_vars      = ["lhf", "shf", "nee"],
        reduction_period = :daily,
    )
    simulation = LandSimulation(
        start_date, stop_date, DT, land;
        set_ic!, updateat = Second(DT), diagnostics = diags,
    )

    solve!(simulation)

    writer = simulation.diagnostics[1].output_writer
    function get_monthly(diag_name, scale = 1.0)
        times, data = ClimaLand.Diagnostics.diagnostic_as_vectors(writer, diag_name)
        months = Dates.month.(date.(times))
        return [mean(data[months .== m]) * scale for m in 1:12]
    end

    nee_m = get_monthly("nee_1d_average", 12.0 * 86400.0)  # mol/m²/s → gC/m²/d
    lhf_m = get_monthly("lhf_1d_average")
    shf_m = get_monthly("shf_1d_average")
    return nee_m, lhf_m, shf_m
end

"""Load monthly obs for a given year from the FLUXNET daily NC file."""
function load_obs_year(year)
    NCDataset(FLUX_NC_PATH, "r") do ds
        dates   = DateTime.(ds["time"][:])
        nee_all = Float64.(coalesce.(ds["NEE_daily"][:], NaN))
        qle_all = Float64.(coalesce.(ds["Qle_daily"][:], NaN))
        qh_all  = Float64.(coalesce.(ds["Qh_daily"][:],  NaN))
        for arr in (nee_all, qle_all, qh_all)
            arr[arr .>= 1.0e19] .= NaN
        end

        nee_m = Float64[]; lhf_m = Float64[]; shf_m = Float64[]
        for mon in 1:12
            mask = (Dates.year.(dates) .== year) .& (Dates.month.(dates) .== mon)
            nee_v = nee_all[mask]; qle_v = qle_all[mask]; qh_v = qh_all[mask]
            nee_f = filter(isfinite, nee_v); push!(nee_m, isempty(nee_f) ? NaN : mean(nee_f))
            lhf_f = filter(isfinite, qle_v); push!(lhf_m, isempty(lhf_f) ? NaN : mean(lhf_f))
            shf_f = filter(isfinite, qh_v);  push!(shf_m, isempty(shf_f) ? NaN : mean(shf_f))
        end
        return nee_m, lhf_m, shf_m
    end
end

function rmse_r2(mod, obs)
    mask = isfinite.(mod) .& isfinite.(obs)
    sum(mask) < 2 && return NaN, NaN
    m = mod[mask]; o = obs[mask]
    rmse = sqrt(mean((m .- o).^2))
    r2   = 1.0 - sum((m .- o).^2) / sum((o .- mean(o)).^2)
    return rmse, r2
end

# ── Run simulations ────────────────────────────────────────────────────────────
mkpath(OUTDIR)

results_mod = Dict{Int, NamedTuple}()
results_obs = Dict{Int, NamedTuple}()

for yr in CHECK_YEARS
    @info "Running year $yr with default parameters..."
    nee_m, lhf_m, shf_m = run_prior_year(yr)
    results_mod[yr] = (; nee = nee_m, lhf = lhf_m, shf = shf_m)
    @info "  Done year $yr"

    nee_o, lhf_o, shf_o = load_obs_year(yr)
    results_obs[yr] = (; nee = nee_o, lhf = lhf_o, shf = shf_o)
end

# ── Aggregate metrics across all 3 years ─────────────────────────────────────
function aggregate_metric(var_key)
    mod_all = vcat([results_mod[yr][var_key] for yr in CHECK_YEARS]...)
    obs_all = vcat([results_obs[yr][var_key] for yr in CHECK_YEARS]...)
    return rmse_r2(mod_all, obs_all)
end

nee_rmse, nee_r2 = aggregate_metric(:nee)
lhf_rmse, lhf_r2 = aggregate_metric(:lhf)
shf_rmse, shf_r2 = aggregate_metric(:shf)

@info "Prior performance (2005/2008/2010 combined):"
@info "  NEE: RMSE=$(round(nee_rmse; digits=2)) gC/m²/d, R²=$(round(nee_r2; digits=3))"
@info "  LHF: RMSE=$(round(lhf_rmse; digits=2)) W/m², R²=$(round(lhf_r2; digits=3))"
@info "  SHF: RMSE=$(round(shf_rmse; digits=2)) W/m², R²=$(round(shf_r2; digits=3))"

# ── Plot ───────────────────────────────────────────────────────────────────────
months = 1:12
var_labels = ["NEE (gC m⁻² d⁻¹)", "LHF (W m⁻²)", "SHF (W m⁻²)"]
var_keys   = [:nee, :lhf, :shf]
metrics    = [(nee_rmse, nee_r2), (lhf_rmse, lhf_r2), (shf_rmse, shf_r2)]

fig = Figure(size = (1000, 800))

for (row, (vk, ylabel, (rmse, r2))) in enumerate(zip(var_keys, var_labels, metrics))
    for (col, yr) in enumerate(CHECK_YEARS)
        ax = Axis(fig[row, col];
            xlabel = col == 2 ? "Month" : "",
            ylabel = col == 1 ? ylabel : "",
            title  = row == 1 ? string(yr) : "",
            xticks = (1:12, string.(1:12)),
        )

        obs_v = results_obs[yr][vk]
        mod_v = results_mod[yr][vk]

        # Obs: black dots for valid months
        valid_months = findall(isfinite, obs_v)
        scatter!(ax, valid_months, obs_v[valid_months];
            color = :black, markersize = 8, label = "FLUXNET")

        # Model: colored line
        lines!(ax, months, mod_v; color = Makie.wong_colors()[row], linewidth = 2,
            label = "ClimaLand prior")

        if col == 1 && row == 1
            axislegend(ax; position = :rt, labelsize = 10)
        end
    end
end

# Overall title with metrics
Label(fig[0, :],
    "DK-Sor prior (default params) vs FLUXNET — 2005/2008/2010\n" *
    "NEE: RMSE=$(round(nee_rmse; digits=2)) gC/m²/d, R²=$(round(nee_r2; digits=3))  |  " *
    "LHF: RMSE=$(round(lhf_rmse; digits=2)) W/m², R²=$(round(lhf_r2; digits=3))  |  " *
    "SHF: RMSE=$(round(shf_rmse; digits=2)) W/m², R²=$(round(shf_r2; digits=3))",
    fontsize = 13, font = :bold,
)

save(joinpath(OUTDIR, "prior_vs_obs.png"), fig; px_per_unit = 2)
@info "Plot saved → $(joinpath(OUTDIR, "prior_vs_obs.png"))"
