"""
Per-year DK-Sor model interface for CalLMIP Phase 1b output generation.

Runs each calendar year 1997–2014 as a separate single-year simulation with a
60-day spinup (PeriodicCalendar recycles the 1997–2014 forcing into the spinup
window), then concatenates the daily output. This mirrors the validated
reference setup (ar/em/calibrate_neon:experiments/integrated/fluxnet/
run_dk_sor_default.jl, which runs single years with a 60-day spinup).

Why per-year, not one continuous 1996→2015 run: the default DK-Sor config at
DT=900s carries a latent numerical instability — on a stiff cold-season /
wet-after-dry transient the coupled soil-energy solve can momentarily produce
T_soil < 0 K, and the soil-CO2 gas-diffusivity term D0 ∝ (T_soil/T_ref)^1.75
then raises a negative base to a fractional power → DomainError
(gas_diffusivity_in_soil, co2_parameterizations.jl). This is NOT fixed by a
smaller timestep, a lower initial saturation, or a different forcing
interpolation (all tested). Single-year runs start clean each year, so a bad
transient in one year cannot propagate, exactly as the reference run avoids it.

Provides:
  run_one_year(year, toml_dict, met_nc_path) -> (surface_data, column_data,
                                                 z_soil, dates)
  save_callmip_diagnostics(all_surface, all_column, z_soil, all_dates, outdir)

Usage (called from run_prior_simulation.jl / run_callmip_simulations.jl):
  include("callmip_model_interface.jl")
"""

import ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.FluxnetSimulations as FluxnetSimulations
import ClimaLand.Simulations: LandSimulation, solve!
using ClimaLand
using ClimaLand.Domains: Column
using ClimaLand.Soil
using ClimaLand.Soil.Biogeochemistry
using ClimaLand.Canopy
using ClimaLand.Snow
using ClimaCore
using ClimaDiagnostics
using ClimaUtilities
import ClimaUtilities.TimeVaryingInputs:
    TimeVaryingInput, LinearInterpolation, PeriodicCalendar, evaluate!
import ClimaUtilities.TimeManager: date
using Dates
using NCDatasets
using Statistics
using JLD2

# ── CalLMIP output period ─────────────────────────────────────────────────────
const SPINUP_DAYS  = 60
const OUTPUT_START = Date(1997, 1, 1)
const OUTPUT_STOP  = Date(2014, 12, 31)
const OUTPUT_YEARS = year(OUTPUT_START):year(OUTPUT_STOP)   # 1997:2014

# CalLMIP surface diagnostic names available for LandModel.
# Note: "soillhf", "soilrn", "soilshf" are NOT in get_possible_diagnostics(LandModel).
const CALLMIP_SURFACE_VARS = [
    "nee", "lhf", "shf", "gpp", "er", "trans", "ct", "lai", "cveg",
]
const CALLMIP_COLUMN_VARS = ["swc", "tsoil", "soc"]

# ─────────────────────────────────────────────────────────────────────────────
# Forcing loader  (sim_start passed in — different per year)
# ─────────────────────────────────────────────────────────────────────────────

function load_forcing_periodic(met_nc_path::String, toml_dict, FT, sim_start::DateTime)
    loc    = FluxnetSimulations.get_location(FT, Val(:DK_Sor))
    height = FluxnetSimulations.get_fluxtower_height(FT, Val(:DK_Sor))
    lat    = loc.lat; long = loc.long; time_offset = loc.time_offset
    atmos_h = height.atmos_h
    method  = LinearInterpolation(PeriodicCalendar())

    earth_param_set = LP.LandParameters(toml_dict)

    function snow_frac(T_K, VPD_hPa)
        T_K - 273.15 < 0.0 ? 1.0 : 0.0
    end

    NCDataset(met_nc_path, "r") do ds
        t_dates = ds["time"][:]
        t_secs  = Float64[
            Second(t - Hour(time_offset) - sim_start).value
            for t in t_dates
        ]

        function tvi(var; scale = 1.0, offset = 0.0)
            raw   = Float64.(coalesce.(ds[var][1, 1, :], NaN))
            valid = .!isnan.(raw) .& (raw .!= -9999.0)
            TimeVaryingInput(t_secs[valid], (raw[valid] .* scale .+ offset); method)
        end

        T_arr    = Float64.(coalesce.(ds["Tair"][1, 1, :],   NaN))
        VPD_arr  = Float64.(coalesce.(ds["VPD"][1, 1, :],    NaN))
        prec_arr = Float64.(coalesce.(ds["Precip"][1, 1, :], NaN))
        valid_p  = .!isnan.(T_arr) .& .!isnan.(prec_arr)

        snow_frac_arr = snow_frac.(T_arr, VPD_arr)
        rain_arr = @. -abs(prec_arr) * (1 - snow_frac_arr)
        snow_arr = @. -abs(prec_arr) * snow_frac_arr

        atmos_T  = tvi("Tair")
        atmos_u  = tvi("Wind")
        atmos_q  = tvi("Qair")
        atmos_P  = tvi("Psurf")
        SW_d     = tvi("SWdown")
        LW_d     = tvi("LWdown")
        c_co2    = tvi("CO2air"; scale = 1e-6)
        valid_idx = findall(valid_p)
        atmos_Pr = TimeVaryingInput(t_secs[valid_idx], rain_arr[valid_idx]; method)
        atmos_Ps = TimeVaryingInput(t_secs[valid_idx], snow_arr[valid_idx]; method)

        atmos = ClimaLand.PrescribedAtmosphere(
            atmos_Pr, atmos_Ps, atmos_T, atmos_u, atmos_q, atmos_P,
            sim_start, atmos_h, toml_dict; c_co2,
        )
        cos_zenith_angle = (t, s) -> ClimaLand.default_cos_zenith_angle(
            t, s;
            insol_params = earth_param_set.insol_params,
            longitude    = long,
            latitude     = lat,
        )
        radiation = ClimaLand.PrescribedRadiativeFluxes(
            FT, SW_d, LW_d, sim_start; cosθs = cos_zenith_angle, toml_dict,
        )
        return (; atmos, radiation)
    end
end

function load_lai_periodic(met_nc_path::String, FT, sim_start::DateTime)
    loc    = FluxnetSimulations.get_location(FT, Val(:DK_Sor))
    method = LinearInterpolation(PeriodicCalendar())
    NCDataset(met_nc_path, "r") do ds
        t_dates = ds["time"][:]
        t_secs  = Float64[
            Second(t - Hour(loc.time_offset) - sim_start).value
            for t in t_dates
        ]
        raw   = Float64.(coalesce.(ds["LAI_alternative"][1, 1, :], NaN))
        valid = .!isnan.(raw)
        maxLAI = maximum(raw[valid])
        TimeVaryingInput(t_secs[valid], raw[valid]; method), maxLAI
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Model builder
# ─────────────────────────────────────────────────────────────────────────────

function build_callmip_model(FT, toml_dict, met_nc_path::String, sim_start::DateTime, DT)
    forcing     = load_forcing_periodic(met_nc_path, toml_dict, FT, sim_start)
    LAI, maxLAI = load_lai_periodic(met_nc_path, FT, sim_start)

    (; dz_tuple, nelements, zmin, zmax) =
        FluxnetSimulations.get_domain_info(FT, Val(:DK_Sor))
    loc  = FluxnetSimulations.get_location(FT, Val(:DK_Sor))
    land_domain = Column(;
        zlim    = (FT(zmin), FT(zmax)),
        nelements, dz_tuple,
        longlat = (loc.long, loc.lat),
    )
    canopy_domain = ClimaLand.Domains.obtain_surface_domain(land_domain)

    (; soil_ν, θ_r, soil_K_sat, soil_vg_n, soil_vg_α,
       ν_ss_om, ν_ss_quartz, ν_ss_gravel,
       z_0m_soil, z_0b_soil, soil_α_PAR, soil_α_NIR, soil_ϵ,
       Ω, χl, α_PAR_leaf, τ_PAR_leaf, α_NIR_leaf, τ_NIR_leaf,
       ϵ_canopy, ac_canopy, Drel, g0,
       SAI, f_root_to_shoot, rooting_depth, h_canopy,
       conductivity_model, retention_model, plant_ν, plant_S_s) =
        FluxnetSimulations.get_parameters(FT, Val(:DK_Sor))

    RAI = maxLAI * f_root_to_shoot
    prognostic_land_components = (:canopy, :snow, :soil, :soilco2)

    soil_albedo = Soil.ConstantTwoBandSoilAlbedo{FT}(;
        PAR_albedo = soil_α_PAR, NIR_albedo = soil_α_NIR)
    retention_parameters = (;
        ν = soil_ν, θ_r, K_sat = soil_K_sat,
        hydrology_cm = vanGenuchten{FT}(; α = soil_vg_α, n = soil_vg_n),
    )
    composition_parameters = (; ν_ss_om, ν_ss_quartz, ν_ss_gravel)

    soil = Soil.EnergyHydrology{FT}(
        land_domain, forcing, toml_dict;
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
        land_domain, co2_drivers, toml_dict,
    )
    radiation_parameters = (;
        Ω, G_Function = CLMGFunction(χl),
        α_PAR_leaf, τ_PAR_leaf, α_NIR_leaf, τ_NIR_leaf,
    )
    radiative_transfer = Canopy.TwoStreamModel{FT}(
        canopy_domain, toml_dict; radiation_parameters, ϵ_canopy,
    )
    biomass = Canopy.PrescribedBiomassModel{FT}(
        land_domain, LAI, toml_dict;
        SAI, RAI, rooting_depth, height = h_canopy,
    )
    canopy = Canopy.CanopyModel{FT}(
        canopy_domain,
        (; atmos = forcing.atmos, radiation = forcing.radiation,
           ground = ClimaLand.PrognosticGroundConditions{FT}()),
        LAI, toml_dict;
        prognostic_land_components,
        radiative_transfer,
        photosynthesis       = Canopy.PModel{FT}(canopy_domain, toml_dict),
        conductance          = Canopy.PModelConductance{FT}(toml_dict; Drel),
        soil_moisture_stress = Canopy.PiecewiseMoistureStressModel{FT}(
            land_domain, toml_dict; soil_params = (; ν = soil_ν, θ_r)),
        hydraulics = Canopy.PlantHydraulicsModel{FT}(
            canopy_domain, toml_dict;
            ν = plant_ν, S_s = plant_S_s, conductivity_model, retention_model),
        energy  = Canopy.BigLeafEnergyModel{FT}(toml_dict; ac_canopy),
        biomass,
    )
    snow = Snow.SnowModel(
        FT, canopy_domain, forcing, toml_dict, DT; prognostic_land_components,
    )
    land = LandModel{FT}(canopy, snow, soil, soilco2, nothing)
    return land, forcing, soil_ν, θ_r
end

function make_callmip_ic(ν, θ_r, atmos)
    function set_ic!(Y, p, t, model)
        FT_l = eltype(Y.soil.ρe_int)
        evaluate!(p.drivers.T, atmos.T, t)

        Y.soil.ϑ_l .= FT_l(θ_r) + (FT_l(ν) - FT_l(θ_r)) * FT_l(0.95)
        Y.soil.θ_i .= FT_l(0)
        ρc_s = ClimaLand.Soil.volumetric_heat_capacity.(
            Y.soil.ϑ_l, Y.soil.θ_i,
            model.soil.parameters.ρc_ds, model.soil.parameters.earth_param_set)
        Y.soil.ρe_int .= ClimaLand.Soil.volumetric_internal_energy.(
            Y.soil.θ_i, ρc_s, p.drivers.T, model.soil.parameters.earth_param_set)

        Y.snow.S .= FT_l(0); Y.snow.S_l .= FT_l(0); Y.snow.U .= FT_l(0)

        if model.canopy.energy isa ClimaLand.Canopy.BigLeafEnergyModel
            Y.canopy.energy.T .= p.drivers.T
        end
        Y.canopy.hydraulics.ϑ_l .= model.canopy.hydraulics.parameters.ν

        if !isnothing(model.soilco2)
            Y.soilco2.CO2 .= FT_l(6e-5)
            Y.soilco2.O2  .= FT_l(0.21)
            SOC_top = FT_l(15.0); SOC_bot = FT_l(0.5)
            τ_soc   = FT_l(1.0 / log(SOC_top / SOC_bot))
            z = ClimaCore.Fields.coordinate_field(axes(Y.soilco2.SOC)).z
            @. Y.soilco2.SOC = SOC_bot + (SOC_top - SOC_bot) * exp(z / τ_soc)
        end
    end
    return set_ic!
end

# ─────────────────────────────────────────────────────────────────────────────
# Diagnostic extraction
# ─────────────────────────────────────────────────────────────────────────────

function extract_daily_diag(simulation, diag_name)
    for d in simulation.diagnostics
        if haskey(d.output_writer.dict, diag_name)
            (times, data) = ClimaLand.Diagnostics.diagnostic_as_vectors(
                d.output_writer, diag_name)
            dates = Date.(times isa Vector{DateTime} ? times : date.(times))
            return dates, Float64.(data)
        end
    end
    error("Diagnostic '$diag_name' not found.")
end

function extract_daily_column_diag(simulation, diag_name)
    for d in simulation.diagnostics
        if haskey(d.output_writer.dict, diag_name)
            raw_dict  = d.output_writer.dict[diag_name]
            times_raw = sort(collect(keys(raw_dict)))
            dates     = Date.(times_raw isa Vector{DateTime} ? times_raw : date.(times_raw))
            n_days    = length(dates)
            first_arr = vec(parent(raw_dict[times_raw[1]]))
            n_z = length(first_arr)
            mat = Matrix{Float64}(undef, n_z, n_days)
            mat[:, 1] = Float64.(first_arr)
            for t in 2:n_days
                mat[:, t] = Float64.(vec(parent(raw_dict[times_raw[t]])))
            end
            return dates, mat
        end
    end
    error("Diagnostic '$diag_name' not found.")
end

# ─────────────────────────────────────────────────────────────────────────────
# Run one calendar year with a 60-day spinup
# ─────────────────────────────────────────────────────────────────────────────

"""
    run_one_year(yr, toml_dict, met_nc_path; FT=Float64, DT=900.0)
        -> (surface_data, column_data, z_soil, dates)

Run calendar year `yr` with a `SPINUP_DAYS`-day spinup preceding it. Returns
daily diagnostics for the output year only (spinup discarded). `surface_data`
maps each surface var to a Vector over the year's days; `column_data` maps each
column var to an (n_z × n_days) Matrix; `dates` are the kept output days.
"""
function run_one_year(yr::Int, toml_dict, met_nc_path::String;
                      FT = Float64, DT = Float64(900))
    sim_start = DateTime(yr, 1, 1) - Day(SPINUP_DAYS)
    sim_stop  = DateTime(yr + 1, 1, 1)
    out_start = Date(yr, 1, 1)
    out_stop  = Date(yr, 12, 31)

    land, forcing, ν, θ_r =
        build_callmip_model(FT, toml_dict, met_nc_path, sim_start, DT)
    set_ic! = make_callmip_ic(ν, θ_r, forcing.atmos)

    output_writer = ClimaDiagnostics.Writers.DictWriter()
    diags_surface = ClimaLand.default_diagnostics(
        land, sim_start, "";
        output_writer, output_vars = CALLMIP_SURFACE_VARS,
        reduction_period = :daily,
    )
    diags_col = ClimaLand.default_diagnostics(
        land, sim_start, "";
        output_writer, output_vars = CALLMIP_COLUMN_VARS,
        reduction_period = :daily,
    )

    simulation = LandSimulation(
        sim_start, sim_stop, DT, land;
        set_ic!, updateat = Second(DT),
        diagnostics = vcat(diags_surface, diags_col),
    )
    solve!(simulation)

    surface_data = Dict{String, Vector{Float64}}()
    ref_dates    = Date[]
    for var in CALLMIP_SURFACE_VARS
        key = "$(var)_1d_average"
        try
            dates, vals = extract_daily_diag(simulation, key)
            keep = (dates .>= out_start) .& (dates .<= out_stop)
            surface_data[var] = vals[keep]
            isempty(ref_dates) && (ref_dates = dates[keep])
        catch
            @warn "  [$yr] Surface diagnostic '$var' not found"
        end
    end

    column_data = Dict{String, Matrix{Float64}}()
    for var in CALLMIP_COLUMN_VARS
        key = "$(var)_1d_average"
        try
            dates, mat = extract_daily_column_diag(simulation, key)
            keep = (dates .>= out_start) .& (dates .<= out_stop)
            column_data[var] = mat[:, keep]
        catch
            @warn "  [$yr] Column diagnostic '$var' not found"
        end
    end

    z_soil = try
        Float64.(vec(parent(
            ClimaCore.Fields.coordinate_field(
                axes(land.soil.domain.fields.z)).z
        )))
    catch
        Float64[]
    end

    return surface_data, column_data, z_soil, ref_dates
end

# ─────────────────────────────────────────────────────────────────────────────
# Save concatenated multi-year diagnostics to JLD2
# ─────────────────────────────────────────────────────────────────────────────

"""
    save_callmip_diagnostics(all_surface, all_column, z_soil, all_dates, outdir)

Save concatenated 1997–2014 diagnostics to `outdir/callmip_diagnostics.jld2`.
`all_surface` maps var → full Vector over all kept days; `all_column` maps
var → (n_z × n_days) Matrix; `all_dates` is the full Vector of days.
"""
function save_callmip_diagnostics(
    all_surface::Dict{String, Vector{Float64}},
    all_column::Dict{String, Matrix{Float64}},
    z_soil::Vector{Float64},
    all_dates::Vector{Date},
    outdir::String,
)
    mkpath(outdir)
    n = length(all_dates)
    @info "Saving $n days of CalLMIP diagnostics " *
          "($(isempty(all_dates) ? "—" : "$(all_dates[1]) – $(all_dates[end]))")"
    jldsave(joinpath(outdir, "callmip_diagnostics.jld2");
        dates        = all_dates,
        surface_data = all_surface,
        column_data  = all_column,
        z_soil,
    )
end
