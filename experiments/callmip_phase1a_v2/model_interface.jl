"""
ClimaCalibrate forward model for CalLMIP Phase 1a DK-Sor EKI calibration.

Runs a single-column ClimaLand simulation for 2 randomly-drawn yearly windows
(selected by EKP minibatcher) and returns NEE + LE + H daily diagnostics.

Current-main API notes:
  - Y.soilco2.O2 (not O2_f); initial value 0.08 kg O2/m³ (not 0.21 vol fraction)
  - Y.soilco2.CO2 initial: 6e-5 kg C/m³
  - Y.soilco2.SOC: prognostic (initialized via SOC profile below)
  - PlantHydraulics merged into Canopy — no separate import
  - PlantHydraulicsModel{FT}(domain, toml_dict): single above-ground compartment
  - PrescribedBiomassModel: height / SAI / RAI via toml_dict or kwargs
  - Timestep: 900 s
  - No make_imp_tendency / make_update_jacobian calls

Constants expected to be defined in the calling script (broadcast via @everywhere):
  OUTPUT_DIR, OBS_FILEPATH, SITE_ID, DT
"""

import ClimaCalibrate
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.Simulations: LandSimulation, solve!
using ClimaLand
using ClimaLand.Domains: Column
using ClimaLand.Soil
using ClimaLand.Soil.Biogeochemistry
using ClimaLand.Canopy
using ClimaCore
using ClimaDiagnostics
using ClimaUtilities
import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput, evaluate!, LinearInterpolation, Flat
import ClimaUtilities.TimeManager: date

import EnsembleKalmanProcesses as EKP
import JLD2
using Dates
using NCDatasets
using Statistics

# ── Site constants ────────────────────────────────────────────────────────────
const _SITE_LAT      = 55.4859   # degrees N
const _SITE_LON      = 11.6446   # degrees E
const _UTC_OFFSET    = 1         # local standard time UTC+1
const _ATMOS_H       = 32.0      # m, sensor height

# ── Soil domain parameters for DK-Sor ────────────────────────────────────────
const _SOIL_ν        = 0.45      # porosity (m³/m³)
const _SOIL_θ_r      = 0.07      # residual water content
const _SOIL_K_sat    = 1e-5      # saturated hydraulic conductivity (m/s)
const _SOIL_VG_N     = 1.6       # van Genuchten n
const _SOIL_VG_α     = 1.6       # van Genuchten α (1/m)
const _SOIL_ν_ss_om  = 0.03
const _SOIL_ν_ss_q   = 0.47
const _SOIL_ν_ss_g   = 0.12

# ── Canopy radiative parameters ───────────────────────────────────────────────
const _χl            = 0.25
const _α_PAR_leaf    = 0.1
const _α_NIR_leaf    = 0.45
const _τ_PAR_leaf    = 0.05
const _τ_NIR_leaf    = 0.25
const _Ω             = 1.0

# ── Helpers ───────────────────────────────────────────────────────────────────

"""Load prescribed forcing from the DK-Sor NetCDF Met file."""
function load_dk_sor_forcing(
    met_nc_path::String,
    sim_start::DateTime,
    toml_dict,
    FT,
)
    lat    = FT(_SITE_LAT)
    lon    = FT(_SITE_LON)
    atmos_h = FT(_ATMOS_H)
    time_offset = _UTC_OFFSET   # hours

    ds = NCDataset(met_nc_path, "r")
    met_times = ds["time"][:]
    close(ds)

    # Convert met times → seconds since sim_start
    t_seconds = [Float64(Second(t - Hour(time_offset) - sim_start).value) for t in met_times]

    function make_tvi(varname; scale = 1.0, offset = 0.0)
        ds = NCDataset(met_nc_path, "r")
        raw = Float64.(coalesce.(ds[varname][1, 1, :], NaN))
        close(ds)
        valid = .!isnan.(raw)
        TimeVaryingInput(t_seconds[valid], (raw[valid] .* scale .+ offset))
    end

    earth_param_set = LP.LandParameters(toml_dict)
    thermo_params   = LP.thermodynamic_parameters(earth_param_set)

    T_air  = make_tvi("Tair"; offset = -273.15)  # K → intermediate for q
    VPD    = make_tvi("Qair")                     # we need to handle humidity specially

    # Read arrays directly for q computation
    ds = NCDataset(met_nc_path, "r")
    T_arr    = Float64.(coalesce.(ds["Tair"][1, 1, :],   NaN))    # K
    q_arr    = Float64.(coalesce.(ds["Qair"][1, 1, :],   NaN))    # kg/kg specific humidity
    P_arr    = Float64.(coalesce.(ds["Psurf"][1, 1, :],  NaN))    # Pa  (note: Psurf not PSurf)
    SW_arr   = Float64.(coalesce.(ds["SWdown"][1, 1, :], NaN))    # W/m²
    LW_arr   = Float64.(coalesce.(ds["LWdown"][1, 1, :], NaN))    # W/m²
    prec_arr = Float64.(coalesce.(ds["Precip"][1, 1, :], NaN))    # kg/m²/s total precip
    wind_arr = Float64.(coalesce.(ds["Wind"][1, 1, :],   NaN))    # m/s
    co2_arr  = Float64.(coalesce.(ds["CO2air"][1, 1, :], 412e-6)) # mol/mol
    close(ds)

    valid = .!isnan.(T_arr) .& .!isnan.(q_arr) .& .!isnan.(P_arr) .&
            .!isnan.(SW_arr) .& .!isnan.(LW_arr) .& .!isnan.(wind_arr)
    tv   = t_seconds[valid]
    flat = LinearInterpolation(Flat())

    # Cyclic 1-year spinup: prepend a shifted copy of the first year of met data
    # so the spinup gap (sim_start → first obs) gets realistic diurnal/seasonal forcing
    t0       = tv[1]                                  # time of first valid obs (s)
    yr1_idx  = tv .<= t0 + 365.25 * 86400            # indices covering first ~year
    tv_pre   = tv[yr1_idx] .- t0                      # shift to start near t=0

    tv_full   = vcat(tv_pre, tv)
    T_full    = vcat(T_arr[valid][yr1_idx],    T_arr[valid])
    q_full    = vcat(q_arr[valid][yr1_idx],    q_arr[valid])
    P_full    = vcat(P_arr[valid][yr1_idx],    P_arr[valid])
    SW_full   = vcat(SW_arr[valid][yr1_idx],   SW_arr[valid])
    LW_full   = vcat(LW_arr[valid][yr1_idx],   LW_arr[valid])
    prec_full = vcat(prec_arr[valid][yr1_idx], prec_arr[valid])
    wind_full = vcat(wind_arr[valid][yr1_idx], wind_arr[valid])
    co2_full  = vcat(co2_arr[valid][yr1_idx],  co2_arr[valid])

    # Split total precip into rain/snow by 0°C threshold (negative = downward flux)
    is_snow   = T_full .< 273.15
    atmos_T   = TimeVaryingInput(tv_full, T_full;    method = flat)
    atmos_q   = TimeVaryingInput(tv_full, q_full;    method = flat)
    atmos_P   = TimeVaryingInput(tv_full, P_full;    method = flat)
    atmos_SW  = TimeVaryingInput(tv_full, SW_full;   method = flat)
    atmos_LW  = TimeVaryingInput(tv_full, LW_full;   method = flat)
    atmos_Pr  = TimeVaryingInput(tv_full, -abs.(prec_full) .* .!is_snow; method = flat)
    atmos_Ps  = TimeVaryingInput(tv_full, -abs.(prec_full) .* is_snow;   method = flat)
    atmos_u   = TimeVaryingInput(tv_full, wind_full; method = flat)
    atmos_co2 = TimeVaryingInput(tv_full, co2_full;  method = flat)

    atmos = ClimaLand.PrescribedAtmosphere(
        atmos_Pr, atmos_Ps, atmos_T, atmos_u, atmos_q, atmos_P,
        sim_start, FT(_ATMOS_H), toml_dict; c_co2 = atmos_co2,
    )

    cos_zenith = (t, s) -> ClimaLand.default_cos_zenith_angle(
        t, s;
        insol_params = earth_param_set.insol_params,
        longitude    = _SITE_LON,
        latitude     = _SITE_LAT,
    )
    radiation = ClimaLand.PrescribedRadiativeFluxes(
        FT, atmos_SW, atmos_LW, sim_start;
        cosθs = cos_zenith, toml_dict,
    )
    return (; atmos, radiation)
end

"""Load LAI from the DK-Sor Met NetCDF file (LAI_alternative column)."""
function load_dk_sor_lai(met_nc_path::String, sim_start::DateTime)
    ds = NCDataset(met_nc_path, "r")
    lai_data = Float64.(coalesce.(ds["LAI_alternative"][1, 1, :], NaN))
    met_times = ds["time"][:]
    close(ds)

    t_seconds = [Float64(Second(t - Hour(_UTC_OFFSET) - sim_start).value)
                 for t in met_times]
    valid = .!isnan.(lai_data)
    flat  = LinearInterpolation(Flat())
    tv    = t_seconds[valid]
    t0    = tv[1]
    yr1_idx  = tv .<= t0 + 365.25 * 86400
    tv_full  = vcat(tv[yr1_idx] .- t0, tv)
    lai_full = vcat(lai_data[valid][yr1_idx], lai_data[valid])
    TimeVaryingInput(tv_full, lai_full; method = flat)
end

"""Build DK-Sor LandModel with site-specific soil/canopy parameters."""
function build_dk_sor_model(
    FT,
    sim_start::DateTime,
    sim_stop::DateTime,
    toml_dict,
    met_nc_path::String,
)
    DT = Float64(900)   # timestep (s)

    forcing = load_dk_sor_forcing(met_nc_path, sim_start, toml_dict, FT)
    LAI     = load_dk_sor_lai(met_nc_path, sim_start)

    # Column domain — 20 layers, 9 m deep, stretched
    land_domain = Column(;
        zlim      = (FT(-9), FT(0)),
        nelements = 20,
        dz_tuple  = (FT(1.5), FT(0.025)),
        longlat   = (FT(_SITE_LON), FT(_SITE_LAT)),
    )
    canopy_domain = ClimaLand.Domains.obtain_surface_domain(land_domain)

    # ── Soil ─────────────────────────────────────────────────────────────────
    ν       = FT(_SOIL_ν)
    θ_r     = FT(_SOIL_θ_r)
    K_sat   = FT(_SOIL_K_sat)
    vg_n    = FT(_SOIL_VG_N)
    vg_α    = FT(_SOIL_VG_α)
    hydrology_cm = ClimaLand.Soil.vanGenuchten{FT}(; α = vg_α, n = vg_n)

    retention_parameters = (;
        ν, hydrology_cm, θ_r, K_sat,
    )
    composition_parameters = (;
        ν_ss_om     = FT(_SOIL_ν_ss_om),
        ν_ss_quartz = FT(_SOIL_ν_ss_q),
        ν_ss_gravel = FT(_SOIL_ν_ss_g),
    )

    prognostic_land_components = (:canopy, :snow, :soil, :soilco2)

    soil = ClimaLand.Soil.EnergyHydrology{FT}(
        land_domain,
        forcing,
        toml_dict;
        prognostic_land_components,
        additional_sources   = (ClimaLand.RootExtraction{FT}(),),
        retention_parameters,
        composition_parameters,
        S_s = FT(1e-3),
    )

    # ── Canopy ────────────────────────────────────────────────────────────────
    radiation_parameters = (;
        Ω            = FT(_Ω),
        G_Function   = CLMGFunction(FT(_χl)),
        α_PAR_leaf   = FT(_α_PAR_leaf),
        τ_PAR_leaf   = FT(_τ_PAR_leaf),
        α_NIR_leaf   = FT(_α_NIR_leaf),
        τ_NIR_leaf   = FT(_τ_NIR_leaf),
    )
    radiative_transfer = Canopy.TwoStreamModel{FT}(canopy_domain, toml_dict;
                                                    radiation_parameters)

    biomass = Canopy.PrescribedBiomassModel{FT}(
        land_domain,
        LAI,
        toml_dict;
        rooting_depth = FT(0.3),
        height        = toml_dict["canopy_height"],
        SAI           = toml_dict["SAI"],
        RAI           = toml_dict["RAI"],
    )
    forcing_nt = (;
        atmos     = forcing.atmos,
        radiation = forcing.radiation,
        ground    = ClimaLand.PrognosticGroundConditions{FT}(),
    )
    canopy = Canopy.CanopyModel{FT}(
        canopy_domain,
        forcing_nt,
        LAI,
        toml_dict;
        prognostic_land_components,
        photosynthesis       = Canopy.PModel{FT}(canopy_domain, toml_dict),
        conductance          = Canopy.PModelConductance{FT}(toml_dict),
        soil_moisture_stress = Canopy.PiecewiseMoistureStressModel{FT}(
            land_domain, toml_dict; soil_params = (; ν, θ_r),
        ),
        biomass,
        radiative_transfer,
        hydraulics           = Canopy.PlantHydraulicsModel{FT}(canopy_domain, toml_dict),
    )

    land = LandModel{FT}(
        forcing,
        LAI,
        toml_dict,
        land_domain,
        DT;
        prognostic_land_components,
        canopy,
        soil,
    )
    return land, forcing, ν, θ_r
end

"""Custom initial-condition setter for the DK-Sor column."""
function make_dk_sor_ic(atmos, ν, θ_r)
    function custom_set_ic!(Y, p, t, model)
        FT = eltype(Y.soil.ϑ_l)
        evaluate!(p.drivers.T, atmos.T, t)

        @. Y.soil.ϑ_l  = FT(θ_r) + (FT(ν) - FT(θ_r)) * FT(0.95)
        Y.soil.θ_i    .= FT(0)
        earth_param_set = ClimaLand.get_earth_param_set(model.soil)
        ρc_s = ClimaLand.Soil.volumetric_heat_capacity.(
            Y.soil.ϑ_l, Y.soil.θ_i,
            model.soil.parameters.ρc_ds,
            earth_param_set,
        )
        Y.soil.ρe_int .= ClimaLand.Soil.volumetric_internal_energy.(
            Y.soil.θ_i, ρc_s, p.drivers.T, earth_param_set,
        )

        Y.snow.S   .= FT(0)
        Y.snow.S_l .= FT(0)
        Y.snow.U   .= FT(0)

        if model.canopy.energy isa ClimaLand.Canopy.BigLeafEnergyModel
            Y.canopy.energy.T .= p.drivers.T
        end

        Y.canopy.hydraulics.ϑ_l .= model.canopy.hydraulics.parameters.ν

        if !isnothing(model.soilco2)
            # CO2: ~412 ppm at 288 K, 1 atm → kg C/m³ air-equivalent
            Y.soilco2.CO2 .= FT(6e-5)
            # O2: ~21% by volume → kg O2/m³ soil
            Y.soilco2.O2  .= FT(0.08)
            # SOC exponential profile (kg C/m³): 15 at surface, 0.5 at -9 m
            SOC_top = FT(15.0)
            SOC_bot = FT(0.5)
            τ_soc   = FT(1.0 / log(SOC_top / SOC_bot))
            z = ClimaCore.Fields.coordinate_field(axes(Y.soilco2.SOC)).z
            @. Y.soilco2.SOC = SOC_bot + (SOC_top - SOC_bot) * exp(z / τ_soc)
        end
    end
    return custom_set_ic!
end

"""
    extract_daily_diag(simulation, diag_name)

Extract a scalar daily diagnostic from the DictWriter.
Returns (dates::Vector{Date}, values::Vector{Float64}).
"""
function extract_daily_diag(simulation, diag_name)
    writer = nothing
    for d in simulation.diagnostics
        if haskey(d.output_writer.dict, diag_name)
            writer = d.output_writer
            break
        end
    end
    if isnothing(writer)
        available = reduce(vcat, [collect(keys(d.output_writer.dict))
                                  for d in simulation.diagnostics])
        error("Diagnostic '$diag_name' not found. Available: $(unique(available))")
    end
    (times, data) = ClimaLand.Diagnostics.diagnostic_as_vectors(writer, diag_name)
    model_dates = Date.(times isa Vector{DateTime} ? times : date.(times))
    return model_dates, Float64.(data)
end

# ── Forward Model ────────────────────────────────────────────────────────────

"""
    ClimaCalibrate.forward_model(iteration, member)

Run the DK-Sor ClimaLand model on the minibatched yearly windows selected for
this iteration.  Saves NEE, LE, H daily diagnostics to JLD2.

Spinup: 2-year cycle of 1997–1998 forcing before the first selected window.
"""
function ClimaCalibrate.forward_model(iteration, member)
    FT = Float64

    climaland_dir = abspath(joinpath(@__DIR__, "..", ".."))
    met_nc_path   = joinpath(climaland_dir, "DK_Sor",
                             "DK-Sor_1997-2014_FLUXNET2015_Met.nc")

    # ── Determine simulation window from minibatch ────────────────────────────
    ekp     = JLD2.load_object(ClimaCalibrate.ekp_path(OUTPUT_DIR, iteration))
    mb_idxs = EKP.get_current_minibatch(ekp)

    # Yearly windows: index i → year (1996 + i), i.e. window 1 = 1997, ..., 16 = 2012
    obs_data = JLD2.load(OBS_FILEPATH)
    window_years = obs_data["window_years"]   # Vector{Int}: [1997, 1998, ..., 2012]

    selected_years  = window_years[mb_idxs]
    year_start      = minimum(selected_years)
    year_end        = maximum(selected_years)

    # 1-year spinup before the first selected year, but never before 1996
    spinup_start = max(year_start - 1, 1996)
    sim_start    = DateTime(spinup_start, 1, 1)
    sim_stop     = DateTime(year_end + 1, 1, 1)   # exclusive end = Jan 1 of year+1

    @info "Member $member (iter $iteration): years $(selected_years), " *
          "sim $sim_start → $sim_stop"

    # ── Load calibration parameters ───────────────────────────────────────────
    calibrate_params_path =
        ClimaCalibrate.parameter_path(OUTPUT_DIR, iteration, member)
    toml_dict = LP.create_toml_dict(FT; override_files = [calibrate_params_path])

    # ── Build model ───────────────────────────────────────────────────────────
    land, forcing, ν, θ_r = build_dk_sor_model(
        FT, sim_start, sim_stop, toml_dict, met_nc_path,
    )
    set_ic! = make_dk_sor_ic(forcing.atmos, ν, θ_r)

    # ── Diagnostics (daily, DictWriter) ──────────────────────────────────────
    output_writer = ClimaDiagnostics.Writers.DictWriter()
    diags = ClimaLand.default_diagnostics(
        land, sim_start, "";          # outdir="" unused (DictWriter handles output)
        output_writer,
        output_vars      = :short,
        reduction_period = :daily,
    )

    # ── Simulation ────────────────────────────────────────────────────────────
    simulation = LandSimulation(
        sim_start, sim_stop, Float64(DT), land;
        set_ic!,
        updateat       = Second(DT),
        user_callbacks = (),
        diagnostics    = diags,
    )
    solve!(simulation)

    # ── Extract and save diagnostics ──────────────────────────────────────────
    member_path = ClimaCalibrate.path_to_ensemble_member(
        OUTPUT_DIR, iteration, member,
    )
    save_eki_diagnostics(simulation, member_path, selected_years)
    return nothing
end

"""Save daily NEE, LE, H for the selected calibration years."""
function save_eki_diagnostics(simulation, member_path, selected_years)
    nee_dates, nee_vals = extract_daily_diag(simulation, "nee_1d_average")
    _,         qle_vals = extract_daily_diag(simulation, "lhf_1d_average")
    _,         qh_vals  = extract_daily_diag(simulation, "shf_1d_average")

    # Keep only days within selected calibration years (drop spinup)
    mask = any(year.(nee_dates) .== yr for yr in selected_years)
    JLD2.jldsave(
        joinpath(member_path, "daily_diagnostics.jld2");
        dates = nee_dates[mask],
        nee   = nee_vals[mask],
        qle   = qle_vals[mask],
        qh    = qh_vals[mask],
    )
end

# ── Observation Map ──────────────────────────────────────────────────────────

"""
    ClimaCalibrate.observation_map(iteration)

Return stacked G ensemble matrix for all minibatch windows.

Layout per member: [NEE_win1..., LE_win1..., H_win1..., NEE_win2..., ...].
Windows are in the order given by the current minibatch index vector.
Missing days are filled with 0.0 (large-variance noise handles these in obs cov).
"""
function ClimaCalibrate.observation_map(iteration)
    ekp          = JLD2.load_object(ClimaCalibrate.ekp_path(OUTPUT_DIR, iteration))
    ensemble_size = EKP.get_N_ens(ekp)
    mb_idxs      = EKP.get_current_minibatch(ekp)

    obs_data     = JLD2.load(OBS_FILEPATH)
    window_years = obs_data["window_years"]
    window_dates = obs_data["window_dates"]   # Vector{Vector{Date}}

    selected_dates = window_dates[mb_idxs]
    n_days_each    = length.(selected_dates)
    g_len          = 3 * sum(n_days_each)

    G_ens = zeros(g_len, ensemble_size)

    for m in 1:ensemble_size
        member_path = ClimaCalibrate.path_to_ensemble_member(
            OUTPUT_DIR, iteration, m,
        )
        diag_path = joinpath(member_path, "daily_diagnostics.jld2")
        try
            d = JLD2.load(diag_path)
            nee_d = Dict(zip(d["dates"], d["nee"]))
            qle_d = Dict(zip(d["dates"], d["qle"]))
            qh_d  = Dict(zip(d["dates"], d["qh"]))

            offset = 0
            for dates in selected_dates
                n = length(dates)
                nee_out = [get(nee_d, dt, 0.0) for dt in dates]
                qle_out = [get(qle_d, dt, 0.0) for dt in dates]
                qh_out  = [get(qh_d,  dt, 0.0) for dt in dates]

                # Convert NEE: mol CO2/m²/s → gC/m²/d
                nee_out .*= 12.0 * 86400.0

                G_ens[(offset + 1):(offset + n),       m] = nee_out
                G_ens[(offset + n + 1):(offset + 2n),  m] = qle_out
                G_ens[(offset + 2n + 1):(offset + 3n), m] = qh_out
                offset += 3n
            end
        catch e
            @error "Error processing member $m" exception = e
            G_ens[:, m] .= NaN
        end
    end
    return G_ens
end
