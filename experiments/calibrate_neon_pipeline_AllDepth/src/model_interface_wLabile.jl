"""
ClimaCalibrate model interface for NEON site soil CO₂ calibration with a
depth-dependent labile carbon / microbial-activity profile.

Extends model_interface.jl by adding a 5th calibrated parameter:
labile_depth_scale (= k, units 1/m) — the exponential decay rate of an extra
depth factor exp(k·z) (z<0, so the factor decays with depth) that is folded
into the prognostic SOC field set in custom_set_ic!.

In the DAMM source term the substrate is `Sx = p_sx · Csom · D_liq · θ_l³`, so
multiplying the SOC field (Csom) by exp(k·z) is mathematically equivalent to a
depth-dependent soluble-carbon fraction p_sx(z) — without modifying ClimaLand's
scalar p_sx. We calibrate the DECAY RATE k (the *shape*), not an amplitude
(which would be degenerate with soilCO2_reference_rate). exp(k·0)=1 keeps the
surface level fixed; k=0 recovers the constant-p_sx behaviour.

Because labile_depth_scale is not a ClimaParams TOML key, it is read directly
via Julia's built-in TOML parser rather than LP.create_toml_dict.

This file is `@everywhere include`d on all workers.
"""

import ClimaCalibrate
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.FluxnetSimulations as FluxnetSimulations
import ClimaLand.Simulations: LandSimulation, solve!
using ClimaLand
using ClimaLand.Domains: Column
using ClimaLand.Soil
using ClimaLand.Soil.Biogeochemistry
using ClimaLand.Canopy
#using ClimaLand.Canopy.PlantHydraulics
using ClimaLand.Snow
using ClimaCore
using ClimaDiagnostics
using ClimaUtilities
import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput, evaluate!
import ClimaUtilities.SpaceVaryingInputs: SpaceVaryingInput
import ClimaUtilities.Utils: searchsortednearest, linear_interpolation
import Interpolations
import ClimaUtilities.TimeManager: date
using Insolation

import EnsembleKalmanProcesses as EKP
import JLD2
import TOML as TOML_pkg
using Dates
using Statistics
using DataFrames
using CSV
using DelimitedFiles

# Read one NEON per-depth retention parameter into a subsurface field.
# Linear-interpolate within the measured range; hold the deepest measured
# value (flat) below it — matching Rosetta's Interpolations.Flat().
function read_neon_profile(csv_path, colname, space, FT)
    df = CSV.read(csv_path, DataFrame)
    valid = .!ismissing.(df[!, colname])
    z_raw = Float64.(df.depth[valid])
    v_raw = Float64.(df[valid, colname])
    si = sortperm(z_raw)                 # ascending z: most-negative → 0
    z = z_raw[si]
    v = v_raw[si]
    z_bot = z[1]                         # deepest (most negative) measurement
    v_bot = v[1]
    zvals = ClimaCore.Fields.coordinate_field(space).z
    return map(zvals) do zc
        val = zc > z_bot ? linear_interpolation(z, v, zc) : v_bot
        FT(val)
    end
end

# ── Model Interface ──────────────────────────────────────────────────────────
# ClimaCalibrate's current API dispatches forward_model / observation_map on a
# user-defined AbstractModelInterface subtype, and calibrate() takes an instance
# as its 3rd argument. This struct carries no state — config still flows through
# module globals / ENV — it exists purely as a dispatch tag (cf. the official
# experiments/calibration/ LandModelInterface, which instead stores a config).
struct NeonLabileModelInterface <: ClimaCalibrate.AbstractModelInterface end

# ── Forward Model ────────────────────────────────────────────────────────────

function ClimaCalibrate.forward_model(::NeonLabileModelInterface, iteration, member)
    FT = Float64
    site_ID_val = FluxnetSimulations.replace_hyphen(SITE_ID)
    climaland_dir = pkgdir(ClimaLand)


    # Get simulation dates from site metadata
    #(; time_offset, lat, long) =
    #    FluxnetSimulations.get_location(FT, Val(site_ID_val))
    time_offset = 0
    metadata = _get_neon_site_metadata(SITE_ID)
    lat = FT(metadata.lat)
    long = FT(metadata.long)
    atmos_h = FT(metadata.atmos_h)

    (start_date, stop_date) =
        FluxnetSimulations.get_data_dates(SITE_ID, time_offset)
    #overwrite with environment variables if provided
    start_date = DateTime(get(ENV, "NEON_START_DATE", string(Date(start_date))))
    stop_date = DateTime(get(ENV, "NEON_STOP_DATE", string(Date(stop_date))))
    @info "Member $member: simulating $start_date to $stop_date"

    # Load calibrated parameters from TOML written by ClimaCalibrate
    calibrate_params_path =
        ClimaCalibrate.parameter_path(OUTPUT_DIR, iteration, member)
    toml_dict =
        LP.create_toml_dict(FT; override_files = [calibrate_params_path])

    # Read labile_depth_scale directly — NOT a ClimaParams key, so LP ignores it.
    # k (1/m): exp(k·z) decays with depth (z<0). k=0 → constant-p_sx behaviour.
    calib_params_raw = TOML_pkg.parsefile(calibrate_params_path)
    labile_depth_scale = FT(calib_params_raw["labile_depth_scale"]["value"])
    @info "Member $member: labile_depth_scale = $labile_depth_scale"

    # Domain
    dz_bottom = FT(2) #FT(1.5),
    dz_top = FT(0.038)
    dz_tuple = (dz_bottom, dz_top)
    nelements = 24#24
    zmin = FT(-6.2)
    zmax = FT(0)

    #(; dz_tuple, nelements, zmin, zmax) =
    #    FluxnetSimulations.get_domain_info(FT, Val(site_ID_val))
    #(; atmos_h) =
    #    FluxnetSimulations.get_fluxtower_height(FT, Val(site_ID_val))

    land_domain = Column(;
        zlim = (zmin, zmax),
        nelements = nelements,
        dz_tuple = dz_tuple,
        longlat = (long, lat),
    )
    canopy_domain = ClimaLand.Domains.obtain_surface_domain(land_domain)
    surface_space = land_domain.space.surface

    # Determine the target model layer for soil CO₂ extraction — one per calibrated
    # depth. The obs JLD2 supplies the per-depth CODES and MEASUREMENT DEPTHS
    # (z_obs_m); the model LAYER is derived here, at runtime, by argmin against the
    # live grid — exactly as the old single-depth path did, just looped. This keeps
    # layer selection tied to the actual domain (independent of the lookup CSV).
    # OBS_FILEPATH is a global set on every process.
    z_field = ClimaCore.Fields.coordinate_field(land_domain.space.subsurface).z
    z_vals = parent(z_field)[:, 1]

    _obs_meta = JLD2.load(OBS_FILEPATH)
    depth_codes = String.(_obs_meta["depth_codes"])
    z_obs_per_code = FT.(_obs_meta["z_obs_m"])                 # ordered like depth_codes
    target_layers = [argmin(abs.(z_vals .- z)) for z in z_obs_per_code]
    @info "Member $member: depths $(depth_codes) (z_obs=$(z_obs_per_code) m) → layers $(target_layers)"

    # Base TOML for non-calibrated parameters (canopy, snow, etc.)
    toml_dict_base = LP.create_toml_dict(FT)

    # ERA5 forcing (global-style)
    (; atmos, radiation) = FluxnetSimulations.prescribed_forcing_fluxnet(
        SITE_ID,
        lat,
        long,
        time_offset,
        atmos_h,
        start_date,
        toml_dict_base,
        FT,
    )

    # Custom canopy aerodynamic coefficients
    toml_dict_base.data["canopy_d_coeff"]["value"] = FT(0.67)
    toml_dict_base.data["canopy_z_0b_coeff"]["value"] = FT(0.013)
    toml_dict_base.data["canopy_z_0m_coeff"]["value"] = FT(0.13)

    # MODIS LAI (global-style)
    #LAI = ClimaLand.Canopy.prescribed_lai_modis(surface_space, start_date, stop_date)
    LAI = ClimaLand.Canopy.prescribed_climatological_lai_modis(surface_space)

    # Build model components
    prognostic_land_components = (:canopy, :snow, :soil, :soilco2)
    forcing = (; atmos, radiation)
    ground = ClimaLand.PrognosticGroundConditions{FT}()
    canopy_forcing = (; atmos, radiation, ground)

    retention_parameters = Soil.soil_vangenuchten_parameters(#rosetta_soil_vangenuchten_parameters(
        land_domain.space.subsurface,
        FT,
    )

    # ── Override retention params with NEON per-site profiles ─────────────────
    # Mutate the fields of the SAME object in place so the shared reference used
    # by the soil model, the canopy soil-moisture-stress model, and LandModel's
    # check_land_equality (θ_high == ν, exact ==) all stay consistent. Do NOT
    # build a second object here.
    neon_dir = "/kiwi-data/Data/groupMembers/evametz/Neon/Neon_data"
    sp = land_domain.space.subsurface

    α_field = read_neon_profile(
        joinpath(neon_dir, "NEON_all_sites_alpha_1_m_2cm_mean.csv"),
        "$(SITE_ID)_alpha_[1/m]", sp, FT)
    n_field = read_neon_profile(
        joinpath(neon_dir, "NEON_all_sites_n_-_2cm_mean.csv"),
        "$(SITE_ID)_n_[-]", sp, FT)

    # hydrology_cm is a field of vanGenuchten structs → rebuild from α, n
    retention_parameters.hydrology_cm .=
        ((α, n) -> ClimaLand.Soil.vanGenuchten{FT}(; α, n)).(α_field, n_field)

    retention_parameters.K_sat .= read_neon_profile(
        joinpath(neon_dir, "NEON_all_sites_Ksat_m_s_2cm_mean.csv"),
        "$(SITE_ID)_Ksat_[m/s]", sp, FT)
    retention_parameters.ν .= read_neon_profile(
        joinpath(neon_dir, "NEON_all_sites_nu_m3_m3_2cm_mean.csv"),
        "$(SITE_ID)_nu_[m3/m3]", sp, FT)
    retention_parameters.θ_r .= read_neon_profile(
        joinpath(neon_dir, "NEON_all_sites_theta_r_m3_m3_2cm_mean.csv"),
        "$(SITE_ID)_theta_r_[m3/m3]", sp, FT)

    # Canopy with PModel (uses base TOML, not calibrated TOML)
    photosynthesis = PModel{FT}(land_domain, toml_dict_base)
    conductance = PModelConductance{FT}(toml_dict_base)
    soil_moisture_stress =
        ClimaLand.Canopy.PiecewiseMoistureStressModel{FT}(
            land_domain, 
            toml_dict_base;
            soil_params = retention_parameters
    )
    biomass = ClimaLand.Canopy.PrescribedBiomassModel{FT}(land_domain, LAI, toml_dict_base)

    canopy = ClimaLand.Canopy.CanopyModel{FT}(
        canopy_domain,
        canopy_forcing,
        LAI,
        toml_dict_base;
        prognostic_land_components,
        photosynthesis,
        conductance,
        soil_moisture_stress,
        biomass,
    )

    # Snow model with zenith-angle-dependent albedo
    α_snow = Snow.ZenithAngleAlbedoModel(toml_dict_base)
    snow = Snow.SnowModel(
        FT,
        canopy_domain,
        forcing,
        toml_dict_base,
        DT;
        prognostic_land_components,
        α_snow,
    )

        # Build the soil model explicitly so we can use the Rosetta (Montzka et al.
    # 2017) van Genuchten retention parameters instead of the default Gupta 2020.
    # rosetta_soil_vangenuchten_parameters returns the same NamedTuple shape
    # (; ν, hydrology_cm, K_sat, θ_r) as the default soil_vangenuchten_parameters.
    soil = Soil.EnergyHydrology{FT}(
        land_domain,
        forcing,
        toml_dict;
        prognostic_land_components,
        additional_sources = (ClimaLand.RootExtraction{FT}(),),
        retention_parameters = retention_parameters
    )


    # Full LandModel — calibrated TOML used for soil/soilCO2 parameters
    land = LandModel{FT}(
        forcing,
        LAI,
        toml_dict,
        land_domain,
        DT;
        prognostic_land_components,
        snow,
        canopy,
        soil,
    )

    # θ_r comes from the NEON CSV (read into retention_parameters above); no
    # per-run uniform override. Mirrors ForwardRun.jl so calibration and the
    # final forward run see identical θ_r.

    # Custom IC with SOC profile from SoilGrids OCD artifact
    base_set_ic! = FluxnetSimulations.make_set_fluxnet_initial_conditions(
        SITE_ID,
        start_date,
        time_offset,
        land,
    )

    ocd_path = ClimaLand.Artifacts.soil_grids_ocd_artifact_path()
    SOC_from_artifact = SpaceVaryingInput(
        ocd_path,
        "ocd",
        land_domain.space.subsurface;
        regridder_type = :InterpolationsRegridder,
        regridder_kwargs = (;
            extrapolation_bc = (
                Interpolations.Periodic(),
                Interpolations.Flat(),
                Interpolations.Flat(),
            ),
            interpolation_method = Interpolations.Linear(),
        ),
    )
    #=
    function custom_set_ic!(Y, p, t, model)
        base_set_ic!(Y, p, t, model)
        Y.soilco2.CO2 .= FT(6e-5)
        Y.soilco2.O2 .= FT(0.08)
        #read csv file with depth and SOC values, then interpolate to model layers
        model_value = ClimaCore.Fields.zeros(land_domain.space.subsurface)
        data = CSV.read("/kiwi-data/Data/groupMembers/evametz/Neon/Neon_data/NEON_all_sites_estimatedOC_2cm_mean_extrapolated.csv", DataFrame)
        valid = .!ismissing.(data[!, "$(SITE_ID)_estimatedOC_kg_m3"])
        #z_bottom = minimum(parent(ClimaCore.Fields.coordinate_field(land_domain.space.subsurface).z))
        raw_z::Vector{Float64} = Float64.(data.depth[valid])
        raw_vals::Vector{Float64} = Float64.(data[valid, "$(SITE_ID)_estimatedOC_kg_m3"])
        sort_idx = sortperm(raw_z)
        # z_bottom is the most negative value → prepend so itp_z stays ascending
        #itp_z::Vector{Float64} = vcat(z_bottom, raw_z[sort_idx])
        #itp_values::Vector{Float64} = vcat(FT(0.5), raw_vals[sort_idx])
        zvalues = ClimaCore.Fields.coordinate_field(axes(Y.soilco2.SOC)).z
        #model_value .= map(zvalues) do z
        #    linear_interpolation(itp_z, itp_values, z)
        model_value .= map(zvalues) do z
            linear_interpolation(raw_z[sort_idx], raw_vals[sort_idx], z)
        end
        Y.soilco2.SOC .= model_value
    end=#

    #= old SOC-only custom_set_ic! (uniform-column SWC from base_set_ic!)
    function custom_set_ic!(Y, p, t, model)
        base_set_ic!(Y, p, t, model)
        #Y.soilco2.CO2 .= FT(6e-5)
        #Y.soilco2.O2 .= FT(0.08)
        model_value = ClimaCore.Fields.zeros(land_domain.space.subsurface)
        data = CSV.read("/kiwi-data/Data/groupMembers/evametz/Neon/Neon_data/NEON_all_sites_estimatedOC_2cm_mean.csv", DataFrame)
        valid = .!ismissing.(data[!, "$(SITE_ID)_estimatedOC_kg_m3"])
        raw_z::Vector{Float64} = Float64.(data.depth[valid])
        sort_idx = sortperm(raw_z)
        raw_vals::Vector{Float64} = Float64.(data[valid, "$(SITE_ID)_estimatedOC_kg_m3"])

        z_extrap_top = (raw_z[sort_idx])[1]
        SOC_extrap_top = (raw_vals[sort_idx])[1]
        SOC_extrap_bot = FT(0.5)
        z_extrap_bot = minimum(parent(ClimaCore.Fields.coordinate_field(land_domain.space.subsurface).z))

        zvalues = ClimaCore.Fields.coordinate_field(axes(Y.soilco2.SOC)).z

        alpha_soc = FT(log(SOC_extrap_top / SOC_extrap_bot) / (z_extrap_bot - z_extrap_top))

        model_value .= map(zvalues) do z
            if z > z_extrap_top
                linear_interpolation(raw_z[sort_idx], raw_vals[sort_idx], z)
            else
                SOC_extrap_top * exp(- alpha_soc * (z - z_extrap_top))
            end
        end
        Y.soilco2.SOC .= model_value
    end
    =#

    # New custom_set_ic! — SOC profile from NEON CSV + soil-moisture profile
    # derived from NEON VSWCMean sensors (time- and plot-mean per depth code).
    # The SOC profile is multiplied by the depth-dependent labile factor
    # exp(labile_depth_scale·z) — the 5th calibrated parameter (see file header).
    function custom_set_ic!(Y, p, t, model)
        base_set_ic!(Y, p, t, model)

        # ── 1. SOC profile (NEON CSV with exponential extrapolation below) ──
        soc_field = ClimaCore.Fields.zeros(land_domain.space.subsurface)
        soc_data = CSV.read(
            "/kiwi-data/Data/groupMembers/evametz/Neon/Neon_data/NEON_all_sites_estimatedOC_2cm_mean.csv",
            DataFrame,
        )
        valid_soc = .!ismissing.(soc_data[!, "$(SITE_ID)_estimatedOC_kg_m3"])
        raw_z::Vector{Float64} = Float64.(soc_data.depth[valid_soc])
        sort_idx_soc = sortperm(raw_z)
        raw_vals::Vector{Float64} =
            Float64.(soc_data[valid_soc, "$(SITE_ID)_estimatedOC_kg_m3"])

        z_extrap_top = (raw_z[sort_idx_soc])[1]
        SOC_extrap_top = (raw_vals[sort_idx_soc])[1]
        SOC_extrap_bot = FT(0.05)
        z_extrap_bot = minimum(parent(
            ClimaCore.Fields.coordinate_field(land_domain.space.subsurface).z,
        ))
        zvalues = ClimaCore.Fields.coordinate_field(axes(Y.soilco2.SOC)).z
        alpha_soc =
            FT(log(SOC_extrap_top / SOC_extrap_bot) / (z_extrap_bot - z_extrap_top))

        # SOC profile × depth-dependent labile factor exp(k·z), in one pass.
        # exp(k·0)=1 keeps the surface fixed, so k controls only the shape
        # (not degenerate with V_ref_sx); k=0 recovers the constant-p_sx case.
        soc_field .= map(zvalues) do z
            soc = if z > z_extrap_top
                linear_interpolation(raw_z[sort_idx_soc], raw_vals[sort_idx_soc], z)
            else
                SOC_extrap_top * exp(-alpha_soc * (z - z_extrap_top))
            end
            soc * exp(labile_depth_scale * z)
        end
        Y.soilco2.SOC .= soc_field

        # ── 2. Soil moisture profile from NEON VSWCMean columns ─────────────
        # Standard NEON soil sensor depths (m, negative = below surface).
        # 501 ≈ 6 cm, 502 ≈ 16 cm, 503 ≈ 26 cm, 504 ≈ 46 cm, 505 ≈ 66 cm,
        # 506 ≈ 86 cm, 507 ≈ 106 cm, 508 ≈ 166 cm
        neon_depths = FT[-0.06, -0.16, -0.26, -0.46, -0.66, -0.86, -1.06, -1.66]
        # NOTE: renamed from `depth_codes` to avoid shadowing forward_model's
        # `depth_codes` (the calibrated CO₂ codes) captured by this closure.
        swc_depth_codes = ["501", "502", "503", "504", "505", "506", "507", "508"]
        n_plots = 5

        csv_path = ClimaLand.Artifacts.experiment_fluxnet_data_path(SITE_ID)
        swc_data = CSV.read(csv_path, DataFrame)
        swc_colnames = names(swc_data)

        swc_per_depth = FT[]
        for code in swc_depth_codes
            vals = Float64[]
            for plot_id in 1:n_plots
                colname = "VSWCMean_$(lpad(plot_id, 3, '0'))_$code"
                colname in swc_colnames || continue
                for v in swc_data[!, colname]
                    (ismissing(v) || isnan(Float64(v))) && continue
                    push!(vals, Float64(v))
                end
            end
            push!(swc_per_depth, isempty(vals) ? FT(NaN) : FT(mean(vals)))
        end

        valid_swc = .!isnan.(swc_per_depth)
        swc_z_valid = neon_depths[valid_swc]
        swc_vals_valid = swc_per_depth[valid_swc]
        sort_idx_swc = sortperm(swc_z_valid)
        swc_z_sorted = swc_z_valid[sort_idx_swc]
        swc_vals_sorted = swc_vals_valid[sort_idx_swc]

        @info "NEON-derived SWC profile" depths = swc_z_sorted theta_l = swc_vals_sorted

        z_top_data = swc_z_sorted[end]
        z_bot_data = swc_z_sorted[1]
        swc_top = swc_vals_sorted[end]
        swc_bot = swc_vals_sorted[1]

        z_soil = ClimaCore.Fields.coordinate_field(axes(Y.soil.ϑ_l)).z
        Y.soil.ϑ_l .= map(z_soil) do z
            if z > z_top_data
                FT(swc_top)
            elseif z < z_bot_data
                FT(swc_bot)
            else
                FT(linear_interpolation(swc_z_sorted, swc_vals_sorted, z))
            end
        end

        ν_field = land.soil.parameters.ν
        θ_r_field = land.soil.parameters.θ_r
        @. Y.soil.ϑ_l =
            clamp(Y.soil.ϑ_l, θ_r_field + FT(1e-4), ν_field - FT(1e-4))
    end
    #=
    function custom_set_ic!(Y, p, t, model)
        base_set_ic!(Y, p, t, model)
        Y.soilco2.CO2 .= FT(6e-5)
        Y.soilco2.O2 .= FT(0.08)
        Y.soilco2.SOC .= SOC_from_artifact
    end=#
    #=
    function custom_set_ic!(Y, p, t, model)
        base_set_ic!(Y, p, t, model)
        Y.soilco2.CO2 .= FT(6e-5)
        Y.soilco2.O2 .= FT(0.08)
        SOC_top = FT(15.0)
        SOC_bot = FT(0.5)
        τ_soc = FT(1.0 / log(SOC_top / SOC_bot))
        z = ClimaCore.Fields.coordinate_field(axes(Y.soilco2.SOC)).z
        @. Y.soilco2.SOC = SOC_bot + (SOC_top - SOC_bot) * exp(z / τ_soc)
    end=#

    # Diagnostics — halfhourly sco2_ppm (and supporting soil variables)
    output_writer = ClimaDiagnostics.Writers.DictWriter()
    output_vars = ["swc", "tsoil", "si", "sco2", "soc", "so2", "sco2_ppm"]
    diags = ClimaLand.default_diagnostics(
        land,
        start_date;
        output_writer = output_writer,
        output_vars,
        reduction_period = :halfhourly,
    )

    simulation = LandSimulation(
        start_date,
        stop_date,
        DT,
        land;
        set_ic! = custom_set_ic!,
        updateat = Second(DT),
        diagnostics = diags,
    )
    solve!(simulation)

    # Extract daily mean sco2_ppm at EACH target layer and save to JLD2
    member_path =
        ClimaCalibrate.path_to_ensemble_member(OUTPUT_DIR, iteration, member)
    save_daily_sco2(simulation, member_path, depth_codes, target_layers)
    return nothing
end

"""
    save_daily_sco2(simulation, member_path, depth_codes, target_layers)

Extract halfhourly sco2_ppm at each `target_layers[i]`, compute daily means, and
save one series per depth code to JLD2. The saved `sco2_ppm_by_code` is a
`Dict(code => daily_mean_vector)` and `dates_by_code` a `Dict(code => dates)`,
both keyed by `depth_codes` — observation_map aligns these to the per-depth obs
dates. Order-matched: `depth_codes[i]` ↔ `target_layers[i]`.
"""
function save_daily_sco2(simulation, member_path, depth_codes, target_layers)
    sco2_by_code = Dict{String, Vector{Float64}}()
    dates_by_code = Dict{String, Vector{Date}}()

    for (code, layer) in zip(depth_codes, target_layers)
        (times, data) = ClimaLand.Diagnostics.diagnostic_as_vectors(
            simulation.diagnostics[1].output_writer,
            "sco2_ppm_30m_average";
            layer = layer,
        )
        model_dates_dt = times isa Vector{DateTime} ? times : date.(times)
        model_df = DataFrame(datetime = model_dates_dt, sco2_ppm = Float64.(data))
        model_df[!, :date] = Date.(model_df.datetime)

        model_daily = combine(
            groupby(model_df, :date),
            :sco2_ppm => mean => :daily_mean,
        )
        sort!(model_daily, :date)

        sco2_by_code[code] = model_daily.daily_mean
        dates_by_code[code] = model_daily.date
    end

    JLD2.jldsave(
        joinpath(member_path, "daily_diagnostics.jld2");
        depth_codes = depth_codes,
        sco2_ppm_by_code = sco2_by_code,
        dates_by_code = dates_by_code,
    )
end

# ── Observation Map ──────────────────────────────────────────────────────────

"""
    ClimaCalibrate.observation_map(::NeonLabileModelInterface, iteration)

Return the G ensemble matrix matching the STACKED (per-depth concatenated)
observation vector. For each ensemble member, the model daily-mean soil CO₂ (ppm)
is aligned to each depth's own observation dates and the per-depth blocks are
concatenated in the SAME order the observation vector was built
(`obs_data["depth_codes"]`), so G row `k` corresponds to y_obs row `k`.
"""
function ClimaCalibrate.observation_map(::NeonLabileModelInterface, iteration)
    ekp = JLD2.load_object(ClimaCalibrate.ekp_path(OUTPUT_DIR, iteration))
    ensemble_size = EKP.get_N_ens(ekp)

    # Stacking order + per-depth obs dates (written by Observation_flag.jl)
    obs_data = JLD2.load(OBS_FILEPATH)
    depth_codes = String.(obs_data["depth_codes"])
    obs_dates_per_code = obs_data["obs_dates_per_code"]
    n_obs = sum(length(obs_dates_per_code[c]) for c in depth_codes)

    @info "Observation map: stacked G length $n_obs over depths $(depth_codes)"

    G_ens = zeros(n_obs, ensemble_size)

    for m in 1:ensemble_size
        member_path = ClimaCalibrate.path_to_ensemble_member(
            OUTPUT_DIR,
            iteration,
            m,
        )
        diag_path = joinpath(member_path, "daily_diagnostics.jld2")

        try
            member_data = JLD2.load(diag_path)
            sco2_by_code = member_data["sco2_ppm_by_code"]
            dates_by_code = member_data["dates_by_code"]

            # Concatenate per-depth blocks in the stacking order; within each depth
            # align model days to that depth's obs dates (missing day → NaN).
            blocks = Vector{Float64}[]
            for code in depth_codes
                model_dict = Dict(zip(dates_by_code[code], sco2_by_code[code]))
                push!(blocks,
                    [get(model_dict, d, NaN) for d in obs_dates_per_code[code]])
            end
            G_ens[:, m] = reduce(vcat, blocks)
        catch e
            @error "Error processing member $m" exception = e
            G_ens[:, m] .= NaN
        end
    end

    return G_ens
end
