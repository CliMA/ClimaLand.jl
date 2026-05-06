"""
ClimaCalibrate model interface for NEON site soil CO₂ calibration.

Defines `forward_model(iteration, member)` and `observation_map(iteration)`.
This file is `@everywhere include`d on all workers.

Calibrates DAMM soil CO₂ parameters against NEON soil CO₂ concentration
observations. Uses ERA5 forcing, MODIS LAI, PModel photosynthesis,
and a prescribed SOC exponential profile.
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
using Dates
using Statistics
using DataFrames
using CSV
using DelimitedFiles

# ── Forward Model ────────────────────────────────────────────────────────────

function ClimaCalibrate.forward_model(iteration, member)
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

    # Domain
    dz_bottom = FT(2) #FT(1.5),
    dz_top = FT(0.038)
    dz_tuple = (dz_bottom, dz_top)
    nelements = 10#24
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

    # Determine target layer for soil CO₂ extraction (~2 cm depth)
    z_field = ClimaCore.Fields.coordinate_field(land_domain.space.subsurface).z
    z_vals = parent(z_field)[:, 1]
    
    target_depth = FT(-1 * parse(Float64, Caldepthnum))
    target_layer = argmin(abs.(z_vals .- target_depth))

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

    # Canopy with PModel (uses base TOML, not calibrated TOML)
    photosynthesis = PModel{FT}(land_domain, toml_dict_base)
    conductance = PModelConductance{FT}(toml_dict_base)
    soil_moisture_stress =
        ClimaLand.Canopy.PiecewiseMoistureStressModel{FT}(land_domain, toml_dict_base)
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
    )

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
        Y.soilco2.CO2 .= FT(0.000412)
        Y.soilco2.O2_f .= FT(0.21)
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
        #Y.soilco2.CO2 .= FT(0.000412)
        #Y.soilco2.O2_f .= FT(0.21)
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

        soc_field .= map(zvalues) do z
            if z > z_extrap_top
                linear_interpolation(raw_z[sort_idx_soc], raw_vals[sort_idx_soc], z)
            else
                FT(0) #SOC_extrap_top * exp(-alpha_soc * (z - z_extrap_top))
            end
        end
        Y.soilco2.SOC .= soc_field

        # ── 2. Soil moisture profile from NEON VSWCMean columns ─────────────
        # Standard NEON soil sensor depths (m, negative = below surface).
        # 501 ≈ 6 cm, 502 ≈ 16 cm, 503 ≈ 26 cm, 504 ≈ 46 cm, 505 ≈ 66 cm,
        # 506 ≈ 86 cm, 507 ≈ 106 cm, 508 ≈ 166 cm
        neon_depths = FT[-0.06, -0.16, -0.26, -0.46, -0.66, -0.86, -1.06, -1.66]
        depth_codes = ["501", "502", "503", "504", "505", "506", "507", "508"]
        n_plots = 5

        csv_path = ClimaLand.Artifacts.experiment_fluxnet_data_path(SITE_ID)
        swc_data = CSV.read(csv_path, DataFrame)
        swc_colnames = names(swc_data)

        swc_per_depth = FT[]
        for code in depth_codes
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
        Y.soilco2.CO2 .= FT(0.000412)
        Y.soilco2.O2_f .= FT(0.21)
        Y.soilco2.SOC .= SOC_from_artifact
    end=#
    #=
    function custom_set_ic!(Y, p, t, model)
        base_set_ic!(Y, p, t, model)
        Y.soilco2.CO2 .= FT(0.000412)
        Y.soilco2.O2_f .= FT(0.21)
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

    # Extract daily mean sco2_ppm at target layer and save to JLD2
    member_path =
        ClimaCalibrate.path_to_ensemble_member(OUTPUT_DIR, iteration, member)
    save_daily_sco2(simulation, member_path, target_layer)
    return nothing
end

"""
    save_daily_sco2(simulation, member_path, target_layer)

Extract halfhourly sco2_ppm at the target layer, compute daily means,
and save to JLD2.
"""
function save_daily_sco2(simulation, member_path, target_layer)
    (times, data) = ClimaLand.Diagnostics.diagnostic_as_vectors(
        simulation.diagnostics[1].output_writer,
        "sco2_ppm_30m_average";
        layer = target_layer,
    )

    model_dates_dt = times isa Vector{DateTime} ? times : date.(times)
    model_df = DataFrame(datetime = model_dates_dt, sco2_ppm = Float64.(data))
    model_df[!, :date] = Date.(model_df.datetime)

    model_daily = combine(
        groupby(model_df, :date),
        :sco2_ppm => mean => :daily_mean,
    )
    sort!(model_daily, :date)

    JLD2.jldsave(
        joinpath(member_path, "daily_diagnostics.jld2");
        dates = model_daily.date,
        sco2_ppm = model_daily.daily_mean,
    )
end

# ── Observation Map ──────────────────────────────────────────────────────────

"""
    ClimaCalibrate.observation_map(iteration)

Return G ensemble matrix matching the filtered observation dates.

The observation vector is daily mean soil CO₂ concentration (ppm) at ~2 cm depth.
"""
function ClimaCalibrate.observation_map(iteration)
    ekp = JLD2.load_object(ClimaCalibrate.ekp_path(OUTPUT_DIR, iteration))
    ensemble_size = EKP.get_N_ens(ekp)

    # Load valid observation dates
    obs_data = JLD2.load(OBS_FILEPATH)
    obs_dates = obs_data["obs_dates"]
    n_obs = length(obs_dates)

    @info "Observation map: $n_obs valid days, G length $n_obs"

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
            model_dates = member_data["dates"]
            sco2_model = member_data["sco2_ppm"]

            # Build date lookup
            model_dict = Dict(zip(model_dates, sco2_model))

            # Extract model values at observation dates
            sco2_out = [get(model_dict, d, NaN) for d in obs_dates]

            G_ens[:, m] = sco2_out
        catch e
            @error "Error processing member $m" exception = e
            G_ens[:, m] .= NaN
        end
    end

    return G_ens
end
