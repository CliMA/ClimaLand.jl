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
using ClimaLand.Canopy.PlantHydraulics
using ClimaLand.Snow
using ClimaCore
using ClimaDiagnostics
using ClimaUtilities
import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput, evaluate!
import ClimaUtilities.TimeManager: date
using Insolation

import EnsembleKalmanProcesses as EKP
import JLD2
using Dates
using Statistics
using DataFrames
using CSV

# ── Forward Model ────────────────────────────────────────────────────────────

function ClimaCalibrate.forward_model(iteration, member)
    FT = Float64
    site_ID_val = FluxnetSimulations.replace_hyphen(SITE_ID)
    climaland_dir = pkgdir(ClimaLand)

    # Get simulation dates from site metadata
    (; time_offset, lat, long) =
        FluxnetSimulations.get_location(FT, Val(site_ID_val))
    (start_date, stop_date) =
        FluxnetSimulations.get_data_dates(SITE_ID, time_offset)

    @info "Member $member: simulating $start_date to $stop_date"

    # Load calibrated parameters from TOML written by ClimaCalibrate
    calibrate_params_path =
        ClimaCalibrate.parameter_path(OUTPUT_DIR, iteration, member)
    toml_dict =
        LP.create_toml_dict(FT; override_files = [calibrate_params_path])

    # Domain
    (; dz_tuple, nelements, zmin, zmax) =
        FluxnetSimulations.get_domain_info(FT, Val(site_ID_val))
    (; atmos_h) =
        FluxnetSimulations.get_fluxtower_height(FT, Val(site_ID_val))

    land_domain = Column(;
        zlim = (zmin, zmax),
        nelements = nelements,
        dz_tuple = dz_tuple,
        longlat = (long, lat),
    )
    canopy_domain = ClimaLand.Domains.obtain_surface_domain(land_domain)
    surface_space = land_domain.space.surface

    # Determine target layer for soil CO₂ extraction (~6 cm depth)
    z_field = ClimaCore.Fields.coordinate_field(land_domain.space.subsurface).z
    z_vals = parent(z_field)[:, 1]
    target_depth = FT(-0.06)
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
    LAI = ClimaLand.Canopy.prescribed_lai_modis(surface_space, start_date, stop_date)

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

    canopy = ClimaLand.Canopy.CanopyModel{FT}(
        canopy_domain,
        canopy_forcing,
        LAI,
        toml_dict_base;
        prognostic_land_components,
        photosynthesis,
        conductance,
        soil_moisture_stress,
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

    # Custom IC with SOC exponential profile
    base_set_ic! = FluxnetSimulations.make_set_fluxnet_initial_conditions(
        SITE_ID,
        start_date,
        time_offset,
        land,
    )

    function custom_set_ic!(Y, p, t, model)
        base_set_ic!(Y, p, t, model)
        # Override soilCO2 ICs with exponential SOC profile
        Y.soilco2.CO2 .= FT(0.000412)
        Y.soilco2.O2_f .= FT(0.21)
        SOC_top = FT(18.0)
        SOC_bot = FT(0.5)
        τ_soc = FT(1.0 / log(SOC_top / SOC_bot))
        z = ClimaCore.Fields.coordinate_field(axes(Y.soilco2.SOC)).z
        @. Y.soilco2.SOC = SOC_bot + (SOC_top - SOC_bot) * exp(z / τ_soc)
    end

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

The observation vector is daily mean soil CO₂ concentration (ppm) at ~6 cm depth.
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
