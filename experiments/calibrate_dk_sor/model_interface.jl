"""
ClimaCalibrate model interface for DK-Sor single-site calibration.

Defines `forward_model(iteration, member)` and `observation_map(iteration)`.
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
using ClimaLand.Canopy.PlantHydraulics
using ClimaLand.Snow
using ClimaCore
using ClimaDiagnostics
using ClimaUtilities
import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput, evaluate!
import ClimaUtilities.TimeManager: date

import EnsembleKalmanProcesses as EKP
import JLD2
using Dates
using NCDatasets
using Statistics

# ── Forward Model ────────────────────────────────────────────────────────────

function ClimaCalibrate.forward_model(iteration, member)
    FT = Float64
    site_ID_val = FluxnetSimulations.replace_hyphen(SITE_ID)
    climaland_dir = pkgdir(ClimaLand)

    # Determine simulation window from current minibatch
    ekp = JLD2.load_object(ClimaCalibrate.ekp_path(OUTPUT_DIR, iteration))
    minibatch_indices = EKP.get_current_minibatch(ekp)

    obs_data = JLD2.load(OBS_FILEPATH)
    cal_years = obs_data["cal_years"]
    minibatch_years = [cal_years[idx] for idx in minibatch_indices]

    # 1-year spinup before first minibatch year, end after last minibatch year
    sim_start = DateTime(minimum(minibatch_years), 1, 1) - Year(1)
    sim_stop = DateTime(maximum(minibatch_years) + 1, 1, 1)
    @info "Member $member: simulating $sim_start to $sim_stop (minibatch years $minibatch_years)"

    # Load calibrated parameters from TOML written by ClimaCalibrate
    calibrate_params_path =
        ClimaCalibrate.parameter_path(OUTPUT_DIR, iteration, member)
    toml_dict =
        LP.create_toml_dict(FT; override_files = [calibrate_params_path])

    # Domain
    (; dz_tuple, nelements, zmin, zmax) =
        FluxnetSimulations.get_domain_info(FT, Val(site_ID_val))
    (; time_offset, lat, long) =
        FluxnetSimulations.get_location(FT, Val(site_ID_val))
    (; atmos_h) =
        FluxnetSimulations.get_fluxtower_height(FT, Val(site_ID_val))

    land_domain = Column(;
        zlim = (zmin, zmax),
        nelements = nelements,
        dz_tuple = dz_tuple,
        longlat = (long, lat),
    )
    canopy_domain = ClimaLand.Domains.obtain_surface_domain(land_domain)

    # Forcing from NetCDF
    met_nc_path = joinpath(
        climaland_dir,
        "DK_Sor",
        "DK-Sor_1997-2014_FLUXNET2015_Met.nc",
    )
    (; atmos, radiation) = FluxnetSimulations.prescribed_forcing_netcdf(
        met_nc_path,
        lat,
        long,
        time_offset,
        atmos_h,
        sim_start,
        toml_dict,
        FT,
    )

    # LAI from met NetCDF
    met_ds = NCDataset(met_nc_path, "r")
    lai_data = Float64.(coalesce.(met_ds["LAI"][1, 1, :], NaN))
    lai_times = met_ds["time"][:]
    close(met_ds)

    lai_seconds = [
        Float64(Second(t - Dates.Hour(time_offset) - sim_start).value)
        for t in lai_times
    ]
    valid_lai = .!isnan.(lai_data)
    LAI = TimeVaryingInput(lai_seconds[valid_lai], lai_data[valid_lai])

    # Build model
    prognostic_land_components = (:canopy, :snow, :soil, :soilco2)
    forcing_nt = (;
        atmos = atmos,
        radiation = radiation,
        ground = ClimaLand.PrognosticGroundConditions{FT}(),
    )

    biomass = Canopy.PrescribedBiomassModel{FT}(
        canopy_domain,
        LAI,
        toml_dict;
        height = FT(25),
    )

    canopy = Canopy.CanopyModel{FT}(
        canopy_domain,
        forcing_nt,
        LAI,
        toml_dict;
        prognostic_land_components,
        photosynthesis = Canopy.PModel{FT}(canopy_domain, toml_dict),
        conductance = Canopy.PModelConductance{FT}(toml_dict),
        soil_moisture_stress = Canopy.PiecewiseMoistureStressModel{FT}(
            land_domain,
            toml_dict,
        ),
        biomass,
    )

    land = LandModel{FT}(
        (; atmos, radiation),
        LAI,
        toml_dict,
        land_domain,
        DT;
        prognostic_land_components,
        canopy,
    )

    # Initial conditions
    function custom_set_ic!(Y, p, t, model)
        earth_param_set = ClimaLand.get_earth_param_set(model.soil)
        evaluate!(p.drivers.T, atmos.T, t)

        (; θ_r, ν, ρc_ds) = model.soil.parameters
        @. Y.soil.ϑ_l = θ_r + (ν - θ_r) / 2
        Y.soil.θ_i .= FT(0.0)
        ρc_s =
            ClimaLand.Soil.volumetric_heat_capacity.(
                Y.soil.ϑ_l,
                Y.soil.θ_i,
                ρc_ds,
                earth_param_set,
            )
        Y.soil.ρe_int .=
            ClimaLand.Soil.volumetric_internal_energy.(
                Y.soil.θ_i,
                ρc_s,
                p.drivers.T,
                earth_param_set,
            )

        Y.snow.S .= FT(0)
        Y.snow.S_l .= FT(0)
        Y.snow.U .= FT(0)
        if model.canopy.energy isa ClimaLand.Canopy.BigLeafEnergyModel
            Y.canopy.energy.T .= p.drivers.T
        end
        n_stem = model.canopy.hydraulics.n_stem
        n_leaf = model.canopy.hydraulics.n_leaf
        for i in 1:(n_stem + n_leaf)
            Y.canopy.hydraulics.ϑ_l.:($i) .=
                model.canopy.hydraulics.parameters.ν
        end

        if !isnothing(model.soilco2)
            Y.soilco2.CO2 .= FT(0.000412)
            Y.soilco2.O2_f .= FT(0.21)
            SOC_top = FT(15.0)
            SOC_bot = FT(0.5)
            τ_soc = FT(1.0 / log(SOC_top / SOC_bot))
            z = ClimaCore.Fields.coordinate_field(axes(Y.soilco2.SOC)).z
            @. Y.soilco2.SOC = SOC_bot + (SOC_top - SOC_bot) * exp(z / τ_soc)
        end
    end

    # Diagnostics — daily short diagnostics (includes nee, lhf, shf)
    output_writer = ClimaDiagnostics.Writers.DictWriter()
    diags = ClimaLand.default_diagnostics(
        land,
        sim_start;
        output_writer = output_writer,
        output_vars = :short,
        reduction_period = :daily,
    )

    simulation = LandSimulation(
        sim_start,
        sim_stop,
        DT,
        land;
        set_ic! = custom_set_ic!,
        updateat = Second(DT),
        diagnostics = diags,
    )
    solve!(simulation)

    # Extract daily diagnostics and save to JLD2
    member_path =
        ClimaCalibrate.path_to_ensemble_member(OUTPUT_DIR, iteration, member)
    save_daily_diagnostics(simulation, member_path)
    return nothing
end

"""
    save_daily_diagnostics(simulation, member_path)

Extract daily NEE, Qle, Qh from the DictWriter and save to JLD2.
"""
function save_daily_diagnostics(simulation, member_path)
    nee_dates, nee_vals = extract_daily_diag(simulation, "nee_1d_average")
    _, qle_vals = extract_daily_diag(simulation, "lhf_1d_average")
    _, qh_vals = extract_daily_diag(simulation, "shf_1d_average")

    JLD2.jldsave(
        joinpath(member_path, "daily_diagnostics.jld2");
        dates = nee_dates,
        nee = nee_vals,
        qle = qle_vals,
        qh = qh_vals,
    )
end

function extract_daily_diag(simulation, diag_name)
    writer = nothing
    for d in simulation.diagnostics
        if haskey(d.output_writer.dict, diag_name)
            writer = d.output_writer
            break
        end
    end
    if isnothing(writer)
        available = String[]
        for d in simulation.diagnostics
            append!(available, collect(keys(d.output_writer.dict)))
        end
        error("Diagnostic '$diag_name' not found. Available: $(unique(available))")
    end

    (times, data) = ClimaLand.Diagnostics.diagnostic_as_vectors(
        writer,
        diag_name,
    )
    model_dates_dt = times isa Vector{DateTime} ? times : date.(times)
    model_dates = Date.(model_dates_dt)
    return model_dates, Float64.(data)
end

# ── Observation Map ──────────────────────────────────────────────────────────

"""
    ClimaCalibrate.observation_map(iteration)

Return G ensemble matrix for the current minibatch's observation years.

The minibatch indices select which yearly Observations are active. We extract
model predictions only for those years' dates, concatenated in the same order
as the stacked observation vector.
"""
function ClimaCalibrate.observation_map(iteration)
    ekp = JLD2.load_object(ClimaCalibrate.ekp_path(OUTPUT_DIR, iteration))
    ensemble_size = EKP.get_N_ens(ekp)

    # Which year indices are in the current minibatch?
    minibatch_indices = EKP.get_current_minibatch(ekp)

    # Load per-year dates
    obs_data = JLD2.load(OBS_FILEPATH)
    year_dates = obs_data["year_dates"]
    cal_years = obs_data["cal_years"]

    # Collect per-year date vectors for minibatch years (preserving order)
    minibatch_years = [cal_years[idx] for idx in minibatch_indices]
    minibatch_year_dates = [year_dates[yr] for yr in minibatch_years]

    # Total G_ens size: each year contributes 3 * n_yr_dates entries
    # EKP stacks by concatenating full year vectors:
    # [yr1_nee, yr1_qle, yr1_qh, yr2_nee, yr2_qle, yr2_qh, ...]
    g_len = sum(3 * length(ds) for ds in minibatch_year_dates)
    @info "Observation map: minibatch years $minibatch_years, G length $g_len"

    G_ens = zeros(g_len, ensemble_size)

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
            nee_model = member_data["nee"]
            qle_model = member_data["qle"]
            qh_model = member_data["qh"]

            # Build date lookup
            model_dict_nee = Dict(zip(model_dates, nee_model))
            model_dict_qle = Dict(zip(model_dates, qle_model))
            model_dict_qh = Dict(zip(model_dates, qh_model))

            # Build G vector year-by-year matching EKP stacking order:
            # [yr1_nee, yr1_qle, yr1_qh, yr2_nee, yr2_qle, yr2_qh, ...]
            result = Float64[]
            for yr_dates in minibatch_year_dates
                nee_yr = [get(model_dict_nee, d, NaN) for d in yr_dates]
                qle_yr = [get(model_dict_qle, d, NaN) for d in yr_dates]
                qh_yr = [get(model_dict_qh, d, NaN) for d in yr_dates]
                # Convert NEE from mol CO2/m2/s to gC/m2/d
                nee_yr .*= 12.0 * 86400.0
                append!(result, nee_yr)
                append!(result, qle_yr)
                append!(result, qh_yr)
            end

            G_ens[:, m] = result
        catch e
            @error "Error processing member $m" exception = e
            G_ens[:, m] .= NaN
        end
    end

    return G_ens
end
