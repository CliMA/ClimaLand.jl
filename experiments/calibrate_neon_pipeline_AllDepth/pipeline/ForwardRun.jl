"""
ForwardRun.jl — `forward_run(run; output_dir, params)` runs the land model once
at a GIVEN parameter set and writes the diagnostic figures.

Function form of the original run_prior_mean_wLabile.jl (frozen in
../calibrate_neon). All configuration is passed as arguments — no ENV, no
top-level globals, no const re-declaration. `params` is a Dict{String,Float64}
of the parameter values to run at (typically the calibrated posterior means);
`labile_depth_scale` defaults to 0.0 when absent (=> exp(0·z)=1, plain model).

Why a forward run is needed at all: the calibration only persists the scalar
soil-CO₂ G vector per ensemble member, NOT the rich diagnostics (SWC, soil T,
O₂, SOC profile, CO₂ budget). So we run the model once at the optimized params
to produce them.

Returns `(; figures_dir)`. This file is `include`d into Main by the driver.
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
import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput, evaluate!
import ClimaUtilities.SpaceVaryingInputs: SpaceVaryingInput
import ClimaUtilities.Utils: searchsortednearest, linear_interpolation
import Interpolations
import ClimaUtilities.TimeManager: date
import ClimaParams as CP
using Insolation
using Dates
using Statistics
using CairoMakie
using CSV
using DataFrames
#import JLD2

# Site metadata helper (_get_neon_site_metadata). Included at TOP LEVEL — not
# inside forward_run — so it compiles once, not on every call. (It is also loaded
# by Calibration.jl; including a second time at top level is harmless but keeps
# ForwardRun.jl self-sufficient when used standalone.)
if !isdefined(@__MODULE__, :_get_neon_site_metadata)
    include(joinpath(pkgdir(ClimaLand), "experiments/calibrate_neon_pipeline_AllDepth/src/site_metadata.jl"))
end

# Default parameter values (used only if a param is missing from `params`).
const _FORWARD_DEFAULTS = Dict{String, Float64}(
    "soilCO2_reference_rate" => 5.1754e-7,
    "soilCO2_activation_energy" => 58230.0,
    "michaelis_constant" => 0.033968,
    "O2_michaelis_constant" => 0.03533,
    "labile_depth_scale" => 0.0,
)

_param(params, name) = Float64(get(params, name, _FORWARD_DEFAULTS[name]))

"""
    forward_run(run; output_dir, params) -> (; figures_dir)

Run the land model for `run`'s site/period at parameter values `params`, writing
all diagnostic figures into `output_dir`.
"""
function forward_run(run; output_dir, params)
    FT = Float64
    site_id = run.site
    spinup_days = run.spinup_days
    DT = run.dt
    climaland_dir = pkgdir(ClimaLand)

    # AllDepth: the calibrated depth codes (obs/model stacking order). Per-depth
    # figures loop over these; the CSV stats pool across them (see end of fn).
    depth_codes = run.cal_depth_codes

    start_date = DateTime(run.start_date)
    stop_date = DateTime(run.stop_date)
    spinup_date = start_date + Day(spinup_days)

    mkpath(output_dir)

    metadata = _get_neon_site_metadata(site_id)
    lat = FT(metadata.lat)
    long = FT(metadata.long)
    atmos_h = FT(metadata.atmos_h)
    time_offset = 0

    # ── Parameters ───────────────────────────────────────────────────────────
    soilCO2_reference_rate = _param(params, "soilCO2_reference_rate")
    soilCO2_activation_energy = _param(params, "soilCO2_activation_energy")
    michaelis_constant = _param(params, "michaelis_constant")
    O2_michaelis_constant = _param(params, "O2_michaelis_constant")
    labile_depth_scale = _param(params, "labile_depth_scale")
    println("forward_run params:")
    for (k, v) in (("soilCO2_reference_rate", soilCO2_reference_rate),
                   ("soilCO2_activation_energy", soilCO2_activation_energy),
                   ("michaelis_constant", michaelis_constant),
                   ("O2_michaelis_constant", O2_michaelis_constant),
                   ("labile_depth_scale", labile_depth_scale))
        println("  $k = $v")
    end

    prior_toml = joinpath(output_dir, "forward_parameters.toml")
    open(prior_toml, "w") do io
        write(io, """
[soilCO2_reference_rate]
value = $(soilCO2_reference_rate)
type = "float"
used_in = ["Land"]

[soilCO2_activation_energy]
value = $(soilCO2_activation_energy)
type = "float"
used_in = ["Land"]

[michaelis_constant]
value = $(michaelis_constant)
type = "float"
used_in = ["Land"]

[O2_michaelis_constant]
value = $(O2_michaelis_constant)
type = "float"
used_in = ["Land"]
""")
    end

    # ── Domain ───────────────────────────────────────────────────────────────
    dz_bottom = FT(2)
    dz_top = FT(0.038)
    dz_tuple = (dz_bottom, dz_top)
    nelements = 24
    zmin = FT(-6.2)
    zmax = FT(0)
    land_domain = Column(; zlim = (zmin, zmax), nelements = nelements,
        dz_tuple = dz_tuple, longlat = (long, lat))
    canopy_domain = ClimaLand.Domains.obtain_surface_domain(land_domain)
    surface_space = land_domain.space.surface

    z_field = ClimaCore.Fields.coordinate_field(land_domain.space.subsurface).z
    z_vals = parent(z_field)[:, 1]

    # Per-depth target layers, derived by argmin against the live grid from each
    # code's measurement depth (z_obs_m) in observations.jld2 — same source and
    # rule as the model interface, so figures/metrics match the calibration.
    obs_jld = joinpath(Config.base_dir_for(run), "observations.jld2")
    z_obs_by_code = Dict{String, Float64}()
    if isfile(obs_jld)
        _meta = JLD2.load(obs_jld)
        _mcodes = String.(_meta["depth_codes"])
        _mz = Float64.(_meta["z_obs_m"])
        for (c, z) in zip(_mcodes, _mz)
            z_obs_by_code[c] = z
        end
    else
        @warn "No observations.jld2 at $obs_jld; falling back to run.cal_depth for all codes"
        for c in depth_codes
            z_obs_by_code[c] = -Float64(run.cal_depth)
        end
    end
    target_layer_of = Dict(c => argmin(abs.(z_vals .- z_obs_by_code[c])) for c in depth_codes)
    for c in depth_codes
        println("Depth $c: z_obs=$(z_obs_by_code[c]) m → layer $(target_layer_of[c]) (z=$(z_vals[target_layer_of[c]]) m)")
    end

    toml_dict_base = LP.create_toml_dict(FT)
    toml_dict = LP.create_toml_dict(FT; override_files = [prior_toml])

    (; atmos, radiation) = FluxnetSimulations.prescribed_forcing_fluxnet(
        site_id, lat, long, time_offset, atmos_h, start_date, toml_dict_base, FT)

    toml_dict_base.data["canopy_d_coeff"]["value"] = FT(0.67)
    toml_dict_base.data["canopy_z_0b_coeff"]["value"] = FT(0.013)
    toml_dict_base.data["canopy_z_0m_coeff"]["value"] = FT(0.13)

    LAI = ClimaLand.Canopy.prescribed_climatological_lai_modis(surface_space)

    # ── Build model ──────────────────────────────────────────────────────────
    prognostic_land_components = (:canopy, :snow, :soil, :soilco2)
    forcing = (; atmos, radiation)
    ground = ClimaLand.PrognosticGroundConditions{FT}()
    canopy_forcing = (; atmos, radiation, ground)

    # Build the Rosetta van Genuchten retention parameters ONCE and share the
    # same object between the soil model (ν) and the canopy soil moisture stress
    # model (θ_high). LandModel's check_land_equality asserts θ_high == ν with
    # exact ==, so they must be the identical field, not two separate regrids.
    retention_parameters = Soil.soil_vangenuchten_parameters(#rosetta_soil_vangenuchten_parameters(
        land_domain.space.subsurface,
        FT,
    )

    photosynthesis = PModel{FT}(land_domain, toml_dict_base)
    conductance = PModelConductance{FT}(toml_dict_base)
    soil_moisture_stress =
        ClimaLand.Canopy.PiecewiseMoistureStressModel{FT}(
            land_domain,
            toml_dict_base;
            soil_params = retention_parameters,
    )
    canopy = ClimaLand.Canopy.CanopyModel{FT}(
        canopy_domain, canopy_forcing, LAI, toml_dict_base;
        prognostic_land_components, photosynthesis, conductance, soil_moisture_stress)

    α_snow = Snow.ZenithAngleAlbedoModel(toml_dict_base)
    snow = Snow.SnowModel(FT, canopy_domain, forcing, toml_dict_base, DT;
        prognostic_land_components, α_snow)

    soil = Soil.EnergyHydrology{FT}(
        land_domain,
        forcing,
        toml_dict;
        prognostic_land_components,
        additional_sources = (ClimaLand.RootExtraction{FT}(),),
        retention_parameters = retention_parameters,
    )

    land = LandModel{FT}(forcing, LAI, toml_dict, land_domain, DT;
        prognostic_land_components, snow, canopy, soil)

    porosity_scale = FT(1)
    land.soil.parameters.ν .*= porosity_scale

    # Optional per-run θ_r override (uniform over depth). Applied before solve!
    # so both the set_ic! clamp (θ_r + 1e-4 floor) and the model dynamics use it.
    if run.theta_r !== nothing
        land.soil.parameters.θ_r .= FT(run.theta_r)
        println("Overriding θ_r with config value: $(run.theta_r)")
    end

    # ── Initial conditions (SOC × labile factor + NEON SWC profile) ──────────
    base_set_ic! = FluxnetSimulations.make_set_fluxnet_initial_conditions(
        site_id, start_date, time_offset, land)

    # Capture labile rate as a concrete-typed local so the map closure stays
    # type-stable (see original note re: ClimaCore broadcast eltype_error).
    k_labile = FT(labile_depth_scale)

    function set_ic!(Y, p, t, model)
        base_set_ic!(Y, p, t, model)

        # Override the base O₂ initial condition (0.08 kg O₂ m⁻³ bulk soil) with
        # zero everywhere; the near-surface O₂ is then set by the diffusive top
        # boundary condition once the run starts.
        #Y.soilco2.O2 .= FT(0)

        soc_field = ClimaCore.Fields.zeros(land_domain.space.subsurface)
        soc_data = CSV.read(
            "/kiwi-data/Data/groupMembers/evametz/Neon/Neon_data/NEON_all_sites_estimatedOC_2cm_mean.csv",
            DataFrame)
        valid_soc = .!ismissing.(soc_data[!, "$(site_id)_estimatedOC_kg_m3"])
        raw_z::Vector{Float64} = Float64.(soc_data.depth[valid_soc])
        sort_idx_soc = sortperm(raw_z)
        raw_vals::Vector{Float64} =
            Float64.(soc_data[valid_soc, "$(site_id)_estimatedOC_kg_m3"])

        z_extrap_top = (raw_z[sort_idx_soc])[1]
        SOC_extrap_top = (raw_vals[sort_idx_soc])[1]
        SOC_extrap_bot = FT(0.05)
        z_extrap_bot = minimum(parent(
            ClimaCore.Fields.coordinate_field(land_domain.space.subsurface).z))
        zvalues = ClimaCore.Fields.coordinate_field(axes(Y.soilco2.SOC)).z
        alpha_soc = FT(log(SOC_extrap_top / SOC_extrap_bot) / (z_extrap_bot - z_extrap_top))

        soc_field .= map(zvalues) do z
            soc = if z > z_extrap_top
                linear_interpolation(raw_z[sort_idx_soc], raw_vals[sort_idx_soc], z)
            else
                SOC_extrap_top * exp(-alpha_soc * (z - z_extrap_top))
            end
            soc * exp(k_labile * z)
        end
        Y.soilco2.SOC .= soc_field

        neon_depths = FT[-0.06, -0.16, -0.26, -0.46, -0.66, -0.86, -1.06, -1.66]
        # NOTE: renamed from `depth_codes` to avoid shadowing forward_run's
        # `depth_codes` (the calibrated CO₂ codes) captured by this closure —
        # otherwise this 8-element list rebinds it and the per-depth loop breaks.
        swc_depth_codes = ["501", "502", "503", "504", "505", "506", "507", "508"]
        n_plots = 5
        csv_path = ClimaLand.Artifacts.experiment_fluxnet_data_path(site_id)
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
        @. Y.soil.ϑ_l = clamp(Y.soil.ϑ_l, θ_r_field + FT(1e-4), ν_field - FT(1e-4))
    end

    # ── Diagnostics ──────────────────────────────────────────────────────────
    output_writer = ClimaDiagnostics.Writers.DictWriter()
    output_vars = ["swc", "tsoil", "si", "sco2", "soc", "hr", "so2", "sco2_ppm", "scd", "scms", "tair", "airp"]
    diags = ClimaLand.default_diagnostics(land, start_date;
        output_writer = output_writer, output_vars, reduction_period = :halfhourly)

    simulation = LandSimulation(start_date, stop_date, DT, land;
        set_ic! = set_ic!, updateat = Second(DT), diagnostics = diags)

    println("Running forward model at given parameters …")
    @time solve!(simulation)
    println("Done.")

    # ── Extract diagnostics ──────────────────────────────────────────────────
    function get_diag_series(sim, name, layer)
        (times, data) = ClimaLand.Diagnostics.diagnostic_as_vectors(
            sim.diagnostics[1].output_writer, name; layer = layer)
        model_dates = times isa Vector{DateTime} ? times : date.(times)
        df = DataFrame(datetime = model_dates, value = Float64.(data))
        df[!, :date] = Date.(df.datetime)
        df = filter(row -> row.date >= Date(spinup_date), df)
        sort!(df, :datetime)
        return df
    end
    function get_diag_layer(sim, name, layer)
        df = get_diag_series(sim, name, layer)
        daily = combine(groupby(df, :date), :value => mean => :daily_mean)
        sort!(daily, :date)
        return daily
    end

    # SOC profile (time-invariant; first timestep)
    soc_vals = Float64[]
    for layer in 1:nelements
        (_, data) = ClimaLand.Diagnostics.diagnostic_as_vectors(
            simulation.diagnostics[1].output_writer, "soc_30m_average"; layer = layer)
        push!(soc_vals, data[1])
    end
    fig_soc = Figure(size = (500, 700))
    ax_soc = Axis(fig_soc[1, 1]; xlabel = "SOC (kg C m⁻³)", ylabel = "Depth (m)",
        title = "$site_id — SOC Profile")
    lines!(ax_soc, soc_vals, z_vals; color = :brown, linewidth = 2.0, label = "Model SOC")
    scatter!(ax_soc, soc_vals, z_vals; color = :brown, markersize = 6)
    axislegend(ax_soc; position = :rb, framevisible = false)
    CairoMakie.save(joinpath(output_dir, "soc_profile_$(site_id).png"), fig_soc)

    # ── NEON observations ────────────────────────────────────────────────────
    # obs_df + its datetime/date columns + rowmean_skipinvalid are loop-invariant
    # (shared across depths). Everything that depends on the depth code is built
    # inside the per-depth loop below.
    csv_path = ClimaLand.Artifacts.experiment_fluxnet_data_path(site_id)
    obs_df = CSV.read(csv_path, DataFrame)
    function rowmean_skipinvalid(row, cols)
        vals = Float64[]
        for c in cols
            v = row[c]
            (!ismissing(v) && !isnan(Float64(v))) && push!(vals, Float64(v))
        end
        return isempty(vals) ? NaN : mean(vals)
    end
    obs_df[!, :datetime] =
        DateTime.(string.(Int.(obs_df.timestamp_fmt)), dateformat"yyyymmddHHMM")
    obs_df[!, :date] = Date.(obs_df.datetime)

    # ── Per-depth loop ───────────────────────────────────────────────────────
    # One iteration per calibrated depth code. Builds that depth's obs+model daily
    # series, writes the 3 per-depth figures (forward_mean, optimized_vs_sensors,
    # co2_budget), and collects per-depth stats. The all-layer profile/O₂ figures
    # (which don't depend on the code) are written once, guarded to the 1st pass.
    # After the loop, per-depth stats are POOLED into the CSV `scatter_stats`.
    per_depth_stats = NamedTuple[]     # one entry per code (for pooling)
    per_depth_porosity = Float64[]     # ν at each depth's target layer
    for code in depth_codes
        obs_depth = code
        target_layer = target_layer_of[code]
        println("── depth $code → layer $target_layer (z=$(z_vals[target_layer]) m) ──")

        sco2_daily = get_diag_layer(simulation, "sco2_ppm_30m_average", target_layer)
        swc_daily = get_diag_layer(simulation, "swc_30m_average", target_layer)
        tsoil_daily = get_diag_layer(simulation, "tsoil_30m_average", target_layer)

        co2_cols = [Symbol("soilCO2concentrationMean_$(lpad(p,3,'0'))_$obs_depth") for p in 1:5]
        obs_df[!, :sco2_mean] = [rowmean_skipinvalid(row, co2_cols) for row in eachrow(obs_df)]
    obs_daily = combine(groupby(obs_df, :date),
        :sco2_mean => (x -> begin
            valid = filter(!isnan, x); length(valid) >= 24 ? mean(valid) : NaN
        end) => :daily_mean)
    obs_daily = filter(row -> row.date >= Date(spinup_date), obs_daily)
    obs_daily = filter(row -> !isnan(row.daily_mean), obs_daily)
    sort!(obs_daily, :date)

    # ── Cap mask: align the obs used for METRICS with the calibration ─────────
    # observations.jld2 stores obs_dates = the soil-CO₂ days that survived the
    # spring-peak mask (Observation_flag.jl). The forward RMSE / correlation /
    # obs-CO₂ means must be computed on the SAME set, otherwise capped spring
    # peaks (which the calibration never saw) inflate the error and bias the
    # means. `obs_daily_cal` is obs_daily restricted to the kept dates; the full
    # `obs_daily` is kept for the time-series figures (so the dropped days stay
    # visible). `capped_co2_dates` are the in-window obs-days that were dropped,
    # used to shade the figures. When the unmasked Observations.jl is used,
    # kept == all obs days, so obs_daily_cal == obs_daily and nothing is shaded.
    obs_jld = joinpath(Config.base_dir_for(run), "observations.jld2")
    if isfile(obs_jld)
        kept_dates = Set(Date.(JLD2.load(obs_jld, "obs_dates")))
        obs_daily_cal = filter(row -> row.date in kept_dates, obs_daily)
        capped_co2_dates = sort!([r.date for r in eachrow(obs_daily)
                                  if r.date <= Date(stop_date) && !(r.date in kept_dates)])
    else
        @warn "No observations.jld2 at $obs_jld; forward metrics use the FULL (unmasked) obs"
        obs_daily_cal = obs_daily
        capped_co2_dates = Date[]
    end

    # Observed SWC daily means at the calibration depth (same depth code + ≥24
    # valid half-hours/day filter as the CO₂ obs above), for the scatter figure.
    swc_obs_cols = [c for c in
        (Symbol("VSWCMean_$(lpad(p,3,'0'))_$obs_depth") for p in 1:5)
        if String(c) in names(obs_df)]
    obs_df[!, :swc_mean] = isempty(swc_obs_cols) ? fill(NaN, nrow(obs_df)) :
        [rowmean_skipinvalid(row, swc_obs_cols) for row in eachrow(obs_df)]
    obs_swc_daily = combine(groupby(obs_df, :date),
        :swc_mean => (x -> begin
            valid = filter(!isnan, x); length(valid) >= 24 ? mean(valid) : NaN
        end) => :daily_mean)
    obs_swc_daily = filter(row -> row.date >= Date(spinup_date), obs_swc_daily)
    obs_swc_daily = filter(row -> !isnan(row.daily_mean), obs_swc_daily)
    sort!(obs_swc_daily, :date)

    # Observed soil temperature daily means at the calibration depth. NEON obs
    # soilTemp is in °C; convert to K (+273.15) to match the model tsoil (K).
    tsoil_obs_cols = [c for c in
        (Symbol("soilTempMean_$(lpad(p,3,'0'))_$obs_depth") for p in 1:5)
        if String(c) in names(obs_df)]
    obs_df[!, :tsoil_mean] = isempty(tsoil_obs_cols) ? fill(NaN, nrow(obs_df)) :
        [rowmean_skipinvalid(row, tsoil_obs_cols) for row in eachrow(obs_df)]
    obs_tsoil_daily = combine(groupby(obs_df, :date),
        :tsoil_mean => (x -> begin
            valid = filter(!isnan, x); length(valid) >= 24 ? mean(valid) : NaN
        end) => :daily_mean)
    obs_tsoil_daily = filter(row -> row.date >= Date(spinup_date), obs_tsoil_daily)
    obs_tsoil_daily = filter(row -> !isnan(row.daily_mean), obs_tsoil_daily)
    obs_tsoil_daily.daily_mean .+= 273.15
    sort!(obs_tsoil_daily, :date)

    # Mark the capped spring-peak days (dropped from the calibration / metrics)
    # with red dots on the obs CO₂ curve. Implemented as scatter! rather than a
    # vspan!/vlines! band because those primitives don't accept Date/DateTime
    # x-coordinates in this Makie version (only scatter!/lines! do). `obs_ref` is
    # the daily obs frame (date => daily_mean) the dots sit on. No-op when empty.
    function _mark_capped!(ax, dates, obs_ref)
        isempty(dates) && return
        capset = Set(dates)
        sub = filter(r -> r.date in capset, obs_ref)
        nrow(sub) == 0 && return
        scatter!(ax, sub.date, Float64.(sub.daily_mean);
            color = (:blue, 0.8), markersize = 10, marker = :xcross,
            label = "capped (excluded from RMSE)")
    end

    # ── Main figure: model vs obs + Tsoil ──────────────────────────────
    fig = Figure(size = (1200, 1000))
    ax1 = Axis(fig[1, 1]; ylabel = "Soil CO₂ (ppm)",
        title = "$site_id — Forward Model vs NEON Obs (depth $obs_depth)")
    lines!(ax1, sco2_daily.date, sco2_daily.daily_mean; color = :blue, linewidth = 1.5,
        label = "Forward model")
    lines!(ax1, obs_daily.date, Float64.(obs_daily.daily_mean); color = :black,
        linewidth = 1.5, label = "NEON obs")
    _mark_capped!(ax1, capped_co2_dates, obs_daily)   # red ×: days excluded from RMSE
    axislegend(ax1; position = :rt, framevisible = false)
    ax2 = Axis(fig[2, 1]; ylabel = "SWC", xlabel = "Date")
    lines!(ax2, swc_daily.date, swc_daily.daily_mean; color = :blue, linewidth = 1.5, label = "SWC (model)")
    axislegend(ax2; position = :rt, framevisible = false)
    ax3 = Axis(fig[3, 1]; xlabel = "Date", ylabel = "T (K)")
    lines!(ax3, tsoil_daily.date, tsoil_daily.daily_mean; color = :red, linewidth = 1.5, label = "T (model)")
    axislegend(ax3; position = :rt, framevisible = false)
    if nrow(swc_daily) > 0
        xl = extrema(swc_daily.date)
        xlims!(ax1, xl...); xlims!(ax2, xl...); xlims!(ax3, xl...)
    end
    main_path = joinpath(output_dir, "forward_mean_$(site_id)_$(code).png")
    CairoMakie.save(main_path, fig)
    println("Saved: $main_path")

    # ── Optimized model vs individual NEON sensor plots ──────────────────────
    function _daily_from_col(df, col)
        valid = .!ismissing.(df[!, col]) .& .!isnan.(Float64.(coalesce.(df[!, col], NaN)))
        isempty(findall(valid)) && return DataFrame(date = Date[], daily_mean = Float64[])
        sub = DataFrame(date = df.date[valid], value = Float64.(df[valid, col]))
        daily = combine(groupby(sub, :date), :value => (x -> begin
            v = filter(!isnan, x); length(v) >= 24 ? mean(v) : NaN
        end) => :daily_mean)
        daily = filter(r -> !isnan(r.daily_mean), daily)
        daily = filter(r -> r.date >= Date(spinup_date), daily)
        sort!(daily, :date)
        return daily
    end

    n_plots_obs = 5
    co2_plot_cols = [Symbol("soilCO2concentrationMean_$(lpad(p,3,'0'))_$obs_depth") for p in 1:n_plots_obs]
    swc_plot_cols = [Symbol("VSWCMean_$(lpad(p,3,'0'))_$obs_depth") for p in 1:n_plots_obs]
    tsoil_plot_cols = [Symbol("soilTempMean_$(lpad(p,3,'0'))_$obs_depth") for p in 1:n_plots_obs]
    obs_colnames = names(obs_df)
    obs_color = (:gray40, 0.6)

    fig_opt = Figure(size = (1200, 1000))
    ax_co2 = Axis(fig_opt[1, 1]; ylabel = "Soil CO₂ (ppm)",
        title = "$site_id — Optimized model vs individual NEON sensors (depth $obs_depth)")
    ax_swc = Axis(fig_opt[2, 1]; ylabel = "SWC (m³/m³)")
    ax_t = Axis(fig_opt[3, 1]; ylabel = "T (K)", xlabel = "Date")
    for (i, col) in enumerate(co2_plot_cols)
        String(col) in obs_colnames || continue
        d = _daily_from_col(obs_df, col); nrow(d) == 0 && continue
        lines!(ax_co2, d.date, d.daily_mean; color = obs_color, linewidth = 1.0,
            label = i == 1 ? "NEON sensor plots (001–$(lpad(n_plots_obs,3,'0')))" : nothing)
    end
    lines!(ax_co2, sco2_daily.date, sco2_daily.daily_mean; color = :firebrick,
        linewidth = 1.8, label = "Optimized model")
    # Mark the capped spring-peak days. The CSV/sensor lines still span the full
    # year (drivers untouched); the red × just flags which days the mask excluded.
    _mark_capped!(ax_co2, capped_co2_dates, obs_daily)
    axislegend(ax_co2; position = :rt, framevisible = false)
    for (i, col) in enumerate(swc_plot_cols)
        String(col) in obs_colnames || continue
        d = _daily_from_col(obs_df, col); nrow(d) == 0 && continue
        lines!(ax_swc, d.date, d.daily_mean; color = obs_color, linewidth = 1.0,
            label = i == 1 ? "NEON sensor plots" : nothing)
    end
    lines!(ax_swc, swc_daily.date, swc_daily.daily_mean; color = :firebrick,
        linewidth = 1.8, label = "Optimized model")
    axislegend(ax_swc; position = :rt, framevisible = false)
    for (i, col) in enumerate(tsoil_plot_cols)
        String(col) in obs_colnames || continue
        d = _daily_from_col(obs_df, col); nrow(d) == 0 && continue
        d.daily_mean .+= 273.15
        lines!(ax_t, d.date, d.daily_mean; color = obs_color, linewidth = 1.0,
            label = i == 1 ? "NEON sensor plots" : nothing)
    end
    lines!(ax_t, tsoil_daily.date, tsoil_daily.daily_mean; color = :firebrick,
        linewidth = 1.8, label = "Optimized model")
    axislegend(ax_t; position = :rt, framevisible = false)
    if nrow(swc_daily) > 0
        xl = extrema(swc_daily.date)
        xlims!(ax_co2, xl...); xlims!(ax_swc, xl...); xlims!(ax_t, xl...)
    end
    CairoMakie.save(joinpath(output_dir, "optimized_vs_sensors_$(site_id)_$(code).png"), fig_opt)

    # ── All-layer profile / O₂ figures: depth-code-independent, write ONCE ────
    if code == depth_codes[1]
    # SoilCO2 O₂ constants for converting fractions ↔ gas-phase mass conc.
    soilco2_params = land.soilco2.parameters
    M_O2 = FT(soilco2_params.M_O2)
    M_C = FT(soilco2_params.M_C)
    O2_f_atm = FT(soilco2_params.O2_f_atm)
    R_gas = FT(LP.gas_constant(soilco2_params.earth_param_set))
    O2_f_lim = FT(1e-4)            # MM attenuation floor (matches RHS, Biogeochemistry.jl)
    stoich_O2_per_C = M_O2 / M_C  # kg O₂ consumed per kg C respired (= 8/3)

    # ── Per-layer profile figure: CO₂ / O₂ / microbial source / SWC by depth ──
    colors = cgrad(:viridis, nelements, categorical = true)
    fig_prof = Figure(size = (1200, 1500))
    axp1 = Axis(fig_prof[1, 1]; ylabel = "Soil CO₂ (ppm)",
        title = "$site_id — Per-layer profiles (forward run)")
    axp2 = Axis(fig_prof[2, 1]; ylabel = "Soil O₂")
    axp3 = Axis(fig_prof[3, 1]; ylabel = "Soil CO₂ Microbial Source")
    axp4 = Axis(fig_prof[4, 1]; xlabel = "Date", ylabel = "Soil Water Content")
    # ── O₂ diagnostic figure (axp5 + axp6 + axp7), separated from fig_prof ──
    fig_o2diag = Figure(size = (1200, 1500))
    # Gas-phase O₂ mass concentration, kg O₂ m⁻³, matching the top-BC variable
    #   O2_f_atm · P_sfc · M_O2 / (R · T_sfc).  Per-layer model value is the
    #   gas-phase concentration recovered from the `so2` volume fraction:
    #   ρ_O2_air = so2 · P · M_O2 / (R · T_soil)  (= Y.soilco2.O2 / θ_eff_o2).
    axp5 = Axis(fig_o2diag[1, 1]; ylabel = "O₂ gas-phase conc (kg m⁻³)")
    # Microbial O₂ tendency (the dY.soilco2.O2 sink term, kg O₂ m⁻³ s⁻¹):
    #   −stoich_O2_per_C · Sm · O2_f/(O2_f + O2_f_lim),  with Sm = `scms`,
    #   O2_f = `so2`. This is the consumption term only (excludes diffusion).
    axp6 = Axis(fig_o2diag[2, 1]; ylabel = "Microbial O₂ tendency (kg m⁻³ s⁻¹)")
    # O₂ volumetric fraction in soil air (the `so2` diagnostic, m³ m⁻³).
    axp7 = Axis(fig_o2diag[3, 1]; xlabel = "Date", ylabel = "Soil O₂ fraction (m³ m⁻³)")
    local last_swc_layer = DataFrame()
    airp_l = get_diag_layer(simulation, "airp_30m_average", 1)  # surface P (z-invariant)
    tair_l = get_diag_layer(simulation, "tair_30m_average", 1)  # surface air T = T_sfc
    for layer in 1:nelements
        sco2_l = get_diag_layer(simulation, "sco2_ppm_30m_average", layer)
        so2_l = get_diag_layer(simulation, "so2_30m_average", layer)
        scms_l = get_diag_layer(simulation, "scms_30m_average", layer)
        swc_l = get_diag_layer(simulation, "swc_30m_average", layer)
        tsoil_l = get_diag_layer(simulation, "tsoil_30m_average", layer)
        last_swc_layer = swc_l
        lines!(axp1, sco2_l.date, sco2_l.daily_mean; linewidth = 1.5, color = colors[layer], label = "layer$layer")
        lines!(axp2, so2_l.date, so2_l.daily_mean; linewidth = 1.5, color = colors[layer], label = "layer$layer")
        lines!(axp3, scms_l.date, scms_l.daily_mean; linewidth = 1.5, color = colors[layer], label = "layer$layer")
        lines!(axp4, swc_l.date, swc_l.daily_mean; linewidth = 1.5, color = colors[layer], label = "layer$layer")
        # ── O₂ diagnostics (axp5 + axp6) → fig_o2diag ──
        # gas-phase O₂ mass conc per layer (joined to airp on date)
        o2c = innerjoin(so2_l, tsoil_l, airp_l, on = :date,
            makeunique = true)
        o2_conc = o2c.daily_mean .* o2c.daily_mean_2 .* M_O2 ./ (R_gas .* o2c.daily_mean_1)
        lines!(axp5, o2c.date, o2_conc; linewidth = 1.5, color = colors[layer], label = "layer$layer")
        # microbial O₂ tendency per layer: −stoich · Sm · O2_f/(O2_f + O2_f_lim)
        # (so2_l = O2_f, scms_l = Sm; daily_mean is Sm after join, daily_mean_1 is O2_f)
        do2 = innerjoin(scms_l, so2_l, on = :date, makeunique = true)
        o2_tend = .-stoich_O2_per_C .* do2.daily_mean .*
            do2.daily_mean_1 ./ (do2.daily_mean_1 .+ O2_f_lim)
        lines!(axp6, do2.date, o2_tend; linewidth = 1.5, color = colors[layer], label = "layer$layer")
        # panel 3: O₂ volumetric fraction (so2) per layer
        lines!(axp7, so2_l.date, so2_l.daily_mean; linewidth = 1.5, color = colors[layer], label = "layer$layer")
    end
    # Atmospheric reference: O2_f_atm · P_sfc · M_O2 / (R · T_sfc).
    atm = innerjoin(airp_l, tair_l, on = :date, makeunique = true)
    o2_atm = O2_f_atm .* atm.daily_mean .* M_O2 ./ (R_gas .* atm.daily_mean_1)
    lines!(axp5, atm.date, o2_atm; linewidth = 2.5, color = :black,
        linestyle = :dash, label = "atmosphere")
    axislegend(axp2; position = :rt, framevisible = false)
    axislegend(axp5; position = :rt, framevisible = false)
    if nrow(last_swc_layer) > 0
        xl = extrema(last_swc_layer.date)
        xlims!(axp1, xl...); xlims!(axp2, xl...); xlims!(axp3, xl...)
        xlims!(axp4, xl...)
        xlims!(axp5, xl...); xlims!(axp6, xl...)
        xlims!(axp7, xl...)
    end
    CairoMakie.save(joinpath(output_dir, "profiles_O2_co2_cms_$(site_id).png"), fig_prof)
    CairoMakie.save(joinpath(output_dir, "o2_diagnostics_$(site_id).png"), fig_o2diag)
    end  # if code == depth_codes[1]  (all-layer profile/O₂ figures, once)

    # ── CO₂ budget at target layer: production / emission / transport ────────
    # 1. Microbial CO₂ production in target layer (kg C m⁻³ s⁻¹), half-hourly.
    scms_hh = get_diag_series(simulation, "scms_30m_average", target_layer)

    # 2. CO₂ emission to atmosphere: surface flux hr (mol CO₂ m⁻² s⁻¹), no layer.
    (hr_times, hr_data) = ClimaLand.Diagnostics.diagnostic_as_vectors(
        simulation.diagnostics[1].output_writer, "hr_30m_average")
    hr_dates = hr_times isa Vector{DateTime} ? hr_times : date.(hr_times)
    hr_df = DataFrame(datetime = hr_dates, value = Float64.(hr_data))
    hr_df[!, :date] = Date.(hr_df.datetime)
    hr_df = filter(row -> row.date >= Date(spinup_date), hr_df)
    sort!(hr_df, :datetime)

    # 3. Diffusive CO₂ transport from below: F = D·Δ(CO2_air_eq)/Δz, with
    #    CO2_air_eq = CO2 / max(ν − swc, ε). Positive = upward into target layer.
    if target_layer > 1
        scd_hh_tl = get_diag_series(simulation, "scd_30m_average", target_layer)
        sco2_hh_tl = get_diag_series(simulation, "sco2_30m_average", target_layer)
        sco2_hh_bl = get_diag_series(simulation, "sco2_30m_average", target_layer - 1)
        swc_hh_tl = get_diag_series(simulation, "swc_30m_average", target_layer)
        swc_hh_bl = get_diag_series(simulation, "swc_30m_average", target_layer - 1)
        ν_tl = parent(land.soil.parameters.ν)[target_layer]
        ν_bl = parent(land.soil.parameters.ν)[target_layer - 1]
        eps_val = 1e-10
        co2_air_eq_tl = sco2_hh_tl.value ./ max.(ν_tl .- swc_hh_tl.value, eps_val)
        co2_air_eq_bl = sco2_hh_bl.value ./ max.(ν_bl .- swc_hh_bl.value, eps_val)
        dz = abs(z_vals[target_layer] - z_vals[target_layer - 1])
        transport_vals = scd_hh_tl.value .* (co2_air_eq_bl .- co2_air_eq_tl) ./ dz
        transport_df = DataFrame(datetime = scd_hh_tl.datetime,
            date = scd_hh_tl.date, value = transport_vals)
    else
        transport_df = DataFrame(datetime = scms_hh.datetime,
            date = scms_hh.date, value = zeros(nrow(scms_hh)))
    end

    fig_budget = Figure(size = (1200, 800))
    ax_prod = Axis(fig_budget[1, 1]; ylabel = "CO₂ Production\n(kg C m⁻³ s⁻¹)",
        title = "$site_id — Layer CO₂ Budget (z = $(round(z_vals[target_layer], digits=3)) m)")
    lines!(ax_prod, scms_hh.datetime, scms_hh.value; color = :darkgreen,
        linewidth = 1.5, label = "Microbial production")
    axislegend(ax_prod; position = :rt, framevisible = false)
    ax_emis = Axis(fig_budget[2, 1]; ylabel = "CO₂ Emission\n(mol CO₂ m⁻² s⁻¹)")
    lines!(ax_emis, hr_df.datetime, hr_df.value; color = :firebrick,
        linewidth = 1.5, label = "Emission to atmosphere (HR)")
    axislegend(ax_emis; position = :rt, framevisible = false)
    ax_trans = Axis(fig_budget[3, 1]; xlabel = "Date",
        ylabel = "CO₂ Transport from below\n(kg CO₂ m⁻³ s⁻¹)")
    lines!(ax_trans, transport_df.datetime, transport_df.value; color = :royalblue,
        linewidth = 1.5, label = "Diffusive transport from below")
    hlines!(ax_trans, [0.0]; color = :gray, linewidth = 0.8, linestyle = :dash)
    axislegend(ax_trans; position = :rt, framevisible = false)
    if nrow(scms_hh) > 0
        t_xlim = extrema(scms_hh.datetime)
        xlims!(ax_prod, t_xlim...); xlims!(ax_emis, t_xlim...); xlims!(ax_trans, t_xlim...)
    end
    CairoMakie.save(joinpath(output_dir, "co2_budget_$(site_id)_$(code).png"), fig_budget)

    # ── Collect this depth's series/porosity for pooled CSV stats (after loop) ──
    push!(per_depth_stats, (;
        code = code,
        sco2_daily = sco2_daily, swc_daily = swc_daily, tsoil_daily = tsoil_daily,
        obs_daily_cal = obs_daily_cal, obs_swc_daily = obs_swc_daily,
        obs_tsoil_daily = obs_tsoil_daily,
    ))
    push!(per_depth_porosity, Float64(parent(land.soil.parameters.ν)[target_layer]))
    end  # for code in depth_codes

    # ── Scatter figure: obs/model CO₂ vs SWC relationships (DISABLED for AllDepth) ──
    # Each panel pairs two daily-mean series on common dates. Title carries the
    # Pearson correlation r and the RMSE of the (x, y) pair.
    #   model CO₂  = sco2_daily,  model SWC = swc_daily,
    #   obs   CO₂  = obs_daily,   obs   SWC = obs_swc_daily.
    function _pair_on_date(dfx, dfy)
        (nrow(dfx) == 0 || nrow(dfy) == 0) && return (Float64[], Float64[])
        j = innerjoin(
            rename(dfx[:, [:date, :daily_mean]], :daily_mean => :x),
            rename(dfy[:, [:date, :daily_mean]], :daily_mean => :y),
            on = :date)
        return (Float64.(j.x), Float64.(j.y))
    end
    _corr(x, y) = length(x) >= 2 && std(x) > 0 && std(y) > 0 ? cor(x, y) : NaN
    _rmse(x, y) = isempty(x) ? NaN : sqrt(mean((x .- y) .^ 2))
    function _scatter_panel!(ax, x, y; with_rmse, base_title)
        r = _corr(x, y)
        t = "$base_title\nr = $(round(r, digits=3)) (n = $(length(x)))"
        if with_rmse
            t *= ", RMSE = $(round(_rmse(x, y), sigdigits=3))"
        end
        ax.title = t
        isempty(x) || scatter!(ax, x, y; color = (:steelblue, 0.6), markersize = 6)
    end

    # Scatter figure disabled for the AllDepth pipeline (per-depth series no longer
    # live at this scope; re-enable per depth inside the loop if needed).
    if false
        obs_mod_co2_x, obs_mod_co2_y = _pair_on_date(obs_daily_cal, sco2_daily)
        obs_mod_swc_x, obs_mod_swc_y = _pair_on_date(obs_swc_daily, swc_daily)
        obs_co2_swc_x, obs_co2_swc_y = _pair_on_date(obs_daily_cal, obs_swc_daily)
        mod_co2_swc_x, mod_co2_swc_y = _pair_on_date(sco2_daily, swc_daily)

        fig_scat = Figure(size = (1100, 950))
        Label(fig_scat[0, 1:2], "$site_id — CO₂ / SWC scatter relationships";
            fontsize = 18, font = :bold)
        sax1 = Axis(fig_scat[1, 1]; xlabel = "Obs soil CO₂ (ppm)", ylabel = "Model soil CO₂ (ppm)")
        _scatter_panel!(sax1, obs_mod_co2_x, obs_mod_co2_y;
            with_rmse = true, base_title = "Obs vs model soil CO₂")
        sax2 = Axis(fig_scat[1, 2]; xlabel = "Obs SWC (m³/m³)", ylabel = "Model SWC (m³/m³)")
        _scatter_panel!(sax2, obs_mod_swc_x, obs_mod_swc_y;
            with_rmse = true, base_title = "Obs vs model soil water content")
        sax3 = Axis(fig_scat[2, 1]; xlabel = "Obs soil CO₂ (ppm)", ylabel = "Obs SWC (m³/m³)")
        _scatter_panel!(sax3, obs_co2_swc_x, obs_co2_swc_y;
            with_rmse = false, base_title = "Obs soil CO₂ vs obs SWC")
        sax4 = Axis(fig_scat[2, 2]; xlabel = "Model soil CO₂ (ppm)", ylabel = "Model SWC (m³/m³)")
        _scatter_panel!(sax4, mod_co2_swc_x, mod_co2_swc_y;
            with_rmse = false, base_title = "Model soil CO₂ vs model SWC")
        CairoMakie.save(joinpath(output_dir, "scatter_co2_swc_$(site_id).png"), fig_scat)
    end

    # mean / min / max over each daily-mean series (NaN when the series is empty).
    # These summarize the obs and model ranges for soil CO₂, SWC, and temperature.
    _dstat(df, f) = nrow(df) == 0 ? NaN : f(Float64.(df.daily_mean))

    # ── Full-period summary stats (INCLUDING spinup) ─────────────────────────
    # These intentionally use the full run window, not the post-spinup series
    # above. Obs come from the model-FORCING columns (TA_F air temp °C, VPD_F in
    # hPa, P_F accumulated precip mm); model air T is the `tair` driver (K).

    # Soil porosity ν — MEAN over the calibrated depths' layers (m³/m³).
    soil_porosity_layer = isempty(per_depth_porosity) ? NaN : mean(per_depth_porosity)

    # Model air temperature mean over the FULL period (K). Pulled directly from
    # the diagnostic writer to bypass get_diag_series' spinup filter; the
    # simulation spans start_date..stop_date, so this is already the full window.
    (_tair_times, _tair_data) = ClimaLand.Diagnostics.diagnostic_as_vectors(
        simulation.diagnostics[1].output_writer, "tair_30m_average"; layer = 1)
    model_tair_mean = isempty(_tair_data) ? NaN : mean(Float64.(_tair_data))

    # Obs bounded to THIS run's full period (obs_df is the whole multi-year CSV).
    obs_full = filter(r -> Date(start_date) <= r.date <= Date(stop_date), obs_df)
    function _col_vals(df, col)
        col in names(df) || return Float64[]
        out = Float64[]
        for v in df[!, col]
            (ismissing(v) || isnan(Float64(v))) && continue
            push!(out, Float64(v))
        end
        return out
    end
    _ta = _col_vals(obs_full, "TA_F")        # °C
    _vpd = _col_vals(obs_full, "VPD_F")      # hPa
    _precip = _col_vals(obs_full, "P_F")     # mm per timestep (accumulated)
    obs_tair_mean = isempty(_ta) ? NaN : mean(_ta) + 273.15   # → K to match model
    obs_vpd_mean = isempty(_vpd) ? NaN : mean(_vpd)
    obs_precip_sum = isempty(_precip) ? NaN : sum(_precip)

    # ── Pool the per-depth series into the master-CSV stats (all depths together) ──
    # RMSE/correlations: concatenate every depth's paired (obs, model) points and
    # compute ONE number over the pool — matches the depth-pooled EKI/Diagnostics
    # RMSE. Range stats (mean/min/max): concatenate each series' daily means over
    # all depths, then reduce. Porosity: mean over depths (computed above).
    _cat(f) = isempty(per_depth_stats) ? Float64[] :
              reduce(vcat, [Float64.(getfield(s, f).daily_mean) for s in per_depth_stats])
    _stat(v, f) = isempty(v) ? NaN : f(v)
    # paired points pooled over depths (obs_daily_cal ↔ model, obs_swc ↔ model swc)
    _pool_pairs(dfx_field, dfy_field) = isempty(per_depth_stats) ? (Float64[], Float64[]) :
        begin
            xs = Float64[]; ys = Float64[]
            for s in per_depth_stats
                px, py = _pair_on_date(getfield(s, dfx_field), getfield(s, dfy_field))
                append!(xs, px); append!(ys, py)
            end
            (xs, ys)
        end
    pooled_obs_mod_co2_x, pooled_obs_mod_co2_y = _pool_pairs(:obs_daily_cal, :sco2_daily)
    pooled_obs_mod_swc_x, pooled_obs_mod_swc_y = _pool_pairs(:obs_swc_daily, :swc_daily)
    pooled_obs_co2_swc_x, pooled_obs_co2_swc_y = _pool_pairs(:obs_daily_cal, :obs_swc_daily)
    pooled_mod_co2_swc_x, pooled_mod_co2_swc_y = _pool_pairs(:sco2_daily, :swc_daily)

    sco2_obs_all = _cat(:obs_daily_cal); sco2_mod_all = _cat(:sco2_daily)
    swc_obs_all = _cat(:obs_swc_daily); swc_mod_all = _cat(:swc_daily)
    tsoil_obs_all = _cat(:obs_tsoil_daily); tsoil_mod_all = _cat(:tsoil_daily)

    # Per-depth CO₂ RMSE (obs_daily_cal ↔ model), keyed by code. NaN for any code
    # not calibrated in this run (so its fixed CSV column reads `missing`).
    rmse_sco2_by_code = Dict{String, Float64}()
    for s in per_depth_stats
        x, y = _pair_on_date(s.obs_daily_cal, s.sco2_daily)
        rmse_sco2_by_code[s.code] = _rmse(x, y)
    end
    _rmse_code(c) = get(rmse_sco2_by_code, c, NaN)

    # Per-depth range stats: reduce one series (by field) for one depth code.
    # Returns NaN (→ missing in CSV) when the code isn't in this run or its series
    # is empty. Keyed by code so each depth's value is independent of ordering.
    _by_code = Dict(s.code => s for s in per_depth_stats)
    function _dstat_code(series_field, reducer, c)
        haskey(_by_code, c) || return NaN
        df = getfield(_by_code[c], series_field)
        nrow(df) == 0 ? NaN : reducer(Float64.(df.daily_mean))
    end

    # Stats handed back to the pipeline for the master CSV (POOLED over depths).
    scatter_stats = (;
        rmse_sco2 = _rmse(pooled_obs_mod_co2_x, pooled_obs_mod_co2_y),
        # per-depth CO₂ RMSE (NaN → missing for codes not in this run)
        rmse_sco2_501 = _rmse_code("501"),
        rmse_sco2_502 = _rmse_code("502"),
        rmse_sco2_503 = _rmse_code("503"),
        rmse_swc = _rmse(pooled_obs_mod_swc_x, pooled_obs_mod_swc_y),
        corr_obs_model_sco2 = _corr(pooled_obs_mod_co2_x, pooled_obs_mod_co2_y),
        corr_obs_model_swc = _corr(pooled_obs_mod_swc_x, pooled_obs_mod_swc_y),
        corr_obs_sco2_swc = _corr(pooled_obs_co2_swc_x, pooled_obs_co2_swc_y),
        corr_model_sco2_swc = _corr(pooled_mod_co2_swc_x, pooled_mod_co2_swc_y),
        # soil CO₂ (ppm) — obs stats on the masked (calibration) set, pooled depths
        obs_sco2_mean = _stat(sco2_obs_all, mean),
        obs_sco2_min = _stat(sco2_obs_all, minimum),
        obs_sco2_max = _stat(sco2_obs_all, maximum),
        model_sco2_mean = _stat(sco2_mod_all, mean),
        model_sco2_min = _stat(sco2_mod_all, minimum),
        model_sco2_max = _stat(sco2_mod_all, maximum),
        # soil water content (m³/m³)
        obs_swc_mean = _stat(swc_obs_all, mean),
        obs_swc_min = _stat(swc_obs_all, minimum),
        obs_swc_max = _stat(swc_obs_all, maximum),
        model_swc_mean = _stat(swc_mod_all, mean),
        model_swc_min = _stat(swc_mod_all, minimum),
        model_swc_max = _stat(swc_mod_all, maximum),
        # soil temperature (K)
        obs_tsoil_mean = _stat(tsoil_obs_all, mean),
        obs_tsoil_min = _stat(tsoil_obs_all, minimum),
        obs_tsoil_max = _stat(tsoil_obs_all, maximum),
        model_tsoil_mean = _stat(tsoil_mod_all, mean),
        model_tsoil_min = _stat(tsoil_mod_all, minimum),
        model_tsoil_max = _stat(tsoil_mod_all, maximum),
        # ── PER-DEPTH range stats (501/502/503); NaN→missing for absent codes ──
        # obs soil CO₂ (ppm) — calibration-masked series
        obs_sco2_mean_501 = _dstat_code(:obs_daily_cal, mean, "501"),
        obs_sco2_min_501  = _dstat_code(:obs_daily_cal, minimum, "501"),
        obs_sco2_max_501  = _dstat_code(:obs_daily_cal, maximum, "501"),
        obs_sco2_mean_502 = _dstat_code(:obs_daily_cal, mean, "502"),
        obs_sco2_min_502  = _dstat_code(:obs_daily_cal, minimum, "502"),
        obs_sco2_max_502  = _dstat_code(:obs_daily_cal, maximum, "502"),
        obs_sco2_mean_503 = _dstat_code(:obs_daily_cal, mean, "503"),
        obs_sco2_min_503  = _dstat_code(:obs_daily_cal, minimum, "503"),
        obs_sco2_max_503  = _dstat_code(:obs_daily_cal, maximum, "503"),
        # model soil CO₂ (ppm)
        model_sco2_mean_501 = _dstat_code(:sco2_daily, mean, "501"),
        model_sco2_min_501  = _dstat_code(:sco2_daily, minimum, "501"),
        model_sco2_max_501  = _dstat_code(:sco2_daily, maximum, "501"),
        model_sco2_mean_502 = _dstat_code(:sco2_daily, mean, "502"),
        model_sco2_min_502  = _dstat_code(:sco2_daily, minimum, "502"),
        model_sco2_max_502  = _dstat_code(:sco2_daily, maximum, "502"),
        model_sco2_mean_503 = _dstat_code(:sco2_daily, mean, "503"),
        model_sco2_min_503  = _dstat_code(:sco2_daily, minimum, "503"),
        model_sco2_max_503  = _dstat_code(:sco2_daily, maximum, "503"),
        # obs soil water content (m³/m³) — calibration-masked series
        obs_swc_mean_501 = _dstat_code(:obs_swc_daily, mean, "501"),
        obs_swc_min_501  = _dstat_code(:obs_swc_daily, minimum, "501"),
        obs_swc_max_501  = _dstat_code(:obs_swc_daily, maximum, "501"),
        obs_swc_mean_502 = _dstat_code(:obs_swc_daily, mean, "502"),
        obs_swc_min_502  = _dstat_code(:obs_swc_daily, minimum, "502"),
        obs_swc_max_502  = _dstat_code(:obs_swc_daily, maximum, "502"),
        obs_swc_mean_503 = _dstat_code(:obs_swc_daily, mean, "503"),
        obs_swc_min_503  = _dstat_code(:obs_swc_daily, minimum, "503"),
        obs_swc_max_503  = _dstat_code(:obs_swc_daily, maximum, "503"),
        # model soil water content (m³/m³)
        model_swc_mean_501 = _dstat_code(:swc_daily, mean, "501"),
        model_swc_min_501  = _dstat_code(:swc_daily, minimum, "501"),
        model_swc_max_501  = _dstat_code(:swc_daily, maximum, "501"),
        model_swc_mean_502 = _dstat_code(:swc_daily, mean, "502"),
        model_swc_min_502  = _dstat_code(:swc_daily, minimum, "502"),
        model_swc_max_502  = _dstat_code(:swc_daily, maximum, "502"),
        model_swc_mean_503 = _dstat_code(:swc_daily, mean, "503"),
        model_swc_min_503  = _dstat_code(:swc_daily, minimum, "503"),
        model_swc_max_503  = _dstat_code(:swc_daily, maximum, "503"),
        # obs soil temperature (K)
        obs_tsoil_mean_501 = _dstat_code(:obs_tsoil_daily, mean, "501"),
        obs_tsoil_min_501  = _dstat_code(:obs_tsoil_daily, minimum, "501"),
        obs_tsoil_max_501  = _dstat_code(:obs_tsoil_daily, maximum, "501"),
        obs_tsoil_mean_502 = _dstat_code(:obs_tsoil_daily, mean, "502"),
        obs_tsoil_min_502  = _dstat_code(:obs_tsoil_daily, minimum, "502"),
        obs_tsoil_max_502  = _dstat_code(:obs_tsoil_daily, maximum, "502"),
        obs_tsoil_mean_503 = _dstat_code(:obs_tsoil_daily, mean, "503"),
        obs_tsoil_min_503  = _dstat_code(:obs_tsoil_daily, minimum, "503"),
        obs_tsoil_max_503  = _dstat_code(:obs_tsoil_daily, maximum, "503"),
        # model soil temperature (K)
        model_tsoil_mean_501 = _dstat_code(:tsoil_daily, mean, "501"),
        model_tsoil_min_501  = _dstat_code(:tsoil_daily, minimum, "501"),
        model_tsoil_max_501  = _dstat_code(:tsoil_daily, maximum, "501"),
        model_tsoil_mean_502 = _dstat_code(:tsoil_daily, mean, "502"),
        model_tsoil_min_502  = _dstat_code(:tsoil_daily, minimum, "502"),
        model_tsoil_max_502  = _dstat_code(:tsoil_daily, maximum, "502"),
        model_tsoil_mean_503 = _dstat_code(:tsoil_daily, mean, "503"),
        model_tsoil_min_503  = _dstat_code(:tsoil_daily, minimum, "503"),
        model_tsoil_max_503  = _dstat_code(:tsoil_daily, maximum, "503"),
        # full-period (incl. spinup) summaries
        soil_porosity_layer = soil_porosity_layer,   # m³/m³, MEAN over cal depths
        obs_tair_mean = obs_tair_mean,               # K (TA_F + 273.15)
        model_tair_mean = model_tair_mean,           # K (tair driver)
        obs_precip_sum_mm = obs_precip_sum,          # mm, summed P_F over period
        obs_vpd_mean_hPa = obs_vpd_mean,             # hPa (VPD_F)
    )

    println("Forward run figures saved to $output_dir")
    return (; figures_dir = output_dir, scatter_stats)
end
