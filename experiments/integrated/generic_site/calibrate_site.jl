# Single-site calibration of ClimaLand parameters at a fluxnet site.
#
# Builds an `EKP.ObservationSeries` from up to `max_n_samples` 1-year windows
# of fluxnet data. Each EKP iteration draws `minibatch_size` samples; the
# forward model runs the corresponding window(s) preceded by a `spinup_days`
# spinup. Daily-mean diagnostics are compared to daily-mean observations, with
# per-variable noise variances given by the config.
#
# Usage:
#   SITE_ID=US-MOz julia --project=.buildkite \
#       experiments/integrated/generic_site/calibrate_site.jl
#
# Optional ENV vars:
#   SITE_ID            (default "US-MOz")
#   CALIBRATION_CONFIG (default "default.jl")
#   N_ITERS            override CALIBRATE_CONFIG.n_iterations
#   LOCAL              (default "true")  — use the per-site Val{} dispatchers
#                                          (the four bundled sites). Anything
#                                          else auto-resolves coordinates from
#                                          FLUXNET2015 metadata.

import ClimaComms
ClimaComms.@import_required_backends
using Dates
using Random
using Statistics
using LinearAlgebra
import JLD2
import ClimaCalibrate
import EnsembleKalmanProcesses as EKP
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.FluxnetSimulations as FluxnetSimulations
using ClimaLand.Simulations: solve!
import ClimaCore
import FileWatching.Pidfile

# ── Configuration ────────────────────────────────────────────────────────────

const CONFIG_FILE = get(ENV, "CALIBRATION_CONFIG", "default.jl")
include(joinpath(@__DIR__, "configs", CONFIG_FILE))
include(joinpath(@__DIR__, "site_drivers.jl"))

const SITE_ID = get(ENV, "SITE_ID", "US-MOz")
const N_ITERATIONS =
    parse(Int, get(ENV, "N_ITERS", string(CALIBRATE_CONFIG.n_iterations)))
const LOCAL_MODE = get(ENV, "LOCAL", "true") == "true"
const FT = Float64
const SHORT_NAMES = CALIBRATE_CONFIG.short_names
const SPINUP_DAYS = CALIBRATE_CONFIG.spinup_days
const CAL_WINDOW_DAYS = CALIBRATE_CONFIG.cal_window_days
const MAX_N_SAMPLES = CALIBRATE_CONFIG.max_n_samples
const MINIBATCH_SIZE = CALIBRATE_CONFIG.minibatch_size

const OUTPUT_DIR = joinpath(
    pkgdir(ClimaLand),
    "experiments",
    "integrated",
    "generic_site",
    "calibration_out",
    "$(SITE_ID)_$(CAL_WINDOW_DAYS)d",
)
mkpath(OUTPUT_DIR)

const SAMPLES_PATH = joinpath(OUTPUT_DIR, "samples.jld2")
const SITE_RESULTS_CSV = joinpath(
    pkgdir(ClimaLand),
    "experiments",
    "integrated",
    "generic_site",
    "calibration_out",
    "site_results.csv",
)

# ── Unit conversions ─────────────────────────────────────────────────────────
# Both observation values (after the per-variable preprocess in
# `get_comparison_data`) and model diagnostic values are converted to a common
# "obs vector" unit before being matched. Noise variances in the config are
# expressed in those units squared.

# `er` model diag is mol CO2 m⁻² s⁻¹; obs after preprocess is mol CO2 m⁻² s⁻¹.
# Both convert to g C m⁻² day⁻¹ via × 12.011 × 86400.
const _GC_PER_MOLCO2_PER_DAY_PER_S = 12.011 * 86400.0

function _obs_unit_factor(short_name)
    short_name == "er" && return _GC_PER_MOLCO2_PER_DAY_PER_S
    return 1.0
end

function _model_unit_factor(short_name)
    short_name == "er" && return _GC_PER_MOLCO2_PER_DAY_PER_S
    return 1.0
end

# ── Coordinate resolution ────────────────────────────────────────────────────

function _resolve_coords(site_ID)
    if LOCAL_MODE
        site_val = FluxnetSimulations.replace_hyphen(site_ID)
        (; time_offset, lat, long) =
            FluxnetSimulations.get_location(FT, Val(site_val))
        (; atmos_h) = FluxnetSimulations.get_fluxtower_height(FT, Val(site_val))
        return (; lat, long, time_offset, atmos_h)
    else
        info = FluxnetSimulations.get_site_info(site_ID)
        return (;
            lat = FT(info.lat),
            long = FT(info.long),
            time_offset = Int(info.time_offset),
            atmos_h = isempty(info.atmospheric_sensor_height) ?
                      error("Missing atmospheric_sensor_height for $site_ID") :
                      FT(first(info.atmospheric_sensor_height)),
        )
    end
end

# ── Sample windows ───────────────────────────────────────────────────────────

"""
    _sample_windows(time_offset)

Enumerate up to `MAX_N_SAMPLES` non-overlapping `CAL_WINDOW_DAYS`-long sample
windows starting at `data_start + spinup` (so the first window's spinup fits
inside the available data). Each entry is `(window_start, window_stop)` in
UTC. Also returns `data_start` so the forward model can compute
`start_offset` against the actual fluxnet record start.
"""
function _sample_windows(time_offset)
    (data_start, data_stop) = FluxnetSimulations.get_data_dates(
        SITE_ID,
        time_offset;
        required_columns = FluxnetSimulations.FLUXNET_FORCING_COLUMNS,
    )
    spinup = Day(SPINUP_DAYS)
    win = Day(CAL_WINDOW_DAYS)
    cursor = data_start + spinup
    samples = NTuple{2,DateTime}[]
    # Allow a 1-day slack on the upper bound — fluxnet records typically end
    # 30 min before the final calendar boundary (last `TIMESTAMP_END` is the
    # start of the last half-hour), so a strict `<=` would reject an exactly
    # 365-day record paired with a 365-day window.
    while cursor + win - Day(1) <= data_stop && length(samples) < MAX_N_SAMPLES
        push!(samples, (cursor, cursor + win))
        cursor += win
    end
    isempty(samples) && error(
        "$SITE_ID has no data window long enough for spinup=$(SPINUP_DAYS)d + window=$(CAL_WINDOW_DAYS)d.",
    )
    return samples, data_start
end

# ── Daily-mean aggregation ───────────────────────────────────────────────────

"""
    daily_means_in_window(comparison_data, varname, window_start, window_stop)

Compute means of `varname` from `comparison_data` over consecutive 24-hour
bins aligned with the simulation's `:daily` diagnostic schedule. Bin `k` covers
`[window_start + Day(k-1), window_start + Day(k))`, so `bin_ends[k] ==
window_start + Day(k)`. Returns `(values, bin_ends)`. Bins with no valid
observations are filled with NaN — the caller should reject the window.
"""
function daily_means_in_window(
    comparison_data,
    varname::Symbol,
    window_start::DateTime,
    window_stop::DateTime,
)
    times = comparison_data.UTC_datetime
    values = getproperty(comparison_data, varname)
    valid = .!ismissing.(values)
    times = times[valid]
    values = Float64.(values[valid])

    bin_ends = DateTime[]
    k = 1
    while window_start + Day(k) <= window_stop
        push!(bin_ends, window_start + Day(k))
        k += 1
    end

    daily = Float64[]
    for (i, bin_end) in enumerate(bin_ends)
        bin_start = i == 1 ? window_start : bin_ends[i - 1]
        mask = (bin_start .<= times) .& (times .< bin_end)
        if !any(mask)
            push!(daily, NaN)
        else
            push!(daily, mean(values[mask]))
        end
    end
    return daily, bin_ends
end

# ── Per-window observation builder ───────────────────────────────────────────

function _build_window_observation(comparison, window_start, window_stop, name)
    per_var = Dict{String,Vector{Float64}}()
    bins = DateTime[]
    for var in SHORT_NAMES
        haskey(NOISE_VARIANCES, var) ||
            error("NOISE_VARIANCES missing entry for short_name '$var'.")
        vals, bin_ends = daily_means_in_window(
            comparison,
            Symbol(var),
            window_start,
            window_stop,
        )
        if any(isnan, vals)
            n_missing = count(isnan, vals)
            error(
                "Missing $(n_missing) daily bins of '$var' at $SITE_ID in $window_start..$window_stop.",
            )
        end
        per_var[var] = vals .* _obs_unit_factor(var)
        isempty(bins) && (bins = bin_ends)
    end

    y_obs = reduce(vcat, per_var[v] for v in SHORT_NAMES)
    n_bins = length(first(values(per_var)))
    noise_diag = reduce(
        vcat,
        fill(NOISE_VARIANCES[v], n_bins) for v in SHORT_NAMES
    )
    return EKP.Observation(Dict(
        "samples" => y_obs,
        "covariances" => Diagonal(noise_diag),
        "names" => name,
    ))
end

"""
    generate_observation_series()

Build an `EKP.ObservationSeries` with one `EKP.Observation` per yearly sample
window. Persists `samples` (the list of `(window_start, window_stop)` tuples)
to `SAMPLES_PATH` so the forward model can dispatch by minibatch index.
"""
function generate_observation_series()
    coords = _resolve_coords(SITE_ID)
    (samples, data_start) = _sample_windows(coords.time_offset)

    comparison = FluxnetSimulations.get_comparison_data(SITE_ID, coords.time_offset)

    obs_vec = EKP.Observation[]
    for (i, (ws, we)) in enumerate(samples)
        push!(
            obs_vec,
            _build_window_observation(
                comparison,
                ws,
                we,
                "$(SITE_ID)_y$(i)_$(SHORT_NAMES |> x -> join(x, "_"))",
            ),
        )
    end

    minibatcher = ClimaCalibrate.minibatcher_over_samples(
        length(obs_vec),
        MINIBATCH_SIZE,
    )
    obs_series = EKP.ObservationSeries(Dict(
        "observations" => obs_vec,
        "minibatcher" => minibatcher,
        "names" => ["$(SITE_ID)_y$i" for i in 1:length(obs_vec)],
    ))

    JLD2.jldsave(
        SAMPLES_PATH;
        samples,
        data_start,
        time_offset = coords.time_offset,
        spinup_days = SPINUP_DAYS,
        cal_window_days = CAL_WINDOW_DAYS,
        short_names = SHORT_NAMES,
        site_ID = SITE_ID,
    )
    @info "Built observation series" SITE_ID n_samples = length(obs_vec) MINIBATCH_SIZE
    return obs_series
end

# ── Forward model ────────────────────────────────────────────────────────────

"""
    _extract_daily_scalar(diags_dict, short_name)

Pull a vector from a DictWriter produced by daily-average diagnostics on a
single-column simulation. The DictWriter stores keys like `"er_1D_average"`
for `short_name="er"`; this finds the matching key and returns the
time-sorted vector of scalar values.
"""
function _extract_daily_scalar(diags_dict, short_name::AbstractString)
    matching =
        [k for k in keys(diags_dict) if startswith(String(k), short_name * "_")]
    isempty(matching) && error(
        "No diagnostic key starting with '$short_name\\_' in DictWriter; available: $(collect(keys(diags_dict)))",
    )
    length(matching) > 1 &&
        @warn "Multiple keys match '$short_name'; using $(matching[1])"
    inner = diags_dict[matching[1]]
    times = sort(collect(keys(inner)))
    return [Float64(parent(inner[t])[1]) for t in times]
end

function _run_window(window_idx, samples, data_start, time_offset, toml_dict)
    (window_start, _) = samples[window_idx]
    spinup = Day(SPINUP_DAYS)
    duration = Day(SPINUP_DAYS + CAL_WINDOW_DAYS)
    sim_start = window_start - spinup
    start_offset = sim_start - data_start
    @assert start_offset >= Dates.Second(0) "Computed start_offset is negative for window $window_idx."

    coords = _resolve_coords(SITE_ID)
    result = FluxnetSimulations.generic_site_simulation(
        SITE_ID;
        coords...,
        toml_dict,
        duration,
        start_offset,
        diagnostic_period = :daily,
        output_vars = SHORT_NAMES,
    )
    solve!(result.simulation)

    writer = first(result.diags).output_writer
    g_window = Float64[]
    for var in SHORT_NAMES
        full = _extract_daily_scalar(writer.dict, var)
        # Drop the spinup prefix; the post-spinup tail must have CAL_WINDOW_DAYS bins.
        post = full[(SPINUP_DAYS + 1):end]
        length(post) == CAL_WINDOW_DAYS || error(
            "Window $window_idx '$var': expected $CAL_WINDOW_DAYS daily bins after spinup, got $(length(post)) (full length $(length(full))).",
        )
        append!(g_window, post .* _model_unit_factor(var))
    end
    return g_window
end

function ClimaCalibrate.forward_model(iteration, member)
    member_path =
        ClimaCalibrate.path_to_ensemble_member(OUTPUT_DIR, iteration, member)
    mkpath(member_path)

    params_path = ClimaCalibrate.parameter_path(OUTPUT_DIR, iteration, member)
    toml_dict = LP.create_toml_dict(FT, override_files = [params_path])

    ekp = JLD2.load_object(ClimaCalibrate.ekp_path(OUTPUT_DIR, iteration))
    minibatch = EKP.get_current_minibatch(ekp)

    state = JLD2.load(SAMPLES_PATH)
    samples = state["samples"]
    data_start = state["data_start"]
    time_offset = state["time_offset"]

    @info "[iter=$iteration mem=$member] running minibatch=$minibatch ($(length(minibatch)) windows)"

    g_member = Float64[]
    for window_idx in minibatch
        g_window = _run_window(
            window_idx,
            samples,
            data_start,
            time_offset,
            toml_dict,
        )
        append!(g_member, g_window)
    end

    JLD2.jldsave(joinpath(member_path, "g.jld2"); g = g_member)
    return nothing
end

# ── Observation map ──────────────────────────────────────────────────────────

function ClimaCalibrate.observation_map(iteration)
    ekp = JLD2.load_object(ClimaCalibrate.ekp_path(OUTPUT_DIR, iteration))
    n_ens = EKP.get_N_ens(ekp)

    first_g_path = joinpath(
        ClimaCalibrate.path_to_ensemble_member(OUTPUT_DIR, iteration, 1),
        "g.jld2",
    )
    isfile(first_g_path) ||
        error("No g.jld2 found for iteration $iteration; forward_model failed.")
    g_dim = length(JLD2.load_object(first_g_path))

    G = fill(NaN, g_dim, n_ens)
    for m in 1:n_ens
        g_path = joinpath(
            ClimaCalibrate.path_to_ensemble_member(OUTPUT_DIR, iteration, m),
            "g.jld2",
        )
        if !isfile(g_path)
            @warn "Missing g.jld2 for iteration $iteration, member $m"
            continue
        end
        G[:, m] = JLD2.load_object(g_path)
    end
    return G
end

# ── Per-site results CSV ─────────────────────────────────────────────────────

function append_site_csv(eki, prior, samples)
    coords = _resolve_coords(SITE_ID)
    drivers = compute_site_drivers(SITE_ID, coords, samples)
    final_mean = EKP.get_ϕ_mean_final(prior, eki)
    pnames = EKP.get_name(prior)
    # Per-parameter SD across the final ensemble in constrained (ϕ) space.
    # EKP doesn't expose a `get_ϕ_cov_final`, so we estimate the spread from
    # the final ensemble points directly. For TransformUnscented these are
    # the unscented sigma points; std across them is a usable posterior SD.
    final_ens = EKP.get_ϕ_final(prior, eki)
    final_sds = vec(std(final_ens; dims = 2))
    final_cost =
        try
            errs = EKP.get_error(eki)
            isempty(errs) ? NaN : last(errs)
        catch
            NaN
        end

    mkpath(dirname(SITE_RESULTS_CSV))

    Pidfile.mkpidlock(SITE_RESULTS_CSV * ".lock"; stale_age = 30) do
        write_header = !isfile(SITE_RESULTS_CSV)
        open(SITE_RESULTS_CSV, "a") do io
            if write_header
                cols = String[
                    "site_id",
                    "n_samples",
                    "n_iterations",
                    "minibatch_size",
                    "spinup_days",
                    "cal_window_days",
                ]
                append!(cols, pnames)
                append!(cols, [n * "_sd" for n in pnames])
                append!(
                    cols,
                    [
                        "mean_TA_K",
                        "mean_VPD_Pa",
                        "annual_cum_P_mm",
                        "mean_SW_in_W_m2",
                        "mean_TS_top_K",
                        "max_LAI",
                        "soil_porosity",
                        "soil_SOC_top_kgC_m3",
                        "IGBP",
                        "lat",
                        "long",
                        "final_eki_cost",
                    ],
                )
                println(io, join(cols, ","))
            end

            row = Any[
                SITE_ID,
                length(samples),
                N_ITERATIONS,
                MINIBATCH_SIZE,
                SPINUP_DAYS,
                CAL_WINDOW_DAYS,
            ]
            append!(row, final_mean)
            append!(row, final_sds)
            push!(row, drivers.mean_TA)
            push!(row, drivers.mean_VPD)
            push!(row, drivers.annual_cum_P)
            push!(row, drivers.mean_SW)
            push!(row, drivers.mean_TS_top)
            push!(row, drivers.max_LAI)
            push!(row, drivers.porosity)
            push!(row, drivers.SOC_top)
            push!(row, drivers.IGBP)
            push!(row, Float64(coords.lat))
            push!(row, Float64(coords.long))
            push!(row, final_cost)

            println(io, join(row, ","))
        end
    end
    @info "Appended site row to $SITE_RESULTS_CSV" SITE_ID
    return nothing
end

# ── Driver ───────────────────────────────────────────────────────────────────

if abspath(PROGRAM_FILE) == @__FILE__
    obs_path = joinpath(OUTPUT_DIR, "observation_series.jld2")
    obs_series = if isfile(obs_path) && isfile(SAMPLES_PATH)
        JLD2.jldopen(obs_path, "r") do f
            f["observation_series"]
        end
    else
        os = generate_observation_series()
        JLD2.jldsave(obs_path; observation_series = os)
        os
    end

    prior = get_calibration_prior()
    samples = JLD2.load(SAMPLES_PATH)["samples"]

    rng = Random.MersenneTwister(CALIBRATE_CONFIG.rng_seed)
    ekp = EKP.EnsembleKalmanProcess(
        obs_series,
        EKP.TransformUnscented(prior; impose_prior = true);
        verbose = true,
        rng,
    )

    @info "Starting calibration" SITE_ID CAL_WINDOW_DAYS SPINUP_DAYS N_ITERATIONS n_samples =
        length(samples) MINIBATCH_SIZE OUTPUT_DIR

    eki = ClimaCalibrate.calibrate(
        ClimaCalibrate.JuliaBackend(),
        ekp,
        N_ITERATIONS,
        prior,
        OUTPUT_DIR,
    )

    @info "Calibration complete." final_mean = EKP.get_ϕ_mean_final(prior, eki)
    append_site_csv(eki, prior, samples)
end
