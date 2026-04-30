# Single-site calibration of PModel parameters at a fluxnet site.
#
# Mirrors the global calibration template (experiments/calibration/configs/gpp_lhf.jl
# + experiments/calibration/run_calibration.jl) but stripped of minibatching,
# multi-backend dispatch, and HPC kwargs. Runs locally with ClimaCalibrate.JuliaBackend.
#
# Usage:
#   SITE_ID=US-MOz julia --project=.buildkite \
#       experiments/integrated/generic_site/calibrate_site.jl
#
# Optional ENV vars:
#   SITE_ID            (default "US-MOz")  — any fluxnet site_ID resolvable via
#                                              ClimaLand.Artifacts.experiment_fluxnet_data_path
#   CALIBRATION_CONFIG (default "default.jl")
#                                          — config file in configs/ that defines
#                                              CALIBRATE_CONFIG, NOISE_VARIANCES,
#                                              and get_calibration_prior(). Use
#                                              "smoke_test.jl" for CI smoke runs.
#   CAL_DURATION_DAYS  (override config)   — length of the calibration window
#                                              starting at the site's data start date
#   N_ITERS            (override config)   — number of EKP iterations
#   LOCAL              (default "true")   — if "true", look up coordinates from the
#                                              per-site Val{} dispatchers (the four
#                                              bundled sites). Set to anything else
#                                              to auto-resolve via FLUXNET2015 metadata.

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

# ── Configuration ────────────────────────────────────────────────────────────

const CONFIG_FILE = get(ENV, "CALIBRATION_CONFIG", "default.jl")
include(joinpath(@__DIR__, "configs", CONFIG_FILE))

const SITE_ID = get(ENV, "SITE_ID", "US-MOz")
const CAL_DURATION_DAYS =
    parse(Int, get(ENV, "CAL_DURATION_DAYS", string(CALIBRATE_CONFIG.cal_duration_days)))
const N_ITERATIONS =
    parse(Int, get(ENV, "N_ITERS", string(CALIBRATE_CONFIG.n_iterations)))
const LOCAL_MODE = get(ENV, "LOCAL", "true") == "true"
const FT = Float64
const SHORT_NAMES = CALIBRATE_CONFIG.short_names

const OUTPUT_DIR = joinpath(
    pkgdir(ClimaLand),
    "experiments",
    "integrated",
    "generic_site",
    "calibration_out",
    "$(SITE_ID)_$(CAL_DURATION_DAYS)d",
)
mkpath(OUTPUT_DIR)

# ── Coordinate resolution ────────────────────────────────────────────────────
# In LOCAL_MODE, resolve from the per-site Val{} dispatchers (works for the four
# bundled sites without the fluxnet2015 artifact). Otherwise let downstream
# code auto-resolve via FLUXNET2015 metadata.

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

# ── Observation generation ───────────────────────────────────────────────────

"""
    monthly_means_in_window(comparison_data, varname, window_start, window_stop)

Compute means of `varname` from `comparison_data` over month-length bins
aligned with the simulation's `:monthly` diagnostic schedule, which uses
`EveryCalendarDtSchedule(Month(1); start_date=window_start)` and therefore
emits at `window_start + Month(k)` for `k = 1, 2, …`. Each bin `k` covers
`[window_start + Month(k-1), window_start + Month(k))`, and only bins whose
end is `<= window_stop` are returned (these are exactly the bins the
diagnostic emits before the simulation stops). Returns `(values, bin_ends)`
where `bin_ends[k] == window_start + Month(k)`.

Errors if any included bin has no valid observations.
"""
function monthly_means_in_window(
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
    while window_start + Month(k) <= window_stop
        push!(bin_ends, window_start + Month(k))
        k += 1
    end

    monthly = Float64[]
    for (i, bin_end) in enumerate(bin_ends)
        bin_start = i == 1 ? window_start : bin_ends[i - 1]
        mask = (bin_start .<= times) .& (times .< bin_end)
        if !any(mask)
            push!(monthly, NaN)
        else
            push!(monthly, mean(values[mask]))
        end
    end
    if any(isnan, monthly)
        bad = findall(isnan, monthly)
        error(
            "Missing months $(bad) of $(varname) at $(SITE_ID) in window $(window_start) — $(window_stop).",
        )
    end
    return monthly, bin_ends
end

"""
    _calibration_window(time_offset)

Resolve the calibration window: starts at the site's data start date and runs
for `CAL_DURATION_DAYS`. Returns `(start_date, stop_date)` as `DateTime`s.

Passes `required_columns = FLUXNET_FORCING_COLUMNS` so the start matches what
`generic_site_simulation` will pick internally — otherwise the observation
window and the simulation window can drift apart at sites with leading
missing forcing data, leading to a g/y length mismatch.
"""
function _calibration_window(time_offset)
    (start_date, _) = FluxnetSimulations.get_data_dates(
        SITE_ID,
        time_offset;
        required_columns = FluxnetSimulations.FLUXNET_FORCING_COLUMNS,
    )
    stop_date = start_date + Day(CAL_DURATION_DAYS)
    return (start_date, stop_date)
end

"""
    generate_observations()

Build an EKP.Observation from the fluxnet site's reported GPP and LHF, reduced
to monthly means over the calibration window. Saves to
`$OUTPUT_DIR/observation.jld2`.
"""
function generate_observations()
    coords = _resolve_coords(SITE_ID)
    (window_start, window_stop) = _calibration_window(coords.time_offset)
    comparison =
        FluxnetSimulations.get_comparison_data(SITE_ID, coords.time_offset)

    monthly_per_var = Dict{String,Vector{Float64}}()
    months = DateTime[]
    for name in SHORT_NAMES
        haskey(NOISE_VARIANCES, name) || error(
            "NOISE_VARIANCES is missing an entry for short_name '$name' (defined in $(CONFIG_FILE)).",
        )
        vals, bin_ends = monthly_means_in_window(
            comparison,
            Symbol(name),
            window_start,
            window_stop,
        )
        monthly_per_var[name] = vals
        isempty(months) && (months = bin_ends)
    end

    # Comparison gpp is in mol/m²/s after the preprocess_func; lhf in W/m².
    y_obs = reduce(vcat, monthly_per_var[name] for name in SHORT_NAMES)

    n_months = length(first(values(monthly_per_var)))
    noise_diag = reduce(
        vcat,
        fill(NOISE_VARIANCES[name], n_months) for name in SHORT_NAMES
    )

    obs = EKP.Observation(
        Dict(
            "samples" => y_obs,
            "covariances" => Diagonal(noise_diag),
            "names" => "$(SITE_ID)_$(CAL_DURATION_DAYS)d_monthly_$(join(SHORT_NAMES, "_"))",
        ),
    )

    obs_path = joinpath(OUTPUT_DIR, "observation.jld2")
    JLD2.jldsave(
        obs_path;
        observation = obs,
        y_obs,
        noise_cov = Diagonal(noise_diag),
        months,
        window_start,
        window_stop,
        site_ID = SITE_ID,
        cal_duration_days = CAL_DURATION_DAYS,
    )
    @info "Saved observation to $obs_path" length_y = length(y_obs) n_months
    return obs
end

# ── Forward model ────────────────────────────────────────────────────────────

"""
    _extract_monthly_scalar(diags_dict, short_name)

Pull a vector from a DictWriter produced by monthly-average diagnostics on a
single-column simulation. The DictWriter stores keys like
`"gpp_1M_average"` for `short_name="gpp"`; this helper finds the matching key
and returns the time-sorted vector of scalar values.
"""
function _extract_monthly_scalar(diags_dict, short_name::AbstractString)
    matching =
        [k for k in keys(diags_dict) if startswith(String(k), short_name * "_")]
    isempty(matching) && error(
        "No diagnostic key starting with '$(short_name)_' in DictWriter; available: $(collect(keys(diags_dict)))",
    )
    length(matching) > 1 &&
        @warn "Multiple keys match '$short_name'; using $(matching[1])"
    inner = diags_dict[matching[1]]
    times = sort(collect(keys(inner)))
    return [Float64(parent(inner[t])[1]) for t in times]
end

function ClimaCalibrate.forward_model(iteration, member)
    member_path =
        ClimaCalibrate.path_to_ensemble_member(OUTPUT_DIR, iteration, member)
    mkpath(member_path)

    params_path = ClimaCalibrate.parameter_path(OUTPUT_DIR, iteration, member)
    toml_dict = LP.create_toml_dict(FT, override_files = [params_path])

    coords = _resolve_coords(SITE_ID)

    @info "[$iteration/$member] Running $SITE_ID for $CAL_DURATION_DAYS days"

    result = FluxnetSimulations.generic_site_simulation(
        SITE_ID;
        coords...,
        toml_dict,
        duration = Day(CAL_DURATION_DAYS),
        diagnostic_period = :monthly,
        output_vars = SHORT_NAMES,
    )
    solve!(result.simulation)

    # Pull g vector from the in-memory DictWriter, in the same variable order
    # as `SHORT_NAMES` so it matches `y_obs` in `generate_observations`.
    writer = first(result.diags).output_writer
    g_member = reduce(
        vcat,
        _extract_monthly_scalar(writer.dict, name) for name in SHORT_NAMES
    )

    JLD2.jldsave(joinpath(member_path, "g.jld2"); g = g_member)
    return nothing
end

# ── Observation map ──────────────────────────────────────────────────────────

function ClimaCalibrate.observation_map(iteration)
    ekp = JLD2.load_object(ClimaCalibrate.ekp_path(OUTPUT_DIR, iteration))
    n_ens = EKP.get_N_ens(ekp)

    # Probe the first available g.jld2 to learn the dimension
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

# ── Driver ───────────────────────────────────────────────────────────────────

if abspath(PROGRAM_FILE) == @__FILE__
    obs_path = joinpath(OUTPUT_DIR, "observation.jld2")
    obs = if isfile(obs_path)
        JLD2.jldopen(obs_path, "r") do f
            f["observation"]
        end
    else
        generate_observations()
    end

    prior = get_calibration_prior()

    rng = Random.MersenneTwister(CALIBRATE_CONFIG.rng_seed)
    ekp = EKP.EnsembleKalmanProcess(
        obs,
        EKP.TransformUnscented(prior; impose_prior = true);
        verbose = true,
        rng,
    )

    @info "Starting calibration" SITE_ID CAL_DURATION_DAYS N_ITERATIONS OUTPUT_DIR
    eki = ClimaCalibrate.calibrate(
        ClimaCalibrate.JuliaBackend(),
        ekp,
        N_ITERATIONS,
        prior,
        OUTPUT_DIR,
    )

    @info "Calibration complete." final_mean = EKP.get_ϕ_mean_final(prior, eki)
end
