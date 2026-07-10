"""
Observation_flag.jl — `generate_observations(run; base_dir)` builds a per-window
`EKP.ObservationSeries` (one `EKP.Observation` per calibration window, each stacked
over depth) + a minibatcher, and saves observations.jld2, WITH a
late-winter→spring "cap" mask applied to the soil-CO₂ observations.

MINIBATCH: instead of one fixed observation for one date range, this builds ONE
`EKP.Observation` per entry in `run.windows`, wraps them with
`ClimaCalibrate.minibatcher_over_samples(n_windows, run.minibatch_size)` into an
`EKP.ObservationSeries`, and saves the per-window `(start, stop)` list plus each
window's per-depth obs-date lists (the `samples` metadata) so the forward model
and observation_map can reconstruct the exact stacking order for whichever
minibatch is live in a given EKP iteration. Windows too thin at any calibrated
depth are SKIPPED (logged), not fatal.

This is a drop-in variant of Observations.jl. The ONLY difference (besides the
window/series structure) is that capped spring-peak days are dropped from the
soil-CO₂ observation vector before writing. Nothing else is touched: this file
only ever reads/writes soil-CO₂ (plus soil-T and air-T, used solely to *decide*
the mask). All other forcing in the CSV (radiation, precipitation, SWC, drivers,
…) is untouched — those are consumed elsewhere (ForwardRun.jl) straight from the
CSV.

Mask definition (the "final setup", matching the Python prototype
plot_neon_springpeak_allyears_uniondepth_batch.py):

  * Operate on DAILY data (plot-mean per half-hour → daily mean requiring ≥24
    valid half-hours), exactly the series that becomes `y_obs`.
  * A day is *capped* when soil CO₂ is ELEVATED (> CO2_FACTOR × baseline) AND the
    soil is COLD (daily soil-T < COLD_SOIL_K; where soil-T is missing, fall back
    to daily air-T < AIRT_FREEZE_K).
  * BASELINE is temperature-conditioned and computed ONCE over ALL YEARS: the
    25% quantile of daily soil CO₂ over all cold-soil days. (Not the fixed-month
    median — a low quantile leaves the cap peaks in the upper tail.)
  * DEPTH-CONSISTENT: the mask is computed per depth (501/502/503), and a date
    capped at ANY depth is dropped at the calibrated depth too (union over
    depths on dates).

All temperatures are converted to K here (the CSV stores soil-T and air-T in °C).

This file is `include`d into Main by the driver, replacing Observations.jl.
"""

using Dates
using Statistics
using LinearAlgebra
import JLD2
import EnsembleKalmanProcesses as EKP
import ClimaCalibrate
import ClimaLand
import ClimaLand.FluxnetSimulations as FluxnetSimulations
using CSV
using DataFrames

# Per-site depth → model-layer lookup (neon_depths_for_site) and its dependency
# _neon_site_key (site_metadata.jl). Included at TOP LEVEL and guarded so a second
# include (Calibration.jl also loads site_metadata.jl) is harmless. Both live in
# this pipeline's self-contained src/ folder.
const _SRC_DIR = normpath(joinpath(@__DIR__, "..", "src"))
if !isdefined(@__MODULE__, :_neon_site_key)
    include(joinpath(_SRC_DIR, "site_metadata.jl"))
end
if !isdefined(@__MODULE__, :neon_depths_for_site)
    include(joinpath(_SRC_DIR, "neon_depth_lookup.jl"))
end

# ── Cap-mask thresholds (mirror plot_neon_springpeak_allyears_uniondepth_batch.py;
#    these would eventually live in Config.jl, see Config.CapMaskConfig sketch). ──
const CAP_COLD_SOIL_K   = 275.15   # soil-T below this ⇒ frozen/thawing cap & baseline pool
const CAP_AIRT_FREEZE_K = 273.15   # air-T below this ⇒ snow/air cap (fallback only)
const CAP_CO2_FACTOR    = 1.3      # × cold-soil baseline to count as "elevated"
const CAP_BASELINE_Q    = 0.25     # quantile of daily CO₂ over cold-soil days
const CAP_MIN_COLD_DAYS = 5        # need this many cold-soil days, else no masking
const CAP_ALL_DEPTHS    = ("501", "502", "503")  # depths unioned for consistency

"""
    capped_dates(obs_df; depth) -> Set{Date}

Daily `is_capped_day` for one depth, returned as the set of capped calendar
dates. Builds the SAME plot-mean → daily(≥24) series as the observation vector,
then flags elevated-and-cold days. Temperatures converted °C → K.
Returns an empty set if the depth's CO₂ columns are absent or there are too few
cold-soil days to define a baseline.
"""
function capped_dates(obs_df, datetimes; depth::AbstractString)
    co2_cols = [Symbol("soilCO2concentrationMean_$(lpad(p, 3, '0'))_$depth") for p in 1:5]
    soilt_cols = [Symbol("soilTempMean_$(lpad(p, 3, '0'))_$depth") for p in 1:5]
    # depth must have CO₂ columns to be maskable; soil-T may be missing (→ air-T fallback)
    any(c -> hasproperty(obs_df, c), co2_cols) || return Set{Date}()

    present(cols) = [c for c in cols if hasproperty(obs_df, c)]
    co2_cols = present(co2_cols)
    soilt_cols = present(soilt_cols)

    function rowmean_skipinvalid(row, cols)
        vals = Float64[]
        for c in cols
            v = row[c]
            (!ismissing(v) && !isnan(Float64(v))) && push!(vals, Float64(v))
        end
        return isempty(vals) ? NaN : mean(vals)
    end

    df = DataFrame(
        date = Date.(datetimes),
        sco2 = [rowmean_skipinvalid(row, co2_cols) for row in eachrow(obs_df)],
        soilt = isempty(soilt_cols) ? fill(NaN, nrow(obs_df)) :
                [rowmean_skipinvalid(row, soilt_cols) for row in eachrow(obs_df)],
        airt = [(v = hasproperty(obs_df, :TA_F) ? row.TA_F : missing;
                 ismissing(v) ? NaN : Float64(v)) for row in eachrow(obs_df)],
    )

    # daily mean: CO₂ requires ≥24 valid half-hours (matches y_obs); T just averages
    daily = combine(
        groupby(df, :date),
        :sco2 => (x -> begin
            valid = filter(!isnan, x)
            length(valid) >= 24 ? mean(valid) : NaN
        end) => :sco2,
        :soilt => (x -> (v = filter(!isnan, x); isempty(v) ? NaN : mean(v))) => :soilt,
        :airt => (x -> (v = filter(!isnan, x); isempty(v) ? NaN : mean(v))) => :airt,
    )
    daily = filter(row -> !isnan(row.sco2), daily)
    isempty(daily) && return Set{Date}()

    soilt_k = daily.soilt .+ 273.15          # °C → K (NaN stays NaN)
    airt_k = daily.airt .+ 273.15
    has_soilt = .!isnan.(soilt_k)
    # cold-soil pool: prefer soil-T, fall back to air-T where soil-T is missing
    cold = [has_soilt[i] ? soilt_k[i] < CAP_COLD_SOIL_K :
                           airt_k[i] < CAP_AIRT_FREEZE_K
            for i in 1:nrow(daily)]

    cold_co2 = daily.sco2[cold]
    length(cold_co2) < CAP_MIN_COLD_DAYS && return Set{Date}()
    baseline = quantile(cold_co2, CAP_BASELINE_Q)

    elevated = daily.sco2 .> CAP_CO2_FACTOR * baseline
    capped = elevated .& cold                # NaN-safe: NaN comparisons are false
    return Set(daily.date[capped])
end

# Row-wise mean / variance across plot columns (skip missing/NaN). Module-level
# so both the per-depth daily builder and the cap mask can reuse them.
function _rowmean_skipinvalid(row, cols)
    vals = Float64[]
    for c in cols
        v = row[c]
        (!ismissing(v) && !isnan(Float64(v))) && push!(vals, Float64(v))
    end
    return isempty(vals) ? NaN : mean(vals)
end
function _rowvar_skipinvalid(row, cols)
    vals = Float64[]
    for c in cols
        v = row[c]
        (!ismissing(v) && !isnan(Float64(v))) && push!(vals, Float64(v))
    end
    return length(vals) >= 2 ? var(vals) : NaN
end

"""
    daily_sco2_for_depth(obs_df; depth, spinup_date, stop_date, capped) -> DataFrame

Plot-mean → daily(≥24 valid half-hours) soil-CO₂ mean + inter-sensor variance for
one `depth` code, trimmed to (spinup, stop], NaN-day-free, and with the shared
`capped` spring-peak dates dropped. Columns: :date, :daily_mean, :daily_var.
Assumes obs_df already carries :date (added once by generate_observations).
"""
function daily_sco2_for_depth(obs_df; depth, spinup_date, stop_date, capped)
    co2_cols = [Symbol("soilCO2concentrationMean_$(lpad(p, 3, '0'))_$depth") for p in 1:5]
    mean_col = [_rowmean_skipinvalid(row, co2_cols) for row in eachrow(obs_df)]
    var_col = [_rowvar_skipinvalid(row, co2_cols) for row in eachrow(obs_df)]
    df = DataFrame(date = obs_df.date, sco2_mean = mean_col, sco2_var = var_col)

    daily = combine(
        groupby(df, :date),
        :sco2_mean => (x -> begin
            valid = filter(!isnan, x)
            length(valid) >= 24 ? mean(valid) : NaN
        end) => :daily_mean,
        :sco2_var => (x -> begin
            valid = filter(!isnan, x)
            isempty(valid) ? NaN : mean(valid)
        end) => :daily_var,
    )
    daily = filter(row -> row.date >= Date(spinup_date), daily)
    daily = filter(row -> row.date <= Date(stop_date), daily)
    daily = filter(row -> !isnan(row.daily_mean), daily)
    daily = filter(row -> !(row.date in capped), daily)
    sort!(daily, :date)
    return daily
end

"""
    build_window_observation(obs_df, depth_lookup, w_start, w_stop, spinup_days,
                             capped, site_id)
        -> (observation, obs_dates_per_code, n_obs_per_code) | nothing

Build ONE stacked-over-depth `EKP.Observation` for the single window
`(w_start, w_stop)` (spinup trimmed off the front, cap-masked dates dropped).

Layout (per-depth concatenation, order = depth_lookup):
    samples     = [y_<code1>; y_<code2>; ...]   (each depth keeps its own valid days)
    covariances = Diagonal([noise_<code1>; noise_<code2>; ...])   (block-diagonal)

Returns `nothing` (→ caller SKIPS the window) if ANY calibrated depth has zero
valid days in the window, so the stacking order stays identical across all kept
windows (see the plan's thin-window decision). Otherwise returns the observation
plus the per-depth obs-date lists + counts the model side needs to realign G.
"""
function build_window_observation(
    obs_df, depth_lookup, w_start::Date, w_stop::Date, spinup_days::Int,
    capped, site_id,
)
    spinup_date = DateTime(w_start) + Day(spinup_days)
    stop_date = DateTime(w_stop)

    y_blocks = Vector{Float64}[]
    noise_blocks = Vector{Float64}[]
    obs_dates_per_code = Dict{String, Vector{Date}}()
    n_obs_per_code = Int[]

    for d in depth_lookup
        code = d.code
        daily = daily_sco2_for_depth(
            obs_df; depth = code, spinup_date = spinup_date,
            stop_date = stop_date, capped = capped,
        )
        n = nrow(daily)
        if n == 0
            @warn "Window $w_start..$w_stop: 0 valid days at depth $code for " *
                  "site $site_id — SKIPPING whole window."
            return nothing
        end
        mean_sensor_var = mean(filter(!isnan, daily.daily_var))
        push!(y_blocks, Float64.(daily.daily_mean))
        push!(noise_blocks, fill(mean_sensor_var, n))
        obs_dates_per_code[code] = daily.date
        push!(n_obs_per_code, n)
    end

    y_obs = reduce(vcat, y_blocks)
    noise_diag = reduce(vcat, noise_blocks)
    observation = EKP.Observation(Dict(
        "samples" => y_obs,
        "covariances" => Diagonal(noise_diag),
        "names" => "neon_sco2_$(site_id)_$(w_start)_$(w_stop)",
    ))
    return (observation, obs_dates_per_code, n_obs_per_code)
end

"""
    generate_observations(run; base_dir) -> obs_filepath

Build <base_dir>/observations.jld2 for `run` as an `EKP.ObservationSeries`: one
per-depth-stacked `EKP.Observation` per window in `run.windows`, wrapped with a
`minibatcher_over_samples(n_kept, run.minibatch_size)` minibatcher.

The cap mask is unchanged and computed ONCE over all years (union over
501/502/503, 25%-quantile cold-soil baseline), then applied within every window.

Saved metadata (the "samples" the model side needs to rebuild the layout for any
minibatch, keeping the function name so the driver is unchanged):
  - `observation_series`       the EKP.ObservationSeries (used by Calibration.jl)
  - `windows`                  ordered kept-window (start, stop) list; index i ↔
                               observation i in the series
  - `depth_codes`, `z_obs_m`   shared per-depth stacking order + measurement depths
  - `obs_dates_per_window`     Vector (per kept window) of Dict(code => dates)
  - `n_obs_per_window_code`    Vector (per kept window) of per-depth day counts
  - `minibatch_size`, `spinup_days`, `site_id`
The model LAYER is still derived at runtime by argmin against the live grid
(forward_model), so it is intentionally NOT stored here.
"""
function generate_observations(run; base_dir)
    FT = Float64
    site_id = run.site
    spinup_days = run.spinup_days
    codes = run.cal_depth_codes

    # per-site depth → model-layer lookup (same order/codes used by the model side)
    depth_lookup = neon_depths_for_site(site_id; codes = codes)

    println("Site: $site_id  (MINIBATCH; $(length(run.windows)) requested windows, ",
        "minibatch_size=$(run.minibatch_size))")
    println("Depth codes: ", join(codes, ", "),
        "  (measurement depths z_obs: ",
        join([string(d.z_obs_m) for d in depth_lookup], ", "), " m)")

    csv_path = ClimaLand.Artifacts.experiment_fluxnet_data_path(site_id)
    println("Loading NEON data from $csv_path")
    obs_df = CSV.read(csv_path, DataFrame)

    obs_df[!, :datetime] =
        DateTime.(string.(Int.(obs_df.timestamp_fmt)), dateformat"yyyymmddHHMM")
    obs_df[!, :date] = Date.(obs_df.datetime)

    # ── soil-CO₂ spring-peak cap mask (unchanged): depth-consistent capped-date set,
    #    union over 501/502/503, baseline over ALL years. Applied within EVERY
    #    window (all-years-baseline decision, see IMPLEMENTATION_PLAN.md). ──
    capped = reduce(
        union,
        (capped_dates(obs_df, obs_df.datetime; depth = d) for d in CAP_ALL_DEPTHS);
        init = Set{Date}(),
    )
    println("Cap mask: $(length(capped)) capped days (union over depths $(CAP_ALL_DEPTHS), all years)")

    # ── Build one Observation per window; skip windows too thin at any depth ──
    obs_vec = EKP.Observation[]
    kept_windows = Tuple{Date, Date}[]
    obs_dates_per_window = Dict{String, Vector{Date}}[]
    n_obs_per_window_code = Vector{Int}[]

    for (w_start, w_stop) in run.windows
        built = build_window_observation(
            obs_df, depth_lookup, w_start, w_stop, spinup_days, capped, site_id,
        )
        built === nothing && continue
        (observation, obs_dates_per_code, n_obs_per_code) = built
        push!(obs_vec, observation)
        push!(kept_windows, (w_start, w_stop))
        push!(obs_dates_per_window, obs_dates_per_code)
        push!(n_obs_per_window_code, n_obs_per_code)
        println("  window $w_start..$w_stop: per-depth days ",
            join(n_obs_per_code, "+"), " (codes ", join(codes, ","), ")")
    end

    n_kept = length(obs_vec)
    n_kept == 0 && error("No usable windows for site $site_id (all skipped)")
    n_kept < length(run.windows) &&
        @warn "Kept $n_kept of $(length(run.windows)) windows for $site_id " *
              "(the rest were too thin)."

    mb_size = min(run.minibatch_size, n_kept)
    mb_size < run.minibatch_size &&
        @warn "minibatch_size ($(run.minibatch_size)) > kept windows ($n_kept); " *
              "using $mb_size."
    minibatcher = ClimaCalibrate.minibatcher_over_samples(n_kept, mb_size)
    window_names =
        ["$(site_id)_$(s)_$(e)" for (s, e) in kept_windows]
    observation_series = EKP.ObservationSeries(obs_vec, minibatcher, window_names)
    println("Built ObservationSeries: $n_kept windows, minibatch_size=$mb_size")

    mkpath(base_dir)
    obs_filepath = joinpath(base_dir, "observations.jld2")
    JLD2.jldsave(
        obs_filepath;
        observation_series = observation_series,
        windows = kept_windows,
        depth_codes = [d.code for d in depth_lookup],
        z_obs_m = [d.z_obs_m for d in depth_lookup],
        obs_dates_per_window = obs_dates_per_window,
        n_obs_per_window_code = n_obs_per_window_code,
        minibatch_size = mb_size,
        site_id = site_id,
        spinup_days = spinup_days,
    )
    println("Saved observation series to $obs_filepath")
    return obs_filepath
end
