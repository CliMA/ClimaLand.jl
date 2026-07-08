"""
Observation_flag.jl — `generate_observations(run; base_dir)` builds the daily
soil-CO₂ observation vector + noise covariance and saves observations.jld2,
WITH a late-winter→spring "cap" mask applied to the soil-CO₂ observations.

This is a drop-in variant of Observations.jl. The ONLY difference is that capped
spring-peak days are dropped from the soil-CO₂ observation vector (`y_obs`,
`obs_dates`, noise) before writing. Nothing else is touched: this file only ever
reads/writes soil-CO₂ (plus soil-T and air-T, used solely to *decide* the mask).
All other forcing in the CSV (radiation, precipitation, SWC, drivers, …) is
untouched — those are consumed elsewhere (ForwardRun.jl) straight from the CSV.

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
    generate_observations(run; base_dir) -> obs_filepath

Build <base_dir>/observations.jld2 for `run`, stacking the daily soil-CO₂
observation over ALL of `run.cal_depth_codes` (per-depth concatenation), with the
soil-CO₂ spring-peak cap mask applied.

Layout (option (a), per-depth concatenation, order = run.cal_depth_codes):
    y_obs     = [y_<code1>; y_<code2>; ...]   (each depth keeps its own valid days)
    noise_cov = Diagonal([noise_<code1>; noise_<code2>; ...])   (block-diagonal)
The per-depth date lists + model layers are saved so observation_map (in the model
interface) rebuilds the same stacking order. Depths→layers come from the shared
lookup (src/neon_depth_lookup.jl → neon_depth_layer_mapping.csv).

The cap mask is unchanged: capped dates (union over 501/502/503, all years,
25%-quantile cold-soil baseline) are dropped at every depth.
"""
function generate_observations(run; base_dir)
    FT = Float64
    site_id = run.site
    spinup_days = run.spinup_days
    codes = run.cal_depth_codes

    # per-site depth → model-layer lookup (same order/codes used by the model side)
    depth_lookup = neon_depths_for_site(site_id; codes = codes)

    start_date = DateTime(run.start_date)
    stop_date = DateTime(run.stop_date)
    spinup_date = start_date + Day(spinup_days)

    println("Site: $site_id")
    println("Data period: $start_date to $stop_date  (spinup until $spinup_date)")
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
    #    union over 501/502/503, baseline over all years. Applied at EVERY depth. ──
    capped = reduce(
        union,
        (capped_dates(obs_df, obs_df.datetime; depth = d) for d in CAP_ALL_DEPTHS);
        init = Set{Date}(),
    )
    println("Cap mask: $(length(capped)) capped days (union over depths $(CAP_ALL_DEPTHS), all years)")

    # ── Build the stacked (per-depth concatenated) observation ───────────────
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
        n == 0 && error("No valid observation days for site $site_id depth code $code")

        mean_sensor_var = mean(filter(!isnan, daily.daily_var))
        push!(y_blocks, Float64.(daily.daily_mean))
        push!(noise_blocks, fill(mean_sensor_var, n))
        obs_dates_per_code[code] = daily.date
        push!(n_obs_per_code, n)
        println("  code $code (z=$(d.z_obs_m) m): ",
            "$n days, noise=$(round(mean_sensor_var, sigdigits = 3)) ppm²")
    end

    y_obs = reduce(vcat, y_blocks)
    noise_diag = reduce(vcat, noise_blocks)
    n_obs = length(y_obs)
    # flat date list in stacking order (back-compat + convenience)
    obs_dates = reduce(vcat, [obs_dates_per_code[d.code] for d in depth_lookup])
    println("Stacked observation vector length: $n_obs ",
        "(per-depth: ", join(n_obs_per_code, "+"), ")")

    observation = EKP.Observation(Dict(
        "samples" => y_obs,
        "covariances" => Diagonal(noise_diag),
        "names" => "neon_sco2_$(site_id)",
    ))

    mkpath(base_dir)
    obs_filepath = joinpath(base_dir, "observations.jld2")
    JLD2.jldsave(
        obs_filepath;
        observation = observation,
        y_obs = y_obs,
        noise_cov = Diagonal(noise_diag),
        obs_dates = obs_dates,
        # multi-depth metadata — obs/model stacking order + per-depth measurement
        # depths. The MODEL LAYER is derived at runtime by argmin against the live
        # grid (forward_model), so it is intentionally NOT stored here.
        depth_codes = [d.code for d in depth_lookup],
        z_obs_m = [d.z_obs_m for d in depth_lookup],
        obs_dates_per_code = obs_dates_per_code,
        n_obs_per_code = n_obs_per_code,
        site_id = site_id,
        spinup_days = spinup_days,
    )
    println("Saved observations to $obs_filepath")
    return obs_filepath
end
