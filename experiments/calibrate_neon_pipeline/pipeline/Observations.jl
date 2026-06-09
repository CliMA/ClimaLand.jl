"""
Observations.jl — `generate_observations(run; base_dir)` builds the daily
soil-CO₂ observation vector + noise covariance and saves observations.jld2.

Function form of the original generate_observations.jl (frozen in
../calibrate_neon). All configuration is passed via the `RunConfig` argument —
no ENV, no top-level globals. Returns the path to the written observations.jld2.

This file is `include`d into Main by the driver.
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

"""
    generate_observations(run; base_dir) -> obs_filepath

Build <base_dir>/observations.jld2 for `run`. Reads NEON soil-CO₂ at the depth
implied by `run.cal_depth` (501 for 0.02 m, 502 for 0.06 m), computes daily means
across plots 001–005 (requiring ≥24 valid half-hours/day), trims to after spinup
and before stop, and estimates noise from inter-sensor variance.
"""
function generate_observations(run; base_dir)
    FT = Float64
    site_id = run.site
    spinup_days = run.spinup_days
    cal_depth_str = string(run.cal_depth)

    obs_depth = Config.obs_depth_code(run.cal_depth)  # validates depth

    start_date = DateTime(run.start_date)
    stop_date = DateTime(run.stop_date)
    spinup_date = start_date + Day(spinup_days)

    println("Site: $site_id")
    println("Data period: $start_date to $stop_date  (spinup until $spinup_date)")

    csv_path = ClimaLand.Artifacts.experiment_fluxnet_data_path(site_id)
    println("Loading NEON data from $csv_path")
    obs_df = CSV.read(csv_path, DataFrame)

    co2_cols = [
        Symbol("soilCO2concentrationMean_$(lpad(p, 3, '0'))_$obs_depth")
        for p in 1:5
    ]

    function rowmean_skipinvalid(row, cols)
        vals = Float64[]
        for c in cols
            v = row[c]
            (!ismissing(v) && !isnan(Float64(v))) && push!(vals, Float64(v))
        end
        return isempty(vals) ? NaN : mean(vals)
    end
    function rowvar_skipinvalid(row, cols)
        vals = Float64[]
        for c in cols
            v = row[c]
            (!ismissing(v) && !isnan(Float64(v))) && push!(vals, Float64(v))
        end
        return length(vals) >= 2 ? var(vals) : NaN
    end

    obs_df[!, :sco2_mean] =
        [rowmean_skipinvalid(row, co2_cols) for row in eachrow(obs_df)]
    obs_df[!, :sco2_var] =
        [rowvar_skipinvalid(row, co2_cols) for row in eachrow(obs_df)]
    obs_df[!, :datetime] =
        DateTime.(string.(Int.(obs_df.timestamp_fmt)), dateformat"yyyymmddHHMM")
    obs_df[!, :date] = Date.(obs_df.datetime)

    daily_df = combine(
        groupby(obs_df, :date),
        :sco2_mean => (x -> begin
            valid = filter(!isnan, x)
            length(valid) >= 24 ? mean(valid) : NaN
        end) => :daily_mean,
        :sco2_var => (x -> begin
            valid = filter(!isnan, x)
            isempty(valid) ? NaN : mean(valid)
        end) => :daily_var,
    )

    daily_df = filter(row -> row.date >= Date(spinup_date), daily_df)
    daily_df = filter(row -> row.date <= Date(stop_date), daily_df)
    daily_df = filter(row -> !isnan(row.daily_mean), daily_df)
    sort!(daily_df, :date)

    y_obs = Float64.(daily_df.daily_mean)
    obs_dates = daily_df.date
    n_obs = length(y_obs)
    println("Valid observation days: $n_obs (after $spinup_days-day spinup)")

    mean_sensor_var = mean(filter(!isnan, daily_df.daily_var))
    noise_diag = fill(mean_sensor_var, n_obs)
    println("Noise variance (inter-sensor): $(round(mean_sensor_var, sigdigits=3)) ppm²")

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
        obs_dates = obs_dates,
        y_obs = y_obs,
        noise_cov = Diagonal(noise_diag),
        site_id = site_id,
        spinup_days = spinup_days,
    )
    println("Saved observations to $obs_filepath")
    return obs_filepath
end
