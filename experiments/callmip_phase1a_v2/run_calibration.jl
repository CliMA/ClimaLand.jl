"""
EKI calibration driver for CalLMIP Phase 1a DK-Sor.

Uses minibatched ObservationSeries with RandomFixedSizeMinibatcher (batch=2).
16 calibration parameters via TransformUnscented (UKI) ensemble.

Usage:
    julia --project=experiments/callmip_phase1a_v2 \\
          experiments/callmip_phase1a_v2/run_calibration.jl
"""

using Dates
using Distributed
import Random
import ClimaCalibrate
import ClimaLand
import EnsembleKalmanProcesses as EKP
import EnsembleKalmanProcesses.ParameterDistributions as PD
import JLD2
using LinearAlgebra

# ── Configuration ────────────────────────────────────────────────────────────
const SITE_ID       = "DK-Sor"
const N_ITERATIONS  = 10
const MINIBATCH_SIZE = 2
const DT            = Float64(900)

const climaland_dir = abspath(joinpath(@__DIR__, "..", ".."))
const OUTPUT_DIR    = joinpath(@__DIR__, "output_eki")
const OBS_FILEPATH  = joinpath(@__DIR__, "observations.jld2")

isdir(OUTPUT_DIR) || mkpath(OUTPUT_DIR)

# ── Priors ────────────────────────────────────────────────────────────────────
include(joinpath(@__DIR__, "priors.jl"))
prior, priors_vec = build_dk_sor_priors()
n_params = length(priors_vec)
println("Prior built for $n_params parameters")

# ── Load observations ─────────────────────────────────────────────────────────
isfile(OBS_FILEPATH) ||
    error("Observations file not found: $OBS_FILEPATH\nRun generate_observations.jl first.")

obs_data  = JLD2.load(OBS_FILEPATH)
observation_vector = obs_data["observation_vector"]
window_years       = obs_data["window_years"]
window_dates       = obs_data["window_dates"]
n_windows          = length(observation_vector)

println("Loaded $n_windows observation windows ($(first(window_years))–$(last(window_years)))")
println("Minibatch size: $MINIBATCH_SIZE → $(n_windows ÷ MINIBATCH_SIZE) batches/iteration")

obs_series = EKP.ObservationSeries(
    Dict(
        "observations" => observation_vector,
        "names"        => string.(window_years),
        "minibatcher"  => ClimaCalibrate.minibatcher_over_samples(
            n_windows, MINIBATCH_SIZE,
        ),
    ),
)

# ── Create EKP ────────────────────────────────────────────────────────────────
rng = Random.MersenneTwister(1234)
ekp = EKP.EnsembleKalmanProcess(
    obs_series,
    EKP.TransformUnscented(prior; impose_prior = true),
    verbose  = true,
    rng      = rng,
    scheduler = EKP.DataMisfitController(terminate_at = 100),
)

N_ens = EKP.get_N_ens(ekp)
println("Ensemble size: $N_ens (= 2 × $n_params + 1 for UKI)")

# ── Add Slurm workers ─────────────────────────────────────────────────────────
curr_backend = ClimaCalibrate.get_backend()
if curr_backend == ClimaCalibrate.CaltechHPCBackend ||
   curr_backend == ClimaCalibrate.GCPBackend
    @info "Slurm detected — launching $N_ens worker jobs."
    addprocs(ClimaCalibrate.SlurmManager())
elseif nworkers() == 1
    @info "No Slurm / only 1 worker — running serially (JuliaBackend fallback)."
end

# ── Broadcast to workers ──────────────────────────────────────────────────────
@everywhere using Distributed
@everywhere import ClimaLand
@everywhere import ClimaLand.Parameters as LP
@everywhere const SITE_ID       = $SITE_ID
@everywhere const OUTPUT_DIR    = $OUTPUT_DIR
@everywhere const OBS_FILEPATH  = $OBS_FILEPATH
@everywhere const DT            = $DT
@everywhere include(joinpath(abspath(joinpath(@__DIR__, "..", "..")),
                              "experiments", "callmip_phase1a_v2", "model_interface.jl"))

# ── Run calibration ───────────────────────────────────────────────────────────
println("\n=== Starting EKI calibration ===")
println("  Backend     : WorkerBackend")
println("  Iterations  : $N_ITERATIONS")
println("  Ensemble sz : $N_ens")
println("  Windows     : $n_windows ($(first(window_years))–$(last(window_years)))")
println("  Minibatch   : $MINIBATCH_SIZE years/iter")
println("  Parameters  : $n_params")
println("  Output dir  : $OUTPUT_DIR")

eki = ClimaCalibrate.calibrate(
    ClimaCalibrate.WorkerBackend,
    ekp,
    N_ITERATIONS,
    prior,
    OUTPUT_DIR,
)

# ── Print final results ───────────────────────────────────────────────────────
final_params = EKP.get_ϕ_mean_final(prior, eki)
param_names  = [only(PD.get_name(d)) for d in priors_vec]

println("\n=== Calibration complete ===")
println("Final parameter means:")
for (name, val) in zip(param_names, final_params)
    println("  $(rpad(name, 40))  $(round(val; sigdigits = 5))")
end

# Save final EKI-optimal params as TOML for use by callmip simulations
eki_optimal_toml = joinpath(OUTPUT_DIR, "eki_optimal_parameters.toml")
open(eki_optimal_toml, "w") do io
    for (name, val) in zip(param_names, final_params)
        used_in = name in ("soilCO2_reference_rate", "michaelis_constant",
                           "O2_michaelis_constant", "soilCO2_activation_energy") ?
                  "[\"Land\"]" : "[\"getindex\"]"
        println(io, "[\"$name\"]")
        println(io, "value = $(Float64(val))")
        println(io, "type  = \"float\"")
        println(io, "used_in = $used_in")
        println(io)
    end
end
println("EKI-optimal parameters saved → $eki_optimal_toml")
