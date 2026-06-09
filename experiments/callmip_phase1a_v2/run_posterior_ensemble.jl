"""
Run 50-member posterior ensemble for CalLMIP Phase 1a DK-Sor.

Draws 50 parameter sets from the MCMC posterior (emulate_sample.jl output),
runs full 1997-2014 simulations, and saves diagnostics for uncertainty analysis.

Usage:
    julia --project=experiments/callmip_phase1a_v2 \\
          experiments/callmip_phase1a_v2/run_posterior_ensemble.jl
"""

using Distributed
import JLD2
import ClimaCalibrate
import ClimaLand
import ClimaLand.Parameters as LP
import EnsembleKalmanProcesses as EKP
import EnsembleKalmanProcesses.ParameterDistributions as PD
using Dates
using Statistics
using Random

const DT              = Float64(900)
const N_ENSEMBLE      = 50
const climaland_dir   = abspath(joinpath(@__DIR__, "..", ".."))
const exp_dir         = @__DIR__
const ces_out_dir     = joinpath(exp_dir, "output_ces")
const ens_output_dir  = joinpath(exp_dir, "output_posterior_ensemble")

isdir(ens_output_dir) || mkpath(ens_output_dir)

include(joinpath(exp_dir, "priors.jl"))
prior, priors_vec = build_dk_sor_priors()
param_names = [only(PD.get_name(d)) for d in priors_vec]

# ── Load MCMC posterior ───────────────────────────────────────────────────────
post_file = joinpath(ces_out_dir, "posterior_fixed_window.jld2")
isfile(post_file) || error("Posterior not found: $post_file\nRun emulate_sample.jl first.")
constrained_posterior = JLD2.load(post_file, "constrained_posterior")
# (n_params × n_samples)

# Draw N_ENSEMBLE samples
rng = Random.MersenneTwister(42)
n_avail = size(constrained_posterior, 2)
idx = rand(rng, 1:n_avail, N_ENSEMBLE)
ensemble_params = constrained_posterior[:, idx]   # (n_params × N_ENSEMBLE)

println("Posterior ensemble: $N_ENSEMBLE members from $n_avail MCMC samples")

# ── Write TOML for each member ────────────────────────────────────────────────
function write_param_toml(path, names, values)
    open(path, "w") do io
        for (name, val) in zip(names, values)
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
end

for m in 1:N_ENSEMBLE
    mpath = ClimaCalibrate.path_to_ensemble_member(ens_output_dir, 0, m)
    isdir(mpath) || mkpath(mpath)
    write_param_toml(
        ClimaCalibrate.parameter_path(ens_output_dir, 0, m),
        param_names,
        ensemble_params[:, m],
    )
end

# ── Broadcast and run ─────────────────────────────────────────────────────────
@everywhere using Distributed
@everywhere import ClimaLand
@everywhere import ClimaLand.Parameters as LP
@everywhere const DT                = $DT
@everywhere const CALLMIP_OUTPUT_DIR = $ens_output_dir
@everywhere include(joinpath(abspath(joinpath(@__DIR__, "..", "..")),
                              "experiments", "callmip_phase1a_v2", "callmip_model_interface.jl"))

println("Running $N_ENSEMBLE posterior ensemble members…")
pmap(m -> ClimaCalibrate.forward_model(0, m), 1:N_ENSEMBLE)

println("\nPosterior ensemble complete. Results in $ens_output_dir")
