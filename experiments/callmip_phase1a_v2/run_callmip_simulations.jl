"""
Run CalLMIP Phase 1a prior and EKI-optimal posterior simulations at DK-Sor.

Produces:
  output_callmip_sims/iteration_000/member_001/callmip_diagnostics.jld2  (prior)
  output_callmip_sims/iteration_000/member_002/callmip_diagnostics.jld2  (posterior)

These are consumed by write_callmip_netcdf.jl to produce the 2 NetCDF files.

Usage:
    julia --project=experiments/callmip_phase1a_v2 \\
          experiments/callmip_phase1a_v2/run_callmip_simulations.jl
"""

using Distributed
import JLD2
import ClimaCalibrate
import ClimaLand
import ClimaLand.Parameters as LP
import EnsembleKalmanProcesses as EKP
import EnsembleKalmanProcesses.ParameterDistributions as PD
using Dates

const DT              = Float64(900)
const climaland_dir   = abspath(joinpath(@__DIR__, "..", ".."))
const exp_dir         = @__DIR__
const CALLMIP_OUTPUT_DIR = joinpath(exp_dir, "output_callmip_sims")
const cal_output_dir  = joinpath(exp_dir, "output_eki")
const ces_out_dir     = joinpath(exp_dir, "output_ces")

isdir(CALLMIP_OUTPUT_DIR) || mkpath(CALLMIP_OUTPUT_DIR)

# ── Load priors ───────────────────────────────────────────────────────────────
include(joinpath(exp_dir, "priors.jl"))
prior, priors_vec = build_dk_sor_priors()
param_names = [only(PD.get_name(d)) for d in priors_vec]

# ── Helper: write TOML ────────────────────────────────────────────────────────
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

# ── Prior parameters ──────────────────────────────────────────────────────────
# Physical prior means = transform_unconstrained_to_constrained(prior, 0s)
prior_params = EKP.transform_unconstrained_to_constrained(prior, zeros(length(param_names)))

println("Prior mean parameters:")
for (n, v) in zip(param_names, prior_params)
    println("  $(rpad(n, 40))  $(round(v; sigdigits=5))")
end

# ── EKI-optimal (posterior) parameters ───────────────────────────────────────
function load_posterior_params()
    # Try fixed_window EKP output first
    fw_cache = joinpath(exp_dir, "output_fixed_window", "G_fixed_window.jld2")
    if isfile(fw_cache)
        u_matrix = JLD2.load(fw_cache, "u_matrix")
        u_mean   = vec(mean(u_matrix; dims=2))
        ϕ = EKP.transform_unconstrained_to_constrained(prior, u_mean)
        @info "EKI-optimal loaded from fixed-window cache"
        return ϕ
    end

    # Fall back to last EKI iteration file
    eki_toml = joinpath(cal_output_dir, "eki_optimal_parameters.toml")
    if isfile(eki_toml)
        @info "Loading EKI-optimal from $eki_toml"
        d = TOML.parsefile(eki_toml)
        return [Float64(d[n]["value"]) for n in param_names]
    end

    # Last resort: load EKP object
    iter = -1
    while isfile(joinpath(cal_output_dir,
                          "iteration_$(lpad(iter+1, 3, '0'))", "eki_file.jld2"))
        iter += 1
    end
    iter == -1 && error("No EKI output found in $cal_output_dir")
    ekp = JLD2.load_object(joinpath(cal_output_dir,
                                    "iteration_$(lpad(iter, 3, '0'))", "eki_file.jld2"))
    ϕ = EKP.get_ϕ_mean_final(prior, ekp)
    @info "EKI-optimal loaded from EKP (iteration $iter)"
    return ϕ
end

using TOML
posterior_params = load_posterior_params()

println("\nPosterior (EKI-optimal) parameters:")
for (n, v) in zip(param_names, posterior_params)
    println("  $(rpad(n, 40))  $(round(v; sigdigits=5))")
end

# ── Prepare member TOML paths ─────────────────────────────────────────────────
for m in (1, 2)
    mpath = ClimaCalibrate.path_to_ensemble_member(CALLMIP_OUTPUT_DIR, 0, m)
    isdir(mpath) || mkpath(mpath)
end

prior_toml = ClimaCalibrate.parameter_path(CALLMIP_OUTPUT_DIR, 0, 1)
post_toml  = ClimaCalibrate.parameter_path(CALLMIP_OUTPUT_DIR, 0, 2)

write_param_toml(prior_toml, param_names, prior_params)
write_param_toml(post_toml,  param_names, posterior_params)

# ── Broadcast and run ─────────────────────────────────────────────────────────
@everywhere using Distributed
@everywhere import ClimaLand
@everywhere import ClimaLand.Parameters as LP
@everywhere const DT               = $DT
@everywhere const CALLMIP_OUTPUT_DIR = $CALLMIP_OUTPUT_DIR
@everywhere include(joinpath(abspath(joinpath(@__DIR__, "..", "..")),
                              "experiments", "callmip_phase1a_v2", "callmip_model_interface.jl"))

println("\nRunning prior simulation (member 1)…")
ClimaCalibrate.forward_model(0, 1)

println("\nRunning posterior simulation (member 2)…")
ClimaCalibrate.forward_model(0, 2)

println("\nCalLMIP simulations complete.")
println("  Prior:     $(ClimaCalibrate.path_to_ensemble_member(CALLMIP_OUTPUT_DIR, 0, 1))")
println("  Posterior: $(ClimaCalibrate.path_to_ensemble_member(CALLMIP_OUTPUT_DIR, 0, 2))")
println("\nNext step: run write_callmip_netcdf.jl")
