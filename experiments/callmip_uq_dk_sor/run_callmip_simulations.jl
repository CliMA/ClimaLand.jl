"""
Run CalLMIP Phase 1 prior and posterior simulations at DK-Sor.

Runs two ClimaLand forward simulations using `callmip_model_interface.jl`:
  1. Prior simulation   — default parameters from prior_mean_parameters.toml
  2. Posterior simulation — EKI-optimal parameters (from emulate_sample.jl output,
                            or directly from run_calibration.jl if CES was not run)

Each simulation saves `callmip_diagnostics.jld2` with all CalLMIP-required
variables to:
   OUTPUT_DIR/iteration_000/member_001/callmip_diagnostics.jld2   (prior)
   OUTPUT_DIR/iteration_000/member_002/callmip_diagnostics.jld2   (posterior)

These files are consumed by `write_callmip_netcdf.jl` to produce the
CalLMIP-compliant NetCDF output files.

Usage
-----
    julia --project=.buildkite \\
          experiments/callmip_uq_dk_sor/run_callmip_simulations.jl

Prerequisite: either
  • emulate_sample.jl has been run → reads posterior from output_posterior_uq/
  • OR run_calibration.jl output exists → uses EKI final mean as posterior
"""

using Distributed
import Random
import JLD2
import ClimaCalibrate
import ClimaLand
import EnsembleKalmanProcesses as EKP
import EnsembleKalmanProcesses.ParameterDistributions as PD
using Dates

# ── Configuration ──────────────────────────────────────────────────────────────
const SITE_ID     = "DK-Sor"
const DT          = Float64(450)        # model time step (s), must match calibration
const MODEL_NAME  = "ClimaLand"         # used in CalLMIP output file names
const MODEL_VER   = "CalLMIP1.0"        # CalLMIP submission version label

const climaland_dir   = abspath(joinpath(@__DIR__, "..", ".."))
const cal_dir         = joinpath(climaland_dir, "experiments", "calibrate_dk_sor")   # calibration artifacts
const exp_dir         = joinpath(climaland_dir, "experiments", "callmip_uq_dk_sor")   # UQ/CalLMIP outputs
const cal_output_dir  = joinpath(cal_dir, "output")              # EKP output
const posterior_dir   = joinpath(exp_dir, "output_posterior_uq")  # CES output
const OUTPUT_DIR      = joinpath(exp_dir, "output_callmip_sims")  # this script's output
const OBS_FILEPATH    = joinpath(cal_dir, "observations.jld2")

isdir(OUTPUT_DIR) || mkpath(OUTPUT_DIR)

# Member indices for this run
const MEMBER_PRIOR     = 1    # prior simulation slot
const MEMBER_POSTERIOR = 2    # posterior simulation slot

# ── Priors (must match Alexis's 12-param backup EXACTLY) ─────────────────────
priors_vec = [
    PD.constrained_gaussian("moisture_stress_c",              0.5,     0.3,     0.01,    5.0),
    PD.constrained_gaussian("pmodel_cstar",                   0.43,    0.15,    0.05,    2.0),
    PD.constrained_gaussian("pmodel_β",                      51.0,    20.0,     5.0,  500.0),
    PD.constrained_gaussian("leaf_Cd",                        0.1,     0.05,    0.005,   1.0),
    PD.constrained_gaussian("canopy_z_0m_coeff",              0.05,    0.03,    0.001,   0.3),
    PD.constrained_gaussian("canopy_z_0b_coeff",              0.001,   0.0005,  1e-5,   0.01),
    PD.constrained_gaussian("canopy_d_coeff",                 0.1,     0.05,    0.001,  0.95),
    PD.constrained_gaussian("canopy_K_lw",                    0.85,    0.25,    0.1,     2.0),
    PD.constrained_gaussian("canopy_emissivity",              0.97,    0.02,    0.9,     1.0),
    PD.constrained_gaussian("soilCO2_pre_exponential_factor", 25000.0, 10000.0, 1000.0, 200000.0),
    PD.constrained_gaussian("michaelis_constant",             0.01,    0.005,   1e-4,    0.1),
    PD.constrained_gaussian("O2_michaelis_constant",          0.01,    0.005,   1e-4,    0.1),
]
prior       = PD.combine_distributions(priors_vec)
param_names = [only(PD.get_name(d)) for d in priors_vec]

# ── Helper: write a parameter TOML file ───────────────────────────────────────
function write_parameter_toml(path, names, values)
    open(path, "w") do io
        for (name, val) in zip(names, values)
            used_in = (name == "soilCO2_pre_exponential_factor" ||
                       name == "michaelis_constant"              ||
                       name == "O2_michaelis_constant") ? "[\"Land\"]" : "[\"getindex\"]"
            println(io, "[\"$name\"]")
            println(io, "value = $(Float64(val))")
            println(io, "type  = \"float\"")
            println(io, "used_in = $used_in")
            println(io)
        end
    end
end

# ── Load parameter sets ────────────────────────────────────────────────────────

# — Prior parameters: stated constrained means from the prior definition above.
# Using transform_unconstrained_to_constrained(prior, zeros(...)) gives values
# far from the physical defaults (e.g. cstar→1.025 instead of 0.43) due to
# the nonlinear logit transform, which kills GPP. Use the μ values directly.
const PRIOR_MEANS = Dict(
    "moisture_stress_c"              => 0.27,
    "pmodel_cstar"                   => 0.43,
    "pmodel_β"                       => 51.0,
    "leaf_Cd"                        => 0.07,
    "canopy_z_0m_coeff"              => 0.02,
    "canopy_z_0b_coeff"              => 0.0007,
    "canopy_d_coeff"                 => 0.007,
    "canopy_K_lw"                    => 0.85,
    "canopy_emissivity"              => 0.98,
    "root_leaf_nitrogen_ratio"       => 1.0,
    "stem_leaf_nitrogen_ratio"       => 0.1,
    "soilCO2_pre_exponential_factor" => 23835.0,
    "michaelis_constant"             => 0.005,
    "O2_michaelis_constant"          => 0.004,
)
prior_params = [PRIOR_MEANS[n] for n in param_names]

println("Prior parameters (from prior distribution means):")
for (n, v) in zip(param_names, prior_params)
    println("  $(rpad(n, 40)) $(round(v; sigdigits = 4))")
end

# — Posterior parameters —
# Prefer output from emulate_sample.jl (MCMC posterior mean).
# Fall back to EKI final mean if CES output is not available.
function load_posterior_params()
    # Try CES output first
    ces_files = isdir(posterior_dir) ?
        filter(f -> startswith(f, "posterior_its") && endswith(f, ".jld2"),
               readdir(posterior_dir)) : String[]
    if !isempty(ces_files)
        path = joinpath(posterior_dir, last(sort(ces_files)))
        @info "Loading posterior from CES output: $path"
        d = JLD2.load(path)
        return d["constrained_ekp_optimal"]   # EKI optimum (best point estimate)
    end

    # Fall back: try latest EKP from calibration output
    iter = -1
    while isfile(joinpath(cal_output_dir,
                          "iteration_$(lpad(iter + 1, 3, '0'))",
                          "eki_file.jld2"))
        iter += 1
    end
    if iter >= 0
        ekp_path = joinpath(cal_output_dir,
                            "iteration_$(lpad(iter, 3, '0'))",
                            "eki_file.jld2")
        @info "Loading posterior from EKI final mean: $ekp_path"
        ekp = JLD2.load_object(ekp_path)
        return EKP.get_ϕ_mean_final(prior, ekp)
    end

    error("Could not find posterior parameters — run emulate_sample.jl or " *
          "run_calibration.jl first.")
end

posterior_params = load_posterior_params()

println("\nPosterior parameters (EKI optimal):")
for (n, v) in zip(param_names, posterior_params)
    println("  $(rpad(n, 40)) $(round(v; sigdigits = 4))")
end

# ── Stage parameter TOMLs in ClimaCalibrate directory layout ──────────────────
for (member, params) in [(MEMBER_PRIOR, prior_params),
                          (MEMBER_POSTERIOR, posterior_params)]
    member_dir = ClimaCalibrate.path_to_ensemble_member(OUTPUT_DIR, 0, member)
    isdir(member_dir) || mkpath(member_dir)
    write_parameter_toml(
        ClimaCalibrate.parameter_path(OUTPUT_DIR, 0, member),
        param_names, params,
    )
end
@info "Parameter TOMLs written for members $MEMBER_PRIOR (prior) " *
      "and $MEMBER_POSTERIOR (posterior)."

# ── Save parameter summary ────────────────────────────────────────────────────
JLD2.jldsave(
    joinpath(OUTPUT_DIR, "callmip_parameters.jld2");
    prior_params     = prior_params,
    posterior_params = posterior_params,
    param_names      = param_names,
    model_name       = MODEL_NAME,
    model_version    = MODEL_VER,
    site_id          = SITE_ID,
)

# ── Broadcast configuration to all workers ────────────────────────────────────
@everywhere using Distributed
@everywhere import ClimaLand
@everywhere const SITE_ID     = $SITE_ID
@everywhere const OUTPUT_DIR  = $OUTPUT_DIR
@everywhere const OBS_FILEPATH = $OBS_FILEPATH
@everywhere const DT          = $DT
@everywhere include(joinpath(abspath(joinpath(@__DIR__, "..", "..")),
    "experiments", "callmip_uq_dk_sor", "callmip_model_interface.jl"))

# ── Run both simulations ──────────────────────────────────────────────────────
@info "Running CalLMIP simulations (2 members)…"
@info "  Member $MEMBER_PRIOR  = Prior (default parameters)"
@info "  Member $MEMBER_POSTERIOR = Posterior (EKI-optimal parameters)"

results = pmap([MEMBER_PRIOR, MEMBER_POSTERIOR]) do m
    label = m == MEMBER_PRIOR ? "prior" : "posterior"
    try
        ClimaCalibrate.forward_model(0, m)
        return (member = m, label = label, status = :ok)
    catch e
        @error "Forward model failed for member $m ($label)" exception = e
        return (member = m, label = label, status = :failed)
    end
end

for r in results
    @info "  Member $(r.member) ($(r.label)): $(r.status)"
end

n_failed = count(r -> r.status == :failed, results)
n_failed == 0 ||
    error("$n_failed simulation(s) failed — check logs before proceeding to " *
          "write_callmip_netcdf.jl")

@info "CalLMIP simulations complete. Run write_callmip_netcdf.jl next."
