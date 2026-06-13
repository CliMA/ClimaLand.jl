"""
Run prior and posterior CalLMIP simulations for DK-Sor Phase 1b.

Runs two full-period ClimaLand simulations (1996-01-01 → 2015-01-01):
  - Prior:     default ClimaParams parameters (CalLMIP Protocol Section 6.1:
               "out-of-the-box simulation using default model parameters")
  - Posterior: posterior mean from emulate_sample.jl
               (falls back to EKI final mean if CES coverage < 0.5)

Outputs (in output_callmip_sims/):
  prior/callmip_diagnostics.jld2
  posterior/callmip_diagnostics.jld2

Run sequentially (each sim takes ~4h on one node).

Usage:
  julia --project=.buildkite \
        experiments/callmip_dksor/run_callmip_simulations.jl
"""

import ClimaComms
ClimaComms.@import_required_backends
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.Simulations: LandSimulation, solve!
using ClimaLand
using ClimaDiagnostics
using ClimaUtilities
import ClimaUtilities.TimeVaryingInputs: evaluate!
import EnsembleKalmanProcesses.ParameterDistributions as PD
using Dates
using JLD2
using Statistics
using Logging

const FT           = Float64
const DT           = Float64(450)
const CLIMALAND_DIR = pkgdir(ClimaLand)
const MET_NC_PATH  = joinpath(CLIMALAND_DIR, "DK_Sor",
    "DK-Sor_1997-2014_FLUXNET2015_Met.nc")
const OUTDIR       = joinpath(CLIMALAND_DIR,
    "experiments/callmip_dksor/output_callmip_sims")
const CAL_OUTDIR   = joinpath(CLIMALAND_DIR,
    "experiments/callmip_dksor/output_calibration")

include(joinpath(@__DIR__, "priors.jl"))
include(joinpath(@__DIR__, "callmip_model_interface.jl"))

# ── Build parameter TOML override file ────────────────────────────────────────
function write_override_toml(path, param_names, param_values)
    open(path, "w") do io
        for (name, val) in zip(param_names, Float64.(param_values))
            println(io, "[\"$(name)\"]")
            println(io, "value = $(val)")
            println(io, "type  = \"float\"")
            println(io)
        end
    end
end

# ── Run one simulation and save diagnostics ────────────────────────────────────
function run_sim(label, param_names, param_values, outdir)
    @info "=== Running $label simulation ==="
    mkpath(outdir)

    # Prior: no overrides → use ClimaParams defaults (CalLMIP Section 6.1)
    if isempty(param_names)
        toml_dict = LP.create_toml_dict(FT)
    else
        tmpfile = tempname() * ".toml"
        write_override_toml(tmpfile, param_names, Float64.(param_values))
        toml_dict = LP.create_toml_dict(FT; override_files = [tmpfile])
        rm(tmpfile; force = true)
    end

    land, forcing, ν, θ_r = build_callmip_model(FT, toml_dict, MET_NC_PATH)
    set_ic! = make_callmip_ic(ν, θ_r, forcing.atmos)

    # Request CalLMIP diagnostic variables available for LandModel.
    # Note: "soillhf", "soilrn", "soilshf" not in get_possible_diagnostics(LandModel).
    output_writer = ClimaDiagnostics.Writers.DictWriter()
    surface_vars = ["nee", "lhf", "shf", "gpp", "er", "trans", "ct", "lai", "cveg"]
    diags_surface = ClimaLand.default_diagnostics(
        land, SIM_START, "";
        output_writer,
        output_vars      = surface_vars,
        reduction_period = :daily,
    )
    diags_col = ClimaLand.default_diagnostics(
        land, SIM_START, "";
        output_writer,
        output_vars      = ["swc", "tsoil", "soc"],
        reduction_period = :daily,
    )

    simulation = LandSimulation(
        SIM_START, SIM_STOP, DT, land;
        set_ic!, updateat = Second(DT),
        diagnostics = vcat(diags_surface, diags_col),
    )

    t_wall = @elapsed begin
        Logging.with_logger(SimpleLogger(stderr, Logging.Info)) do
            solve!(simulation)
        end
    end
    @info "$label simulation done in $(round(t_wall/3600; digits=1)) hours"

    save_callmip_diagnostics(simulation, land, outdir)
    @info "$label diagnostics saved → $outdir/callmip_diagnostics.jld2"
end

# ── Prior simulation (CalLMIP Protocol Section 6.1) ───────────────────────────
# "out-of-the-box simulation using default model parameters in the model
# version being used at each site (i.e. no prior parameter tuning)"
# → empty param lists → run_sim uses LP.create_toml_dict(FT) with no overrides.
@info "Running prior with default ClimaParams parameters (no overrides)"
run_sim("prior", String[], Float64[], joinpath(OUTDIR, "prior"))

# ── Posterior simulation (requires emulate_sample.jl to have run first) ───────
post_path = joinpath(CAL_OUTDIR, "posterior_mean.jld2")
isfile(post_path) || error("posterior_mean.jld2 not found: $post_path\n" *
                            "Run emulate_sample.jl first.")
post_data = JLD2.load(post_path)
posterior_mean   = post_data["posterior_mean"]
param_names_post = post_data["param_names"]
@info "Posterior mean loaded (used_eki_fallback=$(post_data["used_eki_fallback"])):"
for (n, v) in zip(param_names_post, posterior_mean)
    @info "  $(rpad(n, 35)) = $(round(v; sigdigits=4))"
end

run_sim("posterior",
    param_names_post, posterior_mean,
    joinpath(OUTDIR, "posterior"))

@info "Both CalLMIP simulations complete."
@info "Next: run write_callmip_netcdf.jl"
