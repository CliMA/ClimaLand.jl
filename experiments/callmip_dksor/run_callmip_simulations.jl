"""
Run prior and posterior CalLMIP simulations for DK-Sor Phase 1b.

Each simulation covers 1997–2014 run per-year (each year a separate single-year
run with a 60-day spinup; see callmip_model_interface.jl for the rationale),
with the daily output concatenated:
  - Prior:     default ClimaParams parameters (CalLMIP Protocol Section 6.1:
               "out-of-the-box simulation using default model parameters")
  - Posterior: posterior mean from emulate_sample.jl
               (falls back to EKI final mean if CES coverage < 0.5)

Outputs (in output_callmip_sims/):
  prior/callmip_diagnostics.jld2
  posterior/callmip_diagnostics.jld2

Run sequentially (each sim takes several hours on one node).

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
const DT_TRY       = Float64[900, 450, 225]   # per-year retry ladder
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

# ── Run one calendar year, retrying at smaller DT if it crashes early ──────────
function run_year_with_retry(yr, toml_dict)
    for DT in DT_TRY
        @info "  Year $yr at DT=$(Int(DT)) s"
        surf, col, z_soil, dates =
            run_one_year(yr, toml_dict, MET_NC_PATH; FT, DT)
        n = isempty(dates) ? 0 : length(dates)
        n >= 360 && return surf, col, z_soil, dates
        @warn "  Year $yr produced only $n days at DT=$(Int(DT)) — retrying smaller DT"
    end
    error("Year $yr failed at all DT in $(DT_TRY)")
end

# ── Run a full per-year simulation (1997–2014) and save concatenated output ────
function run_sim(label, param_names, param_values, outdir)
    @info "=== Running $label simulation (per-year, 1997–2014) ==="
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

    all_surface = Dict{String, Vector{Float64}}()
    all_column  = Dict{String, Matrix{Float64}}()
    all_dates   = Date[]
    z_soil_ref  = Float64[]

    t_wall = @elapsed begin
        Logging.with_logger(SimpleLogger(stderr, Logging.Info)) do
            for yr in OUTPUT_YEARS
                surf, col, z_soil, dates = run_year_with_retry(yr, toml_dict)
                isempty(z_soil_ref) && (z_soil_ref = z_soil)
                append!(all_dates, dates)
                for (k, v) in surf
                    all_surface[k] = haskey(all_surface, k) ?
                        vcat(all_surface[k], v) : copy(v)
                end
                for (k, v) in col
                    all_column[k] = haskey(all_column, k) ?
                        hcat(all_column[k], v) : copy(v)
                end
                @info "  Year $yr done; total days so far: $(length(all_dates))"
            end
        end
    end
    @info "$label simulation done in $(round(t_wall/3600; digits=1)) hours " *
          "($(length(all_dates)) days)"

    save_callmip_diagnostics(all_surface, all_column, z_soil_ref, all_dates, outdir)
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
