"""
Run the CalLMIP prior simulation for DK-Sor (Phase 1b).

CalLMIP Protocol Section 6.1: "prior simulation, representing an out-of-the-box
simulation using default model parameters in the model version being used at each
site (i.e. no prior parameter tuning or updating of the model to fit the data)."

→ Uses ClimaParams defaults (LP.create_toml_dict(FT) with no overrides).
→ Can be run immediately in parallel with EKI calibration.

Each year 1997–2014 is run as a separate single-year simulation with a 60-day
spinup (see callmip_model_interface.jl for why per-year and not one continuous
run); the daily output is concatenated. If a year crashes at DT=900s, it is
retried at successively smaller DT (450, 225 s).

Output: output_callmip_sims/prior/callmip_diagnostics.jld2

Usage:
  julia --project=.buildkite \
        experiments/callmip_dksor/run_prior_simulation.jl
"""

import ClimaComms
ClimaComms.@import_required_backends
import ClimaLand
import ClimaLand.Parameters as LP
using ClimaLand
using ClimaDiagnostics
using ClimaUtilities
import ClimaUtilities.TimeVaryingInputs: evaluate!
using Dates
using JLD2
using Logging

const FT            = Float64
const DT_TRY        = Float64[900, 450, 225]   # retry ladder per year
const CLIMALAND_DIR = pkgdir(ClimaLand)
const MET_NC_PATH   = joinpath(CLIMALAND_DIR, "DK_Sor",
    "DK-Sor_1997-2014_FLUXNET2015_Met.nc")
const OUTDIR        = joinpath(CLIMALAND_DIR,
    "experiments/callmip_dksor/output_callmip_sims/prior")

include(joinpath(@__DIR__, "callmip_model_interface.jl"))

mkpath(OUTDIR)
@info "=== CalLMIP Prior simulation (per-year, 1997–2014) ==="
@info "Parameters: ClimaParams defaults (no overrides)"
@info "Output: $OUTDIR"

# Default parameters — no TOML overrides
toml_dict = LP.create_toml_dict(FT)

# Run a single year, retrying at smaller DT if a year crashes mid-run.
function run_year_with_retry(yr)
    for DT in DT_TRY
        @info "Year $yr at DT=$(Int(DT)) s"
        surf, col, z_soil, dates = run_one_year(yr, toml_dict, MET_NC_PATH; FT, DT)
        # A full year must yield ~365 days; fewer means it crashed early.
        n = isempty(dates) ? 0 : length(dates)
        if n >= 360
            return surf, col, z_soil, dates
        end
        @warn "Year $yr produced only $n days at DT=$(Int(DT)) — retrying smaller DT"
    end
    error("Year $yr failed at all DT in $(DT_TRY)")
end

all_surface = Dict{String, Vector{Float64}}()
all_column  = Dict{String, Matrix{Float64}}()
all_dates   = Date[]
z_soil_ref  = Float64[]

t_wall = @elapsed begin
    Logging.with_logger(SimpleLogger(stderr, Logging.Info)) do
        for yr in OUTPUT_YEARS
            surf, col, z_soil, dates = run_year_with_retry(yr)
            isempty(z_soil_ref) && (z_soil_ref = z_soil)
            append!(all_dates, dates)
            for (k, v) in surf
                all_surface[k] = isempty(get(all_surface, k, Float64[])) ?
                    copy(v) : vcat(all_surface[k], v)
            end
            for (k, v) in col
                all_column[k] = haskey(all_column, k) ?
                    hcat(all_column[k], v) : copy(v)
            end
            @info "Year $yr done ($(length(dates)) days); total so far: $(length(all_dates))"
        end
    end
end
@info "Prior simulation done in $(round(t_wall/3600; digits=1)) hours " *
      "($(length(all_dates)) days, $(OUTPUT_YEARS))"

save_callmip_diagnostics(all_surface, all_column, z_soil_ref, all_dates, OUTDIR)
@info "Prior diagnostics saved → $OUTDIR/callmip_diagnostics.jld2"
@info "Next: run write_callmip_netcdf.jl for the prior NetCDF"
