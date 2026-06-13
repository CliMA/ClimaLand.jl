"""
Run the CalLMIP prior simulation for DK-Sor (Phase 1b).

CalLMIP Protocol Section 6.1: "prior simulation, representing an out-of-the-box
simulation using default model parameters in the model version being used at each
site (i.e. no prior parameter tuning or updating of the model to fit the data)."

→ Uses ClimaParams defaults (LP.create_toml_dict(FT) with no overrides).
→ Can be run immediately in parallel with EKI calibration.

Output: output_callmip_sims/prior/callmip_diagnostics.jld2

Usage:
  julia --project=.buildkite \
        experiments/callmip_dksor/run_prior_simulation.jl
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
using Dates
using JLD2
using Logging

const FT            = Float64
const DT            = Float64(450)
const CLIMALAND_DIR = pkgdir(ClimaLand)
const MET_NC_PATH   = joinpath(CLIMALAND_DIR, "DK_Sor",
    "DK-Sor_1997-2014_FLUXNET2015_Met.nc")
const OUTDIR        = joinpath(CLIMALAND_DIR,
    "experiments/callmip_dksor/output_callmip_sims/prior")

include(joinpath(@__DIR__, "callmip_model_interface.jl"))

mkpath(OUTDIR)
@info "=== CalLMIP Prior simulation ==="
@info "Parameters: ClimaParams defaults (no overrides)"
@info "Output: $OUTDIR"

# Default parameters — no TOML overrides
toml_dict = LP.create_toml_dict(FT)

land, forcing, ν, θ_r = build_callmip_model(FT, toml_dict, MET_NC_PATH)
set_ic! = make_callmip_ic(ν, θ_r, forcing.atmos)

output_writer = ClimaDiagnostics.Writers.DictWriter()
# Note: "soillhf", "soilrn", "soilshf" are not in get_possible_diagnostics(LandModel)
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
@info "Prior simulation done in $(round(t_wall/3600; digits=1)) hours"

save_callmip_diagnostics(simulation, land, OUTDIR)
@info "Prior diagnostics saved → $OUTDIR/callmip_diagnostics.jld2"
@info "Next: run write_callmip_netcdf.jl for the prior NetCDF"
