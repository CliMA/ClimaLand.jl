"""
Prior mean sanity check for CalLMIP Phase 1a DK-Sor.

Runs a single 1997-01-01 → 2014-12-31 ClimaLand simulation at the prior mean
parameters and checks:
  - No NaN blow-up (>95% of days must be finite)
  - Summer GPP is positive (P-model working)
  - Qle seasonal cycle present
  - Prior NEE RMSE vs FLUXNET is finite and < 10 gC/m²/d mean absolute error

Step 0b gate: this MUST pass before running EKI.

Usage:
    julia --project=experiments/callmip_phase1a_v2 \\
          experiments/callmip_phase1a_v2/run_prior_check.jl
"""

import ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.Simulations: LandSimulation, solve!
using ClimaLand
using ClimaLand.Domains: Column
using ClimaLand.Soil
using ClimaLand.Soil.Biogeochemistry
using ClimaLand.Canopy
using ClimaCore
using ClimaDiagnostics
using ClimaUtilities
import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput, evaluate!
import ClimaUtilities.TimeManager: date

using Dates
using NCDatasets
using Statistics
import JLD2

const FT   = Float64
const DT   = Float64(900)

const climaland_dir  = abspath(joinpath(@__DIR__, "..", ".."))
const met_nc_path    = joinpath(climaland_dir, "DK_Sor",
                                 "DK-Sor_1997-2014_FLUXNET2015_Met.nc")
const flux_nc_path   = joinpath(climaland_dir, "DK_Sor",
                                 "DK-Sor_daily_aggregated_1997-2013_FLUXNET2015_Flux.nc")
const PRIOR_TOML     = joinpath(@__DIR__, "prior_mean_parameters.toml")

include(joinpath(@__DIR__, "model_interface.jl"))

# ── Run ───────────────────────────────────────────────────────────────────────
toml_dict = LP.create_toml_dict(FT; override_files = [PRIOR_TOML])

# DEBUG: 30-day run to identify NaN source in soilco2
sim_start = DateTime(1996, 1, 1)
sim_stop  = DateTime(1996, 2, 1)

println("Prior mean check: $sim_start → $sim_stop")
land, forcing, ν, θ_r = build_dk_sor_model(FT, sim_start, sim_stop, toml_dict, met_nc_path)
set_ic! = make_dk_sor_ic(forcing.atmos, ν, θ_r)

output_writer = ClimaDiagnostics.Writers.DictWriter()
diags = ClimaLand.default_diagnostics(
    land, sim_start, "";   # outdir unused (DictWriter handles output)
    output_writer,
    output_vars      = :short,
    reduction_period = :daily,
)

simulation = LandSimulation(
    sim_start, sim_stop, DT, land;
    set_ic!,
    updateat       = Second(DT),
    user_callbacks = (),
    diagnostics    = diags,
)
solve!(simulation)

# ── Extract daily diagnostics ─────────────────────────────────────────────────
# Print available DictWriter keys for debugging
all_keys = unique(reduce(vcat, [collect(keys(d.output_writer.dict)) for d in simulation.diagnostics]))
println("Available diagnostic keys: ", sort(all_keys))

nee_dates, nee_vals = extract_daily_diag(simulation, "nee_1d_average")
_,         qle_vals = extract_daily_diag(simulation, "lhf_1d_average")
_,         qh_vals  = extract_daily_diag(simulation, "shf_1d_average")
_,         gpp_vals = extract_daily_diag(simulation, "gpp_1d_average")

# Diagnostic summary before spinup filtering
println("NEE first 5: ", nee_vals[1:min(5,end)])
println("GPP first 5: ", gpp_vals[1:min(5,end)])
println("LHF first 5: ", qle_vals[1:min(5,end)])
println("NEE NaN count: $(sum(isnan.(nee_vals))) / $(length(nee_vals))")
println("GPP NaN count: $(sum(isnan.(gpp_vals))) / $(length(gpp_vals))")
println("LHF NaN count: $(sum(isnan.(qle_vals))) / $(length(qle_vals))")

# ── Extra diagnostics to pinpoint NaN source ─────────────────────────────────
_, hr_vals  = extract_daily_diag(simulation, "hr_1d_average")
_, ra_vals  = extract_daily_diag(simulation, "ra_1d_average")
_, er_vals  = extract_daily_diag(simulation, "er_1d_average")
println("HR (heterotrophic resp=soilco2 top_bc*83.26) first 5: ", hr_vals[1:min(5,end)])
println("RA (autotrophic resp) first 5: ", ra_vals[1:min(5,end)])
println("ER (total resp=HR+RA) first 5: ", er_vals[1:min(5,end)])
println("HR NaN count: $(sum(isnan.(hr_vals))) / $(length(hr_vals))")
println("RA NaN count: $(sum(isnan.(ra_vals))) / $(length(ra_vals))")
println("ER NaN count: $(sum(isnan.(er_vals))) / $(length(er_vals))")
# Is Y.soilco2.CO2 itself NaN? (sco2 = soil CO2 diagnostic)
_, sco2_vals = extract_daily_diag(simulation, "sco2_1d_average")
println("SCO2 (soil CO2 conc) first 5: ", sco2_vals[1:min(5,end)])
println("SCO2 NaN count: $(sum(isnan.(sco2_vals))) / $(length(sco2_vals)))")

# Early exit — debug run only
exit(0)

# Discard spinup (keep 1997–2014)
keep = year.(nee_dates) .>= 1997
nee_dates = nee_dates[keep]
nee_vals  = nee_vals[keep]
qle_vals  = qle_vals[keep]
qh_vals   = qh_vals[keep]
gpp_vals  = gpp_vals[keep]

# ── Load FLUXNET observations for comparison ──────────────────────────────────
ds = NCDataset(flux_nc_path, "r")
flux_times = ds["time"][:]
nee_obs_raw = Float64.(coalesce.(ds["NEE_daily"][:], NaN))
qle_obs_raw = Float64.(coalesce.(ds["Qle_daily"][:], NaN))
close(ds)

flux_dates   = Date.(flux_times)
date_to_idx  = Dict{Date, Int}(d => i for (i, d) in enumerate(flux_dates))

nee_obs = [get(date_to_idx, d, 0) == 0 ? NaN :
           nee_obs_raw[date_to_idx[d]] for d in nee_dates]
qle_obs = [get(date_to_idx, d, 0) == 0 ? NaN :
           qle_obs_raw[date_to_idx[d]] for d in nee_dates]

# ── Checks ────────────────────────────────────────────────────────────────────

# Convert model NEE: mol CO2/m²/s → gC/m²/d
nee_model_gC = nee_vals .* 12.0 * 86400.0

# NaN check
n_nan = sum(isnan.(nee_vals))
pct_valid = 1.0 - n_nan / length(nee_vals)
println("\nNaN check: $(round(100*pct_valid; digits=1))% valid days")
@assert pct_valid >= 0.95 "FAIL: >5% NaN in NEE"

# Summer GPP check (JJA: June-August)
summer_mask = month.(nee_dates) .∈ Ref([6, 7, 8])
gpp_summer  = mean(gpp_vals[summer_mask .& .!isnan.(gpp_vals)])
# GPP is in mol CO2/m²/s; convert to gC/m²/d
gpp_summer_gC = gpp_summer * 12.0 * 86400.0
println("Summer (JJA) mean GPP: $(round(gpp_summer_gC; digits=3)) gC/m²/d")
@assert gpp_summer_gC > 0.0 "FAIL: Summer GPP is non-positive (P-model not working)"

# Qle seasonal cycle (summer mean > winter mean)
winter_mask = month.(nee_dates) .∈ Ref([12, 1, 2])
qle_summer  = mean(qle_vals[summer_mask .& .!isnan.(qle_vals)])
qle_winter  = mean(qle_vals[winter_mask .& .!isnan.(qle_vals)])
println("Qle: summer $(round(qle_summer; digits=1)) W/m², winter $(round(qle_winter; digits=1)) W/m²")
@assert qle_summer > qle_winter "FAIL: Qle seasonal cycle absent"

# NEE RMSE vs observations
valid_both  = .!isnan.(nee_obs) .& .!isnan.(nee_model_gC)
nee_rmse    = sqrt(mean((nee_obs[valid_both] .- nee_model_gC[valid_both]) .^ 2))
nee_mae     = mean(abs.(nee_obs[valid_both] .- nee_model_gC[valid_both]))
println("Prior NEE RMSE: $(round(nee_rmse; digits=3)) gC/m²/d  (vs obs, $(sum(valid_both)) days)")
println("Prior NEE MAE : $(round(nee_mae;  digits=3)) gC/m²/d")
@assert isfinite(nee_rmse)    "FAIL: NEE RMSE is Inf/NaN"
@assert nee_mae < 10.0       "FAIL: NEE MAE >= 10 gC/m²/d (model blown up)"

println("\n✓ Prior mean check PASSED — proceed with EKI calibration.")

# Save for record
JLD2.jldsave(
    joinpath(@__DIR__, "prior_check_output.jld2");
    dates = nee_dates, nee_model = nee_model_gC,
    gpp_model = gpp_vals .* 12.0 * 86400.0,
    qle_model = qle_vals, qh_model = qh_vals,
    nee_obs, qle_obs,
    nee_rmse, nee_mae,
)
println("Results saved → $(joinpath(@__DIR__, "prior_check_output.jld2"))")
