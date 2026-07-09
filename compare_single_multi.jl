"""
Compare the single-column (`Column`, `PointSpace` surface) and multi-column
(`MultiColumnFiniteDifferenceSpace`, `PointCloudLevelSpace` surface; N=1 here)
builds of the DK-Sor simulation and report the per-field RMSE between their
states. For the single-point ensemble the two should be near-identical; any
non-zero RMSE flags a difference introduced by the multi-column codepath.

Adapted from `compare_contents` in the `kp/profile` worktree's
`sim_with_new_space.jl`, driven through `run_prior_year(...; return_sim=true)`.

Usage:
    julia --project=.buildkite compare_single_multi.jl
    COMPARE_YEAR=2005 COMPARE_STEPS=500 julia --project=.buildkite compare_single_multi.jl
    COMPARE_FLUXES=1 julia --project=.buildkite compare_single_multi.jl   # also compare monthly NEE/LHF/SHF
"""

include("fluxnet_sim.jl")

import ClimaCore
import ClimaTimeSteppers

const COMPARE_YEAR =
    haskey(ENV, "COMPARE_YEAR") ? parse(Int, ENV["COMPARE_YEAR"]) : first(CHECK_YEARS)
const COMPARE_STEPS =
    haskey(ENV, "COMPARE_STEPS") ? parse(Int, ENV["COMPARE_STEPS"]) : 500

"""
Recursively walk two matching `FieldVector`s and return
`Dict(field_path => rmse)` for every leaf `Field`. `rmse` is `NaN` when the two
`Field`s flatten to different lengths (a genuine structural mismatch).
"""
function state_field_rmses(u1, u2)
    diffs = Dict{String, Float64}()
    _collect_field_rmses!(diffs, u1, u2, String[])
    return diffs
end

function _collect_field_rmses!(diffs, u1, u2, path)
    for name in propertynames(u1)
        push!(path, String(name))
        _collect_field_rmses!(diffs, getproperty(u1, name), getproperty(u2, name), path)
        pop!(path)
    end
    return nothing
end

function _collect_field_rmses!(
    diffs,
    u1::ClimaCore.Fields.Field,
    u2::ClimaCore.Fields.Field,
    path,
)
    # Use parent(), not field2array: the latter only handles scalar-eltype fields
    # and errors on multi-component ones (e.g. two-band PAR/NIR surface fields).
    # parent gives the raw backing array for any eltype/space; vec collapses the
    # singleton node/column dims so the N=1 multi field aligns element-for-element
    # with the single-column field.
    a1 = vec(parent(u1))
    a2 = vec(parent(u2))
    diffs[join(path, ".")] =
        length(a1) == length(a2) ? sqrt(sum((a1 .- a2) .^ 2) / length(a1)) : NaN
    return nothing
end

"""Log every field whose single-vs-multi RMSE is non-zero (or NaN)."""
function report_state_diffs(sim_single, sim_multi; label = "")
    u_rmse = state_field_rmses(sim_single._integrator.u, sim_multi._integrator.u)
    p_rmse = state_field_rmses(sim_single._integrator.p, sim_multi._integrator.p)
    differing(d) = sort(
        filter(kv -> !(last(kv) == 0.0), collect(d));
        by = kv -> (isnan(last(kv)) ? Inf : last(kv)),
        rev = true,
    )
    du, dp = differing(u_rmse), differing(p_rmse)
    @info "State comparison [$label]" prognostic_fields_differing = length(du) cache_fields_differing =
        length(dp)
    for (k, v) in du
        @info @sprintf("  Y.%-45s RMSE = %.3e", k, v)
    end
    for (k, v) in dp
        @info @sprintf("  p.%-45s RMSE = %.3e", k, v)
    end
    return (; prognostic = u_rmse, cache = p_rmse)
end

@info "Building single- and multi-column DK-Sor simulations" COMPARE_YEAR COMPARE_STEPS
sim_single = run_prior_year(COMPARE_YEAR; return_sim = true, multi_col = false)
sim_multi = run_prior_year(COMPARE_YEAR; return_sim = true, multi_col = true)

report_state_diffs(sim_single, sim_multi; label = "initial condition (t=0)")

@info "Stepping both simulations $COMPARE_STEPS steps in lockstep..."
for step in 1:COMPARE_STEPS
    ClimaTimeSteppers.step!(sim_single._integrator)
    ClimaTimeSteppers.step!(sim_multi._integrator)
end

report_state_diffs(sim_single, sim_multi; label = "after $COMPARE_STEPS steps")

# Optional: compare the year's monthly NEE/LHF/SHF end-to-end (reuses rmse_r2).
if haskey(ENV, "COMPARE_FLUXES")
    @info "Running full year $COMPARE_YEAR for both space types (flux comparison)..."
    nee_s, lhf_s, shf_s = run_prior_year(COMPARE_YEAR; multi_col = false)
    nee_m, lhf_m, shf_m = run_prior_year(COMPARE_YEAR; multi_col = true)
    for (name, s, m) in
        (("NEE", nee_s, nee_m), ("LHF", lhf_s, lhf_m), ("SHF", shf_s, shf_m))
        rmse, r2 = rmse_r2(m, s)   # multi vs single
        @info @sprintf("  %s  single-vs-multi RMSE = %.4e  R² = %.5f", name, rmse, r2)
    end
end
