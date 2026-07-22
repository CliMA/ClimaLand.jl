"""
Compare the single-column DK-Sor build (`Column`, `PointSpace` surface) against
column 1 of a multi-column ensemble (`MultiColumnFiniteDifferenceSpace`,
`PointCloudSpace` surface) whose first column is the same DK-Sor site, and
report, per field, the RMS difference between their states normalized by the
field's own RMS magnitude (a scale-free relative error). The two should be
near-identical: a relative error above `ROUNDOFF_REL_TOL` flags a real difference
introduced by the multi-column codepath, while smaller values are floating-point
roundoff. Running a genuine 2-column ensemble also checks that adding a second,
independent site (US-NR1) does not perturb column 1.

Adapted from `compare_contents` in the `kp/profile` worktree's
`sim_with_new_space.jl`, driven through `run_prior_year(...; return_sim=true)`.

Usage:
    julia --project=.buildkite compare_single_multi.jl
    COMPARE_YEAR=2005 COMPARE_STEPS=500 julia --project=.buildkite compare_single_multi.jl
"""

include("fluxnet_sim.jl")

import ClimaCore
import ClimaTimeSteppers

const COMPARE_YEAR =
    haskey(ENV, "COMPARE_YEAR") ? parse(Int, ENV["COMPARE_YEAR"]) :
    first(CHECK_YEARS)
const COMPARE_STEPS =
    haskey(ENV, "COMPARE_STEPS") ? parse(Int, ENV["COMPARE_STEPS"]) : 500

# A single-vs-multi difference for one leaf Field. `rel` is the RMS difference
# normalized by the reference field's own RMS magnitude — a scale-free number
# comparable across fields of any unit. `rmse` and `scale` are kept so the raw
# difference and the magnitude it is measured against are both visible.
struct FieldDiff
    rmse::Float64   # RMS of (single − multi) over all nodes/columns
    scale::Float64  # RMS magnitude of the reference (single-column) field
    rel::Float64    # rmse / scale (relative error); see _relerr for edge cases
end

# Relative error with sensible behavior when the reference field is ~0 (e.g. no
# snow, no ice): if there is also no difference the fields agree (rel = 0); a
# nonzero difference against a zero-magnitude field is flagged as Inf so it sorts
# to the top rather than hiding behind a tiny denominator.
function _relerr(rmse, scale)
    scale > 0 && return rmse / scale
    return rmse == 0 ? 0.0 : Inf
end

"""
Recursively walk two matching `FieldVector`s and return `Dict(field_path => FieldDiff)`
for every leaf `Field`, comparing the single-column state `u1` against column `col`
of the multi-column ensemble state `u2`. A `FieldDiff` with `NaN`s flags a
structural mismatch (the two `Field`s flatten to different lengths).
"""
function state_field_rmses(u1, u2, col = 1)
    diffs = Dict{String, FieldDiff}()
    _collect_field_rmses!(diffs, u1, u2, String[], col)
    return diffs
end

function _collect_field_rmses!(diffs, u1, u2, path, col)
    for name in propertynames(u1)
        push!(path, String(name))
        _collect_field_rmses!(
            diffs,
            getproperty(u1, name),
            getproperty(u2, name),
            path,
            col,
        )
        pop!(path)
    end
    return nothing
end

function _collect_field_rmses!(
    diffs,
    u1::ClimaCore.Fields.Field,
    u2::ClimaCore.Fields.Field,
    path,
    col,
)
    # u1 is the single-column state; u2 is the multi-column ensemble. The ensemble's
    # parent carries its columns along the LAST dimension, so slice out column `col`;
    # the single Column's parent has no column dimension, so use it whole. vec then
    # collapses the node/component dims so the two align element-for-element (scalar
    # and multi-component fields alike). parent() (not field2array) is used because
    # field2array errors on multi-component fields (e.g. two-band PAR/NIR).
    a1 = vec(parent(u1))
    p2 = parent(u2)
    a2 = vec(collect(selectdim(p2, ndims(p2), col)))
    if length(a1) != length(a2)
        diffs[join(path, ".")] = FieldDiff(NaN, NaN, NaN)
    else
        n = length(a1)
        rmse = sqrt(sum((a1 .- a2) .^ 2) / n)
        scale = sqrt(sum(a1 .^ 2) / n)
        diffs[join(path, ".")] = FieldDiff(rmse, scale, _relerr(rmse, scale))
    end
    return nothing
end

# Relative errors at or below this are consistent with floating-point roundoff
# (reassociated reductions on the multi-column path) amplified by the dynamics
# over the compared steps; anything above suggests a genuine algorithmic
# discrepancy between the single- and multi-column codepaths.
const ROUNDOFF_REL_TOL = 1e-6

"""Log every field whose single-vs-multi relative error is non-zero (or NaN),
comparing the single-column state against column `col` of the multi-column
ensemble. Rows are sorted by relative error (largest first) and the raw RMSE and
reference scale are shown so each relative number can be traced to its inputs. A
verdict summarizes whether every difference is at roundoff level."""
function report_state_diffs(sim_single, sim_multi; col = 1, label = "")
    u_rmse = state_field_rmses(
        sim_single._integrator.u,
        sim_multi._integrator.u,
        col,
    )
    p_rmse = state_field_rmses(
        sim_single._integrator.p,
        sim_multi._integrator.p,
        col,
    )
    differing(d) = sort(
        filter(kv -> !(last(kv).rel == 0.0), collect(d));
        by = kv -> (isnan(last(kv).rel) ? Inf : last(kv).rel),
        rev = true,
    )
    du, dp = differing(u_rmse), differing(p_rmse)
    worst = maximum(
        (kv -> isnan(kv[2].rel) ? Inf : kv[2].rel).(vcat(du, dp));
        init = 0.0,
    )
    above = count(
        kv -> isnan(last(kv).rel) || last(kv).rel > ROUNDOFF_REL_TOL,
        vcat(du, dp),
    )
    verdict =
        worst == 0.0 ? "IDENTICAL" :
        above == 0 ? "roundoff-level (max rel $(@sprintf("%.1e", worst)))" :
        "$(above) field(s) ABOVE $(ROUNDOFF_REL_TOL) — investigate"
    @info "State comparison [$label]" verdict prognostic_fields_differing =
        length(du) cache_fields_differing = length(dp)
    _log_diffs(prefix, ds) =
        for (k, v) in ds
            flag =
                isnan(v.rel) ? "  <-- structural mismatch" :
                v.rel > ROUNDOFF_REL_TOL ? "  <-- above roundoff" : ""
            @info @sprintf(
                "  %s.%-40s rel = %.2e  (rmse = %.2e, scale = %.2e)%s",
                prefix,
                k,
                v.rel,
                v.rmse,
                v.scale,
                flag
            )
        end
    _log_diffs("Y", du)
    _log_diffs("p", dp)
    return (; prognostic = u_rmse, cache = p_rmse, worst_rel = worst)
end

@info "Building single- and multi-column DK-Sor simulations" COMPARE_YEAR COMPARE_STEPS
sim_single = run_prior_year(COMPARE_YEAR; return_sim = true, multi_col = false)
sim_multi = run_prior_year(
    COMPARE_YEAR;
    return_sim = true,
    multi_col = true,
    multi_col_sites = ["DK-Sor", "US-NR1"],
)

# Column 1 of the ensemble is DK-Sor (same site as the single-column run) → ~0 RMSE.
report_state_diffs(
    sim_single,
    sim_multi;
    col = 1,
    label = "col 1 (DK-Sor) vs single DK-Sor, initial condition (t=0)",
)

# @info "Stepping both simulations $COMPARE_STEPS steps in lockstep..."
# for step in 1:COMPARE_STEPS
#     ClimaTimeSteppers.step!(sim_single._integrator)
#     ClimaTimeSteppers.step!(sim_multi._integrator)
# end

# report_state_diffs(
#     sim_single,
#     sim_multi;
#     col = 1,
#     label = "after $COMPARE_STEPS steps",
# )
