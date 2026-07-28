"""
Per-field RMSE comparisons between two simulation states, copied from
`compare_single_multi.jl` and `verify_multicol.jl` (kp/multi-cols).

- [`state_field_rmses`](@ref) / [`report_state_diffs`](@ref): single-column
  state (`Column`) vs column `col` of a multi-column ensemble state
  (`MultiColumnFiniteDifferenceSpace`), reported as RMS differences normalized
  by each field's own RMS magnitude (scale-free relative errors).
- [`column_field_rmses`](@ref) / [`report_column_diffs`](@ref): column `col` of
  an N-column state vs column 1 of a 1-column state, reported as raw RMSEs.

Both share one traversal ([`leaf_pairs`](@ref)) and one column extractor
([`flatten`](@ref)); only the metric differs.
"""

import ClimaCore
using Printf

# Root-mean-square of a flat vector.
rms(x) = sqrt(sum(x .^ 2) / length(x))

# Flatten one leaf Field to a plain vector. Without `col`, use the whole field
# (the single Column has no column dimension). With `col`, slice out that
# ensemble column, which is the LAST parent dimension. parent() is used rather
# than field2array so multi-component fields (e.g. two-band PAR/NIR) don't error.
flatten(field::ClimaCore.Fields.Field) = vec(parent(field))
flatten(field::ClimaCore.Fields.Field, col) =
    vec(collect(selectdim(parent(field), ndims(parent(field)), col)))

"""
Walk two matching `FieldVector`s in lockstep and return
`Dict(dotted_path => (leaf_field_1, leaf_field_2))` for every leaf `Field`.
"""
function leaf_pairs(x1, x2)
    pairs =
        Dict{String, Tuple{ClimaCore.Fields.Field, ClimaCore.Fields.Field}}()
    _collect_leaves!(pairs, x1, x2, String[])
    return pairs
end

# A leaf: record the two Fields under the current dotted path.
_collect_leaves!(
    pairs,
    f1::ClimaCore.Fields.Field,
    f2::ClimaCore.Fields.Field,
    path,
) = (pairs[join(path, ".")] = (f1, f2); pairs)

# A branch (FieldVector or named tuple): recurse into each named child.
function _collect_leaves!(pairs, x1, x2, path)
    for name in propertynames(x1)
        push!(path, String(name))
        _collect_leaves!(
            pairs,
            getproperty(x1, name),
            getproperty(x2, name),
            path,
        )
        pop!(path)
    end
    return pairs
end

# A single-vs-multi difference for one leaf Field. `rel` normalizes `rmse` by the
# reference field's own RMS magnitude `scale`, giving a scale-free number
# comparable across fields of any unit.
struct FieldDiff
    rmse::Float64   # RMS of (single − multi)
    scale::Float64  # RMS magnitude of the reference (single-column) field
    rel::Float64    # rmse / scale (see below for the zero-scale cases)
end

# Relative error, handling a ~0 reference field (e.g. no snow/ice): agreeing
# fields give 0; a nonzero difference against a zero-magnitude field gives Inf
# rather than hiding behind a tiny denominator.
_relerr(rmse, scale) = scale > 0 ? rmse / scale : (rmse == 0 ? 0.0 : Inf)

"""
Return `Dict(field_path => FieldDiff)` for every leaf `Field`, comparing the
single-column state `u1` against column `col` of the multi-column ensemble state
`u2`. A `FieldDiff` of `NaN`s flags a structural mismatch (different lengths).
"""
function state_field_rmses(u1, u2, col = 1)
    diffs = Dict{String, FieldDiff}()
    for (path, (f1, f2)) in leaf_pairs(u1, u2)
        a1 = flatten(f1)            # single column, whole
        a2 = flatten(f2, col)       # column `col` of the ensemble
        if length(a1) == length(a2)
            rmse, scale = rms(a1 .- a2), rms(a1)
            diffs[path] = FieldDiff(rmse, scale, _relerr(rmse, scale))
        else
            diffs[path] = FieldDiff(NaN, NaN, NaN)
        end
    end
    return diffs
end

"""Log the relative error (with the raw RMSE and reference scale behind it) of
every leaf field — prognostic (`Y`) and cache (`p`) — comparing the
single-column state against column `col` of the multi-column ensemble. All
fields are listed, sorted by path, so successive reports line up row-for-row."""
function report_state_diffs(sim_single, sim_multi; col = 1, label = "")
    @info "State comparison [$label]"
    for (prefix, u1, u2) in (
        ("Y", sim_single._integrator.u, sim_multi._integrator.u),
        ("p", sim_single._integrator.p, sim_multi._integrator.p),
    )
        for (k, v) in sort!(collect(state_field_rmses(u1, u2, col)); by = first)
            @info @sprintf(
                "  %s.%-40s rel = %.2e  (rmse = %.2e, scale = %.2e)",
                prefix,
                k,
                v.rel,
                v.rmse,
                v.scale
            )
        end
    end
    return nothing
end

"""
Return `Dict(field_path => RMSE)` for every leaf `Field`, comparing column `col`
of the N-column state `uN` against column 1 of the 1-column state `u1`.
"""
function column_field_rmses(uN, u1, col)
    diffs = Dict{String, Float64}()
    for (path, (fN, f1)) in leaf_pairs(uN, u1)
        xN = flatten(fN, col)       # column `col` of the N-column state
        x1 = flatten(f1, 1)         # the only column of the 1-column state
        diffs[path] = length(xN) == length(x1) ? rms(xN .- x1) : NaN
    end
    return diffs
end

"""Log the RMSE of every leaf field between column `col` of an N-column state
and column 1 of a 1-column state. All fields are listed, sorted by path."""
function report_column_diffs(uN, u1, col; label = "")
    diffs = column_field_rmses(uN, u1, col)
    @info "Column-vs-single [$label]"
    for (k, v) in sort!(collect(diffs); by = first)
        @info @sprintf("  Y.%-45s RMSE = %.3e", k, v)
    end
    return diffs
end
