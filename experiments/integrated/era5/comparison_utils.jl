# Utilities for comparing, field by field, two states from simulations on
# column spaces, modeled on ClimaAtmos's `test/restart_utils.jl`. This is
# accomplished by recursively iterating over the fields of the integrator
# while ignoring selected fields.

import ClimaCore
using Printf

"""
    rms(x)

Compute the root mean square of a vector.
"""
rms(x) = sqrt(sum(abs2, x) / length(x))

"""
    mixed_error(x1, x2; abs_floor = 100 * eps(eltype(x1)))

Compute the worst elementwise absolute or relative error of `x1` and `x2`.

For each entry, a relative error is computed when the magnitude of `x1` is
greater than or equal to `abs_floor`. Otherwise, an absolute error is computed.
"""
function mixed_error(x1, x2; abs_floor = 100 * eps(eltype(x1)))
    diff = abs.(x1 .- x2)
    return maximum(ifelse.(abs.(x1) .> abs_floor, diff ./ abs.(x1), diff))
end

"""
    FieldDiff

Store the differences between two fields:
1. The worst elementwise absolute or relative error,
2. The root mean square error,
3. Root mean square (RMS) magnitude of the reference (first) field.
"""
struct FieldDiff
    err::Float64
    rmse::Float64
    scale::Float64
end

function FieldDiff(x1::AbstractVector, x2::AbstractVector)
    # A length mismatch is a structural problem (wrong space/layout), flagged
    # as NaN so it always fails the tolerance test.
    length(x1) == length(x2) || return FieldDiff(NaN, NaN, NaN)
    return FieldDiff(mixed_error(x1, x2), rms(x1 .- x2), rms(x1))
end

"""
    mask_nonfinite(x)

Return zero if `x` is a non-finite value (e.g., `NaN` and `Inf`) and `x`
otherwise.
"""
mask_nonfinite(x) = x .* isfinite.(x)

"""
    flatten(field::ClimaCore.Fields.Field, col::Union{Nothing, Integer})

Flatten a `ClimaCore.Fields.Field` into a vector.

Pass `nothing` for `col` if the space is a single column and pass an integer
to get the `i`th `col`umn of `field`.
"""
# `parent` is used rather than `field2array` so multi-component fields (e.g.
# two-band albedo) don't error; copying to the host with `Array` keeps this
# GPU-compatible.
flatten(field::ClimaCore.Fields.Field, ::Nothing) =
    mask_nonfinite(vec(Array(parent(field))))
# `Fields.column` slices the underlying data by dispatch on its layout.
flatten(field::ClimaCore.Fields.Field, col::Integer) =
    flatten(ClimaCore.Fields.column(field, 1, 1, col), nothing)

"""
    field_diffs(
        v1,
        v2;
        col1 = nothing,
        col2 = nothing,
        name = "",
        ignore = Set{Symbol}(),
    )

Iterate through the fields of `v1` and `v2` while ignoring the fields in
`ignore` and return a dictionary of the `FieldDiff`s between the column `col1`
in `v1` and the column `col2` in `v2`.
"""
function field_diffs(
    v1,
    v2;
    col1 = nothing,
    col2 = nothing,
    name = "",
    ignore = Set{Symbol}(),
)
    diffs = Dict{String, FieldDiff}()
    _field_diffs!(diffs, v1, v2, name, col1, col2, ignore)
    return diffs
end

# Stop recursing when it is a ClimaCore.Fields.Field or a number
function _field_diffs!(
    diffs,
    f1::ClimaCore.Fields.Field,
    f2::ClimaCore.Fields.Field,
    name,
    col1,
    col2,
    ignore,
)
    diffs[name] = FieldDiff(flatten(f1, col1), flatten(f2, col2))
    return diffs
end

function _field_diffs!(diffs, x1::Number, x2::Number, name, col1, col2, ignore)
    diffs[name] =
        FieldDiff(mask_nonfinite([float(x1)]), mask_nonfinite([float(x2)]))
    return diffs
end

# If it is not a ClimaCore.Fields.Field or number, iterate through the fields of
# the struct
function _field_diffs!(diffs, v1, v2, name, col1, col2, ignore)
    props = filter(!in(ignore), propertynames(v1))
    props == filter(!in(ignore), propertynames(v2)) ||
        error("$name: mismatched structure ($(typeof(v1)) vs $(typeof(v2)))")
    for p in props
        _field_diffs!(
            diffs,
            getproperty(v1, p),
            getproperty(v2, p),
            "$name.$p",
            col1,
            col2,
            ignore,
        )
    end
    return diffs
end

"""
    diagnostic_diffs(writer1, writer2; col1 = nothing, col2 = nothing)

Compare every diagnostic stored in two `ClimaDiagnostics.Writers.DictWriter`s.
For each short name, all saved times are compared and the maximum error over the
times in `FieldDiff` is kept.
"""
function diagnostic_diffs(writer1, writer2; col1 = nothing, col2 = nothing)
    diffs = Dict{String, FieldDiff}()
    names = sort!(collect(keys(writer1)))
    names == sort!(collect(keys(writer2))) ||
        error("writers store different diagnostics")
    for diag_name in names
        times = sort!(collect(keys(writer1[diag_name])))
        times == sort!(collect(keys(writer2[diag_name]))) ||
            error("$diag_name: saved at different times")
        worst = FieldDiff(0.0, 0.0, 0.0)
        for t in times
            d = FieldDiff(
                flatten(writer1[diag_name][t], col1),
                flatten(writer2[diag_name][t], col2),
            )
            # `!(<=)` rather than `>` so a NaN (structural mismatch) wins.
            if !(d.err <= worst.err)
                worst = d
            end
        end
        diffs[diag_name] = worst
    end
    return diffs
end

"""
    report_diffs(diffs; label = "")

Print every entry of `diffs` sorted by path.
"""
function report_diffs(diffs; label = "")
    println("Field comparison [$label]")
    for (path, d) in sort(collect(diffs); by = first)
        @printf(
            "  %-55s err = %-12.3g rmse = %-12.3g scale = %-12.3g\n",
            path,
            d.err,
            d.rmse,
            d.scale
        )
    end
    return nothing
end
