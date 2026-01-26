module Isoband

using isoband_jll


export isobands, isolines


struct ReturnValue
    xs::Ptr{Cdouble}
    ys::Ptr{Cdouble}
    ids::Ptr{Cint}
    len::Cint
end

"""
    isobands(xs, ys, zs, low::Real, high::Real)

Create one isoband from the matrix `zs` for the boundaries `low` and `high`.
The rows of `zs` correspond to the linear spaced values in `ys` and the columns to `xs`.
Returns a NamedTuple with two Vector{Float64} fields `x` and `y`, and the
Vector{Int} field `id`. Each unique id marks one polygon, the polygons can be outer polygons or
holes and are given in no particular order. Therefore, they must probably be post-processed to
feed them to plotting packages.
"""
function isobands(xs, ys, zs, low::Real, high::Real)
    results = isobands(xs, ys, zs, Float64[low], Float64[high])
    results[1]
end

function isobands(xs::AbstractVector, ys::AbstractVector, zs::AbstractMatrix, lows::AbstractVector, highs::AbstractVector)
    isobands(Float64.(xs), Float64.(ys), Float64.(zs), Float64.(lows), Float64.(highs))
end

"""
    isobands(xs::Vector{Float64}, ys::Vector{Float64}, zs::Matrix{Float64}, low_values::Vector{Float64}, high_values::Vector{Float64})

Create a vector of isobands from the matrix `zs` for all pairs in `low_values` and `high_values`.
The rows of `zs` correspond to the linear spaced values in `ys` and the columns to `xs`.
Each entry of the return vector is a NamedTuple with two Vector{Float64} fields `x` and `y`, and the
Vector{Int} field `id`. Each unique id marks one polygon, the polygons can be outer polygons or
holes and are given in no particular order. Therefore, they must probably be post-processed to
feed them to plotting packages.
"""
function isobands(
        xs::Vector{Float64},
        ys::Vector{Float64},
        zs::Matrix{Float64},
        low_values::Vector{Float64},
        high_values::Vector{Float64})

    lenx = length(xs)
    leny = length(ys)
    nrow, ncol = size(zs)

    lenx != ncol && error("Length of x ($(length(xs))) must be equal to number of columns in z ($(size(zs, 2)))")
    leny != nrow && error("Length of y $(length(ys)) must be equal to number of rows in z ($(size(zs, 1))")

    length(low_values) != length(high_values) && error("Number of low values ($(length(low_values)))and high values ($(length(high_values))) must be equal.")

    nbands = length(low_values)

    result = ccall((:isobands_impl, libisoband),
        Ptr{ReturnValue},
        (Ptr{Cdouble},
            Cint,
            Ptr{Cdouble},
            Cint,
            Ptr{Cdouble},
            Cint,
            Cint,
            Ptr{Cdouble},
            Ptr{Cdouble},
            Cint),
        xs, length(xs),
        ys, length(ys),
        zs, size(zs, 1), size(zs, 2),
        low_values, high_values, nbands)
    

    returnvalues = unsafe_wrap(Vector{ReturnValue}, result, nbands, own = true)

    groups = map(returnvalues) do rv
        n = rv.len
        xsr = unsafe_wrap(Vector{Cdouble}, rv.xs, n, own = true)
        ysr = unsafe_wrap(Vector{Cdouble}, rv.ys, n, own = true)
        idr = unsafe_wrap(Vector{Cint}, rv.ids, n, own = true)
        (x = xsr, y = ysr, id = idr)
    end

    groups
end

"""
    isolines(xs, ys, zs, value::Real)

Create one isoline from the matrix `zs` for `value`.
The rows of `zs` correspond to the linear spaced values in `ys` and the columns to `xs`.
Returns a NamedTuple with two Vector{Float64} fields `x` and `y`, and the
Vector{Int} field `id`. Each unique id marks one line.
"""
function isolines(xs, ys, zs, value::Real)
    results = isolines(xs, ys, zs, Float64[value])
    results[1]
end

function isolines(xs::AbstractVector, ys::AbstractVector, zs::AbstractMatrix, values::AbstractVector)
    isolines(Float64.(xs), Float64.(ys), Float64.(zs), Float64.(values))
end

"""
    isolines(xs::Vector{Float64}, ys::Vector{Float64}, zs::Matrix{Float64}, values::Vector{Float64})

Create a vector of isolines from the matrix `zs` for all `values`.
The rows of `zs` correspond to the linear spaced values in `ys` and the columns to `xs`.
Each entry of the return vector is a NamedTuple with two Vector{Float64} fields `x` and `y`, and the
Vector{Int} field `id`. Each unique id marks one line.
"""
function isolines(
        xs::Vector{Float64},
        ys::Vector{Float64},
        zs::Matrix{Float64},
        values::Vector{Float64})

    lenx = length(xs)
    leny = length(ys)
    nrow, ncol = size(zs)

    lenx != ncol && error("Length of x ($(length(xs))) must be equal to number of columns in z ($(size(zs, 2)))")
    leny != nrow && error("Length of y $(length(ys)) must be equal to number of rows in z ($(size(zs, 1))")

    nvalues = length(values)

    result = ccall((:isolines_impl, libisoband),
        Ptr{ReturnValue},
        (Ptr{Cdouble},
            Cint,
            Ptr{Cdouble},
            Cint,
            Ptr{Cdouble},
            Cint,
            Cint,
            Ptr{Cdouble},
            Cint),
        xs, length(xs),
        ys, length(ys),
        zs, size(zs, 1), size(zs, 2),
        values, nvalues)
    

    returnvalues = unsafe_wrap(Vector{ReturnValue}, result, nvalues, own = true)

    groups = map(returnvalues) do rv
        n = rv.len
        xsr = unsafe_wrap(Vector{Cdouble}, rv.xs, n, own = true)
        ysr = unsafe_wrap(Vector{Cdouble}, rv.ys, n, own = true)
        idr = unsafe_wrap(Vector{Cint}, rv.ids, n, own = true)
        (x = xsr, y = ysr, id = idr)
    end
    
    groups
end


end
