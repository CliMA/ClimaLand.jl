module Contour

include("interpolate.jl")

export
    ContourLevel,
    Curve2,
    contour,
    contours,
    level,
    levels,
    lines,
    coordinates,
    vertices

import Base: push!, length, eltype, show

"""
    Curve2{T}

Containing an extracted chain of points, where `T` is a vector-like element.
"""
struct Curve2{T}
    vertices::Vector{T}
end
Curve2(::Type{T}) where {T} = Curve2(T[])
show(io::IO, ::MIME"text/plain", c2::Curve2) = write(io, "$(typeof(c2))\n  with $(length(c2.vertices)-1) vertices")
show(io::IO, ::MIME"text/plain", c2s::Vector{Curve2{T}}) where {T} = write(io, "$(typeof(c2s))\n  $(length(c2s)) contour line(s)")

struct ContourLevel{T, L}
    level::L
    lines::Vector{Curve2{T}}
end
ContourLevel(h::T) where {T <: AbstractFloat} = ContourLevel(h, Curve2{NTuple{2,T}}[])
ContourLevel(h::T) where {T} = ContourLevel(Float64(h))
show(io::IO, ::MIME"text/plain", cl::ContourLevel) = write(io, "$(typeof(cl))\n  at $(level(cl)) with $(length(lines(cl))) line(s)")
show(io::IO, ::MIME"text/plain", cls::Vector{ContourLevel{T}}) where {T} = write(io, "$(typeof(cls))\n  $(length(cls)) contour level(s)")
"""
`lines(c)` Extracts an iterable collection of isolines from a contour level.
Use [`coordinates`](@ref) or [`vertices`](@ref) to get the coordinates of a line.
"""
lines(cl::ContourLevel) = cl.lines
"""
`level(c)` Indicates the `z`-value at which the contour level `c` was traced.
"""
level(cl::ContourLevel) = cl.level

struct ContourCollection{Tlevel<:ContourLevel}
    contours::Vector{Tlevel}
end
ContourCollection() = ContourCollection(Float64)
ContourCollection(::Type{Tlevel}) where {Tlevel} = ContourCollection(ContourLevel{Tlevel}[])
show(io::IO, ::MIME"text/plain", cc::ContourCollection) = write(io, "$(typeof(cc))\n with $(length(levels(cc))) level(s).")

"""
Turns the output of [`contours`](@ref) into an iterable with each of the traced
contour levels. Each of the objects support [`level`](@ref),
[`coordinates`](@ref), and [`vertices`](@ref).
"""
levels(cc::ContourCollection) = cc.contours

"""
    contour(x, y, z, level::Number)
    contour(x, y, z, level::Number, VT::Type)

Trace a single contour level, indicated by the argument `level`. The extracted vertex type
maybe be specified by `VT`.

You'll usually call [`lines`](@ref) on the output of `contour`, and then iterate
over the result.
"""
function contour(x, y, z, level::Number; VT=nothing)
    if !(axes(x) == (axes(z,1),) && axes(y) == (axes(z,2),) || axes(x) == axes(y) == axes(z))
        throw(ArgumentError("Incompatible input axes in `Contour.contour`."))
    end
    ET = promote_type(map(eltype, (x, y, z))...)
    ET = ET <: Integer ? Float64 : ET
    VT = VT === nothing ? NTuple{2,ET} : VT
    trace_contour(x, y, z, level, get_level_cells(z, level), VT)
end

"""
`contours` returns a set of isolines.

You'll usually call [`levels`](@ref) on the output of `contours`.
"""
function contours end

"""
    contours(x,y,z,levels;[VT])
Trace the contour levels indicated by the `levels` argument.
The extracted vertex type maybe be specified by the `VT` keyword.
"""
contours(x, y, z, levels; VT=nothing) = ContourCollection([contour(x, y, z, l; VT=VT) for l in levels])

"""
    contours(x,y,z,Nlevels::Integer;[VT])

Trace `Nlevels` contour levels at heights
chosen by the library (using the  [`contourlevels`](@ref) function).
The extracted vertex type maybe be specified by the `VT` keyword.
"""
function contours(x, y, z, Nlevels::Integer;VT=nothing)
    contours(x, y, z, contourlevels(z, Nlevels); VT=VT)
end

"""
`contours(x,y,z;[VT])` Trace 10 automatically chosen contour levels.
"""
contours(x, y, z; VT=nothing) = contours(x, y, z, 10; VT=VT)

"""
`contourlevels(z,n)` Examines the values of `z` and chooses `n` evenly spaced
levels to trace.
"""
function contourlevels(z, n)
    zmin, zmax = extrema(z)
    dz = (zmax - zmin) / (n + 1)
    range(zmin + dz; step = dz, length = n)
end

"""
`coordinates(c)` Returns the coordinates of the vertices of the contour line as
a tuple of lists.
"""
function coordinates(c::Curve2{T}) where {T}
    N = length(c.vertices)
    E = eltype(T)
    xlist = Vector{E}(undef, N)
    ylist = Vector{E}(undef, N)

    for (i, v) in enumerate(c.vertices)
        xlist[i] = v[1]
        ylist[i] = v[2]
    end
    xlist, ylist
end

"""
`vertices(c)`

Returns the vertices of a contour line as a vector of 2-element tuples.
"""
vertices(c::Curve2) = c.vertices

# The marching squares algorithm defines 16 cell types
# based on the edges that a contour line enters and exits
# through. The edges of the cells are identified using
# compass directions, while the vertices are ordered as
# follows:
#
#      N
#  4 +---+ 3
# W  |   |  E
#  1 +---+ 2
#      S
#
# Each cell type is identified with 4 bits, with each
# bit corresponding to a vertex (MSB -> 4, LSB -> 1).
# A bit is set for vertex v_i is set if z(v_i) > h. So a cell
# where a contour line only enters from the W edge and exits
# through the N edge will have the cell type: 0b0111
# Note that there are two cases where there are two
# lines crossing through the same cell: 0b0101, 0b1010.
const N, S, E, W = (UInt8(1)), (UInt8(2)), (UInt8(4)), (UInt8(8))
const NS, NE, NW = N|S, N|E, N|W
const SN, SE, SW = S|N, S|E, S|W
const EN, ES, EW = E|N, E|S, E|W
const WN, WS, WE = W|N, W|S, W|E
const NWSE = NW | 0x10 # special (ambiguous case)
const NESW = NE | 0x10 # special (ambiguous case)

# Maps cell type to crossing types for non-ambiguous cells
const edge_LUT = (SW, SE, EW, NE, 0x0, NS, NW, NW, NS, 0x0, NE, EW, SE, SW)

# The way a contour crossing goes through a cell is labeled
# by combining compass directions (e.g. a NW crossing connects
# the N edge and W edges of the cell).  The Cell type records
# the type of crossing that a cell contains.  While most
# cells will have only one crossing, cell type 5 and 10 will
# have two crossings.
function get_next_edge!(cells::Dict, key, entry_edge::UInt8)
    cell = pop!(cells, key)
    if cell == NWSE
        if entry_edge == N || entry_edge == W
            cells[key] = SE
            cell = NW
        else #SE
            cells[key] = NW
            cell = SE
        end
    elseif cell == NESW
        if entry_edge == N || entry_edge == E
            cells[key] = SW
            cell = NE
        else #SW
            cells[key] = NE
            cell = SW
        end
    end
    return cell ⊻ entry_edge
end

# N, S, E, W
const next_map = ((0,1), (0,-1), (1,0), (-1,0))
const next_edge = (S,N,W,E)

@inline function advance_edge(ind, edge)
    n = trailing_zeros(edge) + 1
    nt = ind .+ next_map[n]
    return nt, next_edge[n]
end

@inline function get_first_crossing(cell)
    if cell == NWSE
        return NW
    elseif cell == NESW
        return NE
    else
        return cell
    end
end

@inline function _get_case(z, h)
    case = z[1] > h ? 0x01 : 0x00
    z[2] > h && (case |= 0x02)
    z[3] > h && (case |= 0x04)
    z[4] > h && (case |= 0x08)
    case
end

function get_level_cells(z, h::Number)
    cells = Dict{Tuple{Int,Int},UInt8}()
    x_ax, y_ax = axes(z)

    @inbounds for xi in first(x_ax):last(x_ax)-1
        for yi in first(y_ax):last(y_ax)-1
            elts = (z[xi, yi], z[xi + 1, yi], z[xi + 1, yi + 1], z[xi, yi + 1])
            case = _get_case(elts, h)

            # Contour does not go through these cells
            if iszero(case) || case == 0x0f
                continue
            end

            # Process ambiguous cells (case 5 and 10) using
            # a bilinear interpolation of the cell-center value.
            if case == 0x05
                cells[(xi, yi)] = 0.25*sum(elts) >= h ? NWSE : NESW
            elseif case == 0x0a
                cells[(xi, yi)] = 0.25*sum(elts) >= h ? NESW : NWSE
            else
                cells[(xi, yi)] = edge_LUT[case]
            end
        end
    end

    return cells
end


# Given a cell and a starting edge, we follow the contour line until we either
# hit the boundary of the input data, or we form a closed contour.
function chase!(cells, curve, x, y, z, h, start, entry_edge, xi_range, yi_range, ::Type{VT}) where VT

    ind = start

    # When the contour loops back to the starting cell, it is possible
    # for it to not intersect with itself.  This happens if the starting
    # cell contains a saddle-point. So a loop is only closed if the
    # contour returns to the starting edge of the starting cell
    loopback_edge = entry_edge

    @inbounds while true
        exit_edge = get_next_edge!(cells, ind, entry_edge)

        push!(curve, interpolate(x, y, z, h, ind, exit_edge, VT))

        ind, entry_edge = advance_edge(ind, exit_edge)

        !((ind[1], ind[2], entry_edge) != (start[1], start[2], loopback_edge) &&
           ind[2] ∈ yi_range && ind[1] ∈ xi_range) && break
    end

    return ind
end


function trace_contour(x, y, z, h::Number, cells::Dict, VT)

    contours = ContourLevel(h, Curve2{VT}[])

    x_ax, y_ax = axes(z)
    xi_range = first(x_ax):last(x_ax)-1
    yi_range = first(y_ax):last(y_ax)-1


    # When tracing out contours, this algorithm picks an arbitrary
    # starting cell, then first follows the contour in one direction
    # until it either ends up where it started # or at one of the boundaries.
    # It then tries to trace the contour in the opposite direction.

    @inbounds while length(cells) > 0
        contour_arr = VT[]

        # Pick initial box
        ind, cell = first(cells)

        # Pick a starting edge
        crossing = get_first_crossing(cell)
        starting_edge = 0x01 << trailing_zeros(crossing)

        # Add the contour entry location for cell (xi_0,yi_0)
        push!(contour_arr, interpolate(x, y, z, h, ind, starting_edge, VT))

        # Start trace in forward direction
        ind_end = chase!(cells, contour_arr, x, y, z, h, ind, starting_edge, xi_range, yi_range, VT)

        if ind == ind_end
            push!(contours.lines, Curve2(contour_arr))
            continue
        end

        ind, starting_edge = advance_edge(ind, starting_edge)

        if ind[2] ∈ yi_range && ind[1] ∈ xi_range
            # Start trace in reverse direction
            chase!(cells, reverse!(contour_arr), x, y, z, h, ind, starting_edge, xi_range, yi_range, VT)
        end

        push!(contours.lines, Curve2(contour_arr))
    end

    return contours

end

end
