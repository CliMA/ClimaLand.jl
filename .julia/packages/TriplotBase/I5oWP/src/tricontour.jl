const XY{T} = Tuple{T, T}
const Polyline{T} = Vector{XY{T}}
const Polylines{T} = Vector{Polyline{T}}
struct Contour{T}
    polylines::Polylines{T}
    level::T
end

Contour(level::T) where T = Contour(Vector{XY{T}}[], level)

"""
    tricontour(x, y, z, t, levels)

Return contour lines of unstructured triangular data. `x`, `y`, and `z` are one-dimensional
input arrays of the same length (the number of points). `t` is an integer array of size
`(3,nt)`, where `nt` is the number of triangles. The coordinates of the `i`th vertex of the
`j`th triangle are given by `(x[t[i,j]], y[t[i,j]])`. Function values at vertices are
specified by the `z` array, and contours levels are defined by `levels.`

If `levels` is an integer, then equally spaced levels between the extrema of `z` will be
generated. Otherwise, `levels` is a vector of increasing values, each value corresponding
to one contour line.
"""
function tricontour(x, y, z, t, levels)
    m = TriMesh(x, y, t)
    if levels isa Integer
        levels = range(minimum(z), maximum(z), length=levels)
    end
    tricontour(m, z, levels)
end

function tricontour(m::TriMesh, z, levels)
    @assert issorted(levels)
    contours = Contour{eltype(levels)}[]
    for level=levels
        push!(contours, generate_unfilled_contours(m, z, level))
    end
    contours
end

function generate_unfilled_contours(m::TriMesh, z, level)
    reset_visited!(m)
    contour = Contour(level)
    add_bdr_lines!(contour.polylines, m, z, level)
    add_interior_lines!(contour.polylines, m, z, level; on_upper=false, filled=false)
    contour
end

function interp(m::TriMesh, p1, p2, z, level)
    t = (z[p2] - level)/(z[p2] - z[p1])
    x = t*m.x[p1] + (1.0-t)*m.x[p2]
    y = t*m.y[p1] + (1.0-t)*m.y[p2]
    x,y
end

function interp(m::TriMesh, edge::TriEdge, z, level)
    p1 = m.t[edge.ie,edge.it]
    p2 = m.t[mod1(edge.ie+1,3),edge.it]
    interp(m, p1, p2, z, level)
end

const exit_edges = [0,3,1,3,2,2,1,0] # Lookup table for exit edges using encoded config

function get_exit_edge(m::TriMesh, it, z, level; on_upper)
    config = (z[m.t[1,it]] >= level) |
             (z[m.t[2,it]] >= level) << 1 |
             (z[m.t[3,it]] >= level) << 2
    on_upper && (config = 7-config)
    exit_edges[config+1]
end

function follow_interior!(polyline::Polyline, m::TriMesh, edge::TriEdge, z, level; end_on_bdr, on_upper=false)
    push!(polyline, interp(m, edge, z, level))
    it = edge.it
    while true
        vi = on_upper ? it + size(m.t,2) : it
        !end_on_bdr && m.visited[vi] && return edge
        m.visited[vi] = true
        ie = get_exit_edge(m, it, z, level; on_upper)
        @assert ie > 0
        edge = TriEdge(ie,it)
        push!(polyline, interp(m, edge, z, level))
        # Move to next triangle
        next_edge = get_neighbor_edge(m,edge)
        it = next_edge.it
        end_on_bdr && it == 0 && return edge
        edge = next_edge
    end
end

function follow_interior(m::TriMesh, edge::TriEdge, z, level::T; end_on_bdr, on_upper=false) where T
    polyline = Polyline{T}()
    edge = follow_interior!(polyline, m, edge, z, level; end_on_bdr, on_upper)
    polyline, edge
end

function add_bdr_lines!(polylines::Polylines, m::TriMesh, z, level)
    for bdr=m.bdr_loops
        nedges = length(bdr)
        end_above = false
        for i=1:nedges
            edge = bdr[i]
            if i == 1
                start_above = z[m.t[edge.ie,edge.it]] >= level
            else
                start_above = end_above
            end
            end_above = z[m.t[mod1(edge.ie+1,3),edge.it]] >= level
            if start_above && !end_above
                polyline,_ = follow_interior(m, edge, z, level; end_on_bdr=true)
                push!(polylines, polyline)
            end
        end
    end
end

function add_interior_lines!(polylines::Polylines, m::TriMesh, z, level; on_upper, filled)
    nt = size(m.t,2)
    for it=1:nt
        vi = on_upper ? it + size(m.t,2) : it
        m.visited[vi] && continue
        m.visited[vi] = true
        ie = get_exit_edge(m, it, z, level; on_upper)
        ie == 0 && continue
        # Found a new outgoing contour line, start from neighbor
        edge = get_neighbor_edge(m, TriEdge(ie,it))
        polyline,_ = follow_interior(m, edge, z, level; on_upper, end_on_bdr=false)
        push!(polyline, first(polyline)) # Make closed loop
        push!(polylines, polyline)
    end
end
