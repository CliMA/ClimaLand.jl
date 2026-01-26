struct FilledContour{T}
    polylines::Polylines{T} # Polylines{T} defined in tricontour.jl
    lower::T
    upper::T
end
FilledContour(lower::T, upper::T) where T = FilledContour(Polylines{T}(), lower, upper)

"""
    tricontourf(x, y, z, t, levels)

Return filled contour lines of unstructured triangular data. `x`, `y`, and `z` are
one-dimensional input arrays of the same length (the number of points). `t` is an integer
array of size `(3,nt)`, where `nt` is the number of triangles. The coordinates of the `i`th
vertex of the `j`th triangle are given by `(x[t[i,j]], y[t[i,j]])`. Function values at
vertices are specified by the `z` array, and contours levels are defined by `levels.`

`tricontourf` results in regions *between* the contour levels, such that specifying `n`
levels will result in `n-1` regions. If `levels` is an array, it will be augmented to
include the minimum and maximum values of the `z` data.
"""
function tricontourf(x, y, z, t, levels)
    m = TriMesh(x, y, t)
    minz,maxz = extrema(z)
    if levels isa Integer
        levels = range(minz-eps(minz), maxz+eps(maxz), length=levels)
    else
        levels = copy(levels)
        minimum(levels) >= minz && pushfirst!(levels, minz-eps(minz))
        maximum(levels) <= maxz && push!(levels, maxz+eps(maxz))
    end
    tricontourf(m, z, levels)
end

function tricontourf(m::TriMesh, z, levels)
    @assert issorted(levels)
    nlevels = length(levels)
    lmin = minimum(levels)
    lmax = (nlevels > 1) ? maximum(levels[1:end-1]) : lmin
    filled_contours = FilledContour{eltype(levels)}[]
    for i=1:nlevels-1
        lower = levels[i]
        upper = levels[i+1]
        push!(filled_contours, generate_filled_contours(m, z, lower, upper))
    end
    filled_contours
end

function generate_filled_contours(m::TriMesh, z, lower, upper)
    reset_visited!(m)
    filled_contour = FilledContour(lower, upper)
    add_filled_bdr_lines!(filled_contour.polylines, m, z, lower, upper)
    add_interior_lines!(filled_contour.polylines, m, z, lower; on_upper=false, filled=true)
    add_interior_lines!(filled_contour.polylines, m, z, upper; on_upper=true,  filled=true)
    filled_contour
end

function follow_boundary!(polyline::Polyline, m::TriMesh, edge::TriEdge, z, lower, upper, on_upper)
    ib,ie = m.bdr_map[edge]
    m.bdr_used[ib] = true

    stop = false
    first_edge = true
    z_end = 0.0
    while !stop
        @assert !m.bdr_visited[ib][ie]
        m.bdr_visited[ib][ie] = true
        if first_edge
            z_start = z[m.t[edge.ie,edge.it]]
        else
            z_start = z_end
        end
        z_end = z[m.t[mod1(edge.ie+1,3),edge.it]]

        if z_end > z_start
            if (on_upper || !first_edge) && (z_end >= lower) && (z_start < lower)
                stop = true
                on_upper = false
            elseif (z_end >= upper) && (z_start < upper)
                stop = true
                on_upper = true
            end
        else
            if (!on_upper || !first_edge) && (z_start >= upper) && (z_end < upper)
                stop = true
                on_upper = true
            elseif (z_start >= lower) && (z_end < lower)
                stop = true
                on_upper = false
            end
        end
        first_edge = false
        if !stop
            ie = mod1(ie+1, length(m.bdr_loops[ib]))
            edge = m.bdr_loops[ib][ie]
            point = m.t[edge.ie,edge.it]
            push!(polyline, (m.x[point], m.y[point]))
        end
    end
    return on_upper,edge
end

function add_filled_bdr_lines!(polylines::Polylines{T}, m::TriMesh, z, lower, upper) where T
    for i=1:length(m.bdr_loops)
        bdr = m.bdr_loops[i]
        for j=1:length(bdr)
            m.bdr_visited[i][j] && continue
            edge = bdr[j]

            z_start = z[m.t[edge.ie,edge.it]]
            z_end = z[m.t[mod1(edge.ie+1,3),edge.it]]

            # Does this boundary edge's z increase through upper level and/or decrease
            # through lower level?
            incr_upper = (z_start < upper && z_end >= upper)
            decr_lower = (z_start >= lower && z_end < lower)

            if decr_lower || incr_upper
                start_edge = edge
                polyline = XY{T}[]
                on_upper = incr_upper
                while true
                    level_int = on_upper ? upper : lower
                    edge = follow_interior!(polyline, m, edge, z, level_int; end_on_bdr=true, on_upper)
                    on_upper,edge = follow_boundary!(polyline, m, edge, z, lower, upper, on_upper)
                    edge == start_edge && break
                end
                # Ensure closed loop
                push!(polyline, first(polyline))
                push!(polylines, polyline)
            end
        end
    end

    # Add full boundaries that lie between the lower and upper levels.  These are boundaries
    # that have not been touched by an internal contour line.
    for i=1:length(m.bdr_loops)
        m.bdr_used[i] && continue
        bdr = m.bdr_loops[i]
        zval = z[m.t[bdr[1].ie,bdr[1].it]]
        if (zval >= lower) && (zval < upper)
            polyline = XY{T}[]
            for j=1:length(bdr)
                edge = bdr[j]
                point = m.t[edge.ie,edge.it]
                push!(polyline, (m.x[point], m.y[point]))
            end
            push!(polyline, first(polyline))
            push!(polylines, polyline)
        end
    end
end
