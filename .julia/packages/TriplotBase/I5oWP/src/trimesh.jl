struct TriEdge
    ie::Int # Edge index in range [1,3]
    it::Int # Triangle index in range [1,nt]
end

struct TriMesh{X,Y,T}
    x::X # x-coordinates of the vertices, vector of length nv
    y::Y # y-coordinates of the verties, vector of length nv
    t::T # Vertex indices of the triangles, array of size (3,nt)
    neighbors::Array{Int,2} # Adjacency list of the mesh, array of size (3,nt)
    # Note: neighbors[ie,it] is index of triangle neighboring edge `ie` of triangle `it`
    bdr_loops::Vector{Vector{TriEdge}}
    bdr_map::Dict{TriEdge,Tuple{Int,Int}} # Map from edge index to boundary indices
    visited::Vector{Bool} # Keep track of visited elements for internal calculations
    bdr_visited::Vector{Vector{Bool}} # Keep track of visited boundary edges
    bdr_used::Vector{Bool} # Keep track of which boundary loops are used for filled contours
end

function TriMesh(x, y, t)
    nv = length(x)
    nt = size(t,2)

    # Compute adjacency of the mesh
    neighbors = zeros(Int,3,nt)
    edgemap = Dict{Tuple{Int,Int},TriEdge}()
    for it=1:nt
        for ie=1:3
            edge = (t[ie,it],t[mod1(ie+1,3),it])
            # Assume mesh is properly oriented, i.e. edges of neighboring triangles are
            # numbered in reversed order
            r_edge = (edge[2], edge[1])
            if haskey(edgemap,r_edge)
                tri_edge = edgemap[r_edge]
                neighbors[ie,it] = tri_edge.it
                neighbors[tri_edge.ie,tri_edge.it] = it
                delete!(edgemap,r_edge)
            else
                edgemap[edge] = TriEdge(ie,it)
            end
        end
    end
    empty!(edgemap)

    # Find all boundary edges, i.e. those with no neighboring triangles
    bdr_i = findall(neighbors .== 0)
    bdr_edges = Set{TriEdge}()
    for i=bdr_i
        push!(bdr_edges, TriEdge(i[1], i[2]))
    end

    # Follow each boundary edge to create boundary loops
    bdr_map = Dict{TriEdge,Tuple{Int,Int}}()
    bdr_loops = Vector{Vector{TriEdge}}()
    while !isempty(bdr_edges)
        bdr = Vector{TriEdge}()
        push!(bdr_loops, bdr)
        first_edge = first(bdr_edges.dict)[1] # same as first(bdr_edges), slightly faster
        edge = first_edge
        while true
            push!(bdr, edge)
            bdr_map[edge] = (length(bdr_loops), length(bdr))
            delete!(bdr_edges, edge)
            it = edge.it
            ie = mod1(edge.ie+1,3)
            point = t[ie,it]
            # Traverse neighbors until we find one with a boundary edge
            while neighbors[ie,it] != 0
                it = neighbors[ie,it]
                ie = findfirst(t[:,it] .== point)
                @assert !isnothing(ie)
            end
            edge = TriEdge(ie,it)
            edge == first_edge && break # Returned to first edge, completed loop
        end
    end

    visited = fill(false,2*nt)
    bdr_visited = Vector{Vector{Bool}}()
    for bdr=bdr_loops
        push!(bdr_visited, fill(false, length(bdr)))
    end
    bdr_used = fill(false, length(bdr_loops))

    TriMesh(x,y,t,neighbors,bdr_loops,bdr_map,visited,bdr_visited,bdr_used)
end

function get_neighbor_edge(m::TriMesh, edge::TriEdge)
    it = m.neighbors[edge.ie,edge.it]
    if it == 0
        return TriEdge(0,0)
    else
        point = m.t[mod1(edge.ie+1,3),edge.it]
        ie = findfirst(m.t[:,it] .== point)
        return TriEdge(ie,it)
    end
end

function reset_visited!(m::TriMesh)
    m.visited .= false
    for bv=m.bdr_visited
        bv .= false
    end
    m.bdr_used .= false
end
