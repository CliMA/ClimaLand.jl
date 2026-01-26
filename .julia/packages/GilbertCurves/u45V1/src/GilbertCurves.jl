module GilbertCurves

export gilbertindices

"""
    gilbertindices(dims::Tuple{Int,Int}; majdim=dims[1] >= dims[2] ? 1 : 2)

Constructs a vector of `CartesianIndex` objects, orderd by a generalized Hilbert
space-filling, starting at `CartesianIndex(1,1)`.  It will end at
`CartesianIndex(dims[1],1)` if `majdim==1`, or `CartesianIndex(1,dims[2])` if `majdim==2`
(or the closest feasible point).
"""
gilbertindices(dims::Tuple{Int,Int}; kwargs...) =
    gilbertorder(CartesianIndices(dims); kwargs...)


"""
    gilbertorder(mat::AbstractMatrix; majdim=size(mat,1) >= size(mat,2) ? 1 : 2)

Constructs a vector of the elements of `mat`, ordered by a generalized Hilbert
space-filling curve. The list will start at `mat[1,1]`, and end at `mat[end,1]` if
`majdim==1` or `mat[1,end]` if `majdim==2`  (or the closest feasible point).
"""
function gilbertorder(mat::AbstractMatrix{T}; majdim=size(mat,1) >= size(mat,2) ? 1 : 2) where {T}
    list = sizehint!(T[], length(mat))
    if majdim == 1
        append_gilbert!(list, mat)
    else
        append_gilbert!(list, permutedims(mat,(2,1)))
    end
end

function append_gilbert!(list, mat::AbstractMatrix)
    # 1 |*    |
    #   | )   |
    # a |v    |
    a,b = size(mat)
    if a == 1 || b == 1
        # single in one dimension
        append!(list, mat)
    elseif 2a > 3b
        # long case: split into two
        #   +-----+
        # 1 |*    |
        #   ||    |
        # a2|v    |
        #   +-----+
        #   |*    |
        #   ||    |
        # a |v    |
        #   +-----+
        a2 = div(a,2)
        if isodd(a2) && a > 2
            a2 += 1
        end
        append_gilbert!(list, mat[1:a2,:])
        append_gilbert!(list, mat[a2+1:a,:])
    else
        # standard case: split into three
        #      b2
        #   +---+---+
        # 1 |*->|*   |
        #   |   ||   |
        # a2|   ||   |
        #   +---+|   |
        #   |   ||   |
        # a |<-*|v   |
        #   +---+----+
        a2 = div(a,2)
        b2 = div(b,2)
        if isodd(b2) && b > 2
            b2 += 1
        end
        append_gilbert!(list, permutedims(mat[1:a2,1:b2],(2,1)))
        append_gilbert!(list, mat[:,b2+1:b])
        append_gilbert!(list, permutedims(mat[a:-1:a2+1,b2:-1:1],(2,1)))
    end
end

"""
    linearindices(list::Vector{CartesianIndex{2}})

Construct an integer matrix `M` containing the integers `1:length(list)` such that
`M[list[i]] == i`.
"""
function linearindices(list::Vector{CartesianIndex{2}})
    cmax = maximum(list)
    L = zeros(Int,cmax.I)
    for (i,c) in enumerate(list)
        L[c] = i
    end
    return L
end

end
