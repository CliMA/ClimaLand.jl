# ported from https://github.com/joelverhagen/2D-Bin-Packing/blob/master/BinPacking/PackingTreeNode.h

mutable struct RectanglePacker{T}
    # the rectangle represented by this node
    area::Rect{2, T}
    filled::Bool
    # a vector of child nodes
    left::Union{RectanglePacker{T}, Nothing}
    right::Union{RectanglePacker{T}, Nothing}
end

function RectanglePacker(area::Rect{2, T}) where {T}
    return RectanglePacker{T}(area, false, nothing, nothing)
end

function Base.push!(node::RectanglePacker{T}, areas::Vector{Rect2{T}}) where T
    sort!(areas, by=GeometryBasics.norm ∘ widths)
    return RectanglePacker{T}[push!(node, area) for area in areas]
end

function Base.push!(node::RectanglePacker{T}, new_rect::Rect) where {T}
    nwidth, nheight = widths(new_rect)
    # if the node is not a leaf (has children)
    if !isnothing(node.left) || !isnothing(node.right)
        # try inserting into the first child
        new_node = push!(node.left, new_rect)
        !isnothing(new_node) && return new_node
        # no room, insert into the second child
        return push!(node.right, new_rect)
    else

        # if there's already an image here, return
        node.filled && return nothing
        area = node.area
        awidth, aheight = widths(area)
        # if the image doesn't fit, return
        if nwidth > awidth || nheight > aheight
            return nothing
        end

        # if the image fits perfectly, accept
        if nwidth == awidth && nheight == aheight
            node.filled = true
            return node
        end

        # otherwise, split the node and create two children

        # decide which way to split

        left, bottom = minimum(area)
        right, top = maximum(area)

        dw = awidth - nwidth
        dh = aheight - nheight

        if dw > dh
            left_rectangle = Rect2{T}(left, bottom, nwidth, aheight)
            right_rectangle = Rect2{T}(left + nwidth, bottom, awidth - nwidth, aheight)
        else
            left_rectangle = Rect2{T}(left, bottom, awidth, nheight)
            right_rectangle = Rect2{T}(left, bottom + nheight, awidth, aheight - nheight)
        end
        node.left = RectanglePacker(left_rectangle)
        node.right = RectanglePacker(right_rectangle)
        # insert into the first child we created
        return push!(node.left, new_rect)
    end
end
