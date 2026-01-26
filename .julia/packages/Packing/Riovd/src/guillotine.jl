# Ported from: https://github.com/juj/RectangleBinPack/blob/master/GuillotineBinPack.cpp

@enum FreeRectChoiceHeuristic begin
    RectBestAreaFit #< -BAF
    RectBestShortSideFit #< -BSSF
    RectBestLongSideFit #< -BLSF
    RectWorstAreaFit #< -WAF
    RectWorstShortSideFit #< -WSSF
    RectWorstLongSideFit #< -WLSF
end

@enum GuillotineSplitHeuristic begin
    SplitShorterLeftoverAxis #< -SLAS
    SplitLongerLeftoverAxis #< -LLAS
    SplitMinimizeArea #< -MINAS, Try to make a single big rectangle at the expense of making the other small.
    SplitMaximizeArea #< -MAXAS, Try to make both remaining rectangles as even-sized as possible.
    SplitShorterAxis #< -SAS
    SplitLongerAxis #< -LAS
end

struct GuillotinePacker{T}
    bin_width::Int
    bin_height::Int
    merge::Bool
    rect_choice::FreeRectChoiceHeuristic
    split_method::GuillotineSplitHeuristic
    free_rectangles::Vector{Rect2{T}}
    used_rectangles::Vector{Rect2{T}}

end

function GuillotinePacker(width::T, height::T; merge::Bool=true,
                rect_choice::FreeRectChoiceHeuristic=RectBestAreaFit,
                split_method::GuillotineSplitHeuristic=SplitMinimizeArea) where T

    packer = GuillotinePacker(width, height, merge, rect_choice, split_method,
                              Rect2{T}[], Rect2{T}[])
    push!(packer.free_rectangles, Rect2(T(0), T(0), width, height))
    return packer
end

function Base.push!(packer::GuillotinePacker, rects::AbstractVector{Rect2{T}}) where T
    rects = copy(rects)
    free_rectangles = packer.free_rectangles
    used_rectangles = packer.used_rectangles
    rect_choice = packer.rect_choice
    split_method = packer.split_method
    # Remember variables about the best packing choice we have made so far during the iteration process.
    best_free_rect = 0
    best_rect = 0
    best_flipped = false
    # Pack rectangles one at a time until we have cleared the rects array of all rectangles.
    # rects will get destroyed in the process.
    while length(rects) > 0
        # Stores the penalty score of the best rectangle placement - bigger=worse, smaller=better.
        best_score = typemax(Int)

        for i in 1:length(free_rectangles)
            force_jump = false
            rect_i = free_rectangles[i]
            for j in 1:length(rects)
                # If this rectangle is a perfect match, we pick it instantly.
                if all(widths(rects[j]) .== widths(rect_i))
                    best_free_rect = i
                    best_rect = j
                    best_flipped = false
                    best_score = typemin(Int)
                    force_jump = true # Force a jump out of the outer loop as well - we got an instant fit.
                    break
                # If flipping this rectangle is a perfect match, pick that then.
                elseif (height(rects[j]) == width(rect_i) && width(rects[j]) == height(rect_i))
                    best_free_rect = i
                    best_rect = j
                    best_flipped = true
                    best_score = typemin(Int)
                    force_jump = true # Force a jump out of the outer loop as well - we got an instant fit.
                    break
                # Try if we can fit the rectangle upright.
                elseif all(widths(rects[j]) .<= widths(rect_i))
                    score = score_by_heuristic(width(rects[j]), height(rects[j]), rect_i, rect_choice)
                    if (score < best_score)
                        best_free_rect = i
                        best_rect = j
                        best_flipped = false
                        best_score = score
                    end
                # If not, then perhaps flipping sideways will make it fit?
                elseif (height(rects[j]) <= width(rect_i) && width(rects[j]) <= height(rect_i))
                    score = score_by_heuristic(height(rects[j]), width(rects[j]), rect_i, rect_choice)
                    if (score < best_score)
                        best_free_rect = i
                        best_rect = j
                        best_flipped = true
                        best_score = score
                    end
                end
            end
            force_jump && break
        end

        # If we didn't manage to find any rectangle to pack, abort.
        (best_score == typemax(Int)) && return nothing

        # Otherwise, we're good to go and do the actual packing.
        xy = minimum(free_rectangles[best_free_rect])
        if best_flipped
            new_node = Rect2(xy, reverse(widths(rects[best_rect])))
        else
            new_node = Rect2(xy, widths(rects[best_rect]))
        end

        # Remove the free space we lost in the bin.
        split_freerect_by_heuristic(packer, free_rectangles[best_free_rect], new_node, split_method)
        splice!(free_rectangles, best_free_rect)

        # Remove the rectangle we just packed from the input list.
        splice!(rects, best_rect)

        # Perform a Rectangle Merge step if desired.
        if packer.merge
            MergeFreeList(packer)
        end

        # Remember the new used rectangle.
        push!(used_rectangles, new_node)
    end
end

function Base.push!(packer::GuillotinePacker, rect::Rect2)
    w, h = widths(rect)
    free_rectangles = packer.free_rectangles
    used_rectangles = packer.used_rectangles
    # Find where to put the new rectangle.
    new_rect, free_node_idx = FindPositionForNewNode(packer, w, h, packer.rect_choice)

    # Abort if we didn't have enough space in the bin.
    if (height(new_rect) == 0)
        return nothing
    end

    # Remove the space that was just consumed by the new rectangle.
    split_freerect_by_heuristic(packer, free_rectangles[free_node_idx], new_rect, packer.split_method)
    splice!(free_rectangles, free_node_idx)

    # Perform a Rectangle Merge step if desired.
    if packer.merge
        MergeFreeList(packer)
    end

    # Remember the new used rectangle.
    push!(used_rectangles, new_rect)

    return new_rect
end

#/ Computes the ratio of used surface area to the total bin area.
function Occupancy(packer::GuillotinePacker)
    used_rectangles = packer.used_rectangles
    #TODO The occupancy rate could be cached/tracked incrementally instead
    # of looping through the list of packed rectangles here.
    surface_area = sum((prod(widths(rect)) for rect in used_rectangles))
    return surface_area / (packer.bin_width * packer.bin_height)

end

#/ Returns the heuristic score value for placing a rectangle of size width*height into free_rect. Does not try to rotate.
function score_by_heuristic(width, height, free_rect::Rect2, rect_choice::FreeRectChoiceHeuristic)
    rect_choice == RectBestAreaFit && return ScoreBestAreaFit(width, height, free_rect)
    rect_choice == RectBestShortSideFit && return ScoreBestShortSideFit(width, height, free_rect)
    rect_choice == RectBestLongSideFit && return ScoreBestLongSideFit(width, height, free_rect)
    rect_choice == RectWorstAreaFit && return ScoreWorstAreaFit(width, height, free_rect)
    rect_choice == RectWorstShortSideFit && return ScoreWorstShortSideFit(width, height, free_rect)
    rect_choice == RectWorstLongSideFit && return ScoreWorstLongSideFit(width, height, free_rect)
end

function ScoreBestAreaFit(w, h, free_rect::Rect2)
    return width(free_rect) * height(free_rect) - w * h
end

function ScoreBestShortSideFit(w, h, free_rect::Rect2)
    leftoverHoriz = abs(width(free_rect) - w)
    leftoverVert = abs(height(free_rect) - h)
    leftover = min(leftoverHoriz, leftoverVert)
    return leftover
end

function ScoreBestLongSideFit(w, h, free_rect::Rect2)
    leftoverHoriz = abs(width(free_rect) - w)
    leftoverVert = abs(height(free_rect) - h)
    leftover = max(leftoverHoriz, leftoverVert)
    return leftover

end

function ScoreWorstAreaFit(width, height, free_rect::Rect2)
    return -ScoreBestAreaFit(width, height, free_rect)
end

function ScoreWorstShortSideFit(width, height, free_rect::Rect2)
    return -ScoreBestShortSideFit(width, height, free_rect)
end

function ScoreWorstLongSideFit(width, height, free_rect::Rect2)
    return -ScoreBestLongSideFit(width, height, free_rect)
end

function FindPositionForNewNode(packer::GuillotinePacker, w, h, rect_choice::FreeRectChoiceHeuristic)
    free_rectangles = packer.free_rectangles
    best_node = Rect2(0,0,0,0)
    free_node_idx = 0
    best_score = typemax(Int)
    #/ Try each free rectangle to find the best one for placement.
    for i in 1:length(free_rectangles)
        rect_i = free_rectangles[i]
        # If this is a perfect fit upright, choose it immediately.
        if (w == width(rect_i) && h == height(rect_i))
            best_node = Rect2(minimum(rect_i), Vec(w, h))
            free_node_idx = i
            break
        # If this is a perfect fit sideways, choose it.
        elseif (h == width(rect_i) && w == height(rect_i))
            best_node = Rect2(minimum(rect_i), Vec(h, w))
            best_score = typemin(Int)
            free_node_idx = i
            break
            # Does the rectangle fit upright?
        elseif (w <= width(rect_i) && h <= height(rect_i))
            score = score_by_heuristic(w, h, rect_i, rect_choice)
            if (score < best_score)
                best_node = Rect2(minimum(rect_i), Vec(w, h))
                best_score = score
                free_node_idx = i
            end
            # Does the rectangle fit sideways?
        elseif (h <= width(rect_i) && w <= height(rect_i))
            score = score_by_heuristic(h, w, rect_i, rect_choice)
            if (score < best_score)
                best_node = Rect2(minimum(rect_i), Vec(h, w))
                best_score = score
                free_node_idx = i
            end
        end
    end
    return best_node, free_node_idx
end

function split_freerect_by_heuristic(packer::GuillotinePacker, free_rect::Rect2, placed_rect::Rect2, method::GuillotineSplitHeuristic)
    # Compute the lengths of the leftover area.
    w, h = widths(free_rect) .- widths(placed_rect)

    # Placing placed_rect into free_rect results in an L-shaped free area, which must be split into
    # two disjorectangles. This can be achieved with by splitting the L-shape using a single line.
    # We have two choices: horizontal or vertical.

    # Use the given heuristic to decide which choice to make.
    split_horizontal = false
    if method == SplitShorterLeftoverAxis
        # Split along the shorter leftover axis.
        split_horizontal = (w <= h)
    elseif method == SplitLongerLeftoverAxis
        # Split along the longer leftover axis.
        split_horizontal = (w > h)
    elseif method == SplitMinimizeArea
        # Maximize the larger area == minimize the smaller area.
        # Tries to make the single bigger rectangle.
        split_horizontal = (width(placed_rect) * h > w * height(placed_rect))
    elseif method == SplitMaximizeArea
        # Maximize the smaller area == minimize the larger area.
        # Tries to make the rectangles more even-sized.
        split_horizontal = (width(placed_rect) * h <= w * height(placed_rect))
    elseif method == SplitShorterAxis
        # Split along the shorter total axis.
        split_horizontal = (width(free_rect) <= height(free_rect))
    elseif method == SplitLongerAxis
        # Split along the longer total axis.
        split_horizontal = (width(free_rect) > height(free_rect))
    else
        error("Method unknown: $(method)")
    end
    # Perform the actual split.
    split_freerect_along_axis(packer, free_rect, placed_rect, split_horizontal)
end

#/ This function will add the two generated rectangles into the free_rectangles array. The caller is expected to
#/ remove the original rectangle from the free_rectangles array after that.
function split_freerect_along_axis(packer::GuillotinePacker, free_rect::Rect2, placed_rect::Rect2, split_horizontal::Bool)
    free_rectangles = packer.free_rectangles
    # Form the two new rectangles.
    x, y = minimum(free_rect)
    if split_horizontal
        bottom = Rect2(x, y + height(placed_rect), width(free_rect), height(free_rect) - height(placed_rect))
        right = Rect2(x + width(placed_rect), y, width(free_rect) - width(placed_rect), height(placed_rect))
    else # Split vertically
        bottom = Rect2(x, y + height(placed_rect), width(placed_rect), height(free_rect) - height(placed_rect))
        right = Rect2(x + width(placed_rect), y, width(free_rect) - width(placed_rect), height(free_rect))
    end

    # Add the new rectangles into the free rectangle pool if they weren't degenerate.
    if (width(bottom) > 0 && height(bottom) > 0)
        push!(free_rectangles, bottom)
    end
    if (width(right) > 0 && height(right) > 0)
        push!(free_rectangles, right)
    end
end

function MergeFreeList(packer::GuillotinePacker)
    free_rectangles = packer.free_rectangles
    # Do a Theta(n^2) loop to see if any pair of free rectangles could me merged into one.
    # Note that we miss any opportunities to merge three rectangles into one. (should call this function again to detect that)
    i = 0; j = 0
    while i < length(free_rectangles)
        i += 1; j = i + 1
        rect_i = free_rectangles[i]
        while j < length(free_rectangles)
            j += 1
            if (width(rect_i) == width(free_rectangles[j]) && minimum(rect_i)[1] == minimum(free_rectangles[j])[1])
                if (minimum(rect_i)[2] == maximum(free_rectangles[j])[2])
                    new_y = minimum(rect_i)[2] - height(free_rectangles[j])
                    new_h = height(rect_i) + height(free_rectangles[j])
                    free_rectangles[i] = Rect2(minimum(rect_i)[1], new_y, Vec(width(rect_i), new_h))
                    splice!(free_rectangles, j)
                    j -= 1
                elseif (maximum(rect_i)[2] == minimum(free_rectangles[j])[2])
                    new_h = height(rect_i) + height(free_rectangles[j])
                    free_rectangles[i] = Rect2(minimum(rect_i), Vec(width(rect_i), new_h))
                    splice!(free_rectangles, j)
                    j -= 1
                end
            elseif (height(rect_i) == height(free_rectangles[j]) && minimum(rect_i)[2] == minimum(free_rectangles[j])[2])
                if (minimum(rect_i)[1] == minimum(free_rectangles[j])[1] + width(free_rectangles[j]))
                    new_x = minimum(rect_i)[1] - width(free_rectangles[j])
                    new_w = width(rect_i) + width(free_rectangles[j])
                    free_rectangles[i] = Rect2(new_x, minimum(rect_i)[2], new_w, height(rect_i))
                    splice!(free_rectangles, j)
                    j -= 1
                elseif maximum(rect_i)[1] == minimum(free_rectangles[j])[1]
                    new_w = Vec(width(rect_i) + width(free_rectangles[j]), height(rect_i))
                    free_rectangles[i] = Rect2(minimum(rect_i), new_w)
                    splice!(free_rectangles, j)
                    j -= 1
                end
            end
        end
    end
end
