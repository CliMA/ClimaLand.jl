const Optional{T} = Union{Nothing, T}

struct RectSides{T}
    left::T
    right::T
    bottom::T
    top::T
end

abstract type Side end

struct Left <: Side end
struct Right <: Side end
struct Top <: Side end
struct Bottom <: Side end
# for protrusion content:
struct TopLeft <: Side end
struct TopRight <: Side end
struct BottomLeft <: Side end
struct BottomRight <: Side end

struct Inner <: Side end
struct Outer <: Side end

abstract type GridDir end
struct Col <: GridDir end
struct Row <: GridDir end

struct RowCols{T <: Union{Number, Vector{Float32}}}
    lefts::T
    rights::T
    tops::T
    bottoms::T
end


"""
    struct Span

Used to specify space that is occupied in a grid. Like 1:1|1:1 for the first square,
or 2:3|1:4 for a rect over the 2nd and 3rd row and the first four columns.
"""
struct Span
    rows::UnitRange{Int64}
    cols::UnitRange{Int64}
end

"""
    mutable struct GridContent{G}

Wraps content elements of a `GridLayout`. It keeps track of the `parent`, the `content` and its position in the grid via `span` and `side`.
"""
mutable struct GridContent{G} # G should be GridLayout but can't be used before definition
    parent::Optional{G}
    content # accessing the content object which can be anything is rare, so avoid overspecialization (type hiding)
    span::Span
    side::Side
    protrusions_handle::Optional{Function}
    reportedsize_handle::Optional{Function}
end

"""
AlignMode that excludes the protrusions from the bounding box. Construct with
`Inside()`.

See also `Outside` and `Mixed`.
"""
struct Inside end

"""
AlignMode that includes the protrusions within the bounding box, plus paddings.

See also `Inside` and `Mixed`.
"""
struct Outside
    padding::RectSides{Float32}
end

"""
    Outside()

Construct an `Outside` AlignMode with no padding.
"""
Outside() = Outside(0f0)
"""
    Outside(padding::Real)

Construct an `Outside` AlignMode with equal padding on all sides.
"""
Outside(padding::Real) = Outside(RectSides{Float32}(padding, padding, padding, padding))
"""
    Outside(left::Real, right::Real, bottom::Real, top::Real)

Construct an `Outside` AlignMode with different paddings on each side.
"""
Outside(left::Real, right::Real, bottom::Real, top::Real) =
    Outside(RectSides{Float32}(left, right, bottom, top))

"""
    Protrusion(p::Float32)

Can be used within a `Mixed` alignmode to override a protrusion manually.
"""
struct Protrusion
    p::Float32
end

"""
AlignMode that is `Inside` where `padding` is `Nothing`, `Outside` where it is
`Real`, and overrides the protrusion with a fixed value where it is a
`Protrusion`.

See also `Inside` and `Outside`.
"""
struct Mixed
    sides::RectSides{Union{Nothing, Float32, Protrusion}}
end

"""
    Mixed(; left = nothing, right = nothing, bottom = nothing, top = nothing)

Construct a `Mixed` AlignMode, which has different behavior on each side.
Arguments that are `nothing` will exclude protrusions from the bounding box on
that side. Those that are real numbers will be padded by that amount and
include protrusions from the bounding box on that side. Arguments that are
`Protrusion` will override the protrusion with a fixed value.
"""
function Mixed(; left = nothing, right = nothing, bottom = nothing, top = nothing)
    sides = map((left, right, bottom, top)) do side
        (side === nothing || side isa Protrusion) ? side : Float32(side)
    end
    Mixed(RectSides{Union{Nothing, Float32, Protrusion}}(sides...))
end

const AlignMode = Union{Inside, Outside, Mixed}


"""
    struct Auto

If used as a `GridLayout`'s row / column size and `trydetermine == true`, signals to the `GridLayout` that the row / column should shrink to match the largest determinable element inside.
If no size of a content element can be determined, the remaining space is split between all `Auto` rows / columns according to their `ratio`.

If used as width / height of a layoutable element and `trydetermine == true`, the element's computed width / height will report the auto width / height if it can be determined.
This enables a parent `GridLayout` to adjust its column / rowsize to the element's width / height.
If `trydetermine == false`, the element's computed width / height will report `nothing` even if an auto width / height can be determined, which will prohibit a parent `GridLayout` from adjusting a row / column to the element's width / height.
This is useful to, e.g., prohibit a `GridLayout` from shrinking a column's width to the width of a super title, even though the title's width can be auto-determined.

The `ratio` is ignored if `Auto` is used as an element size.
"""
struct Auto
    trydetermine::Bool # false for determinable size content that should be ignored
    ratio::Float32 # float ratio in case it's not determinable

    Auto(trydetermine::Bool = true, ratio::Real = 1.0) = new(trydetermine, ratio)
end
Auto(ratio::Real) = Auto(true, ratio)

struct Fixed
    x::Float32
end
struct Relative
    x::Float32
end
struct Aspect
    index::Int
    ratio::Float32
end

const ContentSize = Union{Auto, Fixed, Relative, Aspect}
const GapSize = Union{Fixed, Relative}

struct Dimensions
    inner::Tuple{Optional{Float32}, Optional{Float32}}
    outer::RectSides{Float32}
end

"""
    struct LayoutObservables{G}

`T` is the same type parameter of contained `GridContent`, `G` is `GridLayout` which is defined only after `LayoutObservables`.

A collection of `Observable`s and an optional `GridContent` that are needed to interface with the MakieLayout layouting system.

- `suggestedbbox::Observable{Rect2f}`: The bounding box that an element should place itself in. Depending on the element's `width` and `height` attributes, this is not necessarily equal to the computedbbox.
- `protrusions::Observable{RectSides{Float32}}`: The sizes of content "sticking out" of the main element into the `GridLayout` gaps.
- `reporteddimensions::Observable{Dimensions}`: The dimensions (inner and outer) that the object communicates to the containing `GridLayout`.
- `autosize::Observable{NTuple{2, Optional{Float32}}}`: The width and height that the element reports to its parent `GridLayout`. If the element doesn't want to cause the parent to adjust to its size, autosize can hide the reportedsize from it by being set to `nothing`.
- `computedbbox::Observable{Rect2f}`: The bounding box that the element computes for itself after it has received a suggestedbbox.
- `gridcontent::Optional{GridContent{G}}`: A reference of a `GridContent` if the element is currently placed in a `GridLayout`. This can be used to retrieve the parent layout, remove the element from it or change its position, and assign it to a different layout.
"""
struct LayoutObservables{G} # G again GridLayout
    suggestedbbox::Observable{Rect2f}
    protrusions::Observable{RectSides{Float32}}
    reporteddimensions::Observable{Dimensions}
    autosize::Observable{NTuple{2, Optional{Float32}}}
    computedbbox::Observable{Rect2f}
    gridcontent::Base.RefValue{Optional{GridContent{G}}} # the connecting link to the gridlayout
    block_updates::Base.RefValue{Bool}
end

struct HorizontalAlignment
    x::Float32
end
struct VerticalAlignment
    x::Float32
end

const SizeAttribute = Union{Nothing, Float32, Fixed, Relative, Auto}

mutable struct GridLayout
    parent # this parent is supposed to be any kind of object where it's beneficial
    # to access it through the assigned GridLayout, like a Figure in Makie
    content::Vector{GridContent{GridLayout}}
    size::Tuple{Int, Int}
    offsets::Tuple{Int, Int}
    rowsizes::Vector{ContentSize}
    colsizes::Vector{ContentSize}
    addedrowgaps::Vector{GapSize}
    addedcolgaps::Vector{GapSize}
    alignmode::Observable{AlignMode}
    equalprotrusiongaps::Tuple{Bool, Bool}
    block_updates::Bool
    layoutobservables::LayoutObservables{GridLayout}
    width::Observable{SizeAttribute}
    height::Observable{SizeAttribute}
    tellwidth::Observable{Bool}
    tellheight::Observable{Bool}
    halign::Observable{HorizontalAlignment}
    valign::Observable{VerticalAlignment}
    default_rowgap::GapSize
    default_colgap::GapSize

    function GridLayout(
        parent,
        content, size, offsets, rowsizes, colsizes,
        addedrowgaps, addedcolgaps, alignmode, equalprotrusiongaps,
        layoutobservables, width, height, tellwidth, tellheight, halign, valign, default_rowgap, default_colgap)

        gl = new(parent, content, size, offsets, rowsizes, colsizes,
            addedrowgaps, addedcolgaps, alignmode, equalprotrusiongaps,
            false, layoutobservables, width, height, tellwidth, tellheight,
            halign, valign, default_rowgap, default_colgap)

        validategridlayout(gl)

        gl
    end
end

const Indexables = Union{UnitRange, Int, Colon}
const AutoSize = Union{Nothing, Float32}

struct GridPosition
    layout::GridLayout
    span::Span
    side::Side
end

struct GridSubposition
    parent::Union{GridPosition, GridSubposition}
    rows
    cols
    side::Side
end
