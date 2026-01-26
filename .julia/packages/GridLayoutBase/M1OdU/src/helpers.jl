Base.eltype(::Type{RectSides{T}}) where T = T

"""
Shorthand for `isnothing(optional) ? fallback : optional`
"""
@inline ifnothing(optional, fallback) = isnothing(optional) ? fallback : optional

function Base.foreach(f::Function, contenttype::Type, layout::GridLayout; recursive = true)
    for c in layout.content
        if recursive && c.content isa GridLayout
            foreach(f, contenttype, c.content)
        elseif c.content isa contenttype
            f(c.content)
        end
    end
end

"""
Swaps or rotates the layout positions of the given elements to their neighbor's.
"""
function swap!(layout_elements...)
    gridcontents = gridcontent.(collect(layout_elements))

    # copy relevant fields before gridcontents are mutated
    parents = map(gc -> gc.parent, gridcontents)
    spans = map(gc -> gc.span, gridcontents)
    sides = map(gc -> gc.side, gridcontents)

    for (gc, parent, span, side) in zip(circshift(gridcontents, 1), parents, spans, sides)
        parent[span.rows, span.cols, side] = gc.content
    end
end

function zcumsum(v::AbstractVector{T}) where T
    vpad = [[zero(T)]; v]  # inference-friendly
    cumsum!(vpad, vpad)
end


"""
    grid!(content::Vararg{Pair}; kwargs...)

Creates a GridLayout with all pairs contained in `content`. Each pair consists
of an iterable with row and column spans, and a content object. Each content
object is then placed in the GridLayout at the span from its pair.

Example:

grid!(
    [1, 1] => obj1,
    [1, 2] => obj2,
    [2, :] => obj3,
)
"""
function grid!(content::Vararg{Pair}; kwargs...)
    g = GridLayout(; kwargs...)
    for ((rows, cols), element) in content
        g[rows, cols] = element
    end
    g
end

"""
    hbox!(content::Vararg; kwargs...)

Creates a single-row GridLayout with all elements contained in `content` placed
from left to right.
"""
function hbox!(content::Vararg; kwargs...)
    ncols = length(content)
    g = GridLayout(1, ncols; kwargs...)
    for (i, element) in enumerate(content)
        g[1, i] = element
    end
    g
end

"""
    vbox!(content::Vararg; kwargs...)

Creates a single-column GridLayout with all elements contained in `content` placed
from top to bottom.
"""
function vbox!(content::Vararg; kwargs...)
    nrows = length(content)
    g = GridLayout(nrows, 1; kwargs...)
    for (i, element) in enumerate(content)
        g[i, 1] = element
    end
    g
end



"""
    grid!(content::AbstractMatrix; kwargs...)

Creates a GridLayout filled with matrix-like content. The size of the grid will
be the size of the matrix.
"""
function grid!(content::AbstractMatrix; kwargs...)
    nrows, ncols = size(content)
    g = GridLayout(nrows, ncols; kwargs...)
    for i in 1:nrows, j in 1:ncols
        g[i, j] = content[i, j]
    end
    g
end


left(rect::Rect{2}) = minimum(rect)[1]
right(rect::Rect{2}) = maximum(rect)[1]
bottom(rect::Rect{2}) = minimum(rect)[2]
top(rect::Rect{2}) = maximum(rect)[2]

"""
    BBox(left::Number, right::Number, bottom::Number, top::Number)

Convenience constructor to create a `Rect2` with left, right, bottom and top
extent instead of the usual origin, widths combination.
"""
function BBox(left::T1, right::T2, bottom::T3, top::T4) where {T1 <: Real, T2 <: Real, T3 <: Real, T4 <: Real}
    T = promote_type(T1, T2, T3, T4, Float32) # Float32 to skip Int outputs
    return BBox(T, left, right, bottom, top)
end

function BBox(T::DataType, left::Real, right::Real, bottom::Real, top::Real)
    mini = (left, bottom)
    maxi = (right, top)
    return Rect2{T}(mini, maxi .- mini)
end

function RowCols(ncols::Integer, nrows::Integer)
    return RowCols(
        zeros(Float32, ncols),
        zeros(Float32, ncols),
        zeros(Float32, nrows),
        zeros(Float32, nrows)
    )
end

Base.getindex(rowcols::RowCols, ::Left) = rowcols.lefts
Base.getindex(rowcols::RowCols, ::Right) = rowcols.rights
Base.getindex(rowcols::RowCols, ::Top) = rowcols.tops
Base.getindex(rowcols::RowCols, ::Bottom) = rowcols.bottoms

"""
    eachside(f)
Calls f over all sides (Left, Right, Top, Bottom), and creates a BBox from the result of f(side)
"""
function eachside(f)
    return BBox(f(Left()), f(Right()), f(Bottom()), f(Top()))
end

"""
mapsides(
       f, first::Union{Rect{2}, RowCols}, rest::Union{Rect{2}, RowCols}...
   )::Rect2f
Maps f over all sides of the rectangle like arguments.
e.g.
```
mapsides(BBox(left, right, bottom, top)) do side::Side, side_val::Number
    return ...
end::Rect2f
```
"""
function mapsides(
        f, first::Union{Rect{2}, RowCols}, rest::Union{Rect{2}, RowCols}...
    )
    return eachside() do side
        f(side, getindex.((first, rest...), (side,))...)
    end
end

function set_nrows!(gl, x)
    gl.size = (x, gl.size[2])
end
function set_ncols!(gl, x)
    gl.size = (gl.size[1], x)
end

function set_rowoffset!(gl, x)
    gl.offsets = (x, offset(gl, Col()))
end
function set_coloffset!(gl, x)
    gl.offsets = (offset(gl, Row()), x)
end

offset(gl, ::Row) = offsets(gl)[1]
offset(gl, ::Col) = offsets(gl)[2]

# convert an index into an array from 1:nrow or 1:ncol
# into the respective column / row number that can also be negative if offset
offset(gl, i, ::Row) = i + offset(gl, Row())
offset(gl, i, ::Col) = i + offset(gl, Col())

# convert a column / row number that can also be negative if offset
# to an index from 1:nrow or 1:ncol
unoffset(gl, i, ::Row) = i - offset(gl, Row())
unoffset(gl, i, ::Col) = i - offset(gl, Col())
