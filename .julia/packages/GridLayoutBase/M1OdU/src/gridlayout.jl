"""
    GridLayout(; kwargs...)

Create a `GridLayout` without parent and with size [1, 1].
"""
GridLayout(; kwargs...) = GridLayout(1, 1; kwargs...)

"""
    GridLayout(g::Union{GridPosition, GridSubposition}, args...; kwargs...)

Create a `GridLayout` at position `g` in the parent `GridLayout` of `g` if it is a `GridPosition`
and in a nested child `GridLayout` if it is a `GridSubposition`. The `args` and `kwargs` are passed on to the normal `GridLayout` constructor.
"""
function GridLayout(g::Union{GridPosition, GridSubposition}, args...; kwargs...)
    return g[] = GridLayout(args...; kwargs...)
end

observablify(x::Observable) = x
observablify(x, type=Any) = Observable{type}(x)

get_default_rowgap()::Float64 = DEFAULT_ROWGAP_GETTER[]()
get_default_colgap()::Float64 = DEFAULT_COLGAP_GETTER[]()

Base.convert(::Type{T}, a::Real) where T <: Union{HorizontalAlignment, VerticalAlignment} = T(a)
function Base.convert(::Type{HorizontalAlignment}, s::Symbol)
    if s === :left
        HorizontalAlignment(0.0)
    elseif s === :right
        HorizontalAlignment(1.0)
    elseif s === :center
        HorizontalAlignment(0.5)
    else
        error("Symbol $s is not a horizontal alignment. Possible values are :left, :right, :center.")
    end
end
function Base.convert(::Type{VerticalAlignment}, s::Symbol)
    if s === :bottom
        VerticalAlignment(0.0)
    elseif s === :top
        VerticalAlignment(1.0)
    elseif s === :center
        VerticalAlignment(0.5)
    else
        error("Symbol $s is not a vertical alignment. Possible values are :bottom, :top, :center.")
    end
end

# maybe this is not ideal, converting to an abstract type?
Base.convert(::Type{SizeAttribute}, r::Real) = Float32(r)
Base.convert(::Type{SizeAttribute}, f::Fixed) = f
Base.convert(::Type{SizeAttribute}, r::Relative) = r
Base.convert(::Type{SizeAttribute}, a::Auto) = a

@inline function Base.setproperty!(g::GridLayout, s::Symbol, value)
    T = fieldtype(GridLayout, s)
    if T <: Observable
        setindex!(getfield(g, s), value)
    else
        setfield!(g, s, convert(T, value))
    end
end

"""
    function GridLayout(nrows::Integer, ncols::Integer;
        parent = nothing,
        rowsizes = nothing,
        colsizes = nothing,
        addedrowgaps = nothing,
        addedcolgaps = nothing,
        alignmode = Inside(),
        equalprotrusiongaps = (false, false),
        bbox = nothing,
        width = Auto(),
        height = Auto(),
        tellwidth::Bool = true,
        tellheight::Bool = true,
        halign = :center,
        valign = :center,
        default_rowgap = get_default_rowgap(),
        default_colgap = get_default_colgap(),
        kwargs...)

Create a `GridLayout` with optional parent `parent` with `nrows` rows and `ncols` columns.
"""
function GridLayout(nrows::Integer, ncols::Integer;
        parent = nothing,
        rowsizes = nothing,
        colsizes = nothing,
        addedrowgaps = nothing,
        addedcolgaps = nothing,
        alignmode = Inside(),
        equalprotrusiongaps = (false, false),
        bbox = nothing,
        width = Auto(),
        height = Auto(),
        tellwidth::Bool = true,
        tellheight::Bool = true,
        halign = :center,
        valign = :center,
        default_rowgap = get_default_rowgap(),
        default_colgap = get_default_colgap(),
        kwargs...)

    default_rowgap::GapSize = default_rowgap isa Number ? Fixed(default_rowgap)::Fixed : default_rowgap
    default_colgap::GapSize = default_colgap isa Number ? Fixed(default_colgap)::Fixed : default_colgap
    rowsizes = convert_contentsizes(nrows, rowsizes)
    colsizes = convert_contentsizes(ncols, colsizes)
    addedrowgaps = convert_gapsizes(nrows - 1, addedrowgaps, default_rowgap)
    addedcolgaps = convert_gapsizes(ncols - 1, addedcolgaps, default_colgap)

    content = GridContent[]

    alignmode = observablify(alignmode, AlignMode)
    width_obs = convert(Observable{SizeAttribute}, width)
    height_obs = convert(Observable{SizeAttribute}, height)
    tellwidth_obs = convert(Observable{Bool}, tellwidth)
    tellheight_obs = convert(Observable{Bool}, tellheight)
    halign_obs = convert(Observable{HorizontalAlignment}, halign)
    valign_obs = convert(Observable{VerticalAlignment}, valign)

    layoutobservables = layoutobservables = LayoutObservables(
        width_obs,
        height_obs,
        tellwidth_obs,
        tellheight_obs,
        halign_obs,
        valign_obs;
        suggestedbbox = bbox
    )

    offsets = (0, 0)

    gl = GridLayout(
        parent,
        content,
        (nrows, ncols),
        offsets,
        rowsizes,
        colsizes,
        addedrowgaps,
        addedcolgaps,
        alignmode,
        equalprotrusiongaps,
        layoutobservables,
        width_obs,
        height_obs,
        tellwidth_obs,
        tellheight_obs,
        halign_obs,
        valign_obs,
        default_rowgap,
        default_colgap
    )
    on(computedbboxobservable(gl)) do bbox
        # block_updates can block update! but not setting computedbbox directly through protrusions/etc
        # so this one needs to be blocked as well
        if !gl.block_updates
            align_to_bbox!(gl, bbox)
        end
    end
    gl
end

function update!(gl::GridLayout)
    gl.block_updates && return

    w = determinedirsize(gl, Col())
    h = determinedirsize(gl, Row())

    new_autosize = (w, h)
    new_protrusions = RectSides{Float32}(
        compute_effective_protrusion(gl, Left()),
        compute_effective_protrusion(gl, Right()),
        compute_effective_protrusion(gl, Bottom()),
        compute_effective_protrusion(gl, Top()),
    )

    # if autosize and protrusions didn't change, then we only need to trigger a relayout in this gridlayout
    # but don't retrigger autosize and protrusions which would go up to a possible parent gridlayout

    gc = gridcontent(gl)

    # TODO: if we skip this because protrusions/autosize are unchanged, we get update failures when adding new objects that don't change these values to a gridlayout
    if gc === nothing
        # no parent exists, so protrusions/autosize will trigger computedbbox, which will trigger relayout
        bu = gl.block_updates
        # this will prohibit align_to_bbox! triggering from computedbbox
        gl.block_updates = true
        protrusionsobservable(gl)[] = new_protrusions
        autosizeobservable(gl)[] = new_autosize
        gl.block_updates = bu
        # now we trigger computedbbox ourselves which will call align_to_bbox! and relayout children
        notify(computedbboxobservable(gl))
    else
        # parent exists, so protrusions/autosize will not trigger computedbbox
        # don't update parent by setting protrusions / autosize
        layoutobservables(gl).block_updates[] = true

        protrusionsobservable(gl)[] = new_protrusions
        autosizeobservable(gl)[] = new_autosize

        layoutobservables(gl).block_updates[] = false
        # update parent manually
        update!(gc)
    end

    return
end

function validategridlayout(gl::GridLayout)
    if nrows(gl) < 1
        error("Number of rows can't be smaller than 1")
    end
    if ncols(gl) < 1
        error("Number of columns can't be smaller than 1")
    end

    if length(gl.rowsizes) != nrows(gl)
        error("There are $nrows rows but $(length(gl.rowsizes)) row sizes.")
    end
    if length(gl.colsizes) != ncols(gl)
        error("There are $ncols columns but $(length(gl.colsizes)) column sizes.")
    end
    if length(gl.addedrowgaps) != nrows(gl) - 1
        error("There are $nrows rows but $(length(gl.addedrowgaps)) row gaps.")
    end
    if length(gl.addedcolgaps) != ncols(gl) - 1
        error("There are $ncols columns but $(length(gl.addedcolgaps)) column gaps.")
    end
end

"""
    with_updates_suspended(f::Function, gl::GridLayout; update = true)

Disable layout updates for `gl` and call the function `f`. If `update` is true,
force a layout update after `f` returns.
"""
function with_updates_suspended(f::Function, gl::GridLayout; update = true)
    prev_block_value = gl.block_updates
    gl.block_updates = true
    f()
    gl.block_updates = prev_block_value
    if update
        update!(gl)
    end
    return
end

function connect_layoutobservables!(gc::GridContent)

    disconnect_layoutobservables!(gc::GridContent)

    let content = gc.content
        # gc.protrusions_handle = on(effectiveprotrusionsobservable(content)) do p
        #     if !layoutobservables(content).block_updates[]::Bool
        #         update!(gc)
        #     end
        # end
        gc.reportedsize_handle = on(reporteddimensionsobservable(content)) do c
            if !layoutobservables(content).block_updates[]::Bool
                update!(gc)
            end
        end
    end
end

function disconnect_layoutobservables!(gc::GridContent)
    # if !isnothing(gc.protrusions_handle)
    #     Observables.off(effectiveprotrusionsobservable(gc.content), gc.protrusions_handle)
    #     gc.protrusions_handle = nothing
    # end
    if !isnothing(gc.reportedsize_handle)
        Observables.off(reporteddimensionsobservable(gc.content), gc.reportedsize_handle)
        gc.reportedsize_handle = nothing
    end
end

function add_to_gridlayout!(g::GridLayout, gc::GridContent)
    # to be safe
    remove_from_gridlayout!(gc)

    push!(g.content, gc)

    # let the gridcontent know that it's inside a gridlayout
    gc.parent = g
    # change the parent if the gridcontent contains a gridlayout
    content = gc.content
    if content isa GridLayout
        content.parent = g
    end

    # trigger relayout
    update!(g)
end


function remove_from_gridlayout!(gc::GridContent)
    content = gc.content
    if isnothing(gc.parent)
        if content isa GridLayout
            content.parent = nothing
        end
        return
    end

    i = findfirst(x -> x === gc, gc.parent.content)
    if isnothing(i)
        error("""GridContent had a parent but was not in the content array.
        This must be a bug.""")
    end
    deleteat!(gc.parent.content, i)

    gc.parent = nothing
    # set the parent of a gridlayout content to nothing separately
    # this is mostly for one toplevel parent like a Figure in Makie
    if content isa GridLayout
        content.parent = nothing
    end

    return
end


function convert_contentsizes(n, sizes)::Vector{ContentSize}
    if sizes === nothing
        ContentSize[Auto() for _ in 1:n]
    elseif sizes isa ContentSize
        ContentSize[sizes for _ in 1:n]
    elseif sizes isa Vector{<:ContentSize}
        length(sizes) == n ? sizes : error("$(length(sizes)) sizes instead of $n")
    else
        error("Illegal sizes value $sizes")
    end
end

function convert_gapsizes(n, gaps, defaultsize)::Vector{GapSize}
    if gaps === nothing
        GapSize[defaultsize for _ in 1:n]
    elseif gaps isa GapSize
        GapSize[gaps for _ in 1:n]
    elseif gaps isa Vector{<:GapSize}
        length(gaps) == n ? gaps : error("$(length(gaps)) gaps instead of $n")
    else
        error("Illegal gaps value $gaps")
    end
end

function appendrows!(gl::GridLayout, n::Integer; rowsizes=nothing, addedrowgaps=nothing, update = true)

    rowsizes = convert_contentsizes(n, rowsizes)
    addedrowgaps = convert_gapsizes(n, addedrowgaps, gl.default_rowgap)

    with_updates_suspended(gl, update = update) do
        set_nrows!(gl, nrows(gl) + n)
        append!(gl.rowsizes, rowsizes)
        append!(gl.addedrowgaps, addedrowgaps)
    end
end

function appendcols!(gl::GridLayout, n::Integer; colsizes=nothing, addedcolgaps=nothing, update = true)
    colsizes = convert_contentsizes(n, colsizes)
    addedcolgaps = convert_gapsizes(n, addedcolgaps, gl.default_colgap)

    with_updates_suspended(gl, update = update) do
        set_ncols!(gl, ncols(gl) + n)
        append!(gl.colsizes, colsizes)
        append!(gl.addedcolgaps, addedcolgaps)
    end
end

function prependrows!(gl::GridLayout, n::Integer; rowsizes=nothing, addedrowgaps=nothing, update = true)

    rowsizes = convert_contentsizes(n, rowsizes)
    addedrowgaps = convert_gapsizes(n, addedrowgaps, gl.default_rowgap)

    with_updates_suspended(gl, update = update) do
        set_nrows!(gl, nrows(gl) + n)
        set_rowoffset!(gl, offset(gl, Row()) - n)
        prepend!(gl.rowsizes, rowsizes)
        prepend!(gl.addedrowgaps, addedrowgaps)
    end
end

function prependcols!(gl::GridLayout, n::Integer; colsizes=nothing, addedcolgaps=nothing, update = true)

    colsizes = convert_contentsizes(n, colsizes)
    addedcolgaps = convert_gapsizes(n, addedcolgaps, gl.default_colgap)

    with_updates_suspended(gl, update = update) do
        set_ncols!(gl, ncols(gl) + n)
        set_coloffset!(gl, offset(gl, Col()) - n)
        prepend!(gl.colsizes, colsizes)
        prepend!(gl.addedcolgaps, addedcolgaps)
    end
end


"""
    insertrows!(gl::GridLayout, at::Integer, n::Integer; rowsizes=nothing, addedrowgaps=nothing)

Insert `n` rows at row `at` into `GridLayout` `gl`. The new row sizes and row gaps can be
optionally set with `rowsizes` and `addedrowgaps` keywords.
Objects spanning from at least `at-1` up to or beyond `at` are getting extended to span
over the new rows. Objects from `at` and beyond are pushed back, objects before `at` are
unaffected.
"""
function insertrows!(gl::GridLayout, at::Integer, n::Integer; rowsizes=nothing, addedrowgaps=nothing)

    if !(1 <= at <= nrows(gl))
        error("Invalid row insertion at row $at. GridLayout has $(nrows(gl)) rows.")
    end

    rowsizes = convert_contentsizes(n, rowsizes)
    addedrowgaps = convert_gapsizes(n, addedrowgaps, gl.default_rowgap)

    foreach(gl.content) do gc
        span = gc.span
        rows = span.rows
        newrows = if rows.start < at <= rows.stop
            rows.start : rows.stop + n
        elseif rows.stop < at
            rows
        elseif rows.start >= at
            rows .+ n
        end
        newspan = Span(newrows, span.cols)
        gc.span = newspan
    end

    with_updates_suspended(gl) do
        set_nrows!(gl, nrows(gl) + n)
        splice!(gl.rowsizes, at:at-1, rowsizes)
        splice!(gl.addedrowgaps, at:at-1, addedrowgaps)
    end
end

"""
    insertcols!(gl::GridLayout, at::Integer, n::Integer; colsizes=nothing, addedcolgaps=nothing)

Insert `n` columns at column `at` into `GridLayout` `gl`. The new column sizes and column gaps can be
optionally set with `colsizes` and `addedcolgaps` keywords.
Objects spanning from at least `at-1` up to or beyond `at` are getting extended to span
over the new columns. Objects from `at` and beyond are pushed back, objects before `at` are
unaffected.
"""
function insertcols!(gl::GridLayout, at::Integer, n::Integer; colsizes=nothing, addedcolgaps=nothing)

    if !(1 <= at <= ncols(gl))
        error("Invalid column insertion at column $at. GridLayout has $(ncols(gl)) columns.")
    end

    colsizes = convert_contentsizes(n, colsizes)
    addedcolgaps = convert_gapsizes(n, addedcolgaps, gl.default_colgap)

    foreach(gl.content) do gc
        span = gc.span
        cols = span.cols
        newcols = if cols.start < at <= cols.stop
            cols.start : cols.stop + n
        elseif cols.stop < at
            cols
        elseif cols.start >= at
            cols .+ n
        end
        newspan = Span(span.rows, newcols)
        gc.span = newspan
    end

    with_updates_suspended(gl) do
        set_ncols!(gl, ncols(gl) + n)
        splice!(gl.colsizes, at:at-1, colsizes)
        splice!(gl.addedcolgaps, at:at-1, addedcolgaps)
    end
end

function deleterow!(gl::GridLayout, irow::Integer)
    if !(firstrow(gl) <= irow <= lastrow(gl))
        error("Row $irow does not exist.")
    end

    if nrows(gl) == 1
        error("Can't delete the only row.")
    end

    # new_content = GridContent[]
    to_remove = GridContent[]
    for c in gl.content
        rows = c.span.rows
        newrows = if irow in rows
            # range is one shorter now
            rows.start : rows.stop - 1
        elseif irow > rows.stop
            # content before deleted row stays the same
            rows
        else
            # content completely after is moved forward 1 step
            rows .- 1
        end
        if isempty(newrows)
            # the row span was just one row and now zero, remove the element
            push!(to_remove, c)
        else
            c.span = Span(newrows, c.span.cols)
        end
    end

    for c in to_remove
        remove_from_gridlayout!(c)
    end

    idx = irow - rowoffset(gl)
    deleteat!(gl.rowsizes, idx)
    deleteat!(gl.addedrowgaps, idx == 1 ? 1 : idx - 1)
    set_nrows!(gl, nrows(gl) - 1)
    update!(gl)
end

rowoffset(gl) = offset(gl, Row())
coloffset(gl) = offset(gl, Col())

function deletecol!(gl::GridLayout, icol::Integer)
    if !(firstcol(gl) <= icol <= lastcol(gl))
        error("Col $icol does not exist.")
    end

    if ncols(gl) == 1
        error("Can't delete the only column.")
    end

    to_remove = GridContent[]
    for c in gl.content
        cols = c.span.cols
        newcols = if icol in cols
            # range is one shorter now
            cols.start : cols.stop - 1
        elseif icol > cols.stop
            # content before deleted col stays the same
            cols
        else
            # content completely after is moved forward 1 step
            cols .- 1
        end
        if isempty(newcols)
            # the col span was just one col and now zero, remove the element
            push!(to_remove, c)
        else
            c.span = Span(c.span.rows, newcols)
        end
    end

    for c in to_remove
        remove_from_gridlayout!(c)
    end

    idx = icol - coloffset(gl)
    deleteat!(gl.colsizes, idx)
    deleteat!(gl.addedcolgaps, idx == 1 ? 1 : idx - 1)
    set_ncols!(gl, ncols(gl) - 1)
    update!(gl)
end

function Base.isempty(gl::GridLayout, dir::GridDir, i::Integer)
    !any(gl.content) do c
        span = dir isa Row ? c.span.rows : c.span.cols
        i in span
    end
end

"""
    trim!(gl::GridLayout)

Remove empty rows and columns from `gl`.
"""
function trim!(gl::GridLayout)
    irow = firstrow(gl)
    while irow <= lastrow(gl) && nrows(gl) > 1
        if isempty(gl, Row(), irow)
            deleterow!(gl, irow)
        else
            irow += 1
        end
    end

    icol = firstcol(gl)
    while icol <= lastcol(gl) && ncols(gl) > 1
        if isempty(gl, Col(), icol)
            deletecol!(gl, icol)
        else
            icol += 1
        end
    end
end

function gridnest!(gl::GridLayout, rows::Indexables, cols::Indexables)

    newrows, newcols = adjust_rows_cols!(gl, rows, cols)

    subgl = GridLayout(
        length(newrows), length(newcols);
        parent = nothing,
        colsizes = gl.colsizes[newcols],
        rowsizes = gl.rowsizes[newrows],
        addedrowgaps = gl.addedrowgaps[newrows.start:(newrows.stop-1)],
        addedcolgaps = gl.addedcolgaps[newcols.start:(newcols.stop-1)],
    )

    # remove the content from the parent that is completely inside the replacement grid
    subgl.block_updates = true
    i = 1
    while i <= length(gl.content)
        gc = gl.content[i]

        if (gc.span.rows.start >= newrows.start && gc.span.rows.stop <= newrows.stop &&
            gc.span.cols.start >= newcols.start && gc.span.cols.stop <= newcols.stop)


            subgl[gc.span.rows .- (newrows.start - 1), gc.span.cols .- (newcols.start - 1), gc.side] = gc.content
            continue
            # don't advance i because there's one piece of content less in the queue
            # and the next item is in the same position as the old removed one
        end
        i += 1
    end
    subgl.block_updates = false

    gl[newrows, newcols] = subgl

    subgl
end


function Base.show(io::IO, ::MIME"text/plain", gl::GridLayout)

    function spaceindent(str, n, downconnection)
        joinstr = if downconnection
            "\n" * (" " ^ 1) * "┃" * (" " ^ (n-2))
        else
            "\n" * (" " ^ n)
        end
        join(split(str, "\n"), joinstr)
    end

    rowrange = firstrow(gl):lastrow(gl)
    colrange = firstcol(gl):lastcol(gl)

    println(io, "GridLayout[$rowrange, $colrange] with $(length(gl.content)) children")

    simplespan(span) = span.start == span.stop ? span.start : span

    for (i, c) in enumerate(gl.content)
        rows = c.span.rows
        cols = c.span.cols
        content = c.content

        connector = i == length(gl.content) ? " ┗━ " : " ┣━ "

        if content isa GridLayout
            downconnection = i < length(gl.content)
            str = spaceindent(repr(MIME"text/plain"(), content), 2, downconnection)
            println(io, connector * "[$(simplespan(rows)), $(simplespan(cols))] $str")
        else
            println(io, connector * "[$(simplespan(rows)), $(simplespan(cols))] $(typeof(content))")
        end
    end
end

function Base.show(io::IO, gl::GridLayout)
    print(io, "GridLayout[$(nrows(gl)), $(ncols(gl))] ($(length(gl.content)) children)")
end


"""
    colsize!(gl::GridLayout, i::Integer, s::Union{Aspect, Auto, Fixed, Relative, Real})

Set the size of the `i`th column in `gl`, i.e., `gl[:, i]`.
Passing a real number to `s` has the same behaviour as passing `Fixed(s)`.

See also [Aspect](@ref), [Auto](@ref), [Fixed](@ref), and [Relative](@ref).
"""
function colsize!(gl::GridLayout, i::Integer, s::ContentSize)
    if !(firstcol(gl) <= i <= lastcol(gl))
        error("Can't set size of invalid column $i.")
    end
    i_unoffset = unoffset(gl, i, Col())
    gl.colsizes[i_unoffset] = s
    update!(gl)
end

colsize!(gl::GridLayout, i::Integer, s::Real) = colsize!(gl, i, Fixed(s))

"""
    rowsize!(gl::GridLayout, i::Integer, s::Union{Aspect, Auto, Fixed, Relative, Real})

Set the size of the `i`th row in `gl`, i.e., `gl[i, :]`.
Passing a real number to `s` has the same behaviour as passing `Fixed(s)`.

See also [Aspect](@ref), [Auto](@ref), [Fixed](@ref), and [Relative](@ref).
"""
function rowsize!(gl::GridLayout, i::Integer, s::ContentSize)
    if !(firstrow(gl) <= i <= lastrow(gl))
        error("Can't set size of invalid row $i.")
    end
    i_unoffset = unoffset(gl, i, Row())
    gl.rowsizes[i_unoffset] = s
    update!(gl)
end

rowsize!(gl::GridLayout, i::Integer, s::Real) = rowsize!(gl, i, Fixed(s))

"""
    colgap!(gl::GridLayout, i::Integer, s::Union{Fixed, Relative, Real})
    colgap!(gl::GridLayout, s::Union{Fixed, Relative, Real})

Set the gap between columns in `gl`.  The two-argument version sets all column gaps
in `gl`.  The three-argument version sets the gap between columns `i` and `i+1`.
Passing a real number to `s` has the same behaviour as passing `Fixed(s)`.

See also [Fixed](@ref) and [Relative](@ref).
"""
function colgap!(gl::GridLayout, i::Integer, s::GapSize)
    if !(1 <= i <= (ncols(gl) - 1))
        error("Can't set size of invalid column gap $i.")
    end
    gl.addedcolgaps[i] = s
    update!(gl)
end

colgap!(gl::GridLayout, i::Integer, s::Real) = colgap!(gl, i, Fixed(s))

function colgap!(gl::GridLayout, s::GapSize)
    gl.addedcolgaps .= Ref(s)
    update!(gl)
end

function colgap!(gl::GridLayout, r::Real)
    gl.addedcolgaps .= Ref(Fixed(r))
    update!(gl)
end

"""
    rowgap!(gl::GridLayout, i::Integer, s::Union{Fixed, Relative, Real})
    rowgap!(gl::GridLayout, s::Union{Fixed, Relative, Real})

Set the gap between rows in `gl`.  The two-argument version sets all row gaps
in `gl`.  The three-argument version sets the gap between rows `i` and `i+1`.
Passing a real number to `s` has the same behaviour as passing `Fixed(s)`.

See also [Fixed](@ref) and [Relative](@ref).
"""
function rowgap!(gl::GridLayout, i::Integer, s::GapSize)
    if !(1 <= i <= (nrows(gl) - 1))
        error("Can't set size of invalid row gap $i.")
    end
    gl.addedrowgaps[i] = s
    update!(gl)
end

rowgap!(gl::GridLayout, i::Integer, s::Real) = rowgap!(gl, i, Fixed(s))

function rowgap!(gl::GridLayout, s::GapSize)
    gl.addedrowgaps .= Ref(s)
    update!(gl)
end

function rowgap!(gl::GridLayout, r::Real)
    gl.addedrowgaps .= Ref(Fixed(r))
    update!(gl)
end


function _compute_content_bbox(suggestedbbox, alignmode::Outside)
    pad = alignmode.padding
    BBox(
        left(suggestedbbox) + pad.left,
        right(suggestedbbox) - pad.right,
        bottom(suggestedbbox) + pad.bottom,
        top(suggestedbbox) - pad.top)
end

function _compute_content_bbox(suggestedbbox, ::Inside)
    suggestedbbox
end

function _compute_content_bbox(suggestedbbox, alignmode::Mixed)
    sides = alignmode.sides
    BBox(
        left(suggestedbbox) + (sides.left isa Float32 ? sides.left : 0f0),
        right(suggestedbbox) - (sides.right isa Float32 ? sides.right : 0f0),
        bottom(suggestedbbox) + (sides.bottom isa Float32 ? sides.bottom : 0f0),
        top(suggestedbbox) - (sides.top isa Float32 ? sides.top : 0f0))
end

# compute a grid of maximum protrusions for a gridlayout
function _compute_maxgrid(gl)
    maxgrid = RowCols(ncols(gl), nrows(gl))
    for c in gl.content
        # TODO this RowCols{Int} should actually just be a RectSides and RowCols should always store vectors
        idx_rect = side_indices(gl, c)
        mapsides(idx_rect, maxgrid) do side, idx, grid
            grid[idx] = max(grid[idx], effective_protrusion(c, side, c.side)::Float32)
        end
    end
    maxgrid
end

effective_protrusion(c::GridContent, ::Left) = reporteddimensionsobservable(c.content)[].outer.left
effective_protrusion(c::GridContent, ::Right) = reporteddimensionsobservable(c.content)[].outer.right
effective_protrusion(c::GridContent, ::Top) = reporteddimensionsobservable(c.content)[].outer.top
effective_protrusion(c::GridContent, ::Bottom) = reporteddimensionsobservable(c.content)[].outer.bottom

effective_protrusion(c::GridContent, ::Left, ::Inner) = reporteddimensionsobservable(c.content)[].outer.left
effective_protrusion(c::GridContent, ::Right, ::Inner) = reporteddimensionsobservable(c.content)[].outer.right
effective_protrusion(c::GridContent, ::Top, ::Inner) = reporteddimensionsobservable(c.content)[].outer.top
effective_protrusion(c::GridContent, ::Bottom, ::Inner) = reporteddimensionsobservable(c.content)[].outer.bottom

effective_protrusion(c::GridContent, ::Left, s::Union{Left, BottomLeft, TopLeft}) = something(determinedirsize(c, Col(), s), 0f0)
effective_protrusion(c::GridContent, ::Right, s::Union{Right, BottomRight, TopRight}) = something(determinedirsize(c, Col(), s), 0f0)
effective_protrusion(c::GridContent, ::Bottom, s::Union{Bottom, BottomRight, BottomLeft}) = something(determinedirsize(c, Row(), s), 0f0)
effective_protrusion(c::GridContent, ::Top, s::Union{Top, TopRight, TopLeft}) = something(determinedirsize(c, Row(), s), 0f0)
effective_protrusion(c::GridContent, _, _) = 0f0

sideoffset(gl, ::Union{Right, Left}) = offset(gl, Col())
sideoffset(gl, ::Union{Top, Bottom}) = offset(gl, Row())

function _compute_remaining_horizontal_space(content_bbox, sumcolgaps, leftprot, rightprot, alignmode::Inside)::Float32
    width(content_bbox) - sumcolgaps
end

function _compute_remaining_horizontal_space(content_bbox, sumcolgaps, leftprot, rightprot, alignmode::Outside)::Float32

    width(content_bbox) - sumcolgaps - leftprot - rightprot
end

function _compute_remaining_horizontal_space(content_bbox, sumcolgaps, leftprot, rightprot, alignmode::Mixed)::Float32

    rightal = getside(alignmode, Right())
    leftal = getside(alignmode, Left())
    width(content_bbox) - sumcolgaps -
        (isnothing(leftal) ? zero(leftprot) : isa(leftal, Protrusion) ? leftal.p : leftprot) -
        (isnothing(rightal) ? zero(rightprot) : isa(rightal, Protrusion) ? rightal.p : rightprot)
end

function _compute_remaining_vertical_space(content_bbox, sumrowgaps, topprot, bottomprot, alignmode::Inside)::Float32

    height(content_bbox) - sumrowgaps
end

function _compute_remaining_vertical_space(content_bbox, sumrowgaps, topprot, bottomprot, alignmode::Outside)::Float32

    height(content_bbox) - sumrowgaps - topprot - bottomprot
end

function _compute_remaining_vertical_space(content_bbox, sumrowgaps, topprot, bottomprot, alignmode::Mixed)::Float32

    topal = getside(alignmode, Top())
    bottomal = getside(alignmode, Bottom())
    height(content_bbox) - sumrowgaps -
        (isnothing(bottomal) ? zero(bottomprot) : isa(bottomal, Protrusion) ? bottomal.p : bottomprot) -
        (isnothing(topal) ? zero(topprot) : isa(topal, Protrusion) ? topal.p : topprot)
end

"""
This function solves a grid layout such that the "important lines" fit exactly
into a given bounding box. This means that the protrusions of all objects inside
the grid are not taken into account. This is needed if the grid is itself placed
inside another grid.
"""
function compute_rowcols(gl::GridLayout, suggestedbbox::Rect2f)
    # compute the actual bbox for the content given that there might be outside
    # padding that needs to be removed
    alignmode = gl.alignmode[]
    content_bbox = _compute_content_bbox(suggestedbbox, alignmode)

    # first determine how big the protrusions on each side of all columns and rows are
    maxgrid = _compute_maxgrid(gl)

    # for the outside alignmode
    topprot = maxgrid.tops[1]
    bottomprot = maxgrid.bottoms[end]
    leftprot = maxgrid.lefts[1]
    rightprot = maxgrid.rights[end]

    # compute what size the gaps between rows and columns need to be
    colgaps = maxgrid.lefts[2:end] .+ maxgrid.rights[1:end-1]
    rowgaps = maxgrid.tops[2:end] .+ maxgrid.bottoms[1:end-1]

    # determine the biggest gap
    # using the biggest gap size for all gaps will make the layout more even
    if gl.equalprotrusiongaps[2]
        colgaps = ncols(gl) <= 1 ? [0f0] : fill(maximum(colgaps), ncols(gl) - 1)
    end
    if gl.equalprotrusiongaps[1]
        rowgaps = nrows(gl) <= 1 ? [0f0] : fill(maximum(rowgaps), nrows(gl) - 1)
    end

    # determine the vertical and horizontal space needed just for the gaps
    # again, the gaps are what the protrusions stick into, so they are not actually "empty"
    # depending on what sticks out of the plots
    sumcolgaps = (ncols(gl) <= 1) ? 0f0 : sum(colgaps)
    sumrowgaps = (nrows(gl) <= 1) ? 0f0 : sum(rowgaps)

    # compute what space remains for the inner parts of the plots
    remaininghorizontalspace = _compute_remaining_horizontal_space(content_bbox, sumcolgaps, leftprot, rightprot, alignmode)

    remainingverticalspace = _compute_remaining_vertical_space(content_bbox, sumrowgaps, topprot, bottomprot, alignmode)

    # compute how much gap to add, in case e.g. labels are too close together
    # this is given as a fraction of the space used for the inner parts of the plots
    # so far, but maybe this should just be an absolute pixel value so it doesn't change
    # when resizing the window
    addedcolgaps::Vector{Float32} = map(gl.addedcolgaps) do cg
        if cg isa Fixed
            return cg.x
        elseif cg isa Relative
            return cg.x * remaininghorizontalspace
        else
            return 0f0 # for float type inference
        end
    end
    addedrowgaps::Vector{Float32} = map(gl.addedrowgaps) do rg
        if rg isa Fixed
            return rg.x
        elseif rg isa Relative
            return rg.x * remainingverticalspace
        else
            return 0f0 # for float type inference
        end
    end

    # compute the actual space available for the rows and columns (plots without protrusions)
    spaceforcolumns = remaininghorizontalspace - ((ncols(gl) <= 1) ? 0f0 : sum(addedcolgaps))
    spaceforrows = remainingverticalspace - ((nrows(gl) <= 1) ? 0f0 : sum(addedrowgaps))

    colwidths, rowheights = compute_col_row_sizes(spaceforcolumns, spaceforrows, gl)

    # don't allow smaller widths than 1 px even if it breaks the layout (better than weird glitches)
    colwidths = max.(colwidths, 1f0)
    rowheights = max.(rowheights, 1f0)

    # this is the vertical / horizontal space between the inner lines of all plots
    finalcolgaps = colgaps .+ addedcolgaps
    finalrowgaps = rowgaps .+ addedrowgaps

    # compute the resulting width and height of the gridlayout and compute
    # adjustments for the grid's alignment (this will only matter if the grid is
    # bigger or smaller than the bounding box it occupies)

    gridwidth = sum(colwidths) + sum(finalcolgaps) +
        if alignmode isa Outside
            leftprot + rightprot
        elseif alignmode isa Mixed
            rightal = getside(alignmode, Right())
            leftal = getside(alignmode, Left())
            r = if rightal === nothing
                0f0
            elseif rightal isa Protrusion
                rightal.p
            elseif rightal isa Real
                rightprot
            end
            l = if leftal === nothing
                0f0
            elseif leftal isa Protrusion
                leftal.p
            elseif leftal isa Real
                leftprot
                leftprot
            end
            r + l
        else
            0f0
        end

    gridheight = sum(rowheights) + sum(finalrowgaps) +
        if alignmode isa Outside
            topprot + bottomprot
        elseif alignmode isa Mixed
            bottomal = getside(alignmode, Bottom())
            topal = getside(alignmode, Top())
            b = if bottomal === nothing
                0f0
            elseif bottomal isa Protrusion
                bottomal.p
            elseif bottomal isa Real
                Float32(bottomprot)
            end
            t = if topal === nothing
                0f0
            elseif topal isa Protrusion
                topal.p
            elseif topal isa Real
                Float32(topprot)
            end
            b + t
        else
            0f0
        end
    hal = halign2shift(gl.halign[])
    xadjustment = hal * (width(content_bbox) - gridwidth)

    val = 1 - valign2shift(gl.valign[])
    yadjustment = val * (height(content_bbox) - gridheight)

    # compute the x values for all left and right column boundaries
    xleftcols = if alignmode isa Inside
        xadjustment .+ left(content_bbox) .+ zcumsum(colwidths[1:end-1]) .+
            zcumsum(finalcolgaps)
    elseif alignmode isa Outside
        xadjustment .+ left(content_bbox) .+ zcumsum(colwidths[1:end-1]) .+
            zcumsum(finalcolgaps) .+ leftprot
    elseif alignmode isa Mixed
        leftal = getside(alignmode, Left())
        xadjustment .+ left(content_bbox) .+ zcumsum(colwidths[1:end-1]) .+
            zcumsum(finalcolgaps) .+ (isnothing(leftal) ? zero(leftprot) : isa(leftal, Protrusion) ? leftal.p : leftprot)
    else
        error("Unknown AlignMode of type $(typeof(alignmode))")
    end
    xrightcols = xleftcols .+ colwidths

    # compute the y values for all top and bottom row boundaries
    ytoprows = if alignmode isa Inside
        top(content_bbox) .- yadjustment .- zcumsum(rowheights[1:end-1]) .-
            zcumsum(finalrowgaps)
    elseif alignmode isa Outside
        top(content_bbox) .- yadjustment .- zcumsum(rowheights[1:end-1]) .-
            zcumsum(finalrowgaps) .- topprot
    elseif alignmode isa Mixed
        topal = getside(alignmode, Top())
        top(content_bbox) .- yadjustment .- zcumsum(rowheights[1:end-1]) .-
            zcumsum(finalrowgaps) .- (isnothing(topal) ? zero(topprot) : isa(topal, Protrusion) ? topal.p : topprot)
    else
        error("Unknown AlignMode of type $(typeof(alignmode))")
    end
    ybottomrows = ytoprows .- rowheights

    gridboxes = RowCols(
        xleftcols, xrightcols,
        ytoprows, ybottomrows
    )

    maxgrid, gridboxes
end

"""
    tight_bbox(gl::GridLayout)

Compute the boundingbox that would be the fitting `suggestedbbox`
for the current contents of the `GridLayout`. If a `GridLayout` contains
fixed-size content or aspect-constrained columns, for example, it is
likely that the solved size of the `GridLayout` differs from its
`suggestedbbox`. With the result of `tight_bbox` the `suggestedbbox`
can be adjusted (for example by resizing a figure), so that all
content fits the available space. The size includes possible padding, such
as from an `Outside` alignmode on the `GridLayout`.
"""
function tight_bbox(gl::GridLayout)
    maxgrid, gridboxes = compute_rowcols(gl, suggestedbboxobservable(gl)[])
    base_width = gridboxes.rights[end] - gridboxes.lefts[1]
    base_height = gridboxes.tops[1] - gridboxes.bottoms[end]
    al = gl.alignmode[]
    if al isa Outside
        width = base_width + maxgrid.lefts[1] + maxgrid.rights[end] + al.padding.left + al.padding.right
        height = base_height + maxgrid.tops[1] + maxgrid.bottoms[end] + al.padding.top + al.padding.bottom
        l = gridboxes.lefts[1] - maxgrid.lefts[1] - al.padding.left
        b = gridboxes.bottoms[end] - maxgrid.bottoms[end] - al.padding.bottom
        Rect2f((l, b), (width, height))
    elseif al isa Mixed
        error("Mixed not implemented")
    elseif al isa Inside
        l = gridboxes.lefts[1]
        b = gridboxes.bottoms[end]
        Rect2f((l, b), (base_width, base_height))
    else
        error("Invalid alignmode $al")
    end
end

function align_to_bbox!(gl::GridLayout, suggestedbbox::Rect2f)
    # nesting_level(gl) == 2 && error()
    maxgrid, gridboxes = compute_rowcols(gl, suggestedbbox)

    # now we can solve the content thats inside the grid because we know where each
    # column and row is placed, how wide it is, etc.
    # note that what we did at the top was determine the protrusions of all grid content,
    # but we know the protrusions before we know how much space each plot actually has
    # because the protrusions should be static (like tick labels etc don't change size with the plot)

    for c in gl.content
        idx_rect = side_indices(gl, c)
        bbox_cell = mapsides(idx_rect, gridboxes) do side, idx, gridside
            gridside[idx]
        end

        solving_bbox = bbox_for_solving_from_side(maxgrid, bbox_cell, idx_rect, c.side)

        suggestedbboxobservable(c.content)[] = solving_bbox
    end

    nothing
end

dirlength(gl::GridLayout, c::Col) = ncols(gl)
dirlength(gl::GridLayout, r::Row) = nrows(gl)

function dirgaps(gl::GridLayout, dir::GridDir)
    starts = zeros(Float32, dirlength(gl, dir))
    stops = zeros(Float32, dirlength(gl, dir))
    for c in gl.content
        span = getspan(c, dir)
        start = unoffset(gl, span.start, dir)
        stop = unoffset(gl, span.stop, dir)
        starts[start] = max(starts[start], protrusion(c, startside(dir)))
        stops[stop] = max(stops[stop], protrusion(c, stopside(dir)))
    end
    starts, stops
end

dirsizes(gl::GridLayout, c::Col) = gl.colsizes
dirsizes(gl::GridLayout, r::Row) = gl.rowsizes

"""
Determine the size of a grid layout along one of its dimensions.
`Row` measures from bottom to top and `Col` from left to right.
The size is dependent on the alignmode of the grid, `Outside` includes
protrusions and paddings.
"""
function determinedirsize(gl::GridLayout, gdir::GridDir)::Optional{Float32}
    sum_dirsizes = 0f0

    sizes = dirsizes(gl, gdir)

    for idir in 1:dirlength(gl, gdir)
        # width can only be determined for fixed and auto
        sz = sizes[idir]
        idir_offset = offset(gl, idir, gdir)
        dsize = determinedirsize(idir_offset, gl, gdir)

        if isnothing(dsize)
            # early exit if a colsize can not be determined
            return nothing
        end
        sum_dirsizes += dsize
    end

    dirgapsstart, dirgapsstop = dirgaps(gl, gdir)

    forceequalprotrusiongaps = gl.equalprotrusiongaps[gdir isa Row ? 1 : 2]

    dirgapsizes = if forceequalprotrusiongaps
        innergaps = dirgapsstart[2:end] .+ dirgapsstop[1:end-1]
        m = maximum(innergaps)
        innergaps .= m
    else
        innergaps = dirgapsstart[2:end] .+ dirgapsstop[1:end-1]
    end

    inner_gapsizes = dirlength(gl, gdir) > 1 ? sum(dirgapsizes) : 0

    addeddirgapsizes = gdir isa Row ? gl.addedrowgaps : gl.addedcolgaps

    addeddirgaps = dirlength(gl, gdir) == 1 ? 0 : sum(addeddirgapsizes) do c
        if c isa Fixed
            c.x
        elseif c isa Relative
            error("Auto grid size not implemented with relative gaps")
        end
    end

    inner_size_combined = sum_dirsizes + inner_gapsizes + addeddirgaps
    return if gl.alignmode[] isa Inside
        inner_size_combined
    elseif gl.alignmode[] isa Outside
        paddings = if gdir isa Row
            gl.alignmode[].padding.top + gl.alignmode[].padding.bottom
        else
            gl.alignmode[].padding.left + gl.alignmode[].padding.right
        end
        inner_size_combined + dirgapsstart[1] + dirgapsstop[end] + paddings
    else
        error("Unknown AlignMode of type $(typeof(gl.alignmode[]))")
    end
end


"""
Determine the size of one row or column of a grid layout.
`idir` is the dir index including offset (so can be negative)
"""
function determinedirsize(idir::Integer, gl::GridLayout, dir::GridDir)::Optional{Float32}

    sz = dirsizes(gl, dir)[unoffset(gl, idir, dir)]

    if sz isa Fixed
        # fixed dir size can simply be returned
        return sz.x
    elseif sz isa Relative
        # relative dir size can't be inferred
        return nothing
    elseif sz isa Auto
        # auto dir size can either be determined or not, depending on the
        # trydetermine flag
        !sz.trydetermine && return nothing

        dirsize = nothing
        for c in gl.content
            # content has to be single span to be determinable in size
            if dir isa Row
                singlespanned = c.span.rows.start == c.span.rows.stop == idir
            elseif dir isa Col
                singlespanned = c.span.cols.start == c.span.cols.stop == idir
            end

            # content has to be placed with Inner side, otherwise it's protrusion
            # content
            is_inner = c.side isa Inner

            if singlespanned && is_inner
                s = determinedirsize(c, dir, c.side)
                if !isnothing(s)
                    dirsize = isnothing(dirsize) ? s : max(dirsize, s)
                end
            end
        end
        return dirsize
    end
    nothing
end

# a function that iterates over those sizes that belong to a type T
# while enumerating all indices, so that i can be used to index colwidths / rowheights
# and determinedcols / determinedrows
function filterenum(f, T::Type, iter)
    for (i, value) in enumerate(iter)
        value isa T && f((i, value))
    end
    return
end


function compute_col_row_sizes(spaceforcolumns, spaceforrows, gl)::Tuple{Vector{Float32}, Vector{Float32}}
    # the space for columns and for rows is divided depending on the sizes
    # stored in the grid layout

    # algorithm:

    # 1. get fixed sizes
    # 2. compute relative sizes
    # 3. determine determinable auto sizes
    # 4. compute those aspect sizes that are relative to one of the three above categories
    # 5. at least one side now has to have only undeterminable auto sizes left
    # 6. compute remaining auto sizes for one side
    # 7. compute remaining aspect sizes on other side
    # 8. compute remaining auto sizes on the same side

    colwidths = zeros(ncols(gl))
    rowheights = zeros(nrows(gl))

    determinedcols = zeros(Bool, ncols(gl))
    determinedrows = zeros(Bool, nrows(gl))

    # first fixed sizes
    filterenum(Fixed, gl.colsizes) do (i, fixed)
        colwidths[i] = fixed.x
        determinedcols[i] = true
    end
    filterenum(Fixed, gl.rowsizes) do (i, fixed)
        rowheights[i] = fixed.x
        determinedrows[i] = true
    end

    # then relative sizes
    filterenum(Relative, gl.colsizes) do (i, relative)
        colwidths[i] = relative.x * spaceforcolumns
        determinedcols[i] = true
    end
    filterenum(Relative, gl.rowsizes) do (i, relative)
        rowheights[i] = relative.x * spaceforrows
        determinedrows[i] = true
    end

    # then determinable auto sizes
    filterenum(Auto, gl.colsizes) do (i, auto)
        i_offset = offset(gl, i, Col())
        size = determinedirsize(i_offset, gl, Col())
        if !isnothing(size)
            colwidths[i] = size
            determinedcols[i] = true
        end
    end
    filterenum(Auto, gl.rowsizes) do (i, auto)
        i_offset = offset(gl, i, Row())
        size = determinedirsize(i_offset, gl, Row())
        if !isnothing(size)
            rowheights[i] = size
            determinedrows[i] = true
        end
    end

    # now aspect sizes that refer to already determined counterparts
    filterenum(Aspect, gl.colsizes) do (i, aspect)
        aspectindex_unoffset = unoffset(gl, aspect.index, Row())
        if determinedrows[aspectindex_unoffset]
            colwidths[i] = aspect.ratio * rowheights[aspectindex_unoffset]
            determinedcols[i] = true
        end
    end
    filterenum(Aspect, gl.rowsizes) do (i, aspect)
        aspectindex_unoffset = unoffset(gl, aspect.index, Col())
        if determinedcols[aspectindex_unoffset]
            rowheights[i] = aspect.ratio * colwidths[aspectindex_unoffset]
            determinedrows[i] = true
        end
    end

    remaining_colspace = spaceforcolumns - sum(colwidths)
    remaining_rowspace = spaceforrows - sum(rowheights)

    # if we have aspect sizes left on one side, they can only be determined
    # if the other side has only undeterminable autos left
    n_col_aspects_left = sum(enumerate(gl.colsizes)) do (i, size)
        (size isa Aspect) && (determinedcols[i] == false)
    end
    n_row_aspects_left = sum(enumerate(gl.rowsizes)) do (i, size)
        (size isa Aspect) && (determinedrows[i] == false)
    end

    n_col_autos_left = sum(enumerate(gl.colsizes)) do (i, size)
        (size isa Auto) && (determinedcols[i] == false)
    end
    n_row_autos_left = sum(enumerate(gl.rowsizes)) do (i, size)
        (size isa Auto) && (determinedrows[i] == false)
    end

    if n_col_aspects_left == 0
        let
            indices = Int[]
            ratios = Float64[]
            i_ratios = filterenum(Auto, gl.colsizes) do (i, auto)
                if determinedcols[i] == false
                    push!(indices, i)
                    push!(ratios, auto.ratio)
                end
            end
            sumratios = sum(ratios)
            for (i, ratio) in zip(indices, ratios)
                colwidths[i] = ratio / sumratios * remaining_colspace
                determinedcols[i] = true
            end
        end
    end

    if n_row_aspects_left == 0
        let
            indices = Int[]
            ratios = Float64[]
            i_ratios = filterenum(Auto, gl.rowsizes) do (i, auto)
                if determinedrows[i] == false
                    push!(indices, i)
                    push!(ratios, auto.ratio)
                end
            end
            sumratios = sum(ratios)
            for (i, ratio) in zip(indices, ratios)
                rowheights[i] = ratio / sumratios * remaining_rowspace
                determinedrows[i] = true
            end
        end
    end

    # now if either columns or rows had no aspects left, they should have all sizes determined
    # we run over the aspects again
    filterenum(Aspect, gl.colsizes) do (i, aspect)
        aspectindex_unoffset = unoffset(gl, aspect.index, Row())
        if determinedrows[aspectindex_unoffset]
            colwidths[i] = aspect.ratio * rowheights[aspectindex_unoffset]
            determinedcols[i] = true
        else
            error("Column $i was given an Aspect size relative to row $(aspect.index). This row's size could not be determined in time, therefore the layouting algorithm failed. This probably happened because you used an Aspect row and column size at the same time, which couldn't both be resolved.")
        end
    end
    filterenum(Aspect, gl.rowsizes) do (i, aspect)
        aspectindex_unoffset = unoffset(gl, aspect.index, Col())
        if determinedcols[aspectindex_unoffset]
            rowheights[i] = aspect.ratio * colwidths[aspectindex_unoffset]
            determinedrows[i] = true
        else
            error("Row $i was given an Aspect size relative to column $(aspect.index). This column's size could not be determined in time, therefore the layouting algorithm failed. This probably happened because you used an Aspect row and column size at the same time, which couldn't both be resolved.")
        end
    end

    # if we haven't errored yet, all aspect sizes are done
    # one more pass over the undetermined autos is all that's needed

    remaining_colspace = spaceforcolumns - sum(colwidths)
    remaining_rowspace = spaceforrows - sum(rowheights)

    let
        indices = Int[]
        ratios = Float64[]
        i_ratios = filterenum(Auto, gl.colsizes) do (i, auto)
            if determinedcols[i] == false
                push!(indices, i)
                push!(ratios, auto.ratio)
            end
        end
        sumratios = sum(ratios)
        for (i, ratio) in zip(indices, ratios)
            colwidths[i] = ratio / sumratios * remaining_colspace
            determinedcols[i] = true
        end
    end

    let
        indices = Int[]
        ratios = Float64[]
        i_ratios = filterenum(Auto, gl.rowsizes) do (i, auto)
            if determinedrows[i] == false
                push!(indices, i)
                push!(ratios, auto.ratio)
            end
        end
        sumratios = sum(ratios)
        for (i, ratio) in zip(indices, ratios)
            rowheights[i] = ratio / sumratios * remaining_rowspace
            determinedrows[i] = true
        end
    end


    # now all columns and rows should have their sizes
    ncols_undetermined = sum(.!determinedcols)
    nrows_undetermined = sum(.!determinedrows)

    if ncols_undetermined > 0
        error("After a non-erroring layouting pass, the number of undetermined columns is $ncols_undetermined. This must be a bug.")
    end
    if nrows_undetermined > 0
        error("After a non-erroring layouting pass, the number of undetermined rows is $nrows_undetermined. This must be a bug.")
    end

    colwidths, rowheights
end

function Base.setindex!(g::GridLayout, content, rows::Indexables, cols::Indexables, side::Side = Inner())
    add_content!(g, content, rows, cols, side)
    content
end

function Base.setindex!(g::GridLayout, content_array::AbstractArray{T, 2}) where T
    rowrange = 1:size(content_array, 1)
    colrange = 1:size(content_array, 2)
    g[rowrange, colrange] = content_array
end

function Base.setindex!(g::GridLayout, content_array::AbstractArray{T, 1}) where T
    error("""
        You can only assign a one-dimensional content AbstractArray if you also specify the direction in the layout.
        Valid options are :h for horizontal and :v for vertical.
        Example:
            layout[:h] = contentvector
    """)
end

function Base.setindex!(g::GridLayout, content_array::AbstractArray{T, 1}, h_or_v::Symbol) where T
    if h_or_v == :h
        g[1, 1:length(content_array)] = content_array
    elseif h_or_v == :v
        g[1:length(content_array), 1] = content_array
    else
        error("""
            Invalid direction specifier $h_or_v.
            Valid options are :h for horizontal and :v for vertical.
        """)
    end
end

function Base.setindex!(g::GridLayout, content_array::AbstractArray, rows::Indexables, cols::Indexables)

    rows, cols = to_ranges(g, rows, cols)

    if rows.start < 1
        error("Can't prepend rows using array syntax so far, start row $(rows.start) is smaller than 1.")
    end
    if cols.start < 1
        error("Can't prepend columns using array syntax so far, start column $(cols.start) is smaller than 1.")
    end

    nrows = length(rows)
    ncols = length(cols)
    ncells = nrows * ncols

    if ndims(content_array) == 2
        if size(content_array) != (nrows, ncols)
            error("Content array size is size $(size(content_array)) for $nrows rows and $ncols cols")
        end
        # put the array content into the grid layout in order
        for (i, r) in enumerate(rows), (j, c) in enumerate(cols)
            g[r, c] = content_array[i, j]
        end
    elseif ndims(content_array) == 1
        if length(content_array) != nrows * ncols
            error("Content array size is length $(length(content_array)) for $nrows * $ncols cells")
        end
        # put the content in the layout along columns first, because that is more
        # intuitive
        for (i, (c, r)) in enumerate(Iterators.product(cols, rows))
            g[r, c] = content_array[i]
        end
    else
        error("Can't assign a content array with $(ndims(content_array)) dimensions, only 1 or 2.")
    end
    content_array
end

function GridContent(content, span::Span, side::Side)
    # connect the correct observables
    gc = GridContent{GridLayout}(nothing, content, span, side, nothing, nothing)
    connect_layoutobservables!(gc)
    gc
end

function update!(gc::GridContent)
    p = gc.parent
    if p isa GridLayout
        update!(p)
    end
    return
end

function add_content!(g::GridLayout, content, rows, cols, side::Side)
    # update = false because update is called in add_to_gridlayout! anyway
    rows, cols = adjust_rows_cols!(g, rows, cols; update = false)

    gc = if !isnothing(gridcontent(content))
        # take the existing gridcontent, remove it from its gridlayout if it has one,
        # and modify it with the new span and side
        gridc = gridcontent(content)
        remove_from_gridlayout!(gridc)
        gridc.span = Span(rows, cols)
        gridc.side = side
        gridc
    else
        # make a new one if none existed
        GridContent(content, Span(rows, cols), side)
    end

    layoutobservables(content).gridcontent[] = gc

    connect_layoutobservables!(gc)

    add_to_gridlayout!(g, gc)
end

function Base.lastindex(g::GridLayout, d)
    if d == 1
        lastrow(g)
    elseif d == 2
        lastcol(g)
    else
        error("A grid only has two dimensions, you're indexing dimension $d.")
    end
end

function Base.firstindex(g::GridLayout, d)
    if d == 1
        firstrow(g)
    elseif d == 2
        firstcol(g)
    else
        error("A grid only has two dimensions, you're indexing dimension $d.")
    end
end

function GridPosition(g::GridLayout, rows::Indexables, cols::Indexables, side = Inner())
    span = Span(to_ranges(g, rows, cols)...)
    GridPosition(g, span, side)
end

function Base.getindex(g::GridLayout, rows::Indexables, cols::Indexables, side = Inner())
    GridPosition(g, rows, cols, side)
end

function Base.setindex!(gp::GridPosition, element)
    gp.layout[gp.span.rows, gp.span.cols, gp.side] = element
end

function Base.setindex!(gp::GridPosition, element, rows, cols, side = Inner())
    layout = get_layout_at!(gp, createmissing = true)
    layout[rows, cols, side] = element
    element
end

"""
    nrows(g::GridLayout)

Return the number of rows in `g`.
"""
nrows(g::GridLayout) = size(g)[1]
"""
    ncols(g::GridLayout)

Return the number of columns in `g`.
"""
ncols(g::GridLayout) = size(g)[2]
Base.size(g::GridLayout) = g.size
offsets(g::GridLayout) = g.offsets

Base.in(span1::Span, span2::Span) = span1.rows.start >= span2.rows.start &&
    span1.rows.stop <= span2.rows.stop &&
    span1.cols.start >= span2.cols.start &&
    span1.cols.stop <= span2.cols.stop

"""
    contents(gp::GridPosition; exact::Bool = false)

Retrieve all objects placed in the `GridLayout` at the `Span` and `Side` stored
in the `GridPosition` `gp`. If `exact == true`, elements are only included
if they match the `Span` exactly, otherwise they can also be contained within the spanned layout area.
"""
function contents(gp::GridPosition; exact::Bool = false)
    contents = []
    for c in gp.layout.content
        if exact
            if c.span == gp.span && c.side == gp.side
                push!(contents, c.content)
            end
        else
            if c.span in gp.span && c.side == gp.side
                push!(contents, c.content)
            end
        end
    end
    contents
end

"""
    contents(g::GridLayout)

Retrieve all objects placed in the `GridLayout` `g`, in the order they are stored, extracted from
their containing `GridContent`s.
"""
function contents(g::GridLayout)
    map(g.content) do gc
        gc.content
    end
end


function Base.getindex(gp::Union{GridPosition, GridSubposition}, rows, cols, side = Inner())
    GridSubposition(gp, rows, cols, side)
end

function Base.setindex!(parent::GridSubposition, obj,
    rows, cols, side = GridLayoutBase.Inner())
    layout = get_layout_at!(parent, createmissing = true)
    layout[rows, cols, side] = obj
    obj
end

function Base.setindex!(parent::GridSubposition, obj)
    layout = get_layout_at!(parent.parent, createmissing = true)
    layout[parent.rows, parent.cols, parent.side] = obj
    obj
end

function get_layout_at!(gp::GridPosition; createmissing = false)
    c = contents(gp, exact = true)
    layouts = filter(x -> x isa GridLayoutBase.GridLayout, c)
    if isempty(layouts)
        if createmissing
            return gp[] = GridLayoutBase.GridLayout()
        else
            error("No layout found but `createmissing` is false.")
        end
    elseif length(layouts) == 1
        return first(layouts)
    else
        error("Found more than zero or one GridLayouts at $gp")
    end
end

function get_layout_at!(gsp::GridSubposition; createmissing = false)
    layout = get_layout_at!(gsp.parent; createmissing = createmissing)
    gp = layout[gsp.rows, gsp.cols, gsp.side]
    get_layout_at!(gp, createmissing = createmissing)
end


function contents(g::GridSubposition; exact = false)
    layout = get_layout_at!(g.parent, createmissing = false)
    contents(layout[g.rows, g.cols, g.side], exact = exact)
end

"""
    content(g::Union{GridPosition,GridSubposition})

Return the one object placed in the `GridLayout` at the `Span` and `Side`
stored in the `GridPosition` `g`. If there is more than one object at that
position, throw an error.

See also `contents`.
"""
function content(g::Union{GridPosition,GridSubposition})
    cs = contents(g, exact = true)
    if length(cs) == 1
        return cs[1]
    else
        error("There is not exactly one object at the given GridPosition")
    end
end


function parent(g::GridLayout)
    g.parent
end

function top_parent(g::GridLayout)
    top_parent(parent(g))
end

top_parent(x) = x

function top_parent_grid(g::GridLayout)
    p = parent(g)
    if p isa GridLayout
        top_parent_grid(p)
    else
        g
    end
end

"""
    side_indices(gl::GridLayout, c::GridContent)::RowCols{Int}

Indices of the rows / cols for each side
"""
function side_indices(gl, c::GridContent)
    return RowCols(
        c.span.cols.start - offset(gl, Col()),
        c.span.cols.stop - offset(gl, Col()),
        c.span.rows.start - offset(gl, Row()),
        c.span.rows.stop - offset(gl, Row()),
    )
end

# These functions tell whether an object in a grid touches the left, top, etc. border
# of the grid. This means that it is relevant for the grid's own protrusion on that side.
ismostin(gc::GridContent, grid, ::Left) = gc.span.cols.start == firstcol(grid)
ismostin(gc::GridContent, grid, ::Right) = gc.span.cols.stop == lastcol(grid)
ismostin(gc::GridContent, grid, ::Bottom) = gc.span.rows.stop == lastrow(grid)
ismostin(gc::GridContent, grid, ::Top) = gc.span.rows.start == firstrow(grid)


function protrusion(x::T, side::Side) where T
    protrusions = protrusionsobservable(x)[]::RectSides{Float32}
    if side isa Left
        protrusions.left
    elseif side isa Right
        protrusions.right
    elseif side isa Bottom
        protrusions.bottom
    elseif side isa Top
        protrusions.top
    else
        error("Can't get a protrusion value for side $(typeof(side)), only
        Left, Right, Bottom, or Top.")
    end
end

function protrusion(gc::GridContent, side::Side)
    prot =
        if gc.side isa Inner
            protrusion(gc.content, side)
        elseif gc.side isa Outer
            0.0
        elseif gc.side isa Union{Left, Right}
            if side isa typeof(gc.side)
                determinedirsize(gc.content, Col(), gc.side)
            else
                0.0
            end
        elseif gc.side isa Union{Top, Bottom}
            if side isa typeof(gc.side)
                determinedirsize(gc.content, Row(), gc.side)
            else
                0.0
            end
        elseif gc.side isa TopLeft
            if side isa Top
                determinedirsize(gc.content, Row(), gc.side)
            elseif side isa Left
                determinedirsize(gc.content, Col(), gc.side)
            else
                0.0
            end
        elseif gc.side isa TopRight
            if side isa Top
                determinedirsize(gc.content, Row(), gc.side)
            elseif side isa Right
                determinedirsize(gc.content, Col(), gc.side)
            else
                0.0
            end
        elseif gc.side isa BottomLeft
            if side isa Bottom
                determinedirsize(gc.content, Row(), gc.side)
            elseif side isa Left
                determinedirsize(gc.content, Col(), gc.side)
            else
                0.0
            end
        elseif gc.side isa BottomRight
            if side isa Bottom
                determinedirsize(gc.content, Row(), gc.side)
            elseif side isa Right
                determinedirsize(gc.content, Col(), gc.side)
            else
                0.0
            end
        else
            error("Invalid side $(gc.side)")
        end
    ifnothing(prot, 0.0)
end

getside(m::Mixed, ::Left) = m.sides.left
getside(m::Mixed, ::Right) = m.sides.right
getside(m::Mixed, ::Top) = m.sides.top
getside(m::Mixed, ::Bottom) = m.sides.bottom

function compute_effective_protrusion(gl::GridLayout, side::Side)
    al = gl.alignmode[]
    if al isa Outside
        return 0f0
    elseif al isa Inside
        compute_effective_protrusion_inside(gl, side)
    elseif al isa Mixed
        si = getside(al, side)
        if isnothing(si)
            compute_effective_protrusion_inside(gl, side)
        elseif si isa Protrusion
            si.p
        else
            # Outside alignment
            0f0
        end
    else
        error("Unknown AlignMode of type $(typeof(al))")
    end
end

function compute_effective_protrusion_inside(gl::GridLayout, side::Side)
    prot = 0f0
    for elem in gl.content
        if ismostin(elem, gl, side)
            # take the max protrusion of all elements that are sticking
            # out at this side
            prot = max(effective_protrusion(elem, side, elem.side), prot)
        end
    end
    return prot
end

function inside_protrusion(gl::GridLayout, side::Side)
    prot = 0.0
    for elem in gl.content
        if ismostin(elem, gl, side)
            # take the max protrusion of all elements that are sticking
            # out at this side
            prot = max(protrusion(elem, side), prot)
        end
    end
    return prot
end

function protrusion(gl::GridLayout, side::Side)::Float32
    # when we align with the outside there is by definition no protrusion

    al = gl.alignmode[]
    if al isa Outside
        return 0f0
    elseif al isa Inside
        inside_protrusion(gl, side)
    elseif al isa Mixed
        si = getside(al, side)
        if isnothing(si)
            inside_protrusion(gl, side)
        elseif si isa Protrusion
            si.p
        else
            # Outside alignment
            0f0
        end
    else
        error("Unknown AlignMode of type $(typeof(al))")
    end
end

function bbox_for_solving_from_side(maxgrid::RowCols, bbox_cell::Rect2f, idx_rect::RowCols, side::Side)
    pl = maxgrid.lefts[idx_rect.lefts]
    pr = maxgrid.rights[idx_rect.rights]
    pt = maxgrid.tops[idx_rect.tops]
    pb = maxgrid.bottoms[idx_rect.bottoms]

    l = left(bbox_cell)
    r = right(bbox_cell)
    b = bottom(bbox_cell)
    t = top(bbox_cell)

    if side isa Inner
        bbox_cell
    elseif side isa Outer
        BBox(l - pl, r + pr, b - pb, t + pt)
    elseif side isa Left
        BBox(l - pl, l, b, t)
    elseif side isa Top
        BBox(l, r, t, t + pt)
    elseif side isa Right
        BBox(r, r + pr, b, t)
    elseif side isa Bottom
        BBox(l, r, b - pb, b)
    elseif side isa TopLeft
        BBox(l - pl, l, t, t + pt)
    elseif side isa TopRight
        BBox(r, r + pr, t, t + pt)
    elseif side isa BottomRight
        BBox(r, r + pr, b - pb, b)
    elseif side isa BottomLeft
        BBox(l - pl, l, b - pb, b)
    else
        error("Invalid side $side")
    end
end

startside(c::Col) = Left()
stopside(c::Col) = Right()
startside(r::Row) = Top()
stopside(r::Row) = Bottom()


getspan(gc::GridContent, dir::Col) = gc.span.cols
getspan(gc::GridContent, dir::Row) = gc.span.rows



"""
Determine the size of a protrusion layout along a dimension. This size is dependent
on the `Side` at which the layout is placed in its parent grid. An `Inside` side
means that the protrusion layout reports its width but not its protrusions. `Left`
means that the layout reports only its full width but not its height, because
an element placed in the left protrusion loses its ability to influence height.
"""
function determinedirsize(content, gdir::GridDir, side::Side)
    reportedsize = reporteddimensionsobservable(content)[].inner
    if gdir isa Row
        if side isa Union{Inner, Top, Bottom, TopLeft, TopRight, BottomLeft, BottomRight}
            # TODO: is reportedsize the correct thing to return? or plus protrusions depending on the side
            ifnothing(reportedsize[2], nothing)
        elseif side isa Union{Left, Right}
            nothing
        else
            error("$side not implemented")
        end
    else
        if side isa Union{Inner, Left, Right, TopLeft, TopRight, BottomLeft, BottomRight}
            ifnothing(reportedsize[1], nothing)
        elseif side isa Union{Top, Bottom}
            nothing
        else
            error("$side not implemented")
        end
    end
end

determinedirsize(gc::GridContent, gdir::GridDir, side::Side) = determinedirsize(gc.content, gdir, side)::Optional{Float32}

function to_ranges(g::GridLayout, rows::Indexables, cols::Indexables)
    if rows isa Int
        rows = rows:rows
    elseif rows isa Colon
        rows = firstrow(g):lastrow(g)
    end
    if cols isa Int
        cols = cols:cols
    elseif cols isa Colon
        cols = firstcol(g):lastcol(g)
    end
    rows, cols
end

function adjust_rows_cols!(g::GridLayout, rows, cols; update = true)
    rows, cols = to_ranges(g, rows, cols)

    rowdiff_start = firstrow(g) - rows.start
    if rowdiff_start > 0
        prependrows!(g, rowdiff_start, update = update)
    end
    rowdiff_stop = rows.stop - lastrow(g)
    if rowdiff_stop > 0
        appendrows!(g, rowdiff_stop, update = update)
    end
    coldiff_start = firstcol(g) - cols.start
    if coldiff_start > 0
        prependcols!(g, coldiff_start, update = update)
    end
    coldiff_stop = cols.stop - lastcol(g)

    if coldiff_stop > 0
        appendcols!(g, coldiff_stop, update = update)
    end

    rows, cols
end

firstrow(g) = (1 + offsets(g)[1])
firstcol(g) = (1 + offsets(g)[2])
lastrow(g) = (nrows(g) + offsets(g)[1])
lastcol(g) = (ncols(g) + offsets(g)[2])
