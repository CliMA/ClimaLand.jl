halign2shift(align::HorizontalAlignment)::Float32 = Float32(align.x)
halign2shift(align::Number)::Float32 = Float32(align)
function halign2shift(align::Symbol)::Float32
    align == :left && return 0.0f0
    align == :center && return 0.5f0
    align == :right && return 1.0f0
    error("Invalid horizontal alignment $align (only Real or :left, :center, or :right allowed).")
end

valign2shift(align::VerticalAlignment)::Float32 = Float32(align.x)
valign2shift(align::Number)::Float32 = Float32(align)
function valign2shift(align::Symbol)::Float32
    align == :bottom && return 0.0f0
    align == :center && return 0.5f0
    align == :top && return 1.0f0
    error("Invalid vertical alignment $align (only Real or :bottom, :center, or :top allowed).")
end

function align_shift_tuple(halign, valign)
    return (halign2shift(halign), valign2shift(valign))
end

function LayoutObservables(width::Observable, height::Observable,
        tellwidth::Observable, tellheight::Observable, halign::Observable,
        valign::Observable, alignmode::Observable = Observable{AlignMode}(Inside());
        suggestedbbox = nothing,
        protrusions = nothing,
        gridcontent = nothing,
        block_updates::Bool = false)

    width_obs = convert(Observable{SizeAttribute}, width)
    height_obs = convert(Observable{SizeAttribute}, height)
    sizeobservable = sizeobservable!(width_obs, height_obs)
    alignment = map(align_shift_tuple, halign, valign)

    suggestedbbox_observable = make_suggestedbbox!(suggestedbbox)
    protrusions = make_protrusions!(protrusions)

    tellsizeobservable = map(tuple, tellwidth, tellheight)

    gridcontent_ref = Ref{Optional{GridContent{GridLayout}}}(gridcontent)

    autosizeobservable = Observable{NTuple{2, Optional{Float32}}}((nothing, nothing))
    reporteddimensions = make_reporteddimensions!(sizeobservable, autosizeobservable, tellsizeobservable, protrusions, alignmode)
    computedbbox = make_computedbbox!(suggestedbbox_observable, reporteddimensions, alignment, sizeobservable, autosizeobservable,
        alignmode, protrusions, gridcontent_ref)

    LayoutObservables{GridLayout}(
        suggestedbbox_observable,
        protrusions,
        reporteddimensions,
        autosizeobservable,
        computedbbox,
        gridcontent_ref,
        Ref{Bool}(block_updates),
    )
end

maprectsides(f) = RectSides(map(f, (:left, :right, :bottom, :top))...)

function effective_protrusion(prot::RectSides, @nospecialize(al::AlignMode))
    x = if al isa Inside
        prot
    elseif al isa Outside
        RectSides{Float32}(0, 0, 0, 0)
    elseif al isa Mixed
        maprectsides() do side
            # normal inside mode
            if isnothing(getfield(al.sides, side))
                getfield(prot, side)
            # protrusion override
            elseif getfield(al.sides, side) isa Protrusion
                getfield(al.sides, side).p
            # outside mode
            else
                0f0
            end
        end
    end
end

make_suggestedbbox!(n::Nothing) = Observable(BBox(0, 100, 0, 100))
make_suggestedbbox!(tup::Tuple) = Observable(BBox(tup...))
make_suggestedbbox!(bbox::Rect{2}) = Observable(Rect2f(bbox))
make_suggestedbbox!(observable::Observable{Rect2f}) = observable
function make_suggestedbbox!(observable::Observable{<:Rect{2}})
    bbox = Observable(Rect2f(observable[]))
    on(observable) do o
        bbox[] = Rect2f(o)
    end
    bbox
end

make_protrusions!(p::Nothing) = Observable(RectSides{Float32}(0, 0, 0, 0))
make_protrusions!(p::Observable{RectSides{Float32}}) = p
make_protrusions!(p::RectSides{Float32}) = Observable(p)

function sizeobservable!(widthattr::Observable{SizeAttribute}, heightattr::Observable{SizeAttribute})
    sizeattrs = Observable{Tuple{SizeAttribute, SizeAttribute}}((widthattr[], heightattr[]))
    onany(widthattr, heightattr) do w, h
        sizeattrs[] = (w, h)
    end
    sizeattrs
end

function make_reporteddimensions!(sizeattrs, autosizeobservable::Observable{NTuple{2, Optional{Float32}}}, tellsizeobservable, protrusions, alignmode)

    # set up rsizeobservable with correct type manually
    rdimobservable = Observable{Dimensions}(Dimensions((nothing, nothing), RectSides(0f0, 0f0, 0f0, 0f0)))

    map!(_reporteddimensionsobservable, rdimobservable, sizeattrs, autosizeobservable, tellsizeobservable, protrusions, alignmode)

    # # trigger first value
    # notify(sizeattrs)

    rdimobservable
end

function _reporteddimensionsobservable(
        @nospecialize(sizeattrs::Tuple{SizeAttribute,SizeAttribute}),
        @nospecialize(autosize::Tuple{AutoSize,AutoSize}),
        tellsizeobservable::Tuple{Bool,Bool},
        protrusions,
        alignmode,
    )

    outer = effective_protrusion(protrusions, alignmode)

    wattr, hattr = sizeattrs
    wauto, hauto = autosize
    tellw, tellh = tellsizeobservable

    wsize = computed_size(wattr, wauto, tellw)
    hsize = computed_size(hattr, hauto, tellh)

    inner = if alignmode isa Inside
        (wsize, hsize)
    elseif alignmode isa Outside
        (isnothing(wsize) ? nothing : wsize + protrusions.left + protrusions.right + alignmode.padding.left + alignmode.padding.right,
         isnothing(hsize) ? nothing : hsize + protrusions.top + protrusions.bottom + alignmode.padding.top + alignmode.padding.bottom)
    elseif alignmode isa Mixed
        w = if isnothing(wsize)
            nothing
        else
            w = wsize
            if alignmode.sides.left isa Float32
                w += protrusions.left + alignmode.sides.left
            elseif alignmode.sides.left isa Protrusion
                w
            end
            if alignmode.sides.right isa Float32
                w += protrusions.right + alignmode.sides.right
            elseif alignmode.sides.right isa Protrusion
                w
            end
            w
        end
        h = if isnothing(hsize)
            nothing
        else
            h = hsize
            if alignmode.sides.bottom isa Float32
                h += protrusions.bottom + alignmode.sides.bottom
            elseif alignmode.sides.bottom isa Protrusion
                h
            end
            if alignmode.sides.top isa Float32
                h += protrusions.top + alignmode.sides.top
            elseif alignmode.sides.top isa Protrusion
                h
            end
            h
        end
        (w, h)
    else
        error("Unknown alignmode $alignmode")
    end

    Dimensions(inner, outer)
end

function computed_size(sizeattr, autosize, tellsize)

    if !tellsize
        return nothing
    end

    if sizeattr === nothing
        nothing
    elseif sizeattr isa Real
        Float32(sizeattr)
    elseif sizeattr isa Fixed
        sizeattr.x
    elseif sizeattr isa Relative
        nothing
    elseif sizeattr isa Auto
        autosize
    else
        error("""
            Invalid size attribute $sizeattr.
            Can only be Nothing, Fixed, Relative, Auto or Real""")
    end
end


function make_computedbbox!(
    suggestedbbox::Observable{Rect2f},
    reporteddimensions::Observable{Dimensions},
    alignment::Observable,
    sizeattrs::Observable,
    autosizeobservable::Observable{NTuple{2, Optional{Float32}}},
    alignmode, protrusions, gridcontent_ref)

    computedbbox = Observable(BBox(0, 100, 0, 100))

    # suggestedbbox and alignment don't affect the parent gridlayout
    onany(suggestedbbox, alignment) do sbbox, ali
        update_computedbbox!(computedbbox, sbbox, ali, reporteddimensions[], alignmode[], protrusions[], sizeattrs, autosizeobservable)
    end

    # these will trigger the update of a parent gridlayout if there is one, so an update
    # through suggestedbbox will come back after the parent has updated, so here we only update if
    # the object is standalone
    onany(reporteddimensions) do rdims
        if gridcontent_ref[] === nothing
            update_computedbbox!(computedbbox, suggestedbbox[], alignment[], rdims, alignmode[], protrusions[], sizeattrs, autosizeobservable)
        end
    end

    computedbbox
end

function update_computedbbox!(computedbbox, suggestedbbox, alignment, reporteddimensions, alignmode, protrusions, sizeattrs, autosizeobservable)

    bw = width(suggestedbbox)
    bh = height(suggestedbbox)

    # we only passively retrieve sizeattrs here because if they change
    # they also trigger reportedsize, which triggers this observable, too
    # we only need to know here if there are relative sizes given, because
    # those can only be computed knowing the suggestedbbox
    widthattr, heightattr = sizeattrs[]
    T = eltype(protrusions)

    cwidth, cheight = reporteddimensions.inner
    w_target = T(if isnothing(cwidth)
        if widthattr isa Relative
            widthattr.x * bw
        elseif widthattr isa Nothing
            bw
        elseif widthattr isa Auto
            if isnothing(autosizeobservable[][1])
                # we have no autowidth available anyway
                # take suggested width
                bw
            else
                # use the width that was auto-computed
                autosizeobservable[][1]
            end
        elseif widthattr isa Fixed
            widthattr.x
        elseif widthattr isa Real
            Float32(widthattr)
        else
            error("Unknown width attribute $widthattr")
        end
    else
        cwidth
    end)::T

    h_target = T(if isnothing(cheight)
        if heightattr isa Relative
            heightattr.x * bh
        elseif heightattr isa Nothing
            bh
        elseif heightattr isa Auto
            if isnothing(autosizeobservable[][2])
                # we have no autoheight available anyway
                # take suggested height
                bh
            else
                # use the height that was auto-computed
                autosizeobservable[][2]
            end
        elseif heightattr isa Fixed
            heightattr.x
        elseif heightattr isa Real
            Float32(heightattr)
        else
            error("Unknown height attribute $heightattr")
        end
    else
        cheight
    end)::T

    inner_w, inner_h = if alignmode isa Inside
        (w_target, h_target)
    elseif alignmode isa Outside
        (w_target - protrusions.left - protrusions.right - alignmode.padding.left - alignmode.padding.right,
            h_target - protrusions.top - protrusions.bottom - alignmode.padding.top - alignmode.padding.bottom)
    else
        alignmode = alignmode::Mixed
        let
            w = w_target
            # subtract if outside padding is used via a Float32 value
            # Protrusion and `nothing` are protrusion modes
            if alignmode.sides.left isa Float32
                w -= protrusions.left + alignmode.sides.left
            end
            if alignmode.sides.right isa Float32
                w -= protrusions.right + alignmode.sides.right
            end

            h = h_target
            if alignmode.sides.bottom isa Float32
                h -= protrusions.bottom + alignmode.sides.bottom
            end
            if alignmode.sides.top isa Float32
                h -= protrusions.top + alignmode.sides.top
            end

            w, h
        end
    end

    # how much space is left in the bounding box
    rw = bw - w_target
    rh = bh - h_target

    xshift, yshift = alignment .* (rw, rh)

    if alignmode isa Inside
        # width and height are unaffected
    elseif alignmode isa Outside
        xshift = xshift + protrusions.left + alignmode.padding.left
        yshift = yshift + protrusions.bottom + alignmode.padding.bottom
    else
        am = alignmode::Mixed
        if am.sides.left isa Float32
            xshift += protrusions.left + am.sides.left
        end
        if am.sides.bottom isa Float32
            yshift += protrusions.bottom + am.sides.bottom
        end
    end

    # align the final bounding box in the layout bounding box
    l = left(suggestedbbox) + xshift
    b = bottom(suggestedbbox) + yshift
    r = l + inner_w
    t = b + inner_h
    newbbox = BBox(l, r, b, t)
    # if computedbbox[] != newbbox
    #     computedbbox[] = newbbox
    # end
    computedbbox[] = newbbox

    return
end

"""
    layoutobservables(x::T) where T

Access `x`'s field `:layoutobservables` containing a `LayoutObservables` instance. This should
be overloaded for any type that is layoutable but stores its `LayoutObservables` in
a differently named field.
"""
function layoutobservables(x::T)::LayoutObservables{GridLayout} where T
    if hasfield(T, :layoutobservables) && fieldtype(T, :layoutobservables) === LayoutObservables{GridLayout}
        x.layoutobservables
    else
        error("It's not defined how to get LayoutObservables for type $T, overload this method for layoutable types.")
    end
end

# These are the default API functions to retrieve the layout parts from an object
protrusionsobservable(x) = layoutobservables(x).protrusions
suggestedbboxobservable(x) = layoutobservables(x).suggestedbbox
reporteddimensionsobservable(x) = layoutobservables(x).reporteddimensions
autosizeobservable(x) = layoutobservables(x).autosize
computedbboxobservable(x) = layoutobservables(x).computedbbox
gridcontent(x) = layoutobservables(x).gridcontent[]
