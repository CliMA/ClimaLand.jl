mutable struct DebugRect
    layoutobservables::GridLayoutBase.LayoutObservables{GridLayout}
    width::Observable
    height::Observable
    tellwidth::Observable
    tellheight::Observable
    halign::Observable
    valign::Observable
    leftprot::Observable
    rightprot::Observable
    bottomprot::Observable
    topprot::Observable
    alignmode::Observable
end


observablify(x::Observable) = x
observablify(x, type=Any) = Observable{type}(x)

function DebugRect(; bbox = nothing, width=nothing, height=nothing,
    tellwidth = true, tellheight = true, halign=:center,
    valign=:center, topprot=0.0, leftprot=0.0, rightprot=0.0, bottomprot=0.0,
    alignmode = Inside())

    width = observablify(width)
    height = observablify(height)
    tellwidth = observablify(tellwidth)
    tellheight = observablify(tellheight)
    halign = observablify(halign)
    valign = observablify(valign)
    topprot = observablify(topprot, Float32)
    leftprot = observablify(leftprot, Float32)
    rightprot = observablify(rightprot, Float32)
    bottomprot = observablify(bottomprot, Float32)
    alignmode = observablify(alignmode)

    protrusions::Observable{GridLayoutBase.RectSides{Float32}} = map(leftprot, rightprot, bottomprot, topprot) do l, r, b, t
        GridLayoutBase.RectSides{Float32}(l, r, b, t)
    end

    layoutobservables = GridLayoutBase.LayoutObservables(width,
        height, tellwidth, tellheight, halign, valign, alignmode;
        suggestedbbox = bbox, protrusions = protrusions)

    DebugRect(layoutobservables, height, width, tellwidth, tellheight,
        halign, valign, leftprot, rightprot, bottomprot, topprot, alignmode)
end

Base.show(io::IO, dr::DebugRect) = print(io, "DebugRect()")
