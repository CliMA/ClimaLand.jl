module GridLayoutBase

using GeometryBasics
using Observables
import GeometryBasics: height, width

const DEFAULT_COLGAP = Ref{Float64}(20.0)
const DEFAULT_ROWGAP = Ref{Float64}(20.0)
# These function refs can be mutated by other packages to override the default
# way of retrieving default column and row gap sizes
const DEFAULT_ROWGAP_GETTER = Ref{Function}(() -> DEFAULT_ROWGAP[])
const DEFAULT_COLGAP_GETTER = Ref{Function}(() -> DEFAULT_COLGAP[])

include("types.jl")
include("gridlayout.jl")
include("layoutobservables.jl")
include("helpers.jl")

export GridLayout, GridPosition
export BBox
export LayoutObservables
export Inside, Outside, Mixed, Protrusion
export Fixed, Auto, Relative, Aspect
export width, height, top, bottom, left, right
export with_updates_suspended
export appendcols!, appendrows!, prependcols!, prependrows!, deletecol!, deleterow!, trim!, insertrows!, insertcols!
export gridnest!
export AxisAspect, DataAspect
export colsize!, rowsize!, colgap!, rowgap!
export Left, Right, Top, Bottom, TopLeft, BottomLeft, TopRight, BottomRight
export grid!, hbox!, vbox!
export swap!
export protrusionsobservable, suggestedbboxobservable, reporteddimensionsobservable, autosizeobservable, computedbboxobservable, gridcontent
export ncols, nrows, offsets
export contents, content
export tight_bbox

if Base.VERSION >= v"1.4.2"
    include("precompile.jl")
end

end
