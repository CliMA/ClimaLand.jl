module Packing
using GeometryBasics

include("rectangle.jl")
include("guillotine.jl")

export RectanglePacker, GuillotinePacker

end # module
