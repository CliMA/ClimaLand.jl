module PolygonOps

export inpolygon, HormannAgathos, HaoSun

include("validity_checks.jl")
include("inpolygon.jl")
include("area.jl")

end # module
