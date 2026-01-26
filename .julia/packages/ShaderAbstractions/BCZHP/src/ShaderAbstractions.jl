module ShaderAbstractions

using ColorTypes, FixedPointNumbers, Observables
import GeometryBasics
using GeometryBasics: StaticVector, Mat

include("types.jl")
include("uniforms.jl")
include("context.jl")
include("program.jl")

end # module
