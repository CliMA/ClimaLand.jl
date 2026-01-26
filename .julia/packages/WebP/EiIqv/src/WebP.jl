module WebP

using ColorTypes
using FileIO
using FixedPointNumbers
using ImageCore

include(joinpath(@__DIR__, "Wrapper.jl"))
using .Wrapper

include(joinpath(@__DIR__, "decoding.jl"))
include(joinpath(@__DIR__, "encoding.jl"))
include(joinpath(@__DIR__, "fileio_interface.jl"))

end
