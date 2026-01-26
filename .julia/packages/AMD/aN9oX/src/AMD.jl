module AMD

using LinearAlgebra
using SparseArrays
using SuiteSparse_jll

import Base.show, Base.print

const SS_Int = Base.Sys.WORD_SIZE == 32 ? Int32 : Int64

include("wrappers/amd.jl")
include("wrappers/camd.jl")
include("wrappers/colamd.jl")
include("wrappers/ccolamd.jl")

include("amd_julia.jl")
include("colamd_julia.jl")

end # module
