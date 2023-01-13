module parameterization_viz.jl

using WGLMakie
using JSServe
using SparseArrays

include(joinpath(pkgdir(ClimaLSM), "parameters", "create_parameters.jl"))
export create_lsm_parameters

include("fun_discretisation.jl")
export mat, d1_vec, d2_vec

include("generate_fig.jl")
export param_dashboard

end 
