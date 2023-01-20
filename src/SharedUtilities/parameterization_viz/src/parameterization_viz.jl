module parameterization_viz

FT = Float64 # not sure why this is needed, but bug if I don't

using WGLMakie
using JSServe
using SparseArrays
using Statistics

using ClimaLSM # needed because of create_parameters.jl 
using ClimaLSM.Soil.Biogeochemistry # how can I avoid this deps?
using SurfaceFluxes, Thermodynamics, CLIMAParameters # how can I avoid these deps?
include(joinpath(pkgdir(ClimaLSM), "parameters", "create_parameters.jl")) 
export create_lsm_parameters

include("fun_discretisation.jl")
export mat, d1_vec, d2_vec

include("generate_fig.jl")
export param_dashboard

end 
