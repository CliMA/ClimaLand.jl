module ParamViz

using WGLMakie
using JSServe
using SparseArrays
using Statistics
using Unitful: R, L, mol, K, kJ, °C, m, g, cm, hr, mg, s, μmol
using UnitfulMoles: molC
using Unitful, UnitfulMoles
@compound CO₂

include("struct_and_functions.jl")
export Drivers, Parameters, Constants, Inputs, Output # struct
export mat, d1_vec, d2_vec, parameterisation # functions

include("generate_fig.jl")
export param_dashboard
export webapp

function __init__()
    Unitful.register(ParamViz)
end

end 
