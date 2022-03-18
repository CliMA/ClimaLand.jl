using Test
using Statistics
using DifferentialEquations
using UnPack
using NLsolve
using OrdinaryDiffEq: ODEProblem, solve, Euler, Tsit5
using ClimaCore

if !("." in LOAD_PATH) # for ease of include
    push!(LOAD_PATH, ".")
end
using ClimaLSM
using ClimaLSM.Domains: Column, RootDomain, HybridBox, Plane, Point
using ClimaLSM.Soil
using ClimaLSM.Roots
using ClimaLSM.Pond

FT = Float64
include("./root_test.jl")
include("./soiltest.jl")
include("./lsm_test.jl")
include("./domains.jl")
include("./soil_test_3d.jl")
include("./pond_test.jl")
include("./pond_soil_lsm.jl")
