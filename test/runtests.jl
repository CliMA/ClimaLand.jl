using Test
using Statistics
using DifferentialEquations
using UnPack
using NLsolve
using OrdinaryDiffEq: ODEProblem, solve, Euler
using ClimaCore

if !("." in LOAD_PATH) # for ease of include
    push!(LOAD_PATH, ".")
end
using ClimaLSM
using ClimaLSM.Domains: Column, RootDomain
using ClimaLSM.Soil
using ClimaLSM.Roots


FT = Float64

saved_values = SavedValues(FT, ClimaCore.Fields.FieldVector)

include("./initial_structure_test.jl")
include("./root_test.jl")
include("./soiltest.jl")
