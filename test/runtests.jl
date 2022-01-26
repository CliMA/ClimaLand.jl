using Test

using ClimaLSM
using ClimaLSM.Domains: Column, RootDomain
using ClimaLSM.Soil
using ClimaLSM.Roots
using UnPack
using NLsolve
using OrdinaryDiffEq: ODEProblem, solve, Euler


FT = Float64

include("./initial_structure_test.jl")
include("./root_test.jl")
include("./soiltest.jl")
