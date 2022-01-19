using Test

using ClimaLSM
using ClimaLSM.Domains: Column, RootDomain
include(joinpath("../src", "Soil.jl"))
using .Soil
include(joinpath("../src", "Roots.jl"))
using .Roots
using UnPack
using NLsolve
using OrdinaryDiffEq: ODEProblem, solve, Euler


FT = Float64

include("./initial_structure_test.jl")
include("./root_test.jl")
