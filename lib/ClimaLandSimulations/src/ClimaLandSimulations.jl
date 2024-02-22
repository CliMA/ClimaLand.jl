module ClimaLandSimulations

using CairoMakie
using DataFrames
using LaTeXStrings
using PlotUtils: optimize_ticks
using HTTP
using Statistics
using JSON
using Dates
import ClimaTimeSteppers as CTS
using ClimaLand
using ClimaLand.Domains: Column
import ClimaLand.Parameters as LP

# Directory paths
climalandsimulations_dir = pkgdir(ClimaLandSimulations)
climaland_dir = pkgdir(ClimaLand)

include(joinpath("Fluxnet", "Fluxnet.jl"))
using .Fluxnet

include(joinpath("utilities", "make_domain.jl"))
include(joinpath("utilities", "make_timestepper.jl"))
include(joinpath("utilities", "makie_plots.jl"))
include(joinpath("utilities", "pull_modis.jl"))
include(joinpath("utilities", "climaland_output_dataframe.jl"))

end
