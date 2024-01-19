module ClimaLandSimulations

import SciMLBase
import ClimaTimeSteppers as CTS
using ClimaCore
import CLIMAParameters as CP
using Statistics
using Dates
using Insolation
using StatsBase
using Interpolations
using ClimaLSM
using ClimaLSM.Domains: Column
using ClimaLSM.Soil
using ClimaLSM.Soil.Biogeochemistry
using ClimaLSM.Canopy
using ClimaLSM.Canopy.PlantHydraulics
import ClimaLSM.Parameters as LSMP
using Unitful: R, L, mol, K, kJ, °C, m, g, cm, hr, mg, s, μmol, Pa, W, mm, kPa
using UnitfulMoles: molC
using Unitful, UnitfulMoles
using HTTP
using JSON
using ArtifactWrappers
using DelimitedFiles
using Dierckx
using Thermodynamics
using Formatting
using CairoMakie
using LaTeXStrings
using PlotUtils: optimize_ticks
using DataFrames

@compound CO₂
const FT = Float64

climalsm_dir = pkgdir(ClimaLSM)
savedir =
    joinpath(climalsm_dir, "experiments", "integrated", "fluxnet/figures/")
include(joinpath(climalsm_dir, "parameters", "create_parameters.jl"))
earth_param_set = create_lsm_parameters(FT)

ClimaLandSimulations_dir = pkgdir(ClimaLandSimulations)

include(joinpath("utilities", "makie_plots.jl"))
include(joinpath("utilities", "data_tools.jl"))
include(joinpath("utilities", "pull_MODIS.jl"))
include(joinpath("utilities", "met_drivers_FLUXNET.jl"))
include(joinpath("utilities", "climalsm_output_dataframe.jl"))
include(joinpath("utilities", "inputs_dataframe.jl"))
include("fluxnet_simulation.jl")

function __init__()
    Unitful.register(ClimaLandSimulations)
end

end
