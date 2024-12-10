module Fluxnet

using ClimaLandSimulations
import SciMLBase
import ClimaTimeSteppers as CTS
using ClimaCore
import ClimaParams as CP
using Statistics
using Dates
using Insolation
using StatsBase
using Interpolations
using ClimaLand
using ClimaLand.Domains: Column
using ClimaLand.Soil
using ClimaLand.Soil.Biogeochemistry
using ClimaLand.Canopy
using ClimaLand.Canopy.PlantHydraulics
import ClimaLand.Parameters as LP
import ClimaComms
using Unitful: R, L, mol, K, kJ, °C, m, g, cm, hr, mg, s, μmol, Pa, W, mm, kPa
using UnitfulMoles: molC
using Unitful, UnitfulMoles
using ArtifactWrappers
using DelimitedFiles
using Thermodynamics
using Format
using DataFrames

# Directory paths
climalandsimulations_dir = pkgdir(ClimaLandSimulations)
climaland_dir = pkgdir(ClimaLand)

# ClimaParams
const FT = Float64
earth_param_set = LP.LandParameters(FT)

include("run_fluxnet.jl")

sites = ["US-MOz", "US-Ha1", "US-NR1", "US-Var"]
for site in sites
    include(joinpath("fluxnet_sites", site, "$(site)_simulation.jl"))
    include(joinpath("fluxnet_sites", site, "$(site)_config.jl"))
    include(joinpath("fluxnet_sites", site, "$(site)_parameters.jl"))
end

include(joinpath("fluxnet_utilities", "make_config.jl"))
include(joinpath("fluxnet_utilities", "make_setup.jl"))
include(joinpath("fluxnet_utilities", "make_parameters.jl"))
include(joinpath("fluxnet_utilities", "data_tools.jl"))
include(joinpath("fluxnet_utilities", "make_drivers.jl"))
include(joinpath("fluxnet_utilities", "inputs_dataframe.jl"))

# Register Fluxnet with Unitful.jl
function __init__()
    Unitful.register(Fluxnet)
end

end
