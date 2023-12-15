import SciMLBase
import ClimaTimeSteppers as CTS
using ClimaCore
import CLIMAParameters as CP
using Statistics
using Dates
using Insolation
using StatsBase
using Interpolations
using StatsBase

using ClimaLSM
using ClimaLSM.Domains: Column
using ClimaLSM.Soil
using ClimaLSM.Soil.Biogeochemistry
using ClimaLSM.Canopy
using ClimaLSM.Canopy.PlantHydraulics
import ClimaLSM
import ClimaLSM.Parameters as LSMP
climalsm_dir = pkgdir(ClimaLSM)
include(joinpath(climalsm_dir, "parameters", "create_parameters.jl"))
include(
    joinpath(
        climalsm_dir,
        "experiments",
        "integrated",
        "fluxnet",
        "data_tools.jl",
    ),
)
# include(joinpath(climalsm_dir, "experiments", "integrated", "fluxnet", "plot_utils.jl"))
const FT = Float64
earth_param_set = create_lsm_parameters(FT)

# Read in the site to be run from the command line
if length(ARGS) < 1
    error("Must provide site ID on command line")
end

site_ID = ARGS[1]

# Read all site-specific domain parameters from the simulation file for the site
include(
    joinpath(
        climalsm_dir,
        "experiments/integrated/fluxnet/$site_ID/$(site_ID)_simulation.jl",
    ),
)

include(
    joinpath(climalsm_dir, "experiments/integrated/fluxnet/fluxnet_domain.jl"),
)

# Read all site-specific parameters from the parameter file for the site
include(
    joinpath(
        climalsm_dir,
        "experiments/integrated/fluxnet/$site_ID/$(site_ID)_parameters.jl",
    ),
)

# This reads in the data from the flux tower site and creates
# the atmospheric and radiative driver structs for the model
include(
    joinpath(
        climalsm_dir,
        "experiments/integrated/fluxnet/met_drivers_FLUXNET.jl",
    ),
)
include(
    joinpath(
        climalsm_dir,
        "experiments/integrated/fluxnet/fluxnet_simulation.jl",
    ),
)
