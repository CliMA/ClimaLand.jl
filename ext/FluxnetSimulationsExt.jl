module FluxnetSimulationsExt
import ClimaLand.FluxnetSimulations as FluxnetSimulations
import ClimaUtilities.TimeVaryingInputs:
    TimeVaryingInput, LinearInterpolation, PeriodicCalendar
import ClimaUtilities.TimeManager: ITime, date
using Thermodynamics
using Dates
using DelimitedFiles
using Insolation
import ClimaLand.Parameters as LP
using ClimaLand
using ClimaLand.Soil
using ClimaLand.Canopy
using ClimaLand.PlantHydraulics
using ClimaLand.Domains
using ClimaCore

using NCDatasets
export prescribed_forcing_fluxnet,
    make_set_fluxnet_initial_conditions,
    get_comparison_data,
    get_data_dates,
    get_data_dt,
    replace_hyphen,
    get_parameters,
    get_domain_info,
    get_location

# Get the metadata for the fluxnet sites.
include("fluxnet_simulations/get_fluxnet_metadata.jl")

# Include site-specific configurations, as well as the default generic site.
include("fluxnet_simulations/generic_site.jl")
include("fluxnet_simulations/US-MOz.jl")
include("fluxnet_simulations/US-Ha1.jl")
include("fluxnet_simulations/US-NR1.jl")
include("fluxnet_simulations/US-Var.jl")

# Include the data processing, forcing, and initial conditions utilities.
include("fluxnet_simulations/data_processing.jl")
include("fluxnet_simulations/forcing.jl")
include("fluxnet_simulations/initial_conditions.jl")
end # module FluxnetSimulationsExt
