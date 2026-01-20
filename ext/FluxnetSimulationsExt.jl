module FluxnetSimulationsExt
import ClimaLand.FluxnetSimulations as FluxnetSimulations
import ClimaUtilities.TimeVaryingInputs:
    TimeVaryingInput, LinearInterpolation, PeriodicCalendar
import ClimaUtilities.TimeManager: ITime, date
using Thermodynamics
using Dates
using DelimitedFiles
using DocStringExtensions
using Insolation
import ClimaLand.Parameters as LP
using ClimaLand
using ClimaLand.Canopy
using ClimaLand.PlantHydraulics
export prescribed_forcing_fluxnet,
    make_set_fluxnet_initial_conditions,
    get_comparison_data,
    get_data_dates,
    get_data_dt,
    replace_hyphen,
    get_parameters,
    get_domain_info

# Include site-specific configurations, as well as the default generic site.
include("fluxnet_simulations/generic_site.jl")
include("fluxnet_simulations/US-MOz.jl")
include("fluxnet_simulations/US-Ha1.jl")
include("fluxnet_simulations/US-NR1.jl")
include("fluxnet_simulations/US-Var.jl")
include("fluxnet_simulations/NEON-cper.jl")

# Include the data processing, forcing, and initial conditions utilities.
include("fluxnet_simulations/data_processing.jl")
include("fluxnet_simulations/forcing.jl")
include("fluxnet_simulations/initial_conditions.jl")
end # module FluxnetSimulationsExt
