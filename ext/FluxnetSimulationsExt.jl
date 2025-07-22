module FluxnetSimulationsExt
using ClimaLand
import DelimitedFiles

using ClimaTimeSteppers
using ClimaComms
import ClimaComms: context, device
using SciMLBase
using Dates
import ClimaUtilities.TimeManager: ITime, date
import ClimaDiagnostics

# Include site-specific configurations, as well as the default generic site.
include("fluxnet_simulations/generic_site.jl")
include("fluxnet_simulations/US-MOz.jl")
include("fluxnet_simulations/US-Ha1.jl")
include("fluxnet_simulations/US-NR1.jl")
include("fluxnet_simulations/US-Var.jl")
export get_parameters, get_domain_info

end # module FluxnetSimulationsExt
