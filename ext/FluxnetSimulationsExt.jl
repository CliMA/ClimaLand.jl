module FluxnetSimulationsExt
using ClimaLand
import DelimitedFiles

# Include site-specific configurations, as well as the default generic site.
include("fluxnet_simulations/generic_site.jl")
include("fluxnet_simulations/US-MOz.jl")
include("fluxnet_simulations/US-Ha1.jl")
include("fluxnet_simulations/US-NR1.jl")
include("fluxnet_simulations/US-Var.jl")
# export get_parameters, get_domain_info - TODO uncomment once we implement these

end # module FluxnetSimulationsExt
