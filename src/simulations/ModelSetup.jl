module ModelSetup
using ClimaComms
using ClimaCore
using Interpolations
import ClimaUtilities.Regridders: InterpolationsRegridder
import ClimaUtilities.SpaceVaryingInputs: SpaceVaryingInput
import ClimaUtilities.ClimaArtifacts: @clima_artifact
using ClimaLand

regridder_type = :InterpolationsRegridder
extrapolation_bc =
    (Interpolations.Periodic(), Interpolations.Flat(), Interpolations.Flat())

include("domains.jl")
include("spatial_parameters.jl")
include("model_setup.jl")
end # module
