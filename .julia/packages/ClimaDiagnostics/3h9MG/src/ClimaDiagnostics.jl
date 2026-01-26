module ClimaDiagnostics

include("AbstractTypes.jl")

include("utils.jl")
include("Schedules.jl")
include("DiagnosticVariables.jl")
import .DiagnosticVariables: DiagnosticVariable, average_pre_output_hook!
include("ScheduledDiagnostics.jl")
import .ScheduledDiagnostics: ScheduledDiagnostic
include("Writers.jl")

include("clima_diagnostics.jl")

end
