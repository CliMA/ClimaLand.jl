const FT = Float64;
import ClimaComms
ClimaComms.@import_required_backends

context = ClimaComms.context()
ClimaComms.init(context)
device = ClimaComms.device()

using ClimaUtilities.ClimaArtifacts
import ClimaUtilities.TimeManager: ITime, date

import ClimaDiagnostics
import ClimaUtilities

import ClimaUtilities.TimeVaryingInputs:
    TimeVaryingInput, LinearInterpolation, PeriodicCalendar
import ClimaUtilities.ClimaArtifacts: @clima_artifact
import ClimaParams as CP
using ClimaCore
using ClimaLand
using ClimaLand.Snow
using ClimaLand.Soil
using ClimaLand.Canopy
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.Simulations: LandSimulation, solve!
import ClimaTimeSteppers as CTS

using Dates

import ClimaCalibrate
import JLD2
import EnsembleKalmanProcesses as EKP

interface_path =
    joinpath(pkgdir(ClimaLand), "experiments", "calibration", "models")
include(joinpath(interface_path, "snowy_land.jl"))
include(joinpath(interface_path, "bucket.jl"))

function ClimaCalibrate.forward_model(iteration, member)
    model_type = CALIBRATE_CONFIG.model_type
    ClimaCalibrate.forward_model(iteration, member, model_type)
    return nothing
end
