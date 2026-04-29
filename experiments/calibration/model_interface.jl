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
import TOML

include(joinpath(pkgdir(ClimaLand), "experiments", "calibration", "api.jl"))

"""
    LandModelInterface <: ClimaCalibrate.AbstractModelInterface

A model interface struct for running the land calibration pipeline.

See the ClimaCalibrate.jl documentation for the methods that
`LandModelInterface` should implement.
"""
struct LandModelInterface <: ClimaCalibrate.AbstractModelInterface
    config::CalibrateConfig
end

interface_path =
    joinpath(pkgdir(ClimaLand), "experiments", "calibration", "models")
include(joinpath(interface_path, "snowy_land.jl"))
include(joinpath(interface_path, "bucket.jl"))

"""
    ClimaCalibrate.forward_model(
        model_interface::LandModelInterface,
        iteration,
        member,
    )

Run a land simulation for calibration.

The snowy land model or bucket simulation is ran depending on the the
configuration in the `model_interface`.

This function may be called in parallel depending on the ClimaCalibrate backend
used.
"""
function ClimaCalibrate.forward_model(
    model_interface::LandModelInterface,
    iteration,
    member,
)
    (; config) = model_interface
    (; model_type) = config
    ClimaCalibrate.forward_model(model_interface, iteration, member, model_type)
    return nothing
end

"""
    ClimaCalibrate.model_interface_filepath(::LandModelInterface)

Return a filepath to the definition of the `CouplerModelInterface` struct and
all its associated methods.

This is required to use the `ClimaCalibrate.HPCBackend`s.
"""
function ClimaCalibrate.model_interface_filepath(::LandModelInterface)
    return @__FILE__
end

include(
    joinpath(
        pkgdir(ClimaLand),
        "experiments",
        "calibration",
        "observation_map.jl",
    ),
)
