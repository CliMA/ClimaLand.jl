import SciMLBase
using CUDA
CUDA.precompile_runtime()
import ClimaComms
ClimaComms.@import_required_backends
CUDA.precompile_runtime()
import ClimaTimeSteppers as CTS
using ClimaCore
using ClimaUtilities.ClimaArtifacts
import ClimaUtilities.OnlineLogging: WallTimeInfo, report_walltime

import ClimaDiagnostics
import ClimaAnalysis
import ClimaAnalysis.Visualize as viz
import ClimaUtilities

import ClimaUtilities.TimeVaryingInputs:
    TimeVaryingInput, LinearInterpolation, PeriodicCalendar, Throw
import ClimaUtilities.SpaceVaryingInputs: SpaceVaryingInput
import ClimaUtilities.ClimaArtifacts: @clima_artifact
import ClimaParams as CP

using ClimaLand
using ClimaLand.Snow
using ClimaLand.Soil
using ClimaLand.Canopy
import ClimaLand
import ClimaLand.Parameters as LP
using ClimaLand.Bucket:
    BucketModel, BucketModelParameters, PrescribedBaregroundAlbedo

using Statistics
using Dates
import NCDatasets

import ClimaLand.Simulations: solve!

function timesolve(integrator)
    device = ClimaComms.device()
    return ClimaComms.@elapsed device ClimaLand.Simulations.solve!(integrator)
end

function printstats(time, num_steps, dt; effective_resolution)
    time_per_step = time / num_steps
    np = ClimaComms.nprocs(ClimaComms.context())
    println("($effective_resolution , $np , $time_per_step)")

end

function resolution(; dlat_degrees = 1.0)
    # We estimate 4 * n_horizontal_elements points on the equator
    n_horizontal_elements = convert(Int, ceil(360.0 / (dlat_degrees * 4)))
    effective_resolution = 360.0 / (4 * n_horizontal_elements)
    num_columns = 6 * n_horizontal_elements * n_horizontal_elements
    return n_horizontal_elements, effective_resolution, num_columns
end

function run(func; dt = 450.0, tf = dt * 500, kwargs...)
    # subtract one from total steps because we step once before starting timing
    num_steps = convert(Int, tf / dt) - 1

    np = ClimaComms.nprocs(ClimaComms.context())
    res_array = [0.5, 0.25]
    np > 1 && push!(res_array, res_array[end] / 2)
    np >= 4 && push!(res_array, res_array[end] / 2)
    np >= 8 && push!(res_array, res_array[end] / 2)
    np >= 32 && push!(res_array, res_array[end] / 2)
    for dlat_degrees in res_array
        # solve once to ensure compilation
        timesolve(func(; dlat_degrees, tf = tf / 10, dt, kwargs...))
        @info "\n \n \n \n \n starting timing for $dlat_degrees resolution, $np gpus"
        integrator = func(; dlat_degrees, tf, dt, kwargs...)
        # step once to ensure compilation of anonymous functions
        ClimaLand.Simulations.step!(integrator)
        time = timesolve(integrator)
        _, effective_resolution, num_columns = resolution(; dlat_degrees)
        ClimaComms.iamroot(ClimaComms.context()) &&
            printstats(time, num_steps, dt; effective_resolution)
        @info "\n \n \n \n \n done timing for $dlat_degrees resolution, $np gpus"
    end
end

include("snowy_land.jl")
