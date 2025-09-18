import SciMLBase
import ClimaComms
ClimaComms.@import_required_backends
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
using CairoMakie
import GeoMakie
using Dates
import NCDatasets

import SciMLBase: solve!

function timesolve(integrator)
    device = ClimaComms.device()
    return ClimaComms.@elapsed device solve!(integrator)
end

function printstats(
    time,
    num_steps,
    dt;
    label = "simulation",
    num_columns,
    human_readable = false,
)
    time_per_step = time / num_steps
    steps_in_one_year = 365.25 * 86400 / dt
    steps_completed_in_one_day = 86400 / time_per_step
    sypd = steps_completed_in_one_day / steps_in_one_year
    if human_readable
        @info label num_columns time time_per_step sypd
    else
        println("(", num_columns, ",", time_per_step, ")")
        # println("'$label' ", num_columns, " ", time, " ", time_per_step, " ", sypd)
    end
end

function resolution(; dlat_degrees = 1.0)
    # We estimate 4 * n_horizontal_elements points on the equator
    n_horizontal_elements = convert(Int, ceil(360.0 / (dlat_degrees * 4)))
    effective_resolution = 360.0 / (4 * n_horizontal_elements)
    num_columns = 6 * n_horizontal_elements * n_horizontal_elements
    return n_horizontal_elements, effective_resolution, num_columns
end

function run(
    func;
    label = "",
    plot_attrs = "",
    tf = 86400.0,
    dt = 450.0,
    kwargs...,
)
    num_steps = convert(Int, tf / dt)

    # Warmup
    @time timesolve(func(; tf, dt, kwargs...))
    @time timesolve(func(; tf, dt, kwargs...))

    # Run
    # println("\\addplot[$(plot_attrs)]")
    # println("coordinates {")
    # for dlat_degrees in [8, 4, 2, 1, 0.5, 0.25, 0.125]
    #     time = timesolve(func(; dlat_degrees, tf, dt, kwargs...));
    #     _, _, num_columns = resolution(; dlat_degrees)
    #     printstats(time, num_steps, dt; label = "Full model", num_columns)
    # end
    # println("};")
    # println("\\addlegendentry{$label}")
end

include("snowy_land.jl")
include("bucket.jl")
include("richards.jl")
include("soil_canopy.jl")
include("canopy.jl")
include("soil.jl")
include("snow.jl")
include("lighter_canopy.jl")
