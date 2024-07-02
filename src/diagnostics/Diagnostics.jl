module Diagnostics

import ClimaComms

using ..Bucket: BucketModel

import ..SoilCanopyModel

import ..Domains: top_center_to_surface

import ClimaDiagnostics:
    DiagnosticVariable, ScheduledDiagnostic, average_pre_output_hook!

import ClimaDiagnostics.Schedules: EveryStepSchedule, EveryDtSchedule

import ClimaDiagnostics.Writers: HDF5Writer, NetCDFWriter

include("diagnostic.jl")

end
