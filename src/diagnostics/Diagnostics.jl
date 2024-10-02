module Diagnostics

import Dates: Month, Period

import ClimaComms

using ..Bucket: BucketModel

import ..SoilCanopyModel

import ..LandModel

import ..Soil: EnergyHydrology

import ..Domains: top_center_to_surface

import ClimaDiagnostics:
    DiagnosticVariable, ScheduledDiagnostic, average_pre_output_hook!

import ClimaDiagnostics.Schedules:
    EveryStepSchedule, EveryDtSchedule, EveryCalendarDtSchedule

import ClimaDiagnostics.Writers: HDF5Writer, NetCDFWriter, DictWriter

include("diagnostic.jl")

end
