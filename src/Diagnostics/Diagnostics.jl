module Diagnostics

import ClimaDiagnostics:
    DiagnosticVariable,
    ScheduledDiagnostic,
    average_pre_output_hook!,
    DiagnosticsCallback

import ClimaDiagnostics.DiagnosticVariables: descriptive_short_name

import ClimaDiagnostics.Schedules: EveryStepSchedule, EveryDtSchedule

import ClimaDiagnostics.Writers: HDF5Writer, NetCDFWriter, write_field!

# include("diagnostic.jl")

end
