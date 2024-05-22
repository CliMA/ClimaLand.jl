module Diagnostics

import ClimaComms

import LinearAlgebra: dot

import ClimaCore:
    Fields, Geometry, InputOutput, Meshes, Spaces, Operators, Domains, Grids

import ClimaCore.Utilities: half

import Thermodynamics as TD

import ClimaLand

import ClimaDiagnostics:
    DiagnosticVariable,
    ScheduledDiagnostic,
    average_pre_output_hook!,
    DiagnosticsCallback

import ClimaDiagnostics.DiagnosticVariables: descriptive_short_name

import ClimaDiagnostics.Schedules: EveryStepSchedule, EveryDtSchedule

import ClimaDiagnostics.Writers: HDF5Writer, NetCDFWriter, write_field!

include("diagnostic.jl")

end
