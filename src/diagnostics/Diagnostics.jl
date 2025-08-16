module Diagnostics

import Dates: Month, Period, DateTime

import ClimaComms
using ClimaCore: Spaces, Fields

import ..Parameters as LP

using ..Bucket: BucketModel

import ..SoilCanopyModel

import ..LandModel

import ..Soil: EnergyHydrology

import ..AbstractModel

import ..Canopy:
    CanopyModel,
    medlyn_conductance,
    medlyn_term,
    gs_h2o_pmodel,
    get_Vcmax25,
    MedlynConductanceModel,
    PModelConductance,
    moisture_stress,
    canopy_temperature
import ..Domains:
    top_center_to_surface, AbstractDomain, SphericalShell, HybridBox
import ClimaLand
import ..heaviside

import ClimaDiagnostics:
    DiagnosticVariable, ScheduledDiagnostic, average_pre_output_hook!

import ClimaDiagnostics.Schedules:
    EveryStepSchedule, EveryDtSchedule, EveryCalendarDtSchedule

import ClimaDiagnostics
import ClimaDiagnostics.Writers: HDF5Writer, NetCDFWriter, DictWriter

import ClimaCore.Fields: zeros, field_values
import ClimaCore.Operators: column_integral_definite!
import ClimaUtilities.TimeManager: ITime, date

export close_output_writers

include("diagnostic.jl")

end
