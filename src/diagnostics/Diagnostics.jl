module Diagnostics

import Dates: Month, Period, DateTime

import ClimaComms
using ClimaCore: Spaces, Fields
using LazyBroadcast: lazy
import ..Parameters as LP

import ..AbstractModel, ..AbstractLandModel
using ..Bucket: BucketModel
import ..SoilCanopyModel
import ..LandModel
import ..SoilSnowModel
import ..Soil: EnergyHydrology, Runoff
import ..Soil.Biogeochemistry: SoilCO2Model
import ..Snow: SnowModel
import ..Canopy:
    CanopyModel,
    medlyn_conductance,
    medlyn_term,
    gs_h2o_pmodel,
    get_Vcmax25_leaf,
    get_An_leaf,
    get_Rd_leaf,
    MedlynConductanceModel,
    PModelConductance,
    canopy_temperature,
    ZhouOptimalLAIModel
import ..Domains:
    top_center_to_surface,
    AbstractDomain,
    SphericalShell,
    HybridBox,
    Column,
    SphericalSurface,
    Plane,
    Point

import ClimaLand
import ClimaLand: get_domain
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
