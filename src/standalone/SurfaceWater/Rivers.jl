module Rivers
using ClimaLand
using ClimaCore
using DocStringExtensions
import ClimaCore: Fields

import ClimaLand:
    AbstractExpModel,
    make_compute_exp_tendency,
    prognostic_vars,
    name,
    prognostic_types,
    prognostic_domain_names,
    auxiliary_vars,
    auxiliary_types,
    auxiliary_domain_names,
    make_update_boundary_fluxes,
    FTfromY

using ClimaLand.Domains
export RiverModel

abstract type AbstractSurfaceWaterModel{FT} <: AbstractExpModel{FT} end

"""
    RiverModel{FT, D, R} <: AbstractSurfaceWaterModel{FT}

Interfacer for integrating ClimaLand.jl and ClimaRivers.jl.

$(DocStringExtensions.FIELDS)
"""
struct RiverModel{FT, D,R} <: AbstractSurfaceWaterModel{FT}
    "The domain for the river model (global)"
    domain::D
    "The river model timestep"
    Δt_river::FT
    "time bracket"
    t_bracket::NTuple{2,FT}
    "History length"
    hist_length::Int
    "Runoff model"
    runoff::R
end
ClimaLand.name(model::RiverModel) = :river

ClimaLand.prognostic_vars(model::RiverModel) = (:R_accum,)
ClimaLand.prognostic_types(model::RiverModel{FT}) where {FT} = (FT,)
ClimaLand.prognostic_domain_names(model::RiverModel) = (:surface,)
ClimaLand.auxiliary_vars(model::RiverModel) = (:R_accum_history,:R_inst)
ClimaLand.auxiliary_types(model::RiverModel{FT}) where {FT} = (NTuple{model.hist_length, FT},FT)
ClimaLand.auxiliary_domain_names(model::RiverModel) = (:surface,:surface)

function ClimaLand.make_compute_exp_tendency(model::RiverModel)
    function compute_exp_tendency!(dY, Y, p, t)
        R_inst = p.river.R_inst
        compute_instantaneous_runoff!(R_inst, model.runoff,Y,p,t)
        @. dY.river.R_accum = instantaneous_runoff(model.runoff,Y,p,t)
    end
    return compute_exp_tendency!
end

function compute_instantaneous_runoff!(R_inst, runoff::Standalone,Y,p,t)
    R_inst .= runoff.data(t)
end

function compute_instantaneous_runoff!(R_inst, runoff::Integrated,Y,p,t)
    @. R_inst = p.soil.R_s + p.soil.R_ss # or -> grid_to_basin(R_s+R_ss)
end

function river_callback!(Y,p,t, t_callback)
    # if we have stepped a full Δt_river from the previous step
    if t ≈ t_callback # = t_prev + Δt_river
        p.river.R_accum_history[1:end-1] .= p.river.R_accum_history[2:end]
        p.river.R_accum_history[end] .= Y.river.R_accum
        Y.river.R_accum .= 0
        update_river_state!(p.river.state, p.river.R_accum_history)
    else
        nothing
    end
end
