module Pond
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
    make_update_explicit_boundary_fluxes,
    FTfromY

using ClimaLand.Domains
export PondModel, PrescribedRunoff, surface_runoff

abstract type AbstractSurfaceWaterModel{FT} <: AbstractExpModel{FT} end
abstract type AbstractSurfaceRunoff end
"""
    PondModel{FT, D, R} <: AbstractSurfaceWaterModel{FT}

A stand-in model for models like the snow or river model. In
standalone mode, a prescribed soil infiltration rate
 and precipitation rate
control the rate of change of the pond height variable `η` via an ODE.
In integrated LSM mode, the infiltration into the soil will be computed
via a different method, and also be applied as a flux boundary condition
for the soil model.
$(DocStringExtensions.FIELDS)
"""
struct PondModel{FT, D, R} <: AbstractSurfaceWaterModel{FT}
    "The domain for the pond model"
    domain::D
    "The runoff model for the pond model"
    runoff::R
end
ClimaLand.name(model::AbstractSurfaceWaterModel) = :surface_water
function PondModel{FT}(;
    domain::ClimaLand.Domains.AbstractDomain{FT} = ClimaLand.Domains.Point(
        z_sfc = FT(0),
    ),
    runoff::AbstractSurfaceRunoff,
) where {FT}
    return PondModel{FT, typeof(domain), typeof(runoff)}(domain, runoff)
end


"""
    PrescribedRunoff{F1 <: Function, F2 <: Function} <:  AbstractSurfaceRunoff

The required input for driving the simple pond model: precipitation, as a
function of time, soil effective saturation at a depth `Δz` below the surface,
as a function of time, and soil parameters, which affect infiltration.
"""
struct PrescribedRunoff{F1 <: Function, F2 <: Function} <: AbstractSurfaceRunoff
    "Time dependent precipitation magnitude, given in m/s. Negative is into the soil"
    precip::F1
    "Time dependent infiltration magnitude, given in m/s. Negative is into the soil."
    infil::F2
end

ClimaLand.prognostic_vars(model::PondModel) = (:η,)
ClimaLand.prognostic_types(model::PondModel{FT}) where {FT} = (FT,)
ClimaLand.prognostic_domain_names(model::PondModel) = (:surface,)
ClimaLand.auxiliary_vars(model::PondModel) = (:runoff,)
ClimaLand.auxiliary_types(model::PondModel{FT}) where {FT} = (FT,)
ClimaLand.auxiliary_domain_names(model::PondModel) = (:surface,)

function ClimaLand.make_update_explicit_boundary_fluxes(model::PondModel)
    function update_explicit_boundary_fluxes!(p, Y, t)
        p.surface_water.runoff .= surface_runoff(model.runoff, Y, p, t)
    end
    return update_explicit_boundary_fluxes!
end


function ClimaLand.make_compute_exp_tendency(model::PondModel)
    function compute_exp_tendency!(dY, Y, p, t)
        @. dY.surface_water.η = p.surface_water.runoff
    end
    return compute_exp_tendency!
end

# Runoff > 0 -> into river system.
function surface_runoff(runoff::PrescribedRunoff, Y, p, t)
    FT = FTfromY(Y)
    return @. -FT((runoff.precip(t) - runoff.infil(t)))
end

end
