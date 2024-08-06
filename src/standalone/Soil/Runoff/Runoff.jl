module Runoff
using DocStringExtensions
using ClimaCore
using ClimaCore.Operators: column_integral_definite!
using ClimaLand
import ClimaLand: source!
using ..ClimaLand.Soil:
    AbstractSoilSource,
    AbstractSoilModel,
    RichardsModel,
    EnergyHydrology,
    is_saturated
export TOPMODELRunoff,
    NoRunoff,
    SurfaceRunoff,
    AbstractRunoffModel,
    TOPMODELSubsurfaceRunoff,
    subsurface_runoff_source,
    topmodel_ss_flux,
    update_runoff!


"""
    AbstractRunoffModel

The soil runoff models are only to be used with
the following boundary condition types:
- `ClimaLand.Soil.AtmosDrivenFluxBC`
- `ClimaLand.Soil.RichardsAtmosDrivenFluxBC`.
It must have methods for
- `subsurface_runoff_source` (defined in this module)
- `update_runoff!` (defined in this module)
- `ClimaLand.source!`.
Please see the documentation for these for more details.

Your runoff model must specify the subsurface runoff sink term as well
as the surface runoff implementation.
"""
abstract type AbstractRunoffModel end

"""
    subsurface_runoff_source(runoff::AbstractRunoffModel)

A helper function which returns the subsurface source of the runoff
model `runoff`.
"""
subsurface_runoff_source(runoff::AbstractRunoffModel) = runoff.subsurface_source


"""
    NoRunoff <: AbstractRunoffModel
A concrete type of soil runoff model; the 
default choice, which does not include any
runoff.
"""
struct NoRunoff <: AbstractRunoffModel
    subsurface_source::Nothing
    function NoRunoff()
        return new(nothing)
    end
end

"""
    update_runoff!(p, runoff::NoRunoff, input, _...)

Updates the runoff variables in the cache `p.soil` in place
in the case of NoRunoff: sets infiltration = precipitation.
"""
function update_runoff!(p, runoff::NoRunoff, input, _...)
    p.soil.infiltration .= input
end

runoff_vars(::NoRunoff) = (:infiltration,)
runoff_var_domain_names(::NoRunoff) = (:surface,)
runoff_var_types(::NoRunoff, FT) = (FT,)

"""
    SurfaceRunoff <: AbstractRunoffModel

A simple model for runoff appropriate for single column runs.

Only surface runoff is computed, using a combination of Dunne 
and Hortonian runoff.
"""
struct SurfaceRunoff <: AbstractRunoffModel
    subsurface_source::Nothing
    function SurfaceRunoff()
        return new(nothing)
    end
end

"""
    surface_infiltration(
        f_ic::FT,
        input::FT,
        is_saturated::FT,
    ) where {FT}

Computes the surface infiltration for the simple surface
runoff model. If the soil is saturated at the surface,
all input is converted to runoff (infiltration is zero).

If the soil is not saturated, the maximum of the infiltration
capacity or the input is used as infiltration. Recall that
both are negative (towards the soil).
"""
function surface_infiltration(f_ic::FT, input::FT, is_saturated::FT) where {FT}
    return (1 - is_saturated) * max(f_ic, input)
end

"""
    update_runoff!(
        p,
        runoff::SurfaceRunoff,
        input,
        Y,
        t,
        model::AbstractSoilModel,
)

The update_runoff! function for the SurfaceRunoff model.

Updates the runoff model variables in place in `p.soil` for the SurfaceRunoff 
parameterization:
p.soil.R_s
p.soil.is_saturated
p.soil.infiltration
"""
function update_runoff!(
    p,
    runoff::SurfaceRunoff,
    input,
    Y,
    t,
    model::AbstractSoilModel,
)

    ic = soil_infiltration_capacity(model, Y, p) # should be non-allocating
    ϑ_l = Y.soil.ϑ_l
    FT = eltype(ϑ_l)
    θ_i = model_agnostic_volumetric_ice_content(Y, FT)
    @. p.soil.is_saturated = is_saturated(ϑ_l + θ_i, model.parameters.ν)
    surface_space = model.domain.space.surface
    is_saturated_sfc =
        ClimaLand.Domains.top_center_to_surface(p.soil.is_saturated) # a view

    @. p.soil.infiltration = surface_infiltration(ic, input, is_saturated_sfc)
    @. p.soil.R_s = abs(input .- p.soil.infiltration)
end

runoff_vars(::SurfaceRunoff) =
    (:is_saturated, :R_s, :infiltration, :subsfc_scratch)
runoff_var_domain_names(::SurfaceRunoff) =
    (:subsurface, :surface, :surface, :subsurface)
runoff_var_types(::SurfaceRunoff, FT) = (FT, FT, FT, FT)

# TOPMODEL

"""
    TOPMODELSubsurfaceRunoff{FT} <: AbstractSoilSource{FT}

The TOPMODEL subsurface runoff parameterization, which is implemented
as a sink term in the soil equations.

The runoff flux is given by Equation 12 of Niu et al. (2005),
"A simple TOPMODEL-based runoff parameterization (SIMTOP) for
use in global climate models".

$(DocStringExtensions.FIELDS)
"""
struct TOPMODELSubsurfaceRunoff{FT} <: AbstractSoilSource{FT}
    "The subsurface runoff flux (m/s) when the depth to the water table = 1/f_over; calibrated"
    R_sb::FT
    "A calibrated parameter defining how subsurface runoff decays with depth to water table (1/m ; calibrated)"
    f_over::FT
end

"""
    TOPMODELRunoff{FT <: AbstractFloat, F <: ClimaCore.Fields.Field} <: AbstractRunoffModel

The TOPMODEL surface runoff parameterization, which is affects the surface boundary
condition of the soil model.

The runoff flux is given by Equation 8 of with fsat given
by Equation (11), of Niu et al. (2005),
"A simple TOPMODEL-based runoff parameterization (SIMTOP) for
use in global climate models".

$(DocStringExtensions.FIELDS)
"""
struct TOPMODELRunoff{FT <: AbstractFloat, F <: ClimaCore.Fields.Field} <:
       AbstractRunoffModel
    "A calibrated parameter defining how subsurface runoff decays with depth to water table (1/m ; calibrated)"
    f_over::FT
    "The maximum saturated fraction of a grid cell, computed from the topographic index CDF per grid cell."
    f_max::F
    "The subsurface source term corresponding to this implementation of TOPMODEL."
    subsurface_source::TOPMODELSubsurfaceRunoff{FT}
end

function TOPMODELRunoff{FT}(; f_over::FT, f_max::F, R_sb::FT) where {FT, F}
    subsurface_source = TOPMODELSubsurfaceRunoff{FT}(R_sb, f_over)
    return TOPMODELRunoff{FT, F}(f_over, f_max, subsurface_source)
end


"""
    update_runoff!(p, runoff::TOPMODELRunoff, input, Y,t, model::AbstractSoilModel)

Updates the runoff model variables in place in `p.soil` for the TOPMODELRunoff
parameterization:
p.soil.R_s
p.soil.R_ss
p.soil.h∇
p.soil.infiltration
"""
function update_runoff!(
    p,
    runoff::TOPMODELRunoff,
    input,
    Y,
    t,
    model::AbstractSoilModel,
)
    ϑ_l = Y.soil.ϑ_l
    FT = eltype(ϑ_l)
    θ_i = model_agnostic_volumetric_ice_content(Y, FT)
    @. p.soil.is_saturated = is_saturated(ϑ_l + θ_i, model.parameters.ν)
    column_integral_definite!(p.soil.h∇, p.soil.is_saturated)
    @. p.soil.R_ss = topmodel_ss_flux(
        runoff.subsurface_source.R_sb,
        runoff.f_over,
        model.domain.depth - p.soil.h∇,
    )
    ic = soil_infiltration_capacity(model, Y, p) # should be non-allocating

    precip = p.drivers.P_liq
    @. p.soil.infiltration = topmodel_surface_infiltration(
        p.soil.h∇,
        runoff.f_max,
        runoff.f_over,
        model.domain.depth - p.soil.h∇,
        ic,
        input,
    )
    @. p.soil.R_s = abs(input - p.soil.infiltration)

end

runoff_vars(::TOPMODELRunoff) =
    (:infiltration, :is_saturated, :R_s, :R_ss, :h∇, :subsfc_scratch)
runoff_var_domain_names(::TOPMODELRunoff) =
    (:surface, :subsurface, :surface, :surface, :surface, :subsurface)
runoff_var_types(::TOPMODELRunoff, FT) = (FT, FT, FT, FT, FT, FT)

"""
    model_agnostic_volumetric_ice_content(Y, FT)

Helper function which returns the volumetric ice content stored
in Y.soil.θ_i, if present, as it is for EnergyHydrologyModels, or
else returns a scalar zero of type FT if it is not present,
as it is not for RichardsModel.
"""
model_agnostic_volumetric_ice_content(Y, FT) =
    :θ_i ∈ propertynames(Y.soil) ? Y.soil.θ_i : zero(FT)

"""
    ClimaLand.source!(
        dY::ClimaCore.Fields.FieldVector,
        src::TOPMODELSubsurfaceRunoff,
        Y::ClimaCore.Fields.FieldVector,
        p::NamedTuple,
        model::AbstractSoilModel{FT},
    ) where {FT}

Adjusts dY.soil.ϑ_l in place to account for the loss of
water due to subsurface runoff.

The sink term is given by - R_ss/h∇ H(twc - ν),
where H is the Heaviside function, h∇ is the water table
thickness (defined to be where twc>ν), where twc is the total water content,
 and R_ss is the runoff as a
flux(m/s).
"""
function ClimaLand.source!(
    dY::ClimaCore.Fields.FieldVector,
    src::TOPMODELSubsurfaceRunoff{FT},
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    model::AbstractSoilModel{FT},
) where {FT}
    ϑ_l = Y.soil.ϑ_l
    h∇ = p.soil.h∇
    ϵ = eps(FT)
    @. dY.soil.ϑ_l -= (p.soil.R_ss / max(h∇, ϵ)) * p.soil.is_saturated  # apply only to saturated layers
end

"""
    topmodel_surface_infiltration(h∇, f_max, f_over, z∇, f_ic, precip)

A pointwise function which returns the infiltration into the soil,
given the precipitation flux (m/s), the water table thickness h∇>=0,
the depth to the water table z∇>0, the infiltration capacity flux f_ic,
and the TOPMODEL parameters f_max and f_over.

see: Niu et al. (2005),
"A simple TOPMODEL-based runoff parameterization (SIMTOP) for
use in global climate models", Equations (8) and (11).
"""
function topmodel_surface_infiltration(h∇, f_max, f_over, z∇, f_ic, precip)
    f_sat = f_max * exp(-f_over / 2 * z∇)
    return (1 - f_sat) * max(f_ic, precip)
end

"""
    soil_infiltration_capacity(model::RichardsModel, Y, p)

Computes the soil infiltration capacity on the surface space
 for Richards model.

Currently approximates i_c = -K_sat at the surface.
"""
function soil_infiltration_capacity(model::RichardsModel, Y, p)
    @. p.soil.subsfc_scratch = -1 * model.parameters.K_sat
    return ClimaLand.Domains.top_center_to_surface(p.soil.subsfc_scratch)
end

"""
    soil_infiltration_capacity(model::EnergyHydrology, Y, p)

Computes the soil infiltration capacity on the surface space
 for the full soil model.

Currently approximates i_c = -K_sat*F(θ_i)*G(T) at the surface, where F and G
are the functions which adjust the conductivity for the presence ice and taking into
account the temperature dependence of the viscosity of water.
"""
function soil_infiltration_capacity(model::EnergyHydrology, Y, p)
    (; K_sat, θ_r, Ω, γ, γT_ref) = model.parameters
    surface_space = model.domain.space.surface
    @. p.soil.subsfc_scratch =
        -K_sat *
        ClimaLand.Soil.impedance_factor(
            Y.soil.θ_i / (p.soil.θ_l + Y.soil.θ_i - θ_r),
            Ω,
        ) *
        ClimaLand.Soil.viscosity_factor(p.soil.T, γ, γT_ref)
    return ClimaLand.Domains.top_center_to_surface(p.soil.subsfc_scratch)
end


"""
    topmodel_ss_flux(R_sb::FT, f_over::FT, z∇::FT) where {FT}

A pointwise function which returns the subsurface runoff flux,
given the depth to the water table z∇>0,
and the TOPMODEL parameters R_sb and f_over.

see: Niu et al. (2005),
"A simple TOPMODEL-based runoff parameterization (SIMTOP) for
use in global climate models", Equations (12).
"""
function topmodel_ss_flux(R_sb::FT, f_over::FT, z∇::FT) where {FT}
    return R_sb * exp(-f_over * z∇)
end

end
