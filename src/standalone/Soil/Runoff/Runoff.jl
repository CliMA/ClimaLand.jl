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
    InfiltrationExcess,
    NoRunoff,
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
    update_runoff!(p, runoff::NoRunoff, _...)

Updates the runoff variables in the cache `p.soil` in place
in the case of NoRunoff: sets infiltration = precipitation.
"""
function update_runoff!(p, runoff::NoRunoff, _...)
    p.soil.infiltration .= p.drivers.P_liq
end

"""
TO DO IN THIS PR, docstring
"""
struct InfiltrationExcess <: AbstractRunoffModel
    subsurface_source::Nothing
    function InfiltrationExcess()
        return new(nothing)
    end
end

"""
TO DO IN THIS PR, docstring
"""
function update_runoff!(p, runoff::InfiltrationExcess, Y, t, model)
    # compute infiltration capacity
    precip = p.drivers.P_liq # + melt - interception.... eventually
    flux_ic = soil_infiltration_capacity_flux(model, Y, p) # allocates
   # flux_ic = -i_c
    # If |P| > infiltration capacity i_c (i_c > 0)
    # then infiltration = -i_c, runoff = P-i_c
    # else infiltration = P
    @. p.soil.infiltration = max(flux_ic, precip)
    @. p.soil.R_s = abs(precip - p.soil.infiltration)
# Precip = -20 mm/d, flux_ic = -10 mm/d
# Precip = 0, flux_ic = -10mm/d
end
## Infiltration InfiltrationExcess

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
    update_runoff!(p, runoff::TOPMODELRunoff, Y,t, model::AbstractSoilModel)

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
    Y,
    t,
    model::AbstractSoilModel,
)
    column_integral_definite!(p.soil.h∇, is_saturated(Y, model))
    @. p.soil.R_ss = topmodel_ss_flux(
        runoff.subsurface_source.R_sb,
        runoff.f_over,
        model.domain.depth - p.soil.h∇,
    )
    precip = p.drivers.P_liq
    flux_ic = soil_infiltration_capacity_flux(model, Y, p) # allocates

    @. p.soil.infiltration = topmodel_surface_infiltration(
        p.soil.h∇,
        runoff.f_max,
        runoff.f_over,
        model.domain.depth - p.soil.h∇,
        flux_ic,
        precip,
    )
    @. p.soil.R_s = abs(precip - p.soil.infiltration)

end

"""
    ClimaLand.source!(
        dY::ClimaCore.Fields.FieldVector,
        src::TOPMODELSubsurfaceRunoff,
        Y::ClimaCore.Fields.FieldVector,
        p::NamedTuple,
        model::AbstractSoilModel,
    )

Adjusts dY.soil.ϑ_l in place to account for the loss of
water due to subsurface runoff.

The sink term is given by - R_ss/h∇ H(ϑ_l+θ_i - ν),
where H is the Heaviside function, h∇ is the water table
thickness (defined to be where ϑ_l+θ_i>ν), and R_ss is the runoff as a
flux(m/s).

The θ_i contribution is not include for RichardsModel.
"""
function ClimaLand.source!(
    dY::ClimaCore.Fields.FieldVector,
    src::TOPMODELSubsurfaceRunoff,
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    model::AbstractSoilModel,
)
    FT = eltype(Y.soil.ϑ_l)
    h∇ = p.soil.h∇
    dY.soil.ϑ_l .-= @.(p.soil.R_ss / max(h∇, eps(FT))) .* is_saturated(Y, model) # apply only to saturated layers
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
    soil_infiltration_capacity_flux(model::RichardsModel, Y, p)

Computes the soil infiltration capacity as as flux for Richards model.

Currently approximates i_c = -K_sat at the surface.
"""
function soil_infiltration_capacity_flux(model::RichardsModel, Y, p)
    return ClimaLand.Soil.get_top_surface_field(
        -1 .* model.parameters.K_sat,
        model.domain.space.surface,
    )
end

"""
    soil_infiltration_capacity_flux(model::EnergyHydrology, Y, p)

Computes the soil infiltration capacity as as flux for the full soil model.

Currently approximates i_c = -K_sat*F(θ_i)*G(T) at the surface, where F and G
are the functions which adjust the conductivity for the presence ice and taking into
account the temperature dependence of the viscosity of water.
"""
function soil_infiltration_capacity_flux(model::EnergyHydrology, Y, p)
    (; Ksat, θ_r, Ω, γ, γT_ref) = model.parameters
    K_eff = @. Ksat *
       impedance_factor(Y.soil.θ_i / (p.soil.θ_l + Y.soil.θ_i - θ_r), Ω) *
       viscosity_factor(p.soil.T, γ, γT_ref)
    return ClimaLand.Soil.get_top_surface_field(
        -1 .* Keff,
        model.domain.space.surface,
    )
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
