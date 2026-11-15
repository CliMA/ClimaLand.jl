module Runoff
using DocStringExtensions
using LazyBroadcast: lazy
using ClimaCore
using ClimaCore.Operators: column_integral_definite!
using ClimaLand
import ClimaLand: source!
using ..ClimaLand.Soil:
    AbstractSoilSource,
    AbstractSoilModel,
    RichardsModel,
    EnergyHydrology,
    volumetric_internal_energy_liq
import ClimaLand.Parameters as LP
import ClimaParams as CP
export TOPMODELRunoff,
    NoRunoff,
    SurfaceRunoff,
    AbstractRunoffModel,
    TOPMODELSubsurfaceRunoff,
    subsurface_runoff_source,
    update_infiltration_water_flux!,
    AbstractFrozenInfiltration,
    IceImpedanceInfiltration,
    GradualFrozenInfiltration

"""
    AbstractRunoffModel

The soil runoff models are only to be used with
the following boundary condition types:
- `ClimaLand.Soil.AtmosDrivenFluxBC`
- `ClimaLand.Soil.RichardsAtmosDrivenFluxBC`.
It must have methods for
- `subsurface_runoff_source` (defined in this module)
- `update_infiltration_water_flux!` (defined in this module)
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
    AbstractFrozenInfiltration

Abstract type for the treatment of frozen soil in the infiltration capacity
at the soil surface. The concrete types only change the infiltration capacity
(the surface boundary condition); the conductivity of the soil interior is
not affected.
"""
abstract type AbstractFrozenInfiltration end

"""
    IceImpedanceInfiltration <: AbstractFrozenInfiltration

The default treatment of frozen soil in the infiltration capacity: the
saturated conductivity of the top soil layer is reduced by the ice impedance
factor 10^(-Ω f_i), with f_i the ice fraction of the water in the top layer,
as in the conductivity of the soil interior.
"""
struct IceImpedanceInfiltration <: AbstractFrozenInfiltration end

"""
    GradualFrozenInfiltration{FT} <: AbstractFrozenInfiltration

A treatment of frozen soil in the infiltration capacity following HTESSEL
(ECMWF IFS Documentation CY41R2, Part IV, Eqs. 8.70 and 8.88; Balsamo et al.
2009, doi:10.1175/2008JHM1068.1). A frozen fraction of the top soil layer is
diagnosed from its temperature,

f_fr(T) = 0.5 [1 - sin(π (T - (T_f1 + T_f2)/2) / (T_f1 - T_f2))]

for T_f2 ≤ T ≤ T_f1, with f_fr = 0 above T_f1 and f_fr = 1 below T_f2, where
T_f1 = T_freeze + ΔT_thaw and T_f2 = T_freeze - ΔT_frozen. The infiltration
capacity is the weighted mean of the saturated conductivity and of the
conductivity at the permanent wilting point (matric potential ψ_pwp) of the
top layer,

i_c = -K_sat [(1 - f_fr) + f_fr K_r(ψ_pwp)] G(T),

with G the viscosity factor, in place of the ice impedance factor of
[`IceImpedanceInfiltration`](@ref).

$(DocStringExtensions.FIELDS)
"""
struct GradualFrozenInfiltration{FT} <: AbstractFrozenInfiltration
    "Temperature above freezing at which the frozen fraction is zero, T_f1 - T_freeze (K)"
    ΔT_thaw::FT
    "Temperature below freezing at which the frozen fraction is one, T_freeze - T_f2 (K)"
    ΔT_frozen::FT
    "Matric potential at the permanent wilting point (m), negative"
    ψ_pwp::FT
end

"""
    GradualFrozenInfiltration(
        toml_dict::CP.ParamDict;
        ΔT_thaw = toml_dict["frozen_infiltration_thaw_offset"],
        ΔT_frozen = toml_dict["frozen_infiltration_frozen_offset"],
        ψ_pwp = toml_dict["frozen_infiltration_wilting_point_matric_potential"],
    )

Constructs a `GradualFrozenInfiltration` from the global parameters in
`toml_dict`; defaults are the HTESSEL values (T_f1 = T_freeze + 1 K,
T_f2 = T_freeze - 3 K, ψ_pwp = -153 m, i.e. -15 bar).
"""
function GradualFrozenInfiltration(
    toml_dict::CP.ParamDict;
    ΔT_thaw = toml_dict["frozen_infiltration_thaw_offset"],
    ΔT_frozen = toml_dict["frozen_infiltration_frozen_offset"],
    ψ_pwp = toml_dict["frozen_infiltration_wilting_point_matric_potential"],
)
    FT = CP.float_type(toml_dict)
    return GradualFrozenInfiltration{FT}(ΔT_thaw, ΔT_frozen, ψ_pwp)
end

"""
    frozen_fraction(T::FT, T_freeze::FT, ΔT_thaw::FT, ΔT_frozen::FT) where {FT}

The HTESSEL frozen fraction of soil water diagnosed from temperature (IFS
Documentation CY41R2, Part IV, Eq. 8.70): 0 above T_freeze + ΔT_thaw,
1 below T_freeze - ΔT_frozen, and a sine transition in between.
"""
function frozen_fraction(
    T::FT,
    T_freeze::FT,
    ΔT_thaw::FT,
    ΔT_frozen::FT,
) where {FT}
    T_f1 = T_freeze + ΔT_thaw
    T_f2 = T_freeze - ΔT_frozen
    x = clamp((T - (T_f1 + T_f2) / 2) / (T_f1 - T_f2), -FT(0.5), FT(0.5))
    return (1 - sin(FT(π) * x)) / 2
end

"""
    frozen_infiltration_model(runoff::AbstractRunoffModel)

Returns the treatment of frozen soil in the infiltration capacity used by
the runoff model `runoff`; [`IceImpedanceInfiltration`](@ref) unless the
runoff model carries another one.
"""
frozen_infiltration_model(runoff::AbstractRunoffModel) =
    IceImpedanceInfiltration()

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
    update_infiltration_water_flux!(p, runoff::NoRunoff, input, _...)

Updates the runoff variables in the cache `p.soil` in place
in the case of NoRunoff: sets infiltration = precipitation.
"""
function update_infiltration_water_flux!(p, runoff::NoRunoff, input, _...)
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
    update_infiltration_water_flux!(
        p,
        runoff::SurfaceRunoff,
        input,
        Y,
        t,
        model::AbstractSoilModel,
)

The update_infiltration_water_flux! function for the SurfaceRunoff model.

Updates the runoff model variables in place in `p.soil` for the SurfaceRunoff
parameterization:
p.soil.R_s
p.soil.is_saturated
p.soil.infiltration
"""
function update_infiltration_water_flux!(
    p,
    runoff::SurfaceRunoff,
    input,
    Y,
    t,
    model::AbstractSoilModel,
)

    ic = soil_infiltration_capacity(
        frozen_infiltration_model(runoff),
        model,
        Y,
        p,
    ) # should be non-allocating
    ϑ_l = Y.soil.ϑ_l
    FT = eltype(ϑ_l)
    θ_i = model_agnostic_volumetric_ice_content(Y, FT)
    @. p.soil.is_saturated = is_saturated(ϑ_l + θ_i, model.parameters.ν)
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

This is currently treated implicitly (evaluated at the next value of
ϑ\\_l) but not included in the Jacobian approximation. This is because
`update_infiltration_water_flux` is called as part of the implicit tendency.

$(DocStringExtensions.FIELDS)
"""
struct TOPMODELSubsurfaceRunoff{FT} <: AbstractSoilSource{FT}
    "The subsurface runoff flux (m/s) when the depth to the water table = 1/f_over; calibrated"
    R_sb::FT
    "A calibrated parameter defining how subsurface runoff decays with depth to water table (1/m ; calibrated)"
    f_over::FT
    "Explicit vs implicit stepping"
    explicit::Bool
    TOPMODELSubsurfaceRunoff{FT}(R_sb, f_over) where {FT} =
        new{FT}(R_sb, f_over, false)
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
struct TOPMODELRunoff{
    FT <: AbstractFloat,
    F <: ClimaCore.Fields.Field,
    FI <: AbstractFrozenInfiltration,
} <: AbstractRunoffModel
    "A calibrated parameter defining how subsurface runoff decays with depth to water table (1/m ; calibrated)"
    f_over::FT
    "The maximum saturated fraction of a grid cell, computed from the topographic index CDF per grid cell."
    f_max::F
    "The subsurface source term corresponding to this implementation of TOPMODEL."
    subsurface_source::TOPMODELSubsurfaceRunoff{FT}
    "The treatment of frozen soil in the infiltration capacity"
    frozen_infiltration::FI
end

TOPMODELRunoff{FT, F}(f_over, f_max, subsurface_source) where {FT, F} =
    TOPMODELRunoff{FT, F, IceImpedanceInfiltration}(
        f_over,
        f_max,
        subsurface_source,
        IceImpedanceInfiltration(),
    )

function TOPMODELRunoff{FT}(;
    f_over::FT,
    f_max::F,
    R_sb::FT,
    frozen_infiltration::FI = IceImpedanceInfiltration(),
) where {FT, F, FI}
    subsurface_source = TOPMODELSubsurfaceRunoff{FT}(R_sb, f_over)
    return TOPMODELRunoff{FT, F, FI}(
        f_over,
        f_max,
        subsurface_source,
        frozen_infiltration,
    )
end

function TOPMODELRunoff(
    toml_dict::CP.ParamDict;
    f_max::F,
    f_over = toml_dict["f_over"],
    R_sb = toml_dict["R_sb"],
    frozen_infiltration::FI = IceImpedanceInfiltration(),
) where {F, FI}
    FT = CP.float_type(toml_dict)
    subsurface_source = TOPMODELSubsurfaceRunoff{FT}(R_sb, f_over)
    return TOPMODELRunoff{FT, F, FI}(
        f_over,
        f_max,
        subsurface_source,
        frozen_infiltration,
    )
end

frozen_infiltration_model(runoff::TOPMODELRunoff) = runoff.frozen_infiltration

"""
    update_infiltration_water_flux!(p, runoff::TOPMODELRunoff, input, Y,t, model::AbstractSoilModel)

Updates the runoff model variables in place in `p.soil` for the TOPMODELRunoff
parameterization:
p.soil.R_s
p.soil.R_ss
p.soil.h∇
p.soil.infiltration
"""
function update_infiltration_water_flux!(
    p,
    runoff::TOPMODELRunoff,
    input,
    Y,
    t,
    model::AbstractSoilModel,
)
    ϑ_l = Y.soil.ϑ_l
    ν = model.parameters.ν
    θ_r = model.parameters.θ_r
    FT = eltype(ϑ_l)
    θ_i = model_agnostic_volumetric_ice_content(Y, FT)
    # Surface runoff first
    ic = soil_infiltration_capacity(
        frozen_infiltration_model(runoff),
        model,
        Y,
        p,
    ) # should be non-allocating
    @. p.soil.is_saturated =
        is_saturated(ϑ_l + θ_i - θ_r, ν - θ_r) * (ϑ_l + θ_i - θ_r) / (ν - θ_r) # weighted by how much above saturation, include ice!
    column_integral_definite!(p.soil.h∇, p.soil.is_saturated) # can be bigger than the depth
    @. p.soil.infiltration = topmodel_surface_infiltration(
        runoff.f_max,
        runoff.f_over,
        model.domain.fields.depth - p.soil.h∇,
        ic,
        input,
    )
    @. p.soil.R_s = abs(input - p.soil.infiltration)
    # eek! now do it again but use only liquid water for subsurface runoff
    @. p.soil.is_saturated =
        is_saturated(ϑ_l - θ_r, ν - θ_r) * (ϑ_l - θ_r) / (ν - θ_r) # weighted by how much above saturation
    column_integral_definite!(p.soil.h∇, p.soil.is_saturated) # can be bigger than the depth, only includes liquid water
    @. p.soil.R_ss = topmodel_ss_flux(
        runoff.subsurface_source.R_sb,
        runoff.f_over,
        model.domain.fields.depth - p.soil.h∇,
    )
    update_subsurface_energy_runoff!(p, model)
end

function update_subsurface_energy_runoff!(
    p,
    model::EnergyHydrology{FT},
) where {FT}
    @. p.soil.subsfc_scratch =
        p.soil.is_saturated * volumetric_internal_energy_liq(
            p.soil.T,
            model.parameters.earth_param_set,
        )
    column_integral_definite!(p.soil.R_ess, p.soil.subsfc_scratch) # this actual just computes the average volumetric energy of the liquid in the saturated layers multiplied by the water table height
    @. p.soil.R_ess *= p.soil.R_ss / max(p.soil.h∇, eps(FT)) #this divides by the water table height (to get average energy) and multiplies by volumetric liquid water runoff.
end
update_subsurface_energy_runoff!(p, model::RichardsModel) = nothing

runoff_vars(::TOPMODELRunoff) =
    (:infiltration, :is_saturated, :R_s, :R_ss, :R_ess, :h∇, :subsfc_scratch)
runoff_var_domain_names(::TOPMODELRunoff) =
    (:surface, :subsurface, :surface, :surface, :surface, :surface, :subsurface)
runoff_var_types(::TOPMODELRunoff, FT) = (FT, FT, FT, FT, FT, FT, FT)

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
        model::RichardsModel{FT},
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
    model::RichardsModel{FT},
) where {FT}
    h∇ = p.soil.h∇
    ϵ = eps(FT)
    @. dY.soil.ϑ_l -= (p.soil.R_ss / max(h∇, ϵ)) * p.soil.is_saturated # apply only to saturated layers
    @. dY.soil.∫F_vol_liq_water_dt -= p.soil.R_ss # the integral is designed to be this flux
end

"""
    ClimaLand.source!(
        dY::ClimaCore.Fields.FieldVector,
        src::TOPMODELSubsurfaceRunoff,
        Y::ClimaCore.Fields.FieldVector,
        p::NamedTuple,
        model::EnergyHydrology{FT},
    ) where {FT}

Adjusts dY.soil.ϑ_l and dY.soil.ρe_int in place
 to account for the loss of water and energy due to subsurface runoff.
"""
function ClimaLand.source!(
    dY::ClimaCore.Fields.FieldVector,
    src::TOPMODELSubsurfaceRunoff{FT},
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    model::EnergyHydrology{FT},
) where {FT}
    h∇ = p.soil.h∇
    ϵ = eps(FT)
    @. dY.soil.ϑ_l -= (p.soil.R_ss / max(h∇, ϵ)) * p.soil.is_saturated # apply only to saturated layers
    @. dY.soil.ρe_int -= (p.soil.R_ess / max(h∇, ϵ)) * p.soil.is_saturated # apply only to saturated layers
    @. dY.soil.∫F_vol_liq_water_dt -= p.soil.R_ss # the integral is designed to be this flux
    @. dY.soil.∫F_e_dt -= p.soil.R_ess # the integral is designed to be this flux
end

"""
    topmodel_surface_infiltration(f_max, f_over, z∇, f_ic, precip)

A pointwise function which returns the infiltration into the soil,
given the precipitation flux (m/s),
the depth to the water table z∇>0, the infiltration capacity flux f_ic,
and the TOPMODEL parameters f_max and f_over.

see: Niu et al. (2005),
"A simple TOPMODEL-based runoff parameterization (SIMTOP) for
use in global climate models", Equations (8) and (11).
"""
function topmodel_surface_infiltration(f_max, f_over, z∇, f_ic, precip)
    f_sat = min(f_max * exp(-f_over / 2 * z∇), 1)
    return (1 - f_sat) * max(f_ic, precip)
end

"""
    soil_infiltration_capacity(model::AbstractSoilModel, Y, p)

Computes the soil infiltration capacity on the surface space, using the
treatment of frozen soil of the runoff model in the top boundary condition
of `model` (see [`frozen_infiltration_model`](@ref)).
"""
soil_infiltration_capacity(model::AbstractSoilModel, Y, p) =
    soil_infiltration_capacity(
        frozen_infiltration_model(model.boundary_conditions.top),
        model,
        Y,
        p,
    )

frozen_infiltration_model(bc) =
    hasproperty(bc, :runoff) ? frozen_infiltration_model(bc.runoff) :
    IceImpedanceInfiltration()

"""
    soil_infiltration_capacity(::AbstractFrozenInfiltration, model::RichardsModel, Y, p)

Computes the soil infiltration capacity on the surface space
 for Richards model.

Currently approximates i_c = -K_sat at the surface.
"""
function soil_infiltration_capacity(
    ::AbstractFrozenInfiltration,
    model::RichardsModel,
    Y,
    p,
)
    @. p.soil.subsfc_scratch = -1 * model.parameters.K_sat
    return ClimaLand.Domains.top_center_to_surface(p.soil.subsfc_scratch)
end

"""
    soil_infiltration_capacity(::IceImpedanceInfiltration, model::EnergyHydrology, Y, p)

Computes the soil infiltration capacity on the surface space
 for the full soil model.

Currently approximates i_c = -K_sat*F(θ_i)*G(T) at the surface, where F and G
are the functions which adjust the conductivity for the presence ice and taking into
account the temperature dependence of the viscosity of water.
"""
function soil_infiltration_capacity(
    ::IceImpedanceInfiltration,
    model::EnergyHydrology,
    Y,
    p,
)
    (; K_sat, θ_r, Ω, γ, γT_ref) = model.parameters

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
    soil_infiltration_capacity(fi::GradualFrozenInfiltration, model::EnergyHydrology, Y, p)

Computes the soil infiltration capacity on the surface space for the full
soil model with the HTESSEL treatment of frozen soil (see
[`GradualFrozenInfiltration`](@ref)):
i_c = -K_sat [(1 - f_fr(T)) + f_fr(T) K_r(ψ_pwp)] G(T)
at the surface, with f_fr the frozen fraction diagnosed from the temperature
of the top layer, K_r(ψ_pwp) the relative conductivity at the permanent
wilting point and G the viscosity factor.
"""
function soil_infiltration_capacity(
    fi::GradualFrozenInfiltration,
    model::EnergyHydrology,
    Y,
    p,
)
    (; K_sat, hydrology_cm, γ, γT_ref, earth_param_set) = model.parameters
    (; ΔT_thaw, ΔT_frozen, ψ_pwp) = fi
    FT = eltype(fi.ψ_pwp)
    T_freeze = LP.T_freeze(earth_param_set)
    # Evaluated on the top layer only, into the top level of the scratch field.
    # The top-level views of the other fields are put on the same surface space
    # (scalars pass through).
    ic = ClimaLand.Domains.top_center_to_surface(p.soil.subsfc_scratch)
    sfc_space = axes(ic)
    top(x::ClimaCore.Fields.Field) = ClimaCore.Fields.Field(
        ClimaCore.Fields.field_values(
            ClimaLand.Domains.top_center_to_surface(x),
        ),
        sfc_space,
    )
    top(x) = x
    T_sfc = top(p.soil.T)
    cm_sfc = top(hydrology_cm)
    K_sat_sfc = top(K_sat)
    @. ic =
        -K_sat_sfc *
        gradual_frozen_factor(
            frozen_fraction(T_sfc, T_freeze, ΔT_thaw, ΔT_frozen),
            ClimaLand.Soil.hydraulic_conductivity(
                cm_sfc,
                FT(1),
                ClimaLand.Soil.inverse_matric_potential(cm_sfc, ψ_pwp),
            ),
        ) *
        ClimaLand.Soil.viscosity_factor(T_sfc, γ, γT_ref)
    return ic
end

"""
    gradual_frozen_factor(f_fr::FT, K_r_pwp::FT) where {FT}

The factor (1 - f_fr) + f_fr K_r_pwp multiplying the saturated conductivity in
the infiltration capacity of [`GradualFrozenInfiltration`](@ref), with f_fr
the frozen fraction and K_r_pwp the relative conductivity at the permanent
wilting point.
"""
gradual_frozen_factor(f_fr::FT, K_r_pwp::FT) where {FT} =
    (1 - f_fr) + f_fr * K_r_pwp

"""
    topmodel_ss_flux(R_sb::FT, f_over::FT, z∇::FT) where {FT}

A pointwise function which returns the subsurface runoff flux,
given the depth to the water table z∇>0,
and the TOPMODEL parameters R_sb and f_over.

A negative z∇, which supersaturated layers can give, is treated as a water
table at the surface, so the flux never exceeds R_sb.

see: Niu et al. (2005),
"A simple TOPMODEL-based runoff parameterization (SIMTOP) for
use in global climate models", Equations (12).
"""
function topmodel_ss_flux(R_sb::FT, f_over::FT, z∇::FT) where {FT}
    return R_sb * exp(-f_over * max(z∇, zero(FT)))
end

"""
    is_saturated(twc::FT, ν::FT) where {FT}

A helper function which can be used to indicate whether a layer of soil is
saturated based on if the total volumetric water content, `twc` is greater
than porosity `ν`.
"""
function is_saturated(twc::FT, ν::FT) where {FT}
    return ClimaLand.heaviside(twc - ν)
end

"""
    get_surface_runoff(runoff_model, Y, p)

A helper function which returns the surface runoff variable
based on the type of runoff model being used.
"""
get_surface_runoff(runoff_model::AbstractRunoffModel, Y, p) = nothing
get_surface_runoff(runoff_model::Union{SurfaceRunoff, TOPMODELRunoff}, Y, p) =
    p.soil.R_s

"""
    get_subsurface_runoff(runoff_model, Y, p)

A helper function which returns the subsurface runoff variable
based on the type of runoff model being used.
"""
get_subsurface_runoff(runoff_model::AbstractRunoffModel, Y, p) = nothing
get_subsurface_runoff(runoff_model::TOPMODELRunoff, Y, p) = p.soil.R_ss


"""
    get_saturated_height(runoff_model, Y, p)

A helper function which returns the integrated height of 
saturated soil based on the type of runoff model being used.
"""
get_saturated_height(runoff_model::AbstractRunoffModel, Y, p) = nothing
get_saturated_height(runoff_model::TOPMODELRunoff, Y, p) = p.soil.h∇


"""
    get_soil_fsat(runoff_model, Y, p, depth)

A helper function which returns the surface saturated fraction
 based on the type of runoff model being used.
"""
get_soil_fsat(runoff_model::AbstractRunoffModel, Y, p, depth) = nothing
function get_soil_fsat(runoff_model::TOPMODELRunoff, Y, p, depth)
    @. lazy(
        runoff_model.f_max *
        exp(-runoff_model.f_over / 2 * (depth - p.soil.h∇)),
    )
end

end
