module PlantHydraulics
using ClimaLSM
using ..ClimaLSM.Canopy:
    AbstractCanopyComponent,
    update_canopy_prescribed_field!,
    set_canopy_prescribed_field!,
    AbstractSoilDriver,
    PrescribedSoil
using ClimaCore
using DocStringExtensions

import ClimaLSM:
    initialize_prognostic,
    make_update_aux,
    make_compute_exp_tendency,
    prognostic_vars,
    prognostic_types,
    auxiliary_vars,
    initialize,
    initialize_auxiliary,
    name
export PlantHydraulicsModel,
    AbstractPlantHydraulicsModel,
    flux,
    effective_saturation,
    augmented_liquid_fraction,
    water_retention_curve,
    inverse_water_retention_curve,
    root_flux_per_ground_area!,
    PlantHydraulicsParameters,
    PrescribedTranspiration,
    DiagnosticTranspiration,
    AbstractConductivityModel,
    AbstractRetentionModel,
    LinearRetentionCurve,
    Weibull,
    hydraulic_conductivity,
    PrescribedSiteAreaIndex

"""
    AbstractPlantHydraulicsModel{FT} <: AbstractCanopyComponent{FT}

An abstract type for plant hydraulics models.
"""
abstract type AbstractPlantHydraulicsModel{FT} <: AbstractCanopyComponent{FT} end

ClimaLSM.name(::AbstractPlantHydraulicsModel) = :hydraulics

"""
    AbstractTranspiration{FT <: AbstractFloat}

An abstract type for types representing different models of
transpiration (Prescribed or Diagnostic)
"""
abstract type AbstractTranspiration{FT <: AbstractFloat} end

"""
   PrescribedSiteAreaIndex{FT <:AbstractFloat}

A struct containing the area indices of the plants at a specific site;
LAI varies in time, while SAI and RAI are fixed. No spatial variation is
modeled.

$(DocStringExtensions.FIELDS)
"""
struct PrescribedSiteAreaIndex{FT <: AbstractFloat}
    "A function of simulation time `t` giving the leaf area index (LAI; m2/m2)"
    LAIfunction::Function
    "The constant stem area index (SAI; m2/m2)"
    SAI::FT
    "The constant root area index (RAI; m2/m2)"
    RAI::FT
end

"""
    lai_consistency_check(
        n_stem::Int64,
        n_leaf::Int64,
        area_index::NamedTuple{(:root, :stem, :leaf), Tuple{FT, FT, FT}},
    ) where {FT}

Carries out consistency checks using the area indices supplied and the number of
stem elements being modeled.

Note that it is possible to have a plant with no stem compartments
but with leaf compartments, and that a plant must have leaf compartments
(even if LAI = 0).

Specifically, this checks that:
1. n_leaf > 0
2. if LAI is nonzero or SAI is nonzero, RAI must be nonzero.
3. if SAI > 0, n_stem must be > 0 for consistency. If SAI == 0, n_stem must
be zero.
"""
function lai_consistency_check(
    n_stem::Int64,
    n_leaf::Int64,
    area_index::NamedTuple{(:root, :stem, :leaf), Tuple{FT, FT, FT}},
) where {FT}
    @assert n_leaf > 0
    if area_index.leaf > eps(FT) || area_index.stem > eps(FT)
        @assert area_index.root > eps(FT)
    end
    # If there SAI is zero, there should be no stem compartment
    if area_index.stem < eps(FT)
        @assert n_stem == FT(0)
    else
        # if SAI is > 0, n_stem should be > 0 for consistency
        @assert n_stem > 0
    end

end


"""
    PlantHydraulicsParameters{FT <: AbstractFloat}

A struct for holding parameters of the PlantHydraulics Model.
$(DocStringExtensions.FIELDS)
"""
struct PlantHydraulicsParameters{FT <: AbstractFloat, CP, RP}
    "The area index model for LAI, SAI, RAI"
    ai_parameterization::PrescribedSiteAreaIndex{FT}
    "porosity (m3/m3)"
    ν::FT
    "storativity (m3/m3)"
    S_s::FT
    "Root distribution function P(z)"
    root_distribution::Function
    "Conductivity model and parameters"
    conductivity_model::CP
    "Water retention model and parameters"
    retention_model::RP
end

function PlantHydraulicsParameters(;
    ai_parameterization::PrescribedSiteAreaIndex{FT},
    ν::FT,
    S_s::FT,
    root_distribution::Function,
    conductivity_model,
    retention_model,
) where {FT}
    return PlantHydraulicsParameters{
        FT,
        typeof(conductivity_model),
        typeof(retention_model),
    }(
        ai_parameterization,
        ν,
        S_s,
        root_distribution,
        conductivity_model,
        retention_model,
    )
end

"""
    PlantHydraulicsModel{FT, PS, T} <: AbstractPlantHydraulicsModel{FT}

Defines, and constructs instances of, the PlantHydraulicsModel type, which is used
for simulation flux of water to/from soil, along roots of different depths,
along a stem, to a leaf, and ultimately being lost from the system by
transpiration. Note that the canopy height is specified as part of the 
PlantHydraulicsModel, along with the area indices of the leaves, roots, and
stems.

This model can also be combined with the soil model using ClimaLSM, in which
case the prognostic soil water content is used to determine root extraction, and
the transpiration is also computed diagnostically. In  global run with patches
of bare soil, you can "turn off" the canopy model (to get zero root extraction, zero absorption and
emission, zero transpiration and sensible heat flux from the canopy), by setting:
- n_leaf = 1
- n_stem = 0
- LAI = SAI = RAI = 0.

A plant model can have leaves but no stem, but not vice versa. If n_stem = 0, SAI must be zero.

Finally, the model can be used in Canopy standalone mode by prescribing
the soil matric potential at the root tips or flux in the roots. There is also the
option (intendend only for debugging) to use a prescribed transpiration rate.

$(DocStringExtensions.FIELDS)
"""
struct PlantHydraulicsModel{FT, PS, T} <: AbstractPlantHydraulicsModel{FT}
    "The number of stem compartments for the plant; can be zero"
    n_stem::Int64
    "The number of leaf compartments for the plant; must be >=1"
    n_leaf::Int64
    "The height of the center of each leaf compartment/stem compartment, in meters"
    compartment_midpoints::Vector{FT}
    "The height of the compartments' top faces, in meters. The canopy height is the last element of the vector."
    compartment_surfaces::Vector{FT}
    "The label (:stem or :leaf) of each compartment"
    compartment_labels::Vector{Symbol}
    "Parameters required by the Plant Hydraulics model"
    parameters::PS
    "The transpiration model, of type `AbstractTranspiration`"
    transpiration::T
end

function PlantHydraulicsModel{FT}(;
    n_stem::Int64,
    n_leaf::Int64,
    compartment_midpoints::Vector{FT},
    compartment_surfaces::Vector{FT},
    parameters::PlantHydraulicsParameters{FT},
    transpiration::AbstractTranspiration{FT} = DiagnosticTranspiration{FT}(),
) where {FT}
    args = (parameters, transpiration)
    @assert (n_leaf + n_stem) == length(compartment_midpoints)
    @assert (n_leaf + n_stem) + 1 == length(compartment_surfaces)
    for i in 1:length(compartment_midpoints)
        @assert compartment_midpoints[i] ==
                ((compartment_surfaces[i + 1] - compartment_surfaces[i]) / 2) +
                compartment_surfaces[i]
    end
    compartment_labels = Vector{Symbol}(undef, n_stem + n_leaf)
    for i in 1:(n_stem + n_leaf)
        if i <= n_stem
            compartment_labels[i] = :stem
        else
            compartment_labels[i] = :leaf
        end
    end
    return PlantHydraulicsModel{FT, typeof.(args)...}(
        n_stem,
        n_leaf,
        compartment_midpoints,
        compartment_surfaces,
        compartment_labels,
        args...,
    )
end

"""
    prognostic_vars(model::PlantHydraulicsModel)

A function which returns the names of the prognostic
variables of the `PlantHydraulicsModel`.
"""
prognostic_vars(model::PlantHydraulicsModel) = (:ϑ_l,)

"""
    auxiliary_vars(model::PlantHydraulicsModel)

A function which returns the names of the auxiliary
variables of the `PlantHydraulicsModel`,
the transpiration stress factor `β` (unitless),
the water potential `ψ` (m), the volume flux*cross section `fa` (1/s),
and the volume flux*root cross section in the roots `fa_roots` (1/s),
where the cross section can be represented by an area index.
"""
auxiliary_vars(model::PlantHydraulicsModel) =
    (:β, :ψ, :fa, :fa_roots, :area_index)

"""
    ClimaLSM.prognostic_types(model::PlantHydraulicsModel{FT}) where {FT}

Defines the prognostic types for the PlantHydraulicsModel.
"""
ClimaLSM.prognostic_types(model::PlantHydraulicsModel{FT}) where {FT} =
    (NTuple{model.n_stem + model.n_leaf, FT},)

"""
    ClimaLSM.auxiliary_types(model::PlantHydraulicsModel{FT}) where {FT}

Defines the auxiliary types for the PlantHydraulicsModel.
"""
ClimaLSM.auxiliary_types(model::PlantHydraulicsModel{FT}) where {FT} = (
    FT,
    NTuple{model.n_stem + model.n_leaf, FT},
    NTuple{model.n_stem + model.n_leaf, FT},
    FT,
    NamedTuple{(:root, :stem, :leaf), Tuple{FT, FT, FT}},
)


"""
    set_canopy_prescribed_field!(component::PlantHydraulics{FT},
                                 p,
                                 t0,
                                 ) where {FT}


Sets the canopy prescribed fields pertaining to the PlantHydraulics
component (the area indices) with their initial values at time t0.
"""
function ClimaLSM.Canopy.set_canopy_prescribed_field!(
    component::PlantHydraulicsModel{FT},
    p,
    t0,
) where {FT}
    (; LAIfunction, SAI, RAI) = component.parameters.ai_parameterization
    @. p.canopy.hydraulics.area_index.leaf = LAIfunction(t0)
    @. p.canopy.hydraulics.area_index.stem = SAI
    @. p.canopy.hydraulics.area_index.root = RAI
end


"""
    update_canopy_prescribed_field!(component::PlantHydraulics{FT},
                                    p,
                                    t,
                                    ) where {FT}

Updates the canopy prescribed fields pertaining to the PlantHydraulics
component (the LAI only in this case) with their values at time t.    
"""
function ClimaLSM.Canopy.update_canopy_prescribed_field!(
    component::PlantHydraulicsModel{FT},
    p,
    t,
) where {FT}
    (; LAIfunction) = component.parameters.ai_parameterization
    @. p.canopy.hydraulics.area_index.leaf = LAIfunction(t)
end


"""
    flux(
        z1,
        z2,
        ψ1,
        ψ2,
        K1,
        K2,
    ) where {FT}

Computes the water flux given the absolute potential (pressure/(ρg))
 at the center of the two compartments z1 and z2,
and the conductivity along
the flow path between these two points.

We currently assuming an arithmetic
mean for mean K_sat between the two points (Bonan, 2019; Zhu, 2008)
to take into account the change in K_sat halfway
between z1 and z2; this is incorrect for compartments of differing sizes.
"""
function flux(z1, z2, ψ1, ψ2, K1, K2)
    flux = -(K1 + K2) / 2 * ((ψ2 - ψ1) / (z2 - z1) + 1)
    return flux # (m/s)
end

"""
    AbstractConductivityModel{FT <: AbstractFloat}

An abstract type for the plant hydraulics conductivity model.
"""
abstract type AbstractConductivityModel{FT <: AbstractFloat} end

Base.broadcastable(x::AbstractConductivityModel) = tuple(x)
"""
    AbstractRetentionModel{FT <: AbstractFloat}

An abstract type for the plant retention curve model.
"""
abstract type AbstractRetentionModel{FT <: AbstractFloat} end

Base.broadcastable(x::AbstractRetentionModel) = tuple(x)


"""
    Weibull{FT} <: AbstractConductivityModel{FT}

A concrete type specifying that a Weibull conductivity model is to be used;
the struct contains the require parameters for this model.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct Weibull{FT} <: AbstractConductivityModel{FT}
    "Maximum Water conductivity in the above-ground plant compartments (m/s) at saturation"
    K_sat::FT
    "The absolute water potential in xylem (or xylem water potential) at which ∼63%
    of maximum xylem conductance is lost (Liu, 2020)."
    ψ63::FT
    "Weibull parameter c, which controls shape the shape of the conductance curve (Sperry, 2016)."
    c::FT
end

"""
    hydraulic_conductivity(conductivity_params::Weibull{FT}, ψ::FT) where {FT}

Computes the hydraulic conductivity at a point, using the
Weibull formulation, given the potential ψ.
"""
function hydraulic_conductivity(
    conductivity_params::Weibull{FT},
    ψ::FT,
) where {FT}
    (; K_sat, ψ63, c) = conductivity_params
    if ψ <= FT(0)
        K = exp(-(ψ / ψ63)^c)
    else
        K = FT(1)
    end
    return K * K_sat # (m/s)
end

"""
    LinearRetentionCurve{FT} <: AbstractRetentionModel{FT}

A concrete type specifying that a linear water retention  model is to be used;
the struct contains the require parameters for this model.

When ψ = 0, the effective saturation is one, so the intercept
is not a free parameter, and only the slope must be specified.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct LinearRetentionCurve{FT} <: AbstractRetentionModel{FT}
    "Bulk modulus of elasticity and slope of potential to volume curve. See also Corcuera, 2002, and Christoffersen, 2016."
    a::FT
end

"""
    water_retention_curve(
        S_l::FT,
        b::FT,
        ν::FT,
        S_s::FT) where {FT}

Returns the potential ψ given the effective saturation S at a point, according
to a linear model for the retention curve with parameters specified
by `retention_params`.
"""
function water_retention_curve(
    retention_params::LinearRetentionCurve{FT},
    S_l::FT,
    ν::FT,
    S_s::FT,
) where {FT}
    (; a) = retention_params
    if S_l <= FT(1)
        ψ = 1 / a * (S_l - 1) # ψ(S_l=1)=0.
    else
        ϑ_l = augmented_liquid_fraction(ν, S_l)
        ψ = (ϑ_l - ν) / S_s
    end
    return ψ # (m)
end

"""
    inverse_water_retention_curve(
        ψ::FT,
        b::FT,
        ν::FT,
        S_s::FT) where {FT}

Returns the effective saturation given the potential at a point, according
to the linear retention curve model.
"""
function inverse_water_retention_curve(
    retention_params::LinearRetentionCurve{FT},
    ψ::FT,
    ν::FT,
    S_s::FT,
) where {FT}
    (; a) = retention_params
    if ψ <= FT(0)
        S_l = ψ * a + 1 # ψ(S_l=1)=0.
    else
        ϑ_l = ψ * S_s + ν
        S_l = effective_saturation(ν, ϑ_l)
    end
    return S_l #(m3/m3)
end

"""
    augmented_liquid_fraction(
        ν::FT,
        S_l::FT) where {FT}

Computes the augmented liquid fraction from porosity and
effective saturation.

Augmented liquid fraction allows for
oversaturation: an expansion of the volume of space
available for storage in a plant compartment.
"""
function augmented_liquid_fraction(ν::FT, S_l::FT) where {FT}
    ϑ_l = S_l * ν # ϑ_l can be > ν
    safe_ϑ_l = max(ϑ_l, eps(FT))
    return safe_ϑ_l # (m3 m-3)
end

"""
    effective_saturation(
        ν::FT,
        ϑ_l::FT) where {FT}

Computes the effective saturation given the augmented liquid fraction.
"""
function effective_saturation(ν::FT, ϑ_l::FT) where {FT}
    S_l = ϑ_l / ν # S_l can be > 1
    safe_S_l = max(S_l, eps(FT))
    return safe_S_l # (m3 m-3)
end

"""
    make_compute_exp_tendency(model::PlantHydraulicsModel, _)

A function which creates the compute_exp_tendency! function for the PlantHydraulicsModel.
The compute_exp_tendency! function must comply with a rhs function of OrdinaryDiffEq.jl.

Below, `fa` denotes a flux multiplied by the relevant cross section 
(per unit area ground, or area index, AI). The tendency for the 
ith compartment can be written then as:
∂ϑ[i]/∂t = 1/(AI*dz)[fa[i]-fa[i+1]).

Note that if the area_index is zero because no plant is present, 
AIdz is zero, and the fluxes `fa` appearing in the numerator are 
zero because they are scaled by AI.

To prevent dividing by zero, we change AI/(AI x dz)" to
"AI/max(AI x dz, eps(FT))"
"""
function make_compute_exp_tendency(model::PlantHydraulicsModel, canopy)
    function compute_exp_tendency!(dY, Y, p, t)
        area_index = p.canopy.hydraulics.area_index
        n_stem = model.n_stem
        n_leaf = model.n_leaf
        fa = p.canopy.hydraulics.fa
        fa_roots = p.canopy.hydraulics.fa_roots
        FT = eltype(t)
        @inbounds for i in 1:(n_stem + n_leaf)
            # To prevent dividing by zero, change AI/(AI x dz)" to
            # "AI/max(AI x dz, eps(FT))"
            AIdz =
                max.(
                    getproperty(area_index, model.compartment_labels[i]) .* (
                        model.compartment_surfaces[i + 1] -
                        model.compartment_surfaces[i]
                    ),
                    eps(FT),
                )
            if i == 1
                # All fluxes `fa` are per unit area of ground
                root_flux_per_ground_area!(
                    fa_roots,
                    canopy.soil_driver,
                    model,
                    Y,
                    p,
                    t,
                )
                @inbounds @. dY.canopy.hydraulics.ϑ_l[i] =
                    1 / AIdz * (fa_roots - fa[i])
            else
                @inbounds @. dY.canopy.hydraulics.ϑ_l[i] =
                    1 / AIdz * (fa[i - 1] - fa[i])
            end
        end
    end
    return compute_exp_tendency!
end

"""
    root_flux_per_ground_area!(
        fa::ClimaCore.Fields.Field,
        s::PrescribedSoil{FT},
        model::PlantHydraulicsModel{FT},
        Y::ClimaCore.Fields.FieldVector,
        p::ClimaCore.Fields.FieldVector,
        t::FT,
    )::FT where {FT}

A method which computes the flux between the soil and the stem, via the roots,
and multiplied by the RAI, in the case of a model running without an integrated
soil model.

The returned flux is per unit ground area. This assumes that the stem compartment
is the first element of `Y.canopy.hydraulics.ϑ_l`.
"""
function root_flux_per_ground_area!(
    fa::ClimaCore.Fields.Field,
    s::PrescribedSoil{FT},
    model::PlantHydraulicsModel{FT},
    Y::ClimaCore.Fields.FieldVector,
    p::ClimaCore.Fields.FieldVector,
    t::FT,
) where {FT}

    (; conductivity_model, root_distribution) = model.parameters
    area_index = p.canopy.hydraulics.area_index
    ψ_base = p.canopy.hydraulics.ψ[1]
    root_depths = s.root_depths
    n_root_layers = length(root_depths)
    ψ_soil::FT = s.ψ_soil(t)
    fa .= FT(0.0)
    @inbounds for i in 1:n_root_layers
        if i != n_root_layers
            @. fa +=
                flux(
                    root_depths[i],
                    model.compartment_midpoints[1],
                    ψ_soil,
                    ψ_base,
                    hydraulic_conductivity(conductivity_model, ψ_soil),
                    hydraulic_conductivity(conductivity_model, ψ_base),
                ) *
                root_distribution(root_depths[i]) *
                (root_depths[i + 1] - root_depths[i]) *
                (
                    area_index.root +
                    getproperty(area_index, model.compartment_labels[1])
                ) / 2
        else
            @. fa +=
                flux(
                    root_depths[i],
                    model.compartment_midpoints[1],
                    ψ_soil,
                    ψ_base,
                    hydraulic_conductivity(conductivity_model, ψ_soil),
                    hydraulic_conductivity(conductivity_model, ψ_base),
                ) *
                root_distribution(root_depths[i]) *
                (model.compartment_surfaces[1] - root_depths[n_root_layers]) *
                (
                    area_index.root +
                    getproperty(area_index, model.compartment_labels[1])
                ) / 2
        end
    end
end

"""
    PrescribedTranspiration{FT} <: AbstractTranspiration{FT}

A concrete type used for dispatch when computing the transpiration
from the leaves, in the case where transpiration is prescribed.
"""
struct PrescribedTranspiration{FT} <: AbstractTranspiration{FT}
    T::Function
end

"""
    transpiration_per_ground_area(
        transpiration::PrescribedTranspiration{FT},
        Y,
        p,
        t::FT,
    )::FT where {FT}

A method which computes the transpiration in meters/sec between the leaf
and the atmosphere,
in the case of a standalone plant hydraulics model with prescribed
transpiration rate.

Transpiration should be per unit ground area, not per leaf area.
"""
function transpiration_per_ground_area(
    transpiration::PrescribedTranspiration{FT},
    _,
    _,
    t::FT,
)::FT where {FT}
    return transpiration.T(t) # (m/s)
end

"""
    DiagnosticTranspiration{FT} <: AbstractTranspiration{FT}

A concrete type used for dispatch in the case where transpiration is computed
diagnostically, as a function of prognostic variables and parameters,
and stored in `p` during the `update_aux!` step.
"""
struct DiagnosticTranspiration{FT} <: AbstractTranspiration{FT} end

"""
    transpiration_per_ground_area(transpiration::DiagnosticTranspiration, Y, p, t)

Returns the transpiration computed diagnostically using local conditions.
In this case, it just returns the value which was computed and stored in
the `aux` state during the update_aux! step.

Transpiration should be per unit ground area, not per leaf area.
"""
function transpiration_per_ground_area(
    transpiration::DiagnosticTranspiration,
    Y,
    p,
    t,
)
    @inbounds return p.canopy.conductance.transpiration
end

end
