module PlantHydraulics
using ClimaLSM
using ..ClimaLSM.Canopy: AbstractCanopyComponent
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
    PrescribedSoilPressure,
    PrescribedTranspiration,
    DiagnosticTranspiration,
    AbstractRootExtraction,
    Layers

"""
    AbstractPlantHydraulicsModel{FT} <: AbstractCanopyComponent{FT}

An abstract type for plant hydraulics models.
"""
abstract type AbstractPlantHydraulicsModel{FT} <: AbstractCanopyComponent{FT} end

ClimaLSM.name(::AbstractPlantHydraulicsModel) = :hydraulics

"""
    AbstractRootExtraction{FT <: AbstractFloat}

An abstract type for types representing different models of
water exchange between soil and plants.
Currently, prescribed soil matric potential and prescribed flux models
are supported for standalone plant hydraulics.
Use within an LSM requires types defined within ClimaLSM,
and include a prognostic soil pressure for models with
both soil and roots.
"""
abstract type AbstractRootExtraction{FT <: AbstractFloat} end

"""
    AbstractTranspiration{FT <: AbstractFloat}

An abstract type for types representing different models of
transpiration (Prescribed or Diagnostic)
"""
abstract type AbstractTranspiration{FT <: AbstractFloat} end


"""
    PlantHydraulicsParameters{FT <: AbstractFloat}

A struct for holding parameters of the PlantHydraulics Model.
$(DocStringExtensions.FIELDS)
"""
struct PlantHydraulicsParameters{FT <: AbstractFloat}
    "RAI, SAI, LAI in [m2 m-2]"
    area_index::NamedTuple{(:root, :stem, :leaf), Tuple{FT, FT, FT}}
    "water conductivity in the above-ground plant compartments (m/s) when pressure is zero, a maximum"
    K_sat::NamedTuple{(:root, :stem, :leaf), Tuple{FT, FT, FT}}
    "van Genuchten parameter (1/m)"
    vg_α::FT
    "van Genuchten parameter (unitless)"
    vg_n::FT
    "van Genuchten parameter (unitless)"
    vg_m::FT
    "porosity (m3/m3)"
    ν::FT
    "storativity (m3/m3)"
    S_s::FT
    "Root distribution function P(z)"
    root_distribution::Function
end

"""
    PlantHydraulicsModel{FT, PS, RE, T, B} <: AbstractPlantHydraulicsModel{FT}

Defines, and constructs instances of, the PlantHydraulicsModel type, which is used
for simulation flux of water to/from soil, along roots of different depths,
along a stem, to a leaf, and ultimately being lost from the system by
transpiration.

This model can be used in standalone mode by prescribing the transpiration rate
and soil matric potential at the root tips or flux in the roots, or with a
dynamic soil model using `ClimaLSM`.
$(DocStringExtensions.FIELDS)
"""
struct PlantHydraulicsModel{FT, PS, RE, T} <: AbstractPlantHydraulicsModel{FT}
    "The number of stem compartments for the plant"
    n_stem::Int64
    "The number of leaf compartments for the plant"
    n_leaf::Int64
    "The depth of the root tips, in meters"
    root_depths::Vector{FT}
    "The height of the center of each leaf compartment/stem compartment, in meters"
    compartment_midpoints::Vector{FT}
    "The height of the compartments' top faces, in meters"
    compartment_surfaces::Vector{FT}
    "The label (:stem or :leaf) of each compartment"
    compartment_labels::Vector{Symbol}
    "Parameters required by the Plant Hydraulics model"
    parameters::PS
    "The root extraction model, of type `AbstractRootExtraction`"
    root_extraction::RE
    "The transpiration model, of type `AbstractTranspiration`"
    transpiration::T
end

function PlantHydraulicsModel{FT}(;
    root_depths::Vector{FT},
    n_stem::Int64,
    n_leaf::Int64,
    compartment_midpoints::Vector{FT},
    compartment_surfaces::Vector{FT},
    parameters::PlantHydraulicsParameters{FT},
    root_extraction::AbstractRootExtraction{FT},
    transpiration::AbstractTranspiration{FT} = DiagnosticTranspiration{FT}(),
) where {FT}
    args = (parameters, root_extraction, transpiration)
    @assert n_leaf != 0
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
        root_depths,
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
auxiliary_vars(model::PlantHydraulicsModel) = (:β, :ψ, :fa, :fa_roots)

"""
    Base.zero(x::Type{NTuple{N, FT}}) where {N, FT}

A `zero` method for NTuples of floats.

Used when initializing the state vector, when it consists
of an NTuple at each coordinate point.
"""
function Base.zero(x::Type{NTuple{N, FT}}) where {N, FT}
    ntuple(i -> FT(0), N)
end

"""

    Copied over from EDMF. This lets us index into ClimaCore
Fields of NTuples.
"""
Base.@propagate_inbounds Base.getindex(
    field::ClimaCore.Fields.Field,
    i::Integer,
) = Base.getproperty(field, i)

"""
    Copied over from EDMF. This is a helper struct which
lets us set indices of ClimaCore
Fields of Ntuples by indexing.
"""
struct Cent{I <: Integer}
    i::I
end

"""
    Copied over from EDMF. This lets us set components of ClimaCore
Fields of NTuples by indexing.
"""
Base.@propagate_inbounds Base.setindex!(
    field::ClimaCore.Fields.Field,
    v,
    i::Cent,
) = Base.setindex!(ClimaCore.Fields.field_values(field), v, i.i)

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
)

"""
    flux(
        z1,
        z2,
        ψ1,
        ψ2,
        vg_α,
        vg_n,
        vg_m,
        ν,
        S_s,
        K_sat1,
        K_sat2,
) where {FT}

Computes the water flux given 1) the absolute pressure (represented by corresponding
water column height) at two points located at z1 and z2,
corresponding here to the middle of each compartment
(for eg. a stem or leaf compartment), and 2) the maximum conductivity along
the flow path between these two points, assuming an arithmetic
mean for mean K_sat between the two points (Bonan, 2019; Zhu, 2008)
to take into account the change in K_sat halfway
between z1 and z2. The expression is given in full in the Clima docs
in section 6.4 "Bulk plant hydraulics model".
"""
function flux(z1, z2, ψ1, ψ2, vg_α, vg_n, vg_m, ν, S_s, K_sat1, K_sat2)

    S_l1 = inverse_water_retention_curve(vg_α, vg_n, vg_m, ψ1, ν, S_s)
    S_l2 = inverse_water_retention_curve(vg_α, vg_n, vg_m, ψ2, ν, S_s)

    K1_ψ1 = hydraulic_conductivity(K_sat1, vg_m, S_l1)
    K2_ψ2 = hydraulic_conductivity(K_sat2, vg_m, S_l2)

    flux = -(K1_ψ1 + K2_ψ2) / 2 * ((ψ2 - ψ1) / (z2 - z1) + 1)

    return flux # (m/s)
end

"""
     hydraulic_conductivity(K_sat::FT, vg_m::FT, S::FT) where {FT}

A point-wise function returning the hydraulic conductivity, using the
van Genuchten formulation.
"""
function hydraulic_conductivity(K_sat::FT, vg_m::FT, S::FT) where {FT}
    if S < FT(1)
        K = sqrt(S) * (FT(1) - (FT(1) - S^(FT(1) / vg_m))^vg_m)^FT(2)
    else
        K = FT(1)
    end
    return K * K_sat # (m/s)
end

"""
    augmented_liquid_fraction(
        ν::FT,
        S_l::FT) where {FT}

Computes the augmented liquid fraction from porosity and
effective saturation. Augmented liquid fraction allows for
oversaturation: an expansion of the volume of space
available for storage in a plant compartment. It is analogous
to the augmented liquid fraction state variable in the soil model.
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

Computes the effective saturation (volumetric water content relative to
porosity).
"""
function effective_saturation(ν::FT, ϑ_l::FT) where {FT}
    S_l = ϑ_l / ν # S_l can be > 1
    safe_S_l = max(S_l, eps(FT))
    return safe_S_l # (m3 m-3)
end

"""
    water_retention_curve(
        vg_α::FT,
        vg_n::FT,
        vg_m::FT,
        S_l::FT,
        ν::FT,
        S_s::FT) where {FT}

Converts augmented liquid fraction (ϑ_l) to effective saturation
(S_l), and then effective saturation to absolute pressure, represented by
the height (ψ) of the water column that would give rise to this pressure.
Pressure for both the unsaturated and saturated regimes are calculated.
Units are in meters.
"""
function water_retention_curve(
    vg_α::FT,
    vg_n::FT,
    vg_m::FT,
    S_l::FT,
    ν::FT,
    S_s::FT,
) where {FT}
    if S_l <= FT(1.0)
        ψ = -((S_l^(-FT(1) / vg_m) - FT(1)) * vg_α^(-vg_n))^(FT(1) / vg_n)
    else
        ϑ_l = augmented_liquid_fraction(ν, S_l)
        ψ = (ϑ_l - ν) / S_s
    end
    return ψ # (m)
end

"""
    inverse_water_retention_curve(
        vg_α::FT,
        vg_n::FT,
        vg_m::FT,
        ψ::FT,
        ν::FT,
        S_s::FT) where {FT}

Converts absolute pressure, represented by the height (ψ)
of the water column that would give rise to this pressure,
to effective saturation (S_l).
"""
function inverse_water_retention_curve(
    vg_α::FT,
    vg_n::FT,
    vg_m::FT,
    ψ::FT,
    ν::FT,
    S_s::FT,
) where {FT}
    if ψ <= FT(0.0)
        S_l = ((-vg_α * ψ)^vg_n + FT(1.0))^(-vg_m)
    else
        ϑ_l = ψ * S_s + ν
        S_l = effective_saturation(ν, ϑ_l)
    end
    return S_l
end

"""
    make_compute_exp_tendency(model::PlantHydraulicsModel, _)

A function which creates the compute_exp_tendency! function for the PlantHydraulicsModel.
The compute_exp_tendency! function must comply with a rhs function of OrdinaryDiffEq.jl.

Below, `fa` denotes a flux multiplied by the relevant cross section (per unit area ground).
"""
function make_compute_exp_tendency(model::PlantHydraulicsModel, _)
    function compute_exp_tendency!(dY, Y, p, t)
        (; vg_α, vg_n, vg_m, S_s, ν, K_sat, area_index) = model.parameters
        n_stem = model.n_stem
        n_leaf = model.n_leaf
        fa = p.canopy.hydraulics.fa
        fa_roots = p.canopy.hydraulics.fa_roots

        @inbounds for i in 1:(n_stem + n_leaf)
            AIdz =
                area_index[model.compartment_labels[i]] * (
                    model.compartment_surfaces[i + 1] -
                    model.compartment_surfaces[i]
                )
            if i == 1
                # All fluxes `fa` are per unit area of ground
                root_flux_per_ground_area!(
                    fa_roots,
                    model.root_extraction,
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
    PrescribedSoilPressure{FT} <: AbstractRootExtraction{FT}

A concrete type used for dispatch when computing the `root_flux_per_ground_area!`,
in the case where the soil matric potential at each root layer is prescribed.
"""
struct PrescribedSoilPressure{FT} <: AbstractRootExtraction{FT}
    ψ_soil::Function
end

"""
    root_flux_per_ground_area!(
        fa::ClimaCore.Fields.Field,
        re::PrescribedSoilPressure{FT},
        model::PlantHydraulicsModel{FT},
        Y::ClimaCore.Fields.FieldVector,
        p::ClimaCore.Fields.FieldVector,
        t::FT,
    )::FT where {FT}

A method which computes the flux between the soil and the stem, via the roots,
and multiplied by the RAI, 
in the case of a standalone plant hydraulics model with prescribed soil pressure (in m)
at the root tips.

The returned flux is per unit ground area. This assumes that the stem compartment 
is the first element of `Y.canopy.hydraulics.ϑ_l`.
"""
function root_flux_per_ground_area!(
    fa::ClimaCore.Fields.Field,
    re::PrescribedSoilPressure{FT},
    model::PlantHydraulicsModel{FT},
    Y::ClimaCore.Fields.FieldVector,
    p::ClimaCore.Fields.FieldVector,
    t::FT,
) where {FT}
    (; vg_α, vg_n, vg_m, ν, S_s, K_sat, area_index, root_distribution) =
        model.parameters
    ψ_base = p.canopy.hydraulics.ψ[1]
    n_root_layers = length(model.root_depths)
    ψ_soil::FT = re.ψ_soil(t)
    @inbounds for i in 1:n_root_layers
        if i != n_root_layers
            @. fa +=
                flux(
                    model.root_depths[i],
                    model.compartment_midpoints[1],
                    ψ_soil,
                    ψ_base,
                    vg_α,
                    vg_n,
                    vg_m,
                    ν,
                    S_s,
                    K_sat[:root],
                    K_sat[model.compartment_labels[1]],
                ) *
                root_distribution(model.root_depths[i]) *
                (model.root_depths[i + 1] - model.root_depths[i]) *
                area_index[:root]
        else
            @. fa +=
                flux(
                    model.root_depths[i],
                    model.compartment_midpoints[1],
                    ψ_soil,
                    ψ_base,
                    vg_α,
                    vg_n,
                    vg_m,
                    ν,
                    S_s,
                    K_sat[:root],
                    K_sat[model.compartment_labels[1]],
                ) *
                root_distribution(model.root_depths[i]) *
                (FT(0) - model.root_depths[n_root_layers]) *
                area_index[:root]

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
