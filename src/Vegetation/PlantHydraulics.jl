module PlantHydraulics
#=
    PlantHydraulics
This module contains everything needed to run a vegetation model
in standalone mode.
The vegetation model is assumed to have a set of prognostic `Y` and
auxiliary `p` variables, which describe the state of the
vegetation system. The system is evolved in time by solving 
equations of the form
```
\frac{d Y}{d t} = f(Y, t, p(Y, t; \ldots);\ldots),
```
i.e. ordinary differential equations depending on state `Y`,
auxiliary functions of the state `p`, and other parameters 
represented by the ellipses. For example, `p` may represent
the transpiration rate, which must be computed each time step
based on the vegetation prognostic state, functions of time, 
like the atmosphere state, and other parameters.
Currently, only a simple plant hydraulics model is supported,
but our plan is to include much more complex representations
of the vegetation.
Addition of additional versions of vegetation
models requires defining a model type (of super type 
`AbstractVegetationModel`), and extending the methods
imported by Models.jl, as needed, for computing the
right hand side functions of the ordinary differential equations
and the functions which update auxiliary variables whenever
the right hand side is evaluated.
This code base assumes that DifferentialEquations.jl
will be used for evolving the system in time,
and that the array-like object being stepped forward
is a `ClimaCore.Fields.FieldVector`. 
While a simple array may be sufficient for vegetation models,
the `FieldVector` type is used in order to make use of 
`ClimaCore` functionality when solving PDEs (required for
other components of Land Surface Models) and for ease of
handling multi-column models.
To simulate land surfaces with multiple components (vegetation,
soil, rivers, etc), the ClimaLSM.jl package should be used.
That package will use the methods of this function for advancing
the system forward in time, extending methods as needed to account
for interactions between components.
=#
using ClimaLSM
using ClimaCore
using UnPack
using DocStringExtensions
import ClimaCore: Fields

using ClimaLSM.Domains: AbstractDomain
import ClimaLSM:
    AbstractModel,
    initialize_prognostic,
    make_update_aux,
    make_rhs,
    make_ode_function,
    prognostic_vars,
    prognostic_types,
    auxiliary_vars,
    initialize,
    initialize_auxiliary,
    name
export PlantHydraulicsModel,
    AbstractVegetationModel,
    flux,
    effective_saturation,
    augmented_liquid_fraction,
    water_retention_curve,
    inverse_water_retention_curve,
    flux_out_roots,
    PlantHydraulicsParameters,
    PrescribedSoilPressure,
    PrescribedTranspiration,
    AbstractRootExtraction,
    Layers

"""
    AbstractVegetationModel{FT} <: AbstractModel{FT}

An abstract type for vegetation models.
Concrete types include a plant hydraulics model, but future types will
include multi-layer canopy models and possibly a big leaf model.
"""
abstract type AbstractVegetationModel{FT} <: AbstractModel{FT} end

ClimaLSM.name(::AbstractVegetationModel) = :vegetation
ClimaLSM.domain(::AbstractVegetationModel) = :surface

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
transpiration.
Currently, only a PrescribedTranspiration is supported.
"""
abstract type AbstractTranspiration{FT <: AbstractFloat} end


"""
    PlantHydraulicsParameters{FT <: AbstractFloat}

A struct for holding parameters of the PlantHydraulics Model.
$(DocStringExtensions.FIELDS)
"""
struct PlantHydraulicsParameters{FT <: AbstractFloat, PSE}
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
    "Physical Constants and other Clima-wide parameters"
    earth_param_set::PSE
end

"""
    PlantHydraulicsModel{FT, PS, D, RE, T, B} <: AbstractVegetationModel{FT}

Defines, and constructs instances of, the PlantHydraulicsModel type, which is used
for simulation flux of water to/from soil, along roots of different depths,
along a stem, to a leaf, and ultimately being lost from the system by
transpiration. 

This model can be used in standalone mode by prescribing the transpiration rate
and soil matric potential at the root tips or flux in the roots, or with a 
dynamic soil model using `ClimaLSM`.
$(DocStringExtensions.FIELDS)
"""
struct PlantHydraulicsModel{FT, PS, D, RE, T} <: AbstractVegetationModel{FT}
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
    "The plant hydraulics model domain, of type `AbstractDomain`"
    domain::D
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
    parameters::PlantHydraulicsParameters{FT, PSE},
    domain::AbstractDomain{FT},
    root_extraction::AbstractRootExtraction{FT},
    transpiration::AbstractTranspiration{FT},
) where {FT, PSE}
    args = (parameters, domain, root_extraction, transpiration)
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
the water potential (m) and volume flux*cross section `fa` (1/s),
where the cross section can be represented by an area index.
"""
auxiliary_vars(model::PlantHydraulicsModel) = (:ψ, :fa)

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
    NTuple{model.n_stem + model.n_leaf, FT},
    NTuple{model.n_stem + model.n_leaf, FT},
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
    safe_ϑ_l = max(ϑ_l, 0)
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
    safe_S_l = max(S_l, 0)
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
    make_rhs(model::VegetationModel)

A function which creates the rhs! function for the PlantHydraulicsModel.
The rhs! function must comply with a rhs function of OrdinaryDiffEq.jl.

Note that this rhs function updates the auxiliary variables as well. We chose to
do so here because it reduces to a single loop over the layers, instead of looping
over them in the auxiliary function as well. 

Below, `fa` denotes a flux multiplied by the relevant cross section (per unit area ground).
"""
function make_rhs(model::PlantHydraulicsModel)
    function rhs!(dY, Y, p, t)
        @unpack vg_α, vg_n, vg_m, S_s, ν, K_sat, area_index = model.parameters
        n_stem = model.n_stem
        n_leaf = model.n_leaf

        ψ = p.vegetation.ψ
        ϑ_l = Y.vegetation.ϑ_l
        fa = p.vegetation.fa

        # initialize first index
        @inbounds @. ψ[1] = water_retention_curve(
            vg_α,
            vg_n,
            vg_m,
            effective_saturation(ν, ϑ_l[1]),
            ν,
            S_s,
        )
        # boundary flux * cross section at the bottom
        fa0 =
            flux_out_roots(model.root_extraction, model, Y, p, t) .*
            area_index[:root]

        @inbounds for i in 1:(n_stem + n_leaf)
            AIdz =
                area_index[model.compartment_labels[i]] * (
                    model.compartment_surfaces[i + 1] -
                    model.compartment_surfaces[i]
                )

            # If we arent at the top compartment, compute the flux*area
            # between the current compartment and the one above it
            if i != (n_stem + n_leaf)
                # Compute potential in compartment above - the potential in current
                # compartment was already computed during the last iteration.
                @. ψ[i + 1] = water_retention_curve(
                    vg_α,
                    vg_n,
                    vg_m,
                    effective_saturation(ν, ϑ_l[i + 1]),
                    ν,
                    S_s,
                )
                # Compute the flux*area between the current compartment `i`
                # and the compartment above.
                @. fa[i] =
                    flux(
                        model.compartment_midpoints[i],
                        model.compartment_midpoints[i + 1],
                        ψ[i],
                        ψ[i + 1],
                        vg_α,
                        vg_n,
                        vg_m,
                        ν,
                        S_s,
                        K_sat[model.compartment_labels[i]],
                        K_sat[model.compartment_labels[i + 1]],
                    ) * (
                        area_index[model.compartment_labels[i]] +
                        area_index[model.compartment_labels[i + 1]]
                    ) / 2
            else  # Apply upper boundary condition * area index
                fa[i] .=
                    transpiration(model.transpiration, t) .* area_index[:leaf]
            end

            if i == 1
                @. dY.vegetation.ϑ_l[i] = 1 / AIdz * (fa0 - fa[i])
            else
                @. dY.vegetation.ϑ_l[i] = 1 / AIdz * (fa[i - 1] - fa[i])
            end
        end
    end

    return rhs!
end

"""
    PrescribedSoilPressure{FT} <: AbstractRootExtraction{FT}

A concrete type used for dispatch when computing the `flux_out_roots`,
in the case where the soil matric potential at each root layer is prescribed.
"""
struct PrescribedSoilPressure{FT} <: AbstractRootExtraction{FT}
    ψ_soil::Function
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
    flux_out_roots(
        re::PrescribedSoilPressure{FT},
        model::PlantHydraulicsModel{FT},
        Y::ClimaCore.Fields.FieldVector,
        p::ClimaCore.Fields.FieldVector,
        t::FT,
    )::FT where {FT}

A method which computes the flux between the soil and the stem, via the roots,
in the case of a standalone plant hydraulics model with prescribed soil pressure (in m)
at the root tips.
This assumes that the stem compartment is the first element of `Y.vegetation.ϑ_l`.
"""
function flux_out_roots(
    re::PrescribedSoilPressure{FT},
    model::PlantHydraulicsModel{FT},
    Y::ClimaCore.Fields.FieldVector,
    p::ClimaCore.Fields.FieldVector,
    t::FT,
)::ClimaCore.Fields.Field where {FT}
    @unpack vg_α, vg_n, vg_m, ν, S_s, K_sat, root_distribution =
        model.parameters
    ψ_stem = p.vegetation.ψ[1]
    n_root_layers = length(model.root_depths)
    sum_flux_out_roots = similar(ψ_stem)
    dz_roots = (
        vcat(model.root_depths, [0.0])[2:end] -
        vcat(model.root_depths, [0.0])[1:(end - 1)]
    )
    ψ_soil = re.ψ_soil(t)
    @inbounds for i in 1:n_root_layers
        @. sum_flux_out_roots +=
            flux(
                model.root_depths[i],
                model.compartment_midpoints[1],
                ψ_soil,
                ψ_stem,
                vg_α,
                vg_n,
                vg_m,
                ν,
                S_s,
                K_sat[:root],
                K_sat[:stem],
            ) *
            root_distribution(model.root_depths[i]) *
            dz_roots[i]
    end

    return sum_flux_out_roots
end

"""
    transpiration(
        transpiration::PrescribedTranspiration{FT},
        t::FT,
    )::FT where {FT}

A method which computes the transpiration in meters/sec between the leaf
and the atmosphere,
in the case of a standalone plant hydraulics model with prescribed
transpiration rate.
"""
function transpiration(
    transpiration::PrescribedTranspiration{FT},
    t::FT,
)::FT where {FT}
    return transpiration.T(t) # (m/s)
end

end
