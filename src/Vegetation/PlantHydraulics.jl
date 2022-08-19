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

using ClimaLSM.Domains: AbstractVegetationDomain, PlantHydraulicsDomain
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
    PrescribedRootFlux,
    PrescribedTranspiration,
    AbstractRootExtraction

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
Currently, prescribed soil pressure and prescribed flux models 
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
    "water conductivity in the different plant compartments (m/s) when pressure is zero, a maximum"
    K_sat::NamedTuple{(:root, :stem, :leaf), Tuple{FT, FT, FT}}
    "van Genuchten parameter"
    vg_α::FT
    "van Genuchten parameter"
    vg_n::FT
    "van Genuchten parameter"
    vg_m::FT
    "porosity"
    ν::FT
    "storativity"
    S_s::FT
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
and soil pressure at the root tips or flux in the roots, or with a 
dynamic soil model using `ClimaLSM`.
$(DocStringExtensions.FIELDS)
"""
struct PlantHydraulicsModel{FT, PS, D, RE, T} <: AbstractVegetationModel{FT}
    "Parameters required by the Plant Hydraulics model"
    parameters::PS
    "The plant hydraulics model domain, of type `AbstractVegetationDomain`"
    domain::D
    "The root extraction model, of type `AbstractRootExtraction`"
    root_extraction::RE
    "The transpiration model, of type `AbstractTranspiration`"
    transpiration::T
end

function PlantHydraulicsModel{FT}(;
    parameters::PlantHydraulicsParameters{FT, PSE},
    domain::AbstractVegetationDomain{FT},
    root_extraction::AbstractRootExtraction{FT},
    transpiration::AbstractTranspiration{FT},
) where {FT, PSE}
    args = (parameters, domain, root_extraction, transpiration)
    return PlantHydraulicsModel{FT, typeof.(args)...}(args...)
end

"""
    prognostic_vars(model::PlantHydraulicsModel)
A function which returns the names of the prognostic 
variables of the `PlantHydraulicsModel`.
"""
prognostic_vars(model::PlantHydraulicsModel) = (:ϑ_l,)
prognostic_types(model::PlantHydraulicsModel{FT}) where {FT} = (FT,)

"""
    flux(
        z1,
        z2,
        p1,
        p2,
        vg_α,
        vg_n,
        vg_m,
        ν,
        S_s,
        K_sat1,
        K_sat2,
) where {FT} 
Computes the water flux given 1) the absolute pressure at two points located
at z1 and z2, corresponding here to the middle of each compartment 
(for eg. a stem or leaf compartment), and 2) the maximum conductivity along 
the flow path between these two points, assuming an arithmetic 
mean for mean K_sat between the two points (Bonan, 2019; Zhu, 2008) 
to take into account the change in K_sat halfway 
between z1 and z2. The expression is given in full in the Clima docs
in section 6.4 "Bulk plant hydraulics model".
"""
function flux(
    z1,
    z2,
    p1,
    p2,
    vg_α,
    vg_n,
    vg_m,
    ν,
    S_s,
    K_sat1,
    K_sat2,
) where {FT}

    S_l1 = inverse_water_retention_curve(vg_α, vg_n, vg_m, p1, ν, S_s)
    S_l2 = inverse_water_retention_curve(vg_α, vg_n, vg_m, p2, ν, S_s)

    K1_p1 = hydraulic_conductivity(K_sat1, vg_m, S_l1)
    K2_p2 = hydraulic_conductivity(K_sat2, vg_m, S_l2)

    flux = -(K1_p1 + K2_p2) / 2 * ((p2 - p1) / (z2 - z1) + 1)

    return flux # units of [m s-1]
end

"""
     hydraulic_conductivity(K_sat::FT, m::FT, S::FT) where {FT}

A point-wise function returning the hydraulic conductivity, using the
van Genuchten formulation.
"""
function hydraulic_conductivity(K_sat::FT, m::FT, S::FT) where {FT}
    if S < FT(1)
        K = sqrt(S) * (FT(1) - (FT(1) - S^(FT(1) / m))^m)^FT(2)
    else
        K = FT(1)
    end
    return K * K_sat
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
    return safe_ϑ_l # units of [m3 m-3]
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
    return safe_S_l # units of [m3 m-3]
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
(S_l), and then effective saturation to absolute pressure (p), 
for both the unsaturated and saturated regimes. Note pressure 
is in SI units [m].
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
        p = -((S_l^(-FT(1) / vg_m) - FT(1)) * vg_α^(-vg_n))^(FT(1) / vg_n)
    else
        ϑ_l = augmented_liquid_fraction(ν, S_l)
        p = (ϑ_l - ν) / S_s
    end
    return p # units of (m)
end

"""
    inverse_water_retention_curve(
        vg_α::FT,
        vg_n::FT,
        vg_m::FT,
        p::FT,
        ν::FT,
        S_s::FT) where {FT}
Converts absolute pressure (p) to effective saturation (S_l).
"""
function inverse_water_retention_curve(
    vg_α::FT,
    vg_n::FT,
    vg_m::FT,
    p::FT,
    ν::FT,
    S_s::FT,
) where {FT}
    if p <= FT(0.0)
        S_l = ((-vg_α * p)^vg_n + FT(1.0))^(-vg_m)
    else
        ϑ_l = p * S_s + ν
        S_l = effective_saturation(ν, ϑ_l)
    end
    return S_l
end

"""
    make_rhs(model::VegetationModel)
A function which creates the rhs! function for the PlantHydraulicsModel.
The rhs! function must comply with a rhs function of OrdinaryDiffEq.jl.
"""
function make_rhs(model::PlantHydraulicsModel)
    function rhs!(dY, Y, p, t)
        @unpack vg_α, vg_n, vg_m, S_s, ν, K_sat = model.parameters

        n_stem = model.domain.n_stem
        n_leaf = model.domain.n_leaf

        S_l = effective_saturation.(ν, Y.vegetation.ϑ_l)
        pressure = water_retention_curve.(vg_α, vg_n, vg_m, S_l, ν, S_s)
        for i in 1:(n_stem + n_leaf)
            if i == 1
                flux_in = flux_out_roots(model.root_extraction, model, Y, p, t)
            else
                flux_in = flux(
                    model.domain.compartment_midpoints[i - 1],
                    model.domain.compartment_midpoints[i],
                    pressure[i - 1],
                    pressure[i],
                    vg_α,
                    vg_n,
                    vg_m,
                    ν,
                    S_s,
                    K_sat[model.domain.compartment_labels[i - 1]],
                    K_sat[model.domain.compartment_labels[i]],
                )
            end

            if i == (n_stem + n_leaf)
                flux_out = transpiration(model.transpiration, t)
            else
                flux_out = flux(
                    model.domain.compartment_midpoints[i],
                    model.domain.compartment_midpoints[i + 1],
                    pressure[i],
                    pressure[i + 1],
                    vg_α,
                    vg_n,
                    vg_m,
                    ν,
                    S_s,
                    K_sat[model.domain.compartment_labels[i]],
                    K_sat[model.domain.compartment_labels[i + 1]],
                )
            end
            dY.vegetation.ϑ_l[i] =
                1 / (
                    model.domain.compartment_surfaces[i + 1] -
                    model.domain.compartment_surfaces[i]
                ) * (flux_in - flux_out)
        end
    end

    return rhs!
end

"""
    PrescribedSoilPressure{FT} <: AbstractRootExtraction{FT}
A concrete type used for dispatch when computing the `flux_out_roots`,
in the case where the soil pressure at each root layer is prescribed.
"""
struct PrescribedSoilPressure{FT} <: AbstractRootExtraction{FT}
    p_soil::Function
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
)::FT where {FT}
    @unpack vg_α, vg_n, vg_m, ν, S_s, K_sat = model.parameters

    S_l_stem = effective_saturation(ν, Y.vegetation.ϑ_l[1])

    #How to ungeneralize

    p_stem = water_retention_curve.(vg_α, vg_n, vg_m, S_l_stem, ν, S_s)
    return sum(
        flux.(
            model.domain.root_depths,
            model.domain.compartment_midpoints[1],
            re.p_soil(t),
            p_stem,
            vg_α,
            vg_n,
            vg_m,
            ν,
            S_s,
            K_sat[:root],
            K_sat[model.domain.compartment_labels[1]],
        ),
    )
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
    return transpiration.T(t)
end

end
