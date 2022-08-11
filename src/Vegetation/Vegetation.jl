# To do: move constants out of src/ code and pass through,via structs. 
# Convert from a single plant to a bulk plant (units change, required
# input changes). 
module Vegetation
#=
    Vegetation
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
and that the array-like objected being stepped forward
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

using ClimaLSM.Domains: AbstractVegetationDomain, VegetationDomain
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
export VegetationModel,
    AbstractVegetationModel,
    flux,
    effective_saturation,
    augmented_liquid_fraction,
    van_Genuchten_volume_to_pressure,
    van_Genuchten_pressure_to_volume,
    volume_to_pressure,
    pressure_to_volume,
    flux_out_roots,
    VegetationParameters,
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
Currently, only a prescribed soil pressure is supported for
standalone plant hydraulics. Use within an LSM requires
types defined within ClimaLSM, and include a prognostic
soil pressure for models with both soil and roots.
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
    VegetationParameters{FT <: AbstractFloat}
A struct for holding parameters of the Vegetation Model.
$(DocStringExtensions.FIELDS)
"""
struct VegetationParameters{FT <: AbstractFloat, PSE}
    "controls the shape and steepness of conductivity vs. absolute pressure curve, for roots: unitless"
    a_root::FT
    "controls the steepness of the relative conductivity vs. absolute pressure curve, for roots: inverse m"
    b_root::FT
    "controls the shape and steepness of relative conductivity vs. absolute pressure curve, for stems: unitless"
    a_stem::FT
    "controls the steepness of the conductivity vs. pressure curve, for stems: inverse m"
    b_stem::FT
    "controls the shape and steepness of relative conductivity vs. pressure curve, for leaves: unitless"
    a_leaf::FT
    "controls the steepness of the conductivity vs. pressure curve, for leaves: inverse m"
    b_leaf::FT
    "water conductivity in roots (m/s) when pressure is zero, a maximum"
    K_sat_root::FT
    "water conductivity in stems (m/s) when pressure is zero, a maximum"
    K_sat_stem::FT
    "water conductivity in leaves (m/s) when pressure is zero, a maximum"
    K_sat_leaf::FT
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
    VegetationModel{FT, PS, D, RE, T, B} <: AbstractVegetationModel{FT}
Defines, and constructs instances of, the VegetationModel type, which is used
for simulation flux of water to/from soil, along roots of different depths,
along a stem, to a leaf, and ultimately being lost from the system by
transpiration. 
This model can be used in standalone mode by prescribing the transpiration rate
and soil pressure at the root tips or flux in the roots, or with a 
dynamic soil model using `ClimaLSM`.
$(DocStringExtensions.FIELDS)
"""
struct VegetationModel{FT, PS, D, RE, T} <: AbstractVegetationModel{FT}
    "Parameters required by the vegetation model"
    parameters::PS
    "The vegetation model domain, of type `AbstractVegetationDomain`"
    domain::D
    "The root extraction model, of type `AbstractRootExtraction`"
    root_extraction::RE
    "The transpiration model, of type `AbstractTranspiration`"
    transpiration::T
end

function VegetationModel{FT}(;
    parameters::VegetationParameters{FT, PSE},
    domain::AbstractVegetationDomain{FT},
    root_extraction::AbstractRootExtraction{FT},
    transpiration::AbstractTranspiration{FT},
) where {FT, PSE}
    args = (parameters, domain, root_extraction, transpiration)
    return VegetationModel{FT, typeof.(args)...}(args...)
end

"""
    prognostic_vars(model::VegetationModel)
A function which returns the names of the prognostic 
variables of the `VegetationModel`.
"""
prognostic_vars(model::VegetationModel) = (:ϑ_l,)
prognostic_types(model::VegetationModel{FT}) where {FT} = (FT,)

"""
    flux(
        z1::FT,
        z2::FT, 
        p1::FT, 
        p2::FT, 
        a1::FT, 
        a2::FT, 
        b1::FT, 
        b2::FT, 
        K_sat1::FT, 
        K_sat2::FT) where {FT} 
Computes the water flux given 1) the absolute pressure at two points located
at z1 and z2, corresponding here to the middle of each compartment 
(for eg. a stem or leaf compartment), and 2) the maximum conductivity along 
the flow path between these two points, assuming an arithmetic 
mean for mean K_sat between the two points (Bonan, 2019; Zhu, 2008) 
to take into account the change in K_sat halfway 
betwwen z1 and z2. The expression is given in full in the Clima docs
in section 6.4 "Bulk plant hydraulics model".
"""
function flux(
    z1::FT,
    z2::FT, 
    p1::FT, 
    p2::FT, 
    a1::FT, 
    a2::FT, 
    b1::FT, 
    b2::FT, 
    K_sat1::FT, 
    K_sat2::FT) where {FT}  
        u1 = a1 * exp(b1 * p1) 
        u2 = a1 * exp(b1 * p2) 
        num1 = log(u1 + FT(1))
        num2 = log(u2 + FT(1))
        c1 = K_sat1 * (a1 + FT(1)) / a1
        term1 = -c1 / b1 * (num2 - num1) /(z2 - z1)

        c2 = K_sat2 * (a2 + FT(1)) / a2
        term2_up = c2 * (FT(1) - FT(1) / (FT(1) + a2*exp(b2 * p2)))
        term2_do = c1 * (FT(1) - FT(1) / (FT(1) + a1*exp(b1 * p1)))
        term2 = -(term2_up - term2_do)/2
        flux = term1 + term2  
    return flux # units of [m s-1]
end

"""
    augmented_liquid_fraction(
        ν::FT, 
        S_l::FT) where {FT}
Computes the augmented liquid fraction from porosity and
effective saturation. Augmented liquid fraction allows for
oversaturation: an expansion of the volume of space
available for storage in a plant compartment.
"""
function augmented_liquid_fraction(
    ν::FT, 
    S_l::FT) where {FT}
    ϑ_l = S_l * ν 
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
function effective_saturation(
    ν::FT, 
    ϑ_l::FT) where {FT}
    S_l = ϑ_l / ν # S_l can be > 1
    safe_S_l = max(S_l, 0)
    return safe_S_l # units of [m3 m-3]
end

"""
    van_Genuchten_volume_to_pressure(
        α::FT, 
        n::FT, 
        m::FT, 
        S_l::FT) where {FT}
Converts effective saturation to absolute pressure, in unsaturated case.
"""
function van_Genuchten_volume_to_pressure(
        α::FT, 
        n::FT, 
        m::FT, 
        S_l::FT) where {FT}
        p = -((S_l^(-FT(1) / m) - FT(1)) * α^(-n))^(FT(1) / n) 
    return p # units of [m]
end

"""
    van_Genuchten_pressure_to_volume(
        α::FT, 
        n::FT, 
        m::FT, 
        p::FT) where {FT}
Converts absolute pressure to effective saturation, in unsaturated case.
"""
function van_Genuchten_pressure_to_volume(
        α::FT, 
        n::FT, 
        m::FT,
        p::FT) where {FT}
        S_l = ((-α * p)^n + FT(1.0))^(-m)
    return S_l # units of [m3 m-3]
end

"""
    volume_to_pressure(
        vg_α::FT,
        vg_n::FT,
        vg_m::FT,
        ϑ_l::FT,
        ν::FT,
        S_s::FT) where {FT}
Converts augmented liquid fraction to effective saturation, 
and then effective saturation to pressure, for both the 
the unsaturated and saturated regimes. Currently this is using 
vG parameters for loamy type soil (change in future pr).
"""
function volume_to_pressure(
            vg_α::FT,
            vg_n::FT,
            vg_m::FT,
            ϑ_l::FT,
            ν::FT,
            S_s::FT) where {FT}
    S_l = effective_saturation(ν, ϑ_l)
    if S_l <= FT(1.0)
        p = van_Genuchten_volume_to_pressure(vg_α, vg_n, vg_m, S_l)
    else
        p = (ϑ_l - ν) / S_s
    end
    return p # units of (m)
end

"""
    pressure_to_volume(p::FT) where {FT}
Computes the pressure (p)  given the augmented liquid fraction (ϑ_l).
Currently this is using vG parameters for loamy type soil 
(change in future pr).
"""
function pressure_to_volume(
    vg_α::FT,
    vg_n::FT,
    vg_m::FT,
    p::FT,
    ν::FT,
    S_s::FT) where {FT}    
    if p <= FT(0.0)
        S_l = van_Genuchten_pressure_to_volume(vg_α, vg_n, vg_m, p)
    else
        S_l = p * S_s + ν 
    end
        ϑ_l = augmented_liquid_fraction(ν, S_l)   
    return ϑ_l
end

"""
    make_rhs(model::VegetationModel)
A function which creates the rhs! function for the VegetationModel.
The rhs! function must comply with a rhs function of OrdinaryDiffEq.jl.
"""
function make_rhs(model::VegetationModel)
    function rhs!(dY, Y, p, t)
        @unpack a_root,
        a_stem,
        b_root,
        b_stem,
        K_sat_root,
        K_sat_stem, 
        vg_α,
        vg_n,
        vg_m,
        ν,
        S_s = model.parameters

        z_ground, z_stem_top, z_leaf_top = model.domain.compartment_surfaces
        z_stem_midpoint, z_leaf_midpoint = model.domain.compartment_midpoints

        p_stem = volume_to_pressure.(vg_α,vg_n,vg_m,Y.vegetation.ϑ_l[1],ν,S_s)
        p_leaf = volume_to_pressure.(vg_α,vg_n,vg_m,Y.vegetation.ϑ_l[2],ν,S_s)

        # Fluxes are in meters/second
        flux_in_stem = flux_out_roots(model.root_extraction, model, Y, p, t)

        flux_out_stem = flux(
            z_stem_midpoint,
            z_leaf_midpoint, 
            p_stem, 
            p_leaf, 
            a_root, 
            a_stem, 
            b_root, 
            b_stem, 
            K_sat_root, 
            K_sat_stem)

            dY.vegetation.ϑ_l[1] = 1/(z_stem_top-z_ground)*(flux_in_stem - flux_out_stem)            
            dY.vegetation.ϑ_l[2] =
            1/(z_leaf_top-z_stem_top)*(flux_out_stem - transpiration(model.transpiration, t))
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
    PrescribedRootFlux{FT} <: AbstractRootExtraction{FT}
A concrete type used for dispatch when computing the `flux_out_roots`,
in the case where the total flux from all the roots prescribed.
"""
struct PrescribedRootFlux{FT} <: AbstractRootExtraction{FT}
    RF::Function
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
        model::VegetationModel{FT},
        Y::ClimaCore.Fields.FieldVector,
        p::ClimaCore.Fields.FieldVector,
        t::FT,
    )::FT where {FT}
A method which computes the flux between the soil and the stem, via the roots,
in the case of a standalone vegetation model with prescribed soil pressure (in m)
at the root tips.
This assumes that the stem compartment is the first element of `Y.vegetation.ϑ_l`.
"""
function flux_out_roots(
    re::PrescribedSoilPressure{FT},
    model::VegetationModel{FT},
    Y::ClimaCore.Fields.FieldVector,
    p::ClimaCore.Fields.FieldVector,
    t::FT,
)::FT where {FT}
    @unpack vg_α, 
    vg_n, 
    vg_m, 
    ν, 
    S_s, 
    a_root, 
    a_stem, 
    b_root, 
    b_stem, 
    K_sat_root, 
    K_sat_stem = model.parameters
    p_stem = volume_to_pressure.(vg_α,vg_n,vg_m,Y.vegetation.ϑ_l[1],ν,S_s)
    return sum(
        flux.(model.domain.root_depths,
        model.domain.compartment_midpoints[1], 
        re.p_soil(t), 
        p_stem, 
        a_root, 
        a_stem, 
        b_root, 
        b_stem, 
        K_sat_root, 
        K_sat_stem),
    )
end

"""
    flux_out_roots(
        flux_out_roots::PrescribedRootFlux{FT},
        t::FT,
    )::FT where {FT}
A method which prescribes the flux between the soil and the stem, via the roots,
in the case of a standalone vegetation model. This assumes that the 
stem compartment is the first element of `Y.vegetation.ϑ_l`.
"""
function flux_out_roots(
    flux_out_roots::PrescribedRootFlux{FT},
    t::FT,
)::FT where {FT}
    return flux_out_roots.RF(t)
end

"""
    transpiration(
        transpiration::PrescribedTranspiration{FT},
        t::FT,
    )::FT where {FT}
A method which computes the transpiration in meters/sec between the leaf
and the atmosphere,
in the case of a standalone vegetation model with prescribed
transpiration rate.
"""
function transpiration(
    transpiration::PrescribedTranspiration{FT},
    t::FT,
)::FT where {FT}
    return transpiration.T(t)
end

end
