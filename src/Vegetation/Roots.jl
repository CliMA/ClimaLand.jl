# To do: convert all units to SI, move constants out of code and use ClimaParameters,
# convert from a single plant to a bulk plant (units change, required input changes).
# we should also change the name to Vegetation or biophysics as appropriate.
module Roots
#=
    Roots
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
To simulated land surfaces with multiple components (vegetation,
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
using CLIMAParameters: AbstractEarthParameterSet

using ClimaLSM.Domains: AbstractVegetationDomain, RootDomain
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
export RootsModel,
    AbstractVegetationModel,
    flow,
    theta_to_p,
    p_to_theta,
    flow_out_roots,
    RootsParameters,
    PrescribedSoilPressure,
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
    RootsParameters{FT <: AbstractFloat}
A struct for holding parameters of the Root Model. Eventually to be used with ClimaParameters.
$(DocStringExtensions.FIELDS)
"""
struct RootsParameters{FT <: AbstractFloat, PSE <: AbstractEarthParameterSet}
    "controls the shape and steepness of conductance vs. pressure curve, for roots: unitless"
    a_root::FT
    "controls the steepness of the relative conductance vs. pressure curve, for roots: inverse MPa"
    b_root::FT
    "controls the shape and steepness of relative conductance vs. pressure curve, for stems: unitless"
    a_stem::FT
    "controls the steepness of the conductance vs. pressure curve, for stems: inverse MPa"
    b_stem::FT
    "the physical size of the stem, in moles"
    size_reservoir_stem_moles::FT
    "the physical size of the leaves, in moles"
    size_reservoir_leaf_moles::FT
    "water conductance in roots (moles/s/MPa) when pressure is zero, a maximum"
    K_max_root_moles::FT
    "water conductance in stems (moles/s/MPa) when pressure is zero, a maximum"
    K_max_stem_moles::FT
    "Physical Constants and other clima-wide parameters"
    earth_param_set::PSE
end

"""
    RootsModel{FT, PS, D, RE, T, B} <: AbstractVegetationModel{FT}
Defines, and constructs instances of, the RootsModel type, which is used
for simulation flow of water to/from soil, along roots of different depths,
along a stem, to a leaf, and ultimately being lost from the system by
transpiration. 
This model can be used in standalone mode by prescribing the transpiration rate
and soil pressure at the root tips, or with a dynamic soil model using `ClimaLSM`.
$(DocStringExtensions.FIELDS)
"""
struct RootsModel{FT, PS, D, RE, T} <: AbstractVegetationModel{FT}
    "Parameters required by the root model"
    parameters::PS
    "The root model domain, of type `AbstractVegetationDomain`"
    domain::D
    "The root extraction model, of type `AbstractRootExtraction`"
    root_extraction::RE
    "The transpiration model, of type `AbstractTranspiration`"
    transpiration::T
end

function RootsModel{FT}(;
    parameters::RootsParameters{FT, PSE},
    domain::AbstractVegetationDomain{FT},
    root_extraction::AbstractRootExtraction{FT},
    transpiration::AbstractTranspiration{FT},
) where {FT, PSE}
    args = (parameters, domain, root_extraction, transpiration)
    return RootsModel{FT, typeof.(args)...}(args...)
end

"""
    prognostic_vars(model::RootsModel)
A function which returns the names of the prognostic 
variables of the `RootsModel`.
"""
prognostic_vars(model::RootsModel) = (:rwc,)
prognostic_types(model::RootsModel{FT}) where {FT} = (FT,)
"""
    function flow(
        z1::FT,
        z2::FT,
        p1::FT,
        p2::FT,
        a::FT,
        b::FT,
        Kmax::FT,
    ) where {FT}
Computes the flow of water (moles/sec)  given the height and pressures
at two points. Here, `a`, `b, and `Kmax` are parameters
which parameterize the hydraulic conductance of the pathway along which
the flow occurs.
"""
function flow(
    z1::FT,
    z2::FT,
    p1::FT,
    p2::FT,
    a::FT,
    b::FT,
    Kmax::FT,
)::FT where {FT}
    u1, u2, A, B, flow_approx = vc_integral_approx(z1, z2, p1, p2, a, b, Kmax)
    flow = vc_integral(u1, u2, A, B, flow_approx)
    return flow
end

"""
    vc_integral_approx(
        z1::FT,
        z2::FT,
        p1::FT,
        p2::FT,
        a::FT,
        b::FT,
        Kmax::FT,
    ) where {FT}
Approximates the vc integral given the height and pressures
at two points. Here, `a`, `b, and `Kmax` are parameters
which parameterize the hydraulic conductance of the pathway along which
the flow occurs.
"""
function vc_integral_approx(
    z1::FT,
    z2::FT,
    p1::FT,
    p2::FT,
    a::FT,
    b::FT,
    Kmax::FT,
) where {FT}
    rhog_MPa = FT(0.0098)
    u1 = a * exp(b * p1)
    u2 = a * exp(b * p2)
    num1 = log(u1 + FT(1))
    num2 = log(u2 + FT(1))
    c = Kmax * (a + FT(1)) / a
    d = rhog_MPa * (z2 - z1)
    flow_approx = -c / b * (num2 - num1) * (p2 - p1 + d) / (p2 - p1) # this is NaN if p2 = p1
    A = c * d + flow_approx
    B = -c * flow_approx / (b * A)
    return u1, u2, A, B, flow_approx
end

"""
    vc_integral(u1::FT, u2::FT, A::FT, B::FT, flow_approx::FT) where {FT}
Computes the vc integral given the approximate flow.
"""
function vc_integral(u1::FT, u2::FT, A::FT, B::FT, flow_approx::FT) where {FT}
    flow = B * log((u2 * A + flow_approx) / (u1 * A + flow_approx))
    return flow
end

"""
    theta_to_p(theta::FT) where {FT}
Computes the volumetric water content (moles/moles) given pressure (p).
Currently this is using appropriate vG parameters for loamy type soil.
First the head (m) is computed, and then converted to a pressure in MPa.
"""
function theta_to_p(theta::FT) where {FT}
    α = FT(5.0) # inverse meters
    n = FT(2.0)
    m = FT(0.5)
    rhog_MPa = FT(0.0098) #MPa/m
    p = -((theta^(-FT(1) / m) - FT(1)) * α^(-n))^(FT(1) / n) * rhog_MPa
    return p
end


"""
    p_to_theta(p::FT) where {FT}
Computes the pressure (p)  given the volumetric water content (theta).
Currently this is using appropriate vG parameters for loamy type soil.
The pressure (MPa)  must be converted to meters (head) for use in the
van Genuchten formula.
"""
function p_to_theta(p::FT) where {FT}
    α = FT(5.0) # inverse meters
    n = FT(2.0)
    m = FT(0.5)
    rhog_MPa = FT(0.0098) #MPa/m
    theta = ((-α * (p / rhog_MPa))^n + FT(1.0))^(-m)
    return theta
end


"""
    make_rhs(model::RootsModel)
A function which creates the rhs! function for the RootsModel.
The rhs! function must comply with a rhs function of OrdinaryDiffEq.jl.
"""
function make_rhs(model::RootsModel)
    function rhs!(dY, Y, p, t)
        @unpack a_root,
        b_root,
        K_max_root_moles,
        size_reservoir_stem_moles,
        a_stem,
        b_stem,
        K_max_stem_moles,
        size_reservoir_leaf_moles = model.parameters

        z_stem, z_leaf = model.domain.compartment_heights

        p_stem = theta_to_p(Y.vegetation.rwc[1] / size_reservoir_stem_moles)
        p_leaf = theta_to_p(Y.vegetation.rwc[2] / size_reservoir_leaf_moles)

        # Flows are in moles/second
        flow_in_stem = flow_out_roots(model.root_extraction, model, Y, p, t)

        flow_out_stem = flow(
            z_stem,
            z_leaf,
            p_stem,
            p_leaf,
            a_stem,
            b_stem,
            K_max_stem_moles,
        )

        dY.vegetation.rwc[1] = flow_in_stem - flow_out_stem
        dY.vegetation.rwc[2] =
            flow_out_stem - transpiration(model.transpiration, t)
    end
    return rhs!
end

"""
    PrescribedSoilPressure{FT} <: AbstractRootExtraction{FT}
A concrete type used for dispatch when computing the `flow_out_roots`,
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
    flow_out_roots(
        re::PrescribedSoilPressure{FT},
        model::RootsModel{FT},
        Y::ClimaCore.Fields.FieldVector,
        p::ClimaCore.Fields.FieldVector,
        t::FT,
    )::FT where {FT}
A method which computes the flow between the soil and the stem, via the roots,
in the case of a standalone root model with prescribed soil pressure (in MPa)
at the root tips.
This assumes that the stem compartment is the first element of `Y.roots.rwc`.
"""
function flow_out_roots(
    re::PrescribedSoilPressure{FT},
    model::RootsModel{FT},
    Y::ClimaCore.Fields.FieldVector,
    p::ClimaCore.Fields.FieldVector,
    t::FT,
)::FT where {FT}
    @unpack a_root, b_root, K_max_root_moles, size_reservoir_stem_moles =
        model.parameters
    p_stem = theta_to_p(Y.vegetation.rwc[1] / size_reservoir_stem_moles)
    return sum(
        flow.(
            model.domain.root_depths,
            model.domain.compartment_heights[1],
            re.p_soil(t),
            p_stem,
            a_root,
            b_root,
            K_max_root_moles,
        ),
    )
end

"""
    transpiration(
        transpiration::PrescribedTranspiration{FT},
        t::FT,
    )::FT where {FT}
A method which computes the transpiration in moles/sec between the leaf
and the atmosphere,
in the case of a standalone root model with prescribed transpiration rate.
"""
function transpiration(
    transpiration::PrescribedTranspiration{FT},
    t::FT,
)::FT where {FT}
    return transpiration.T(t)
end

end
