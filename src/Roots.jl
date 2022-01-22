module Roots
using ClimaLSM
using ClimaCore
using UnPack
import ClimaCore: Fields
using ClimaLSM.ComponentExchanges: AbstractComponentExchange
using ClimaLSM.Domains: AbstractPlantDomain, RootDomain

import ClimaLSM:
    AbstractModel,
    initialize_prognostic,
    make_update_aux,
    make_rhs,
    make_ode_function,
    prognostic_vars,
    auxiliary_vars,
    initialize,
    initialize_auxiliary
export RootsModel,
    compute_flow,
    theta_to_p,
    p_to_theta,
    RootsPrescribedExchange,
    RootsParameters

"""
    RootsParameters{FT <: AbstractFloat}


A struct for holding parameters of the Root Model. Eventually to be used with ClimaParameters.
"""
struct RootsParameters{FT <: AbstractFloat}
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
end


"""
    RootsModel

Defines, and constructs instances of, the RootsModel type, which is used
for simulation flow of water to/from soil, along roots of different depths,
along a stem, to a leaf, and ultimately being lost from the system by
transpiration. 

This model can be used in standalone mode by prescribing the transpiration rate
and soil pressure at the root tips (boundary_exchanges of type `RootsPrescribedExchange`),
or with a dynamic soil model (boundary exchanges of type TBD).

The RootModel domain must be of type `AbstractPlantDomain`.
"""
struct RootsModel{FT, PS, D, B} <: AbstractModel{FT}
    param_set::PS
    domain::D
    boundary_exchanges::B
    model_name::Symbol
end

function RootsModel{FT}(;
    param_set,
    domain::AbstractPlantDomain{FT},
    boundary_exchanges::AbstractComponentExchange{FT},
) where {FT}
    return RootsModel{
        FT,
        typeof(param_set),
        typeof(domain),
        typeof(boundary_exchanges),
    }(
        param_set,
        domain,
        boundary_exchanges,
        :roots,
    )
end




"""
    RootsPrescribedExchange{FT} <: AbstractComponentExchange{FT}

The component exchange type to be used for the Root model in standalone mode.

The user provides a function of time returning the transpiration rate (moles/sec),
 and a function of time returning a 
vector of soil pressures (MPa). The different elements of
this vector hold the soil pressure at the root tips, with depths indicated by
the different elements of the `root_depths` vector of `RootDomain`. 
"""
struct RootsPrescribedExchange{FT} <: AbstractComponentExchange{FT}
    "Time dependent transpiration, given in moles/sec"
    T::Function
    "Time dependent soil pressure at root tips, given in MPa"
    p_soil::Function
end


prognostic_vars(model::RootsModel) = (:rwc,)

"""
    function compute_flow(
        z_do::FT,
        z_up::FT,
        p_do::FT,
        p_up::FT,
        a::FT,
        b::FT,
        Kmax::FT,
    ) where {FT}

Computes the flow of water (moles/sec)  given the height and pressures
at two points, `do` and `up`. Here, `a`, `b, and `Kmax` are parameters
which parameterize the hydraulic conductance of the pathway along which
the flow occurs.
"""
function compute_flow(
    z_do::FT,
    z_up::FT,
    p_do::FT,
    p_up::FT,
    a::FT,
    b::FT,
    Kmax::FT,
) where {FT}
    u_do, u_up, A, B, flow_approx =
        vc_integral_approx(z_do, z_up, p_do, p_up, a, b, Kmax)
    flow::FT = vc_integral(u_do, u_up, A, B, flow_approx)
    return flow
end

"""
    vc_integral(
        u_do::FT,
        u_up::FT,
        A::FT,
        B::FT,
        flow_approx::FT,
    ) where {FT}

Approximates the vc integral given the height and pressures
at two points, `do` and `up`. Here, `a`, `b, and `Kmax` are parameters
which parameterize the hydraulic conductance of the pathway along which
the flow occurs.
"""
function vc_integral_approx(
    z_do::FT,
    z_up::FT,
    p_do::FT,
    p_up::FT,
    a::FT,
    b::FT,
    Kmax::FT,
) where {FT}
    rhog_MPa = FT(0.0098)
    u_do = a * exp(b * p_do)
    u_up = a * exp(b * p_up)
    num_do = log(u_do + FT(10))
    num_up = log(u_up + FT(10))
    c = Kmax * (a + FT(10)) / a
    d = rhog_MPa * (z_up - z_do)
    flow_approx = -c / b * (num_up - num_do) * (p_up - p_do + d) / (p_up - p_do) # this is NaN if p_up = p_do
    A = c * d + flow_approx
    B = -c * flow_approx / (b * A)
    return u_do, u_up, A, B, flow_approx
end

"""
    vc_integral(
        u_do::FT,
        u_up::FT,
        A::FT,
        B::FT,
        flow_approx::FT,
    ) where {FT}

Computes the vc integral given the approximate flow.
"""
function vc_integral(
    u_do::FT,
    u_up::FT,
    A::FT,
    B::FT,
    flow_approx::FT,
) where {FT}
    flow = B * log((u_up * A + flow_approx) / (u_do * A + flow_approx))
    return flow
end

"""
    theta_to_p(theta::FT) where {FT}

Computes the volumetric water content (theta) given pressure (p)
"""
function theta_to_p(theta::FT) where {FT}
    p = (theta - FT(1)) * FT(5)
    return p
end


"""
    p_to_theta(p::FT) where {FT}

Computes the pressure (p)  given the volumetric water content (theta).
"""
function p_to_theta(p::FT) where {FT}
    theta = p / FT(5) + FT(1)
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
        size_reservoir_leaf_moles = model.param_set

        z_stem, z_leaf = model.domain.compartment_heights

        p_stem = theta_to_p(Y.roots.rwc[1] / size_reservoir_stem_moles)
        p_leaf = theta_to_p(Y.roots.rwc[2] / size_reservoir_leaf_moles)

        # Flows are in moles/second
        flow_out_roots =
            compute_flow_out_roots(model.boundary_exchanges, model, Y, p, t)

        flow_in_stem = sum(flow_out_roots)
        flow_out_stem = compute_flow(
            z_stem,
            z_leaf,
            p_stem,
            p_leaf,
            a_stem,
            b_stem,
            K_max_stem_moles,
        )

        dY.roots.rwc[1] = flow_in_stem - flow_out_stem
        dY.roots.rwc[2] =
            flow_out_stem -
            compute_transpiration(model.boundary_exchanges, Y, p, t)
    end
    return rhs!
end

"""
    function compute_flow_out_roots(
        boundary_exchanges::RootsPrescribedExchange{FT},
        model::RootsModel{FT},
        Y::ClimaCore.Fields.FieldVector,
        p::ClimaCore.Fields.FieldVector,
        t::FT,
    )::Vector{FT} where {FT}

A method which computes the flow between the soil and the stem, via the roots,
in the case of a standalone root model with prescribed soil pressure (in MPa)
at the root tips.

This assumes that the stem compartment is the first element of `Y.roots.rwc`.
"""
function compute_flow_out_roots(
    boundary_exchanges::RootsPrescribedExchange{FT},
    model::RootsModel{FT},
    Y::ClimaCore.Fields.FieldVector,
    p::ClimaCore.Fields.FieldVector,
    t::FT,
)::Vector{FT} where {FT}
    @unpack a_root, b_root, K_max_root_moles, size_reservoir_stem_moles =
        model.param_set
    p_stem = theta_to_p(Y.roots.rwc[1] / size_reservoir_stem_moles)

    flow =
        compute_flow.(
            model.domain.root_depths,
            model.domain.compartment_heights[1],
            boundary_exchanges.p_soil(t),
            p_stem,
            a_root,
            b_root,
            K_max_root_moles,
        )
    return flow
end


"""
    compute_transpiration(boundary_exchanges::RootsPrescribedExchange{FT},
        Y::ClimaCore.Fields.FieldVector,
        p::ClimaCore.Fields.FieldVector,
        t::FT)::FT where {FT}

A method which computes the transpiration in moles/sec between the leaf
and the atmosphere,
in the case of a standalone root model with prescribed transpiration rate.
"""
function compute_transpiration(
    boundary_exchanges::RootsPrescribedExchange{FT},
    Y::ClimaCore.Fields.FieldVector,
    p::ClimaCore.Fields.FieldVector,
    t::FT,
)::FT where {FT}
    return boundary_exchanges.T(t)
end






end
