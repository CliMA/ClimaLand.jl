module Soil
using ClimaLSM
using ClimaCore
import ClimaCore: Fields, Operators, Geometry
using ClimaLSM.Domains
import ClimaLSM:
    AbstractModel,
    initialize,
    initialize_auxiliary,
    initialize_prognostic,
    make_update_aux,
    make_rhs,
    make_ode_function,
    prognostic_vars,
    auxiliary_vars
export RichardsModel

"""
    RichardsModel

A standin for a model used to simulate water flow in soil via Richards equation.
"""
struct RichardsModel{FT, PS, D} <: AbstractModel{FT}
    param_set::PS
    domain::D
    model_name::Symbol
end

function RichardsModel{FT}(; param_set, domain::AbstractDomain{FT}) where {FT}
    return RichardsModel{FT, typeof(param_set), typeof(domain)}(
        param_set,
        domain,
        :soil,
    )
end

prognostic_vars(model::RichardsModel) = (:ϑ_l,)

auxiliary_vars(model::RichardsModel) = (:z,)

end
