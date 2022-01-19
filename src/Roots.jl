module Roots
using ClimaLSM
using ClimaCore
import ClimaCore: Fields
using ClimaCore.Domains
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
export RootsModel

"""
    RootsModel

A standin for a model used to simulate water content in roots.
"""
struct RootsModel{FT, ps, d} <: AbstractModel{FT}
    param_set::ps
    domain::d
    model_name::Symbol
end

function RootsModel{FT}(;
    param_set = param_set,
    domain = domain <: AbstractDomain{FT},
) where {FT}
    return RootsModel{FT, typeof(param_set), typeof(domain)}(
        param_set,
        domain,
        :roots,
    )
end

prognostic_vars(model::RootsModel) = (:rwc,)

end
