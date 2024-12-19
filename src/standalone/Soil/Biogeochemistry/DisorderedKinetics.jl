module DisorderedKinetics
using Distributions
using ClimaLand
using DocStringExtensions
using ClimaCore
import ...Parameters as LP
import ClimaCore: Fields, Operators, Geometry, Spaces

import ClimaLand.Domains: AbstractDomain
import ClimaLand:
    AbstractExpModel,
    make_update_aux,
    make_compute_exp_tendency,
    make_update_boundary_fluxes,
    prognostic_vars,
    auxiliary_vars,
    name,
    prognostic_types,
    auxiliary_types,
    prognostic_domain_names,
    auxiliary_domain_names,
    TopBoundary,
    BottomBoundary,
    AbstractBC,
    boundary_flux!,
    AbstractSource,
    source!

"""
DisorderedKineticsModelParameters{FT <: AbstractFloat, PSE}

A struct for storing parameters of the `DisorderedKineticsModel`.

All of these parameters are currently treated as global constants.
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct DisorderedKineticsModelParameters{FT <: AbstractFloat}
    "the mean of the log-normal distribution of the soil C decay rate []"
    mu::FT
    "the standard deviation of the log-normal distribution of the soil C decay rate []"
    sigma::FT
end


"""
AbstractSoilCarbonModel{FT} <: ClimaLand.AbstractExpModel{FT}

An abstract model type for soil biogeochemistry models.
"""
abstract type AbstractSoilCarbonModel{FT} <:
        ClimaLand.AbstractExpModel{FT} end

"""
DisorderedKineticsModel

A model for simulating soil organic carbon dynamics using a disordered kinetics model.

$(DocStringExtensions.FIELDS)
"""
struct DisorderedKineticsModel{FT, D, BC, S} <:
    AbstractSoilCarbonModel{FT}
"the parameter set"
parameters::PS
"the number of soil carbon pools for discretization"
npools::Int64
"the soil domain, using ClimaCore.Domains"
domain::D
"the boundary conditions, of type NamedTuple"
boundary_conditions::BC
"A tuple of sources, each of type AbstractSource"
sources::S
end


"""
DisorderedKineticsModel{FT}(;
        parameters::DisorderedKineticsModelParameters{FT},
        npools::Int64,
        domain::ClimaLand.AbstractDomain,
        boundary_conditions::NamedTuple,
        sources::Tuple,
    ) where {FT, BC}

A constructor for `DisorderedKineticsModel`.
"""
function DisorderedKineticsModel{FT}(;
    parameters::DisorderedKineticsModelParameters{FT},
    npools::Int64,
    domain::ClimaLand.AbstractDomain,
    boundary_conditions::BC,
    sources::Tuple,
) where {FT, BC}
    args = (parameters, npools, domain, boundary_conditions, sources)
    DisorderedKineticsModel{FT, typeof.(args)...}(args...)
end

ClimaLand.name(model::DisorderedKineticsModel) = :soilC
ClimaLand.prognostic_vars(::DisorderedKineticsModel) = (:C,)
ClimaLand.prognostic_types(::DisorderedKineticsModel{FT}) where {FT} = (NTuple{model.npools, FT},)
ClimaLand.prognostic_domain_names(::DisorderedKineticsModel) = (:surface,)

ClimaLand.auxiliary_vars(model::DisorderedKineticsModel) = (:ks,)
ClimaLand.auxiliary_types(model::DisorderedKineticsModel{FT}) where {FT} = (NTuple{model.npools, FT},)
ClimaLand.auxiliary_domain_names(model::DisorderedKineticsModel) = (:surface)

"""
    make_compute_exp_tendency(model::DisorderedKineticsModel)

An extension of the function `make_compute_exp_tendency`, for the disordered kinetics model.
This function creates and returns a function which computes the entire
right hand side of the ODE for `C`, and updates `dY.soilC.C` in place
with that value. These quantities will be stepped explicitly.

This has been written so as to work with Differential Equations.jl.
"""
function ClimaLand.make_compute_exp_tendency(model::DisorderedKineticsModel)
    function compute_exp_tendency!(dY, Y, p, t)
        
        @. dY.soilC.C = -p.ks * Y.soilC.C

        # Source terms are added in here
        for src in model.sources
            source!(dY, src, Y, p, model.parameters)
        end
    end
    return compute_exp_tendency!
end

"""
    AbstractCarbonSource{FT} <: ClimaLand.AbstractSource{FT}

An abstract type for soil CO2 sources. There are two sources:
roots and microbes, in struct RootProduction and MicrobeProduction.
"""
abstract type AbstractSoilCarbonSource{FT} <: ClimaLand.AbstractSource{FT} end

"""
    ClimaLand.source!(dY::ClimaCore.Fields.FieldVector,
                          src::AbstractSoilCarbonSource,
                          Y::ClimaCore.Fields.FieldVector,
                          p::NamedTuple,
                          params)

A method which extends the ClimaLand source! function for the
case of microbe production of CO2 in soil.
"""
function ClimaLand.source!(
    dY::ClimaCore.Fields.FieldVector,
    src::AbstractSoilCarbonSource,
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    params,
)
    @. dY.soilco2.C += src * pdf(LogNormal(params.mu, params.sigma), p.ks)
end


"""
    make_update_aux(model::DisorderedKineticsModel)

An extension of the function `make_update_aux`, for the disordered kinetics model.
This function creates and returns a function which updates the auxiliary
variables `p.soilC.variable` in place.
This has been written so as to work with Differential Equations.jl.

This function is empty because the auxiliary variable `ks` is not updated in time.
"""
function ClimaLand.make_update_aux(model::DisorderedKineticsModel)
    function update_aux!(p, Y, t) end
    return update_aux!
end

Base.broadcastable(ps::DisorderedKineticsModelParameters) = tuple(ps)

end