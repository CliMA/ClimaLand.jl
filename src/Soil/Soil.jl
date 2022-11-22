module Soil
#=
    Soil

This module contains everything needed to run a soil model
in standalone mode.

The soil model is assumed to have a set of prognostic `Y` and
auxiliary `p` variables, which describe the state of the
soil system. The system is evolved in time by solving 
equations of the form

```

\frac{\partial Y}{\partial t} = D(Y, t, p(Y, t; \ldots);\ldots),

```

i.e. partial (or ordinary) differential equations depending 
on state `Y`,
auxiliary functions of the state `p`, and other parameters 
represented by the ellipses. The operator `D` indicates a generic
nonlinear differential operator.  Not every model
requires auxilary variables, but these are useful for storing
quantities that are needed multiple times per right hand side
evaluation. An example would be the temperature, which is computed
from the prognostic variables of water content and internal
energy.

Currently, both the Richardson Richards Equation (RRE; hydrology alone)
and an integrated soil energy and hydrology model are supported.

Addition of additional versions of soil
models requires defining a model type (of super type 
`AbstractSoilModel`), and extending the methods
imported by Models.jl, as needed, for computing the
right hand side functions of the ordinary differential equations
and the functions which update auxiliary variables whenever
the right hand side is evaluated.

This code base assumes that DifferentialEquations.jl
will be used for evolving the system in time,
and that the array-like objected being stepped forward
is a `ClimaCore.Fields.FieldVector`. 

The `FieldVector` type is used in order to make use of 
`ClimaCore` functionality when solving PDEs (`ClimaCore` handles
all of the spatial discretization and operator functionality)
 and for ease of handling multi-column models.

To simulate land surfaces with multiple components (vegetation,
soil, rivers, etc), the ClimaLSM.jl package should be used.
That package will use the methods of this function for advancing
the system forward in time, extending methods as needed to account
for interactions between components.
=#

using ClimaLSM
using UnPack
using DocStringExtensions
using ClimaCore
import ..Parameters as LSMP
import ClimaCore: Fields, Operators, Geometry, Spaces

import ClimaLSM.Domains: Column, HybridBox, SphericalShell
import ClimaLSM:
    AbstractModel,
    make_update_aux,
    make_rhs,
    prognostic_vars,
    auxiliary_vars,
    name,
    prognostic_types,
    auxiliary_types,
    AbstractSource,
    AbstractBC,
    source!,
    boundary_flux
export RichardsModel,
    RichardsParameters,
    EnergyHydrology,
    EnergyHydrologyParameters,
    FluxBC,
    StateBC,
    FreeDrainage,
    AbstractSoilModel,
    AbstractSoilSource,
    FluxBC,
    FreeDrainage,
    PhaseChange

"""
    AbstractSoilBC <: ClimaLSM. AbstractBC

An abstract type for soil-specific types of boundary conditions, like free drainage.
"""
abstract type AbstractSoilBC <: ClimaLSM.AbstractBC end

"""
    AbstractSoilSource{FT} <:  ClimaLSM.AbstractSource{FT}

An abstract type for types of source terms for the soil equations.

In standalone mode, the only supported source type is freezing and 
thawing. ClimaLSM.jl creates additional sources to include as
necessary e.g. root extraction (not available in stand alone mode).
"""
abstract type AbstractSoilSource{FT} <: ClimaLSM.AbstractSource{FT} end

"""
    AbstractSoilModel{FT} <: AbstractModel{FT} 

The abstract type for all soil models.

Currently, we only have plans to support a RichardsModel, simulating
the flow of liquid water through soil via the Richardson-Richards equation,
and a fully integrated soil heat and water model, with phase change.
"""
abstract type AbstractSoilModel{FT} <: ClimaLSM.AbstractModel{FT} end

ClimaLSM.name(::AbstractSoilModel) = :soil
ClimaLSM.domain(::AbstractSoilModel) = :subsurface

"""
   horizontal_components!(dY::ClimaCore.Fields.FieldVector,
                          domain::Column, _...)
Updates dY in place by adding in the tendency terms resulting from
horizontal derivative operators.

In the case of a column domain, there are no horizontal
contributions to the right hand side.
"""
function horizontal_components!(
    dY::ClimaCore.Fields.FieldVector,
    domain::Column,
    _...,
) end

"""
   dss!(dY::ClimaCore.Fields.FieldVector,domain::Column)

Computes the appropriate weighted direct stiffness summation based on
the domain type, updates `dY` in place.

For column domains, no dss is needed.
"""
function dss!(dY::ClimaCore.Fields.FieldVector, domain::Column) end

"""
   dss!(dY::ClimaCore.Fields.FieldVector,domain::Union{HybridBox, SphericalShell})

Computes the appropriate weighted direct stiffness summation based on
the domain type, updates `dY` in place.
"""
function dss!(
    dY::ClimaCore.Fields.FieldVector,
    domain::Union{HybridBox, SphericalShell},
)
    for key in propertynames(dY.soil)
        Spaces.weighted_dss!(getproperty(dY.soil, key))
    end
end

"""
   StateBC <: AbstractSoilBC

A simple concrete type of boundary condition, which enforces a
state boundary condition f(p,t) at either the top or bottom of the domain.
"""
struct StateBC <: AbstractSoilBC
    bc::Function
end

"""
   FluxBC <: AbstractSoilBC

A simple concrete type of boundary condition, which enforces a
normal flux boundary condition f(p,t) at either the top or bottom of the domain.
"""
struct FluxBC <: AbstractSoilBC
    bc::Function
end

"""
    FreeDrainage{FT} <: AbstractSoilBC
A concrete type of soil boundary condition, for use at 
the BottomBoundary only, where the flux is set to be
`F = -K∇h = -K`.
"""
struct FreeDrainage{FT} <: AbstractSoilBC end

"""
    ClimaLSM.boundary_flux(bc::FluxBC, _, Δz, _...)::ClimaCore.Fields.Field

A method of boundary fluxes which returns the desired flux.

We add a field of zeros in order to convert the bc (float) into
a field.
"""
function ClimaLSM.boundary_flux(
    bc::FluxBC,
    boundary::ClimaLSM.AbstractBoundary,
    Δz::ClimaCore.Fields.Field,
    p::ClimaCore.Fields.FieldVector,
    t::FT,
    _...,
)::ClimaCore.Fields.Field where {FT}
    return bc.bc(p, t) .+ ClimaCore.Fields.zeros(axes(Δz))
end


"""
    ClimaLSM.boundary_flux(rre_bc::StateBC, ::ClimaLSM.TopBoundary, Δz, p, t, params)::ClimaCore.Fields.Field

A method of boundary fluxes which converts a state boundary condition on θ_l at the top of the
domain into a flux of liquid water.
"""
function ClimaLSM.boundary_flux(
    rre_bc::StateBC,
    ::ClimaLSM.TopBoundary,
    Δz,
    p,
    t,
    params,
)::ClimaCore.Fields.Field
    # Approximate K_bc ≈ K_c, ψ_bc ≈ ψ_c (center closest to the boundary)
    p_len = Spaces.nlevels(axes(p.soil.K))
    K_c = Fields.level(p.soil.K, p_len)
    ψ_c = Fields.level(p.soil.ψ, p_len)

    # Calculate pressure head using boundary condition
    @unpack vg_α, vg_n, vg_m, θ_r, ν, S_s = params
    θ_bc = rre_bc.bc(p, t)
    ψ_bc = @. pressure_head(vg_α, vg_n, vg_m, θ_r, θ_bc, ν, S_s)

    # Pass in (ψ_bc .+ Δz) as x_2 to account for contribution of gravity in RRE
    return ClimaLSM.diffusive_flux(K_c, ψ_bc .+ Δz, ψ_c, Δz)
end

"""
    ClimaLSM.boundary_flux(rre_bc::StateBC, ::ClimaLSM.BottomBoundary, Δz, p, t, params)::ClimaCore.Fields.Field

A method of boundary fluxes which converts a state boundary condition on θ_l at the bottom of the
domain into a flux of liquid water.
"""
function ClimaLSM.boundary_flux(
    rre_bc::StateBC,
    ::ClimaLSM.BottomBoundary,
    Δz,
    p,
    t,
    params,
)::ClimaCore.Fields.Field
    # Approximate K_bc ≈ K_c, ψ_bc ≈ ψ_c (center closest to the boundary)
    K_c = Fields.level(p.soil.K, 1)
    ψ_c = Fields.level(p.soil.ψ, 1)

    # Calculate pressure head using boundary condition
    @unpack vg_α, vg_n, vg_m, θ_r, ν, S_s = params
    θ_bc = rre_bc.bc(p, t)
    ψ_bc = @. pressure_head(vg_α, vg_n, vg_m, θ_r, θ_bc, ν, S_s)

    # At the bottom boundary, ψ_c is at larger z than ψ_bc
    #  so we swap their order in the derivative calc
    # Pass in (ψ_c .+ Δz) as x_2 to account for contribution of gravity in RRE
    return ClimaLSM.diffusive_flux(K_c, ψ_c .+ Δz, ψ_bc, Δz)
end


"""
    ClimaLSM.boundary_flux(heat_bc::StateBC, ::ClimaLSM.TopBoundary, Δz, p, t)::ClimaCore.Fields.Field

A method of boundary fluxes which converts a state boundary condition on temperature at the top of the
domain into a flux of energy.
"""
function ClimaLSM.boundary_flux(
    heat_bc::StateBC,
    ::ClimaLSM.TopBoundary,
    Δz,
    p,
    t,
)::ClimaCore.Fields.Field
    # Approximate κ_bc ≈ κ_c (center closest to the boundary)
    p_len = Spaces.nlevels(axes(p.soil.T))
    T_c = Fields.level(p.soil.T, p_len)
    κ_c = Fields.level(p.soil.κ, p_len)

    T_bc = heat_bc.bc(p, t)
    return ClimaLSM.diffusive_flux(κ_c, T_bc, T_c, Δz)
end

"""
    ClimaLSM.boundary_flux(heat_bc::StateBC, ::ClimaLSM.BottomBoundary, Δz, p, t)::ClimaCore.Fields.Field

A method of boundary fluxes which converts a state boundary condition on temperature at the bottom of the
domain into a flux of energy.
"""
function ClimaLSM.boundary_flux(
    heat_bc::StateBC,
    ::ClimaLSM.BottomBoundary,
    Δz,
    p,
    t,
)::ClimaCore.Fields.Field
    # Approximate κ_bc ≈ κ_c (center closest to the boundary)
    T_c = Fields.level(p.soil.T, 1)
    κ_c = Fields.level(p.soil.κ, 1)

    T_bc = heat_bc.bc(p, t)
    return ClimaLSM.diffusive_flux(κ_c, T_c, T_bc, Δz)
end


"""
    ClimaLSM.boundary_flux(bc::FreeDrainage{FT}, ::ClimaLSM.BottomBoundary, Δz, p, t)::ClimaCore.Fields.Field

A method of boundary fluxes which enforces free drainage at the bottom
of the domain.
"""
function ClimaLSM.boundary_flux(
    bc::FreeDrainage{FT},
    boundary::ClimaLSM.BottomBoundary,
    Δz,
    p,
    t,
    params,
)::ClimaCore.Fields.Field where {FT}
    K_c = Fields.level(p.soil.K, 1)
    return -1 .* K_c
end

"""
     heaviside(x::FT)::FT where {FT}

Computes the heaviside function.
"""
function heaviside(x::FT)::FT where {FT}
    if x > eps(FT)
        return FT(1.0)
    else
        return FT(0.0)
    end
end
include("./rre.jl")
include("./energy_hydrology.jl")
include("./soil_hydrology_parameterizations.jl")
include("./soil_heat_parameterizations.jl")
include("Biogeochemistry/Biogeochemistry.jl")
using .Biogeochemistry
end
