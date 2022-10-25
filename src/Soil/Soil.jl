# To do: add in SoilHeatWaterModel, different boundary conditions, freeze thaw.
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
    auxiliary_types
export RichardsModel,
    RichardsParameters,
    EnergyHydrology,
    EnergyHydrologyParameters,
    boundary_flux,
    AbstractBC,
    FluxBC,
    StateBC,
    RootExtraction,
    AbstractSoilModel,
    AbstractSoilSource,
    source!,
    AbstractBoundary,
    TopBoundary,
    BottomBoundary,
    diffusive_flux,
    get_Δz,
    PhaseChange

"""
    AbstractBC{FT <: AbstractFloat}

An abstract type for types of boundary conditions, which will include
prescribed functions of space and time as Dirichlet conditions or
Neumann conditions, in addition to other 
convenient soil-specific conditions, like free drainage.
"""
abstract type AbstractBC{FT <: AbstractFloat} end

"""
    AbstractSoilSource{FT <: AbstractFloat}

An abstract type for types of source terms for the soil equations.

In standalone mode, the only supported source type is freezing and 
thawing. ClimaLSM.jl creates additional sources to include as
necessary e.g. root extraction (not available in stand alone mode).
"""
abstract type AbstractSoilSource{FT <: AbstractFloat} end

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
    AbstractBoundary{}

An abstract type to indicate which boundary we are doing calculations for.
Currently, we support the top boundary (TopBoundary)
and bottom boundary (BottomBoundary).
"""
abstract type AbstractBoundary{} end

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
# Todo: consider if we should introduce an AbstractHybridDomain so as 
# to avoid listing with Union. Possibly not worth it at this point, but if
# other domains
# are added or if we end up doing this often, it might be nice. But it may
# also be nice to be explicit.
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
   StateBC{FT} <: AbstractBC{FT}

A simple concrete type of boundary condition, which enforces a
constant normal state at either the top or bottom of the domain.
"""
struct StateBC{FT} <: AbstractBC{FT}
    bc::FT
end

"""
   FluxBC{FT} <: AbstractBC{FT}

A simple concrete type of boundary condition, which enforces a
constant normal flux at either the top or bottom of the domain.
"""
struct FluxBC{FT} <: AbstractBC{FT}
    bc::FT
end

"""
    TopBoundary{} <: AbstractBoundary{}

A simple object which should be passed into a function to
indicate that we are considering the top boundary of the soil.
"""
struct TopBoundary <: AbstractBoundary end

"""
    BottomBoundary{} <: AbstractBoundary{}

A simple object which should be passed into a function to
indicate that we are considering the bottom boundary of the soil.
"""
struct BottomBoundary <: AbstractBoundary end

"""
    get_Δz(z)

A function to return a tuple containing the distance between the top boundary
and its closest center, and the bottom boundary and its closest center, 
both as Fields.
"""
function get_Δz(z)
    # Extract the differences between levels of the face space
    fs = ClimaLSM.Domains.obtain_face_space(axes(z))
    z_face = ClimaCore.Fields.coordinate_field(fs).z
    Δz = ClimaCore.Fields.Δz_field(z_face)

    Δz_top = Fields.level(
        Δz,
        ClimaCore.Utilities.PlusHalf(ClimaCore.Spaces.nlevels(fs) - 1),
    )
    Δz_bottom = Fields.level(Δz, ClimaCore.Utilities.PlusHalf(0))
    return Δz_top ./ 2, Δz_bottom ./ 2
end

"""
    diffusive_flux(K, x_2, x_1, Δz)

Calculates the diffusive flux of a quantity x (water content, temp, etc).
Here, x_2 = x(z + Δz) and x_1 = x(z), so x_2 is at a larger z by convention.
"""
function diffusive_flux(K, x_2, x_1, Δz)
    return @. -K * (x_2 - x_1) / Δz
end

"""
    boundary_flux(bc::AbstractBC, bound_type::AbstractBoundary, Δz::FT, _...)

A function which returns the correct boundary flux (of type FluxBC) given
    any boundary condition (BC) of type `FluxBC` or `StateBC`.
A StateBC must be converted into a flux value before being returned.
"""
function boundary_flux end

# A flux BC requires no conversion calculations
function boundary_flux(bc::FluxBC, _, Δz, _...)
    return bc.bc .+ ClimaCore.Fields.zeros(axes(Δz))
end

# Convert water state to flux at top boundary
function boundary_flux(rre_bc::StateBC, ::TopBoundary, Δz, p, params)
    # Approximate K_bc ≈ K_c, ψ_bc ≈ ψ_c (center closest to the boundary)
    p_len = Spaces.nlevels(axes(p.soil.K))
    K_c = Fields.level(p.soil.K, p_len)
    ψ_c = Fields.level(p.soil.ψ, p_len)

    # Calculate pressure head using boundary condition
    @unpack vg_α, vg_n, vg_m, θ_r, ν, S_s = params
    ψ_bc = @. pressure_head(vg_α, vg_n, vg_m, θ_r, rre_bc.bc, ν, S_s)

    # Pass in (ψ_bc .+ Δz) as x_2 to account for contribution of gravity in RRE
    return diffusive_flux(K_c, ψ_bc .+ Δz, ψ_c, Δz)
end

# Convert water state to flux at bottom boundary
function boundary_flux(rre_bc::StateBC, ::BottomBoundary, Δz, p, params)
    # Approximate K_bc ≈ K_c, ψ_bc ≈ ψ_c (center closest to the boundary)
    K_c = Fields.level(p.soil.K, 1)
    ψ_c = Fields.level(p.soil.ψ, 1)

    # Calculate pressure head using boundary condition
    @unpack vg_α, vg_n, vg_m, θ_r, ν, S_s = params
    ψ_bc = @. pressure_head(vg_α, vg_n, vg_m, θ_r, rre_bc.bc, ν, S_s)

    # At the bottom boundary, ψ_c is at larger z than ψ_bc
    #  so we swap their order in the derivative calc
    # Pass in (ψ_c .+ Δz) as x_2 to account for contribution of gravity in RRE
    return diffusive_flux(K_c, ψ_c .+ Δz, ψ_bc, Δz)
end

# Convert heat state to flux at top boundary
function boundary_flux(heat_bc::StateBC, ::TopBoundary, Δz, p)
    # Approximate κ_bc ≈ κ_c (center closest to the boundary)
    p_len = Spaces.nlevels(axes(p.soil.T))
    T_c = Fields.level(p.soil.T, p_len)
    κ_c = Fields.level(p.soil.κ, p_len)

    return diffusive_flux(κ_c, heat_bc.bc, T_c, Δz)
end

# Convert heat state to flux at bottom boundary
function boundary_flux(heat_bc::StateBC, ::BottomBoundary, Δz, p)
    # Approximate κ_bc ≈ κ_c (center closest to the boundary)
    T_c = Fields.level(p.soil.T, 1)
    κ_c = Fields.level(p.soil.κ, 1)

    # At the bottom boundary, T_c is at larger z than T_bc
    #  so we swap their order in the derivative calc
    return diffusive_flux(κ_c, T_c, heat_bc.bc, Δz)
end

"""
     source!(dY::ClimaCore.Fields.FieldVector,
             src::AbstractSoilSource,
             Y::ClimaCore.Fields.FieldVector,
             p::ClimaCore.Fields.FieldVector
             )

A stub function, which is extended by ClimaLSM.

"""
function source!(
    dY::ClimaCore.Fields.FieldVector,
    src::AbstractSoilSource,
    Y::ClimaCore.Fields.FieldVector,
    p::ClimaCore.Fields.FieldVector,
) end

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

end
