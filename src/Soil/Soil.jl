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

Currently, only Richards equation is supported, but there
are immediate plans to add in a model with coupled water and
energy, including phase changes.

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
import ClimaCore: Fields, Operators, Geometry

import ClimaLSM.Domains: coordinates
import ClimaLSM:
    AbstractModel,
    make_update_aux,
    make_rhs,
    prognostic_vars,
    auxiliary_vars,
    name
export RichardsModel,
    RichardsParameters,
    FluxBC,
    RootExtraction,
    AbstractSoilModel,
    AbstractSoilSource,
    source

"""
    AbstractSoilBoundaryConditions{FT <: AbstractFloat}

An abstract type for types of boundary conditions, which will include
prescribed functions of space and time as Dirichlet conditions or
Neumann conditions, in addition to other 
convenient soil-specific conditions, like free drainage.
"""
abstract type AbstractSoilBoundaryConditions{FT <: AbstractFloat} end

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
the flow of liquid water through soil, and a fully integrated soil heat
and water model, with phase change.
"""
abstract type AbstractSoilModel{FT} <: AbstractModel{FT} end

ClimaLSM.name(::AbstractSoilModel) = :soil

"""
    RichardsParameters{FT <: AbstractFloat}

A struct for storing parameters of the `RichardModel`.
$(DocStringExtensions.FIELDS)
"""
struct RichardsParameters{FT <: AbstractFloat}
    "The porosity of the soil (m^3/m^3)"
    ν::FT
    "The van Genuchten parameter α (1/m)"
    vg_α::FT
    "The van Genuchten parameter n"
    vg_n::FT
    "The van Genuchten parameter m"
    vg_m::FT
    "The saturated hydraulic conductivity (m/s)"
    Ksat::FT
    "The specific storativity (1/m)"
    S_s::FT
    "The residual water fraction (m^3/m^3"
    θ_r::FT
end

"""
    RichardsModel

A model for simulating the flow of water in a porous medium
by solving the Richardson-Richards Equation.

$(DocStringExtensions.FIELDS)
"""
struct RichardsModel{FT, PS, D, C, BC, S} <: AbstractSoilModel{FT}
    "the parameter set"
    param_set::PS
    "the soil domain, using ClimaCore.Domains"
    domain::D
    "the domain coordinates"
    coordinates::C
    "the boundary conditions, of type AbstractSoilBoundaryConditions"
    boundary_conditions::BC
    "A tuple of sources, each of type AbstractSoilSource"
    sources::S
end

"""
    RichardsModel{FT}(;
        param_set::RichardsParameters{FT},
        domain::D,
        boundary_conditions::AbstractSoilBoundaryConditions{FT},
        sources::Tuple,
    ) where {FT, D}

A constructor for a `RichardsModel`.
"""
function RichardsModel{FT}(;
    param_set::RichardsParameters{FT},
    domain::D,
    boundary_conditions::AbstractSoilBoundaryConditions{FT},
    sources::Tuple,
) where {FT, D}
    coords = coordinates(domain)
    args = (param_set, domain, coords, boundary_conditions, sources)
    RichardsModel{FT, typeof.(args)...}(args...)
end

"""
    coordinates(model::RichardsModel)

A extension of the `coordinates` function, which returns the coordinates
of a model domain. 

The coordinates are stored in the model because they are required in
computing the right hand side. 
"""
coordinates(model::RichardsModel) = model.coordinates

"""
    volumetric_liquid_fraction(ϑ_l::FT, ν_eff::FT) where {FT}

A pointwise function returning the volumetric liquid fraction
given the augmented liquid fraction and the effective porosity.
"""
function volumetric_liquid_fraction(ϑ_l::FT, ν_eff::FT) where {FT}
    if ϑ_l < ν_eff
        θ_l = ϑ_l
    else
        θ_l = ν_eff
    end
    return θ_l
end


"""
     matric_potential(α::FT, n::FT, m::FT, S::FT) where {FT}

A point-wise function returning the matric potential, using the
van Genuchten formulation.
"""
function matric_potential(α::FT, n::FT, m::FT, S::FT) where {FT}
    ψ_m = -((S^(-FT(1) / m) - FT(1)) * α^(-n))^(FT(1) / n)
    return ψ_m
end

"""
    effective_saturation(porosity::FT, ϑ_l::FT, θr::FT) where {FT}

A point-wise function computing the effective saturation.
"""
function effective_saturation(porosity::FT, ϑ_l::FT, θr::FT) where {FT}
    ϑ_l_safe = max(ϑ_l, θr + eps(FT))
    S_l = (ϑ_l_safe - θr) / (porosity - θr)
    return S_l
end

"""
    pressure_head(
        α::FT,
        n::FT,
        m::FT,
        θ_r::FT,
        ϑ_l::FT,
        ν_eff::FT,
        S_s::FT,
    ) where {FT}

A point-wise function returning the pressure head in
variably saturated soil, using the van Genuchten formulation
for matric potential if the soil is not saturated, and
an approximation of the positive pressure in the soil if the
soil is saturated.
"""
function pressure_head(
    α::FT,
    n::FT,
    m::FT,
    θ_r::FT,
    ϑ_l::FT,
    ν_eff::FT,
    S_s::FT,
) where {FT}
    S_l_eff = effective_saturation(ν_eff, ϑ_l, θ_r)
    if S_l_eff <= FT(1.0)
        ψ = matric_potential(α, n, m, S_l_eff)
    else
        ψ = (ϑ_l - ν_eff) / S_s
    end
    return ψ
end

"""
     hydraulic_conductivity(Ksat::FT, m::FT, S::FT) where {FT}

A point-wise function returning the hydraulic conductivity, using the
van Genuchten formulation.
"""
function hydraulic_conductivity(Ksat::FT, m::FT, S::FT) where {FT}
    if S < FT(1)
        K = sqrt(S) * (FT(1) - (FT(1) - S^(FT(1) / m))^m)^FT(2)
    else
        K = FT(1)
    end
    return K * Ksat
end

"""
    make_rhs(model::RichardsModel)

An extension of the function `make_rhs`, for the Richardson-
Richards equation. 

This function creates and returns a function which computes the entire
right hand side of the PDE for `ϑ_l`, and updates `dY.soil.ϑ_l` in place
with that value.

This has been written so as to work with Differential Equations.jl.
"""
function make_rhs(model::RichardsModel)
    function rhs!(dY, Y, p, t)
        @unpack ν, vg_α, vg_n, vg_m, Ksat, S_s, θ_r = model.param_set
        top_flux_bc, bot_flux_bc =
            boundary_fluxes(model.boundary_conditions, p, t)
        z = model.coordinates
        interpc2f = Operators.InterpolateC2F()
        gradc2f_water = Operators.GradientC2F()
        divf2c_water = Operators.DivergenceF2C(
            top = Operators.SetValue(Geometry.WVector(top_flux_bc)),
            bottom = Operators.SetValue(Geometry.WVector(bot_flux_bc)),
        )
        @. dY.soil.ϑ_l =
            -(divf2c_water(-interpc2f(p.soil.K) * gradc2f_water(p.soil.ψ + z)))
        for src in model.sources
            dY.soil.ϑ_l .+= source(src, Y, p)
        end

    end
    return rhs!
end

"""
    prognostic_vars(soil::RichardsModel)

A function which returns the names of the prognostic variables
of `RichardsModel`.
"""
prognostic_vars(soil::RichardsModel) = (:ϑ_l,)

"""
    auxiliary_vars(soil::RichardsModel)

A function which returns the names of the auxiliary variables 
of `RichardsModel`.

Note that auxiliary variables are not needed for such a simple model.
We could instead compute the conductivity and matric potential within the
rhs function explicitly, rather than compute and store them in the 
auxiliary vector `p`. We did so in this case as a demonstration.
"""
auxiliary_vars(soil::RichardsModel) = (:K, :ψ)

"""
    make_update_aux(model::RichardsModel)

An extension of the function `make_update_aux`, for the Richardson-
Richards equation. 

This function creates and returns a function which updates the auxiliary
variables `p.soil.variable` in place.

This has been written so as to work with Differential Equations.jl.
"""
function make_update_aux(model::RichardsModel)
    function update_aux!(p, Y, t)
        @unpack ν, vg_α, vg_n, vg_m, Ksat, S_s, θ_r = model.param_set
        @. p.soil.K = hydraulic_conductivity(
            Ksat,
            vg_m,
            effective_saturation(ν, Y.soil.ϑ_l, θ_r),
        )
        @. p.soil.ψ = pressure_head(vg_α, vg_n, vg_m, θ_r, Y.soil.ϑ_l, ν, S_s)
    end
    return update_aux!
end

"""
   FluxBC{FT} <: AbstractSoilBoundaryConditions{FT}

A simple concrete type of boundary condition, which enforces
constant fluxes at the top and bottom of the domain.
"""
struct FluxBC{FT} <: AbstractSoilBoundaryConditions{FT}
    top_flux_bc::FT
    bot_flux_bc::FT
end

"""
    boundary_fluxes(bc::FluxBC, _...)

A function which returns the correct boundary flux
given the boundary condition type `FluxBC`.

This is a trivial example, but a more complex one would be e.g.
Dirichlet conditions on the state, which then must be converted into
a flux before being applied as a boundary condition.
"""
function boundary_fluxes(bc::FluxBC, _...)
    return bc.top_flux_bc, bc.bot_flux_bc
end

"""
    source(src::AbstractSoilSource, _...) end

A stub function, which is extended by ClimaLSM.

Once we have the freeze thaw source function, we do not need this
stub anymore
"""
function source(src::AbstractSoilSource, _...) end

end
