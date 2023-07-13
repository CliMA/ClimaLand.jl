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
using DocStringExtensions
using LinearAlgebra
using ClimaCore
import ..Parameters as LSMP
import ClimaCore: Fields, Operators, Geometry, Spaces
using Thermodynamics

import ClimaLSM.Domains: Column, HybridBox, SphericalShell
using ClimaLSM: AbstractTridiagonalW

import ClimaLSM:
    AbstractImExModel,
    make_update_aux,
    make_compute_exp_tendency,
    make_compute_imp_tendency,
    make_set_initial_aux_state,
    make_update_jacobian,
    prognostic_vars,
    auxiliary_vars,
    prognostic_domain_names,
    auxiliary_domain_names,
    name,
    prognostic_types,
    auxiliary_types,
    AbstractSource,
    source!,
    heaviside,
    surface_temperature,
    surface_specific_humidity,
    surface_albedo,
    surface_emissivity,
    surface_air_density,
    surface_height
export RichardsModel,
    RichardsParameters,
    RichardsTridiagonalW,
    EnergyHydrology,
    EnergyHydrologyParameters,
    AbstractSoilModel,
    AbstractSoilSource,
    PhaseChange
"""
    AbstractSoilModel{FT} <: ClimaLSM.AbstractImExModel{FT}

The abstract type for all soil models.

Currently, we only have plans to support a RichardsModel, simulating
the flow of liquid water through soil via the Richardson-Richards equation,
and a fully integrated soil heat and water model, with phase change.
"""
abstract type AbstractSoilModel{FT} <: ClimaLSM.AbstractImExModel{FT} end

ClimaLSM.name(::AbstractSoilModel) = :soil
ClimaLSM.domain_name(::AbstractSoilModel) = :subsurface

"""
    AbstractSoilSource{FT} <:  ClimaLSM.AbstractSource{FT}

An abstract type for types of source terms for the soil equations.

In standalone mode, the only supported source type is freezing and
thawing. ClimaLSM.jl creates additional sources to include as
necessary e.g. root extraction (not available in stand alone mode).
"""
abstract type AbstractSoilSource{FT} <: ClimaLSM.AbstractSource{FT} end


function ClimaLSM.make_set_initial_aux_state(model::AbstractSoilModel)
    update_aux! = make_update_aux(model)
    function set_initial_aux_state!(p, Y0, t0)
        set_initial_parameter_field!(
            model.boundary_conditions,
            p,
            model.domain.surface.space,
        )

        update_aux!(p, Y0, t0)
    end
    return set_initial_aux_state!
end

ClimaLSM.name(::AbstractSoilModel) = :soil
    #=
function ClimaLSM.initialize_auxiliary(
    model::AbstractSoilModel{FT},
    _,
) where {FT}
    model_name = name(model)
    vars = auxiliary_vars(model)
    types = auxiliary_types(model)
    domain_names = auxiliary_domain_names(model)

    auxiliary_state = map(zip(types, domain_names)) do (T, domain_name)
        zero_instance = zero(T)
        space =
            domain_name == :surface ?
            ClimaLSM.Domains.obtain_surface_space(model.domain.space) :
            model.domain.space
        cds = ClimaCore.Fields.coordinate_field(space)
        map(_ -> zero_instance, cds)
    end
    return NamedTuple{vars}(auxiliary_state)
end
=#

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
   horizontal_components!(dY::ClimaCore.Fields.FieldVector,
                          domain::Union{HybridBox, SphericalShell},
                          lateral_flow::Val{false},
                          _...)
Updates dY in place by adding in the tendency terms resulting from
horizontal derivative operators.

In the case of a 3D domain, for either the `RichardsModel` or the
`EnergyHydrology` model, if the `lateral_flow` flag is set to false,
there are no horizontal contributions to the right hand side.
"""
function horizontal_components!(
    dY::ClimaCore.Fields.FieldVector,
    domain::Union{HybridBox, SphericalShell},
    lateral_flow::Val{false},
    _...,
) end
include("Runoff/Runoff.jl")
using .Runoff
include("./retention_models.jl")
include("./rre.jl")
include("./energy_hydrology.jl")
include("./boundary_conditions.jl")
include("./soil_hydrology_parameterizations.jl")
include("./soil_heat_parameterizations.jl")
include("Biogeochemistry/Biogeochemistry.jl")
using .Biogeochemistry
end
