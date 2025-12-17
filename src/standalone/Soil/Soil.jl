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
soil, rivers, etc), the ClimaLand.jl package should be used.
That package will use the methods of this function for advancing
the system forward in time, extending methods as needed to account
for interactions between components.
=#

using ClimaLand
using LazyBroadcast: lazy
using DocStringExtensions
using LinearAlgebra
using ClimaCore
using ClimaCore.MatrixFields
import ClimaCore.MatrixFields: @name, ⋅
import ..Parameters as LP
import ClimaCore: Fields, Operators, Geometry, Spaces
using Thermodynamics
import ClimaParams as CP
using SurfaceFluxes
using StaticArrays
import SurfaceFluxes.Parameters as SFP
import ClimaUtilities.SpaceVaryingInputs: SpaceVaryingInput
import ClimaLand.Domains: Column, HybridBox, SphericalShell
import ClimaLand:
    AbstractImExModel,
    make_update_aux,
    make_update_imp_aux,
    make_compute_exp_tendency,
    make_compute_imp_tendency,
    make_update_boundary_fluxes,
    make_compute_jacobian,
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
    surface_albedo,
    surface_emissivity,
    surface_height,
    surface_resistance,
    turbulent_fluxes!,
    get_drivers,
    total_liq_water_vol_per_area!,
    total_energy_per_area!,
    return_momentum_fluxes
export RichardsModel,
    RichardsParameters,
    EnergyHydrology,
    EnergyHydrologyParameters,
    AbstractSoilModel,
    AbstractSoilSource,
    PhaseChange



"""
    AbstractSoilSource{FT} <:  ClimaLand.AbstractSource{FT}

An abstract type for types of source terms for the soil equations.

In standalone mode, the only supported source type is freezing and
thawing. ClimaLand.jl creates additional sources to include as
necessary e.g. root extraction (not available in stand alone mode).
"""
abstract type AbstractSoilSource{FT} <: ClimaLand.AbstractSource{FT} end

"""
    AbstractSoilModel{FT} <: ClimaLand.AbstractImExModel{FT}

The abstract type for all soil models.

Currently, we only have plans to support a RichardsModel, simulating
the flow of liquid water through soil via the Richardson-Richards equation,
and a fully integrated soil heat and water model, with phase change.
"""
abstract type AbstractSoilModel{FT} <: ClimaLand.AbstractImExModel{FT} end

ClimaLand.name(::AbstractSoilModel) = :soil
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

"""
    append_source(src::AbstractSoilSource, srcs::Tuple)::Tuple
Appends `src` to the tuple of sources `srcs` if `src` is of type `AbstractSoilSource`.
"""
append_source(src::AbstractSoilSource, srcs::Tuple)::Tuple = (srcs..., src)

"""
    append_source(src::Nothing , srcs::Tuple)::Tuple
Appends `src` to the tuple of sources `srcs` if `src` is of type `AbstractSoilSource`.
"""
append_source(src::Nothing, srcs::Tuple)::Tuple = srcs

include("./retention_models.jl")
include("./rre.jl")
include("./energy_hydrology.jl")
include("./soil_heat_parameterizations.jl")
include("Runoff/Runoff.jl")
using .Runoff
include("./boundary_conditions.jl")
include("./soil_albedo.jl")
include("./soil_hydrology_parameterizations.jl")
include("./spatially_varying_parameters.jl")
include("Biogeochemistry/Biogeochemistry.jl")
using .Biogeochemistry


# Soil model constructor useful for working with simulations forced by
# the atmosphere
"""
    EnergyHydrology{FT}(domain, forcing, toml_dict::CP.ParamDict;
                         prognostic_land_components = (:soil),
                         albedo = CLMTwoBandSoilAlbedo{FT}(; clm_soil_albedo_parameters(domain.space.surface)...),
                         runoff::Runoff.AbstractRunoffModel = Runoff.TOPMODELRunoff(toml_dict,
                                                                f_max = topmodel_fmax(domain.space.surface, FT),
                                                            ),
                         retention_parameters = soil_vangenuchten_parameters(domain.space.subsurface, FT),
                         composition_parameters = soil_composition_parameters(domain.space.subsurface, FT),
                         S_s = ClimaCore.Fields.zeros(domain.space.subsurface) .+ 1e-3,
                         z_0m = toml_dict["soil_momentum_roughness_length"],
                         z_0b = toml_dict["soil_scalar_roughness_length"],
                         emissivity = toml_dict["emissivity_bare_soil"],
                         additional_sources = (),
                         ) where {FT <: AbstractFloat}

Creates a EnergyHydrology model with the given float type `FT`, `domain`,
`toml_dict`, `forcing`, and prognostic land components.

The argument `forcing` should be a NamedTuple containing two fields: `atmos` and
`radiation`.

When running the soil model in standalone mode, `prognostic_land_components =
(:soil,)`, while for running integrated land models, this should be a list of
the component models. This value of this argument must be the same across all
components in the land model.

Default spatially varying parameters (for retention curve parameters,
composition, and specific storativity) are provided but can be changed with
keyword arguments. Note that these parameters must all be of the same type:
either `FT` or ClimaCore Fields. By default they are Fields read in from data,
so in practice this means if some values are provided as Floats, all of these
parameter defaults must be overwritten as Floats.

`retention_parameters` should be a NamedTuple with the following fields:
- `ν`: soil porosity
- `hydrology_cm`: the hydrology closure model being used
- `K_sat`: saturated hydraulic conductivity
- `θ_r`: residual water content
`composition_parameters` should be a NamedTuple with the following fields:
- `ν_ss_om`: soil organic matter volume fraction
- `ν_ss_quartz`: soil quartz volume fraction
- `ν_ss_gravel`: soil gravel volume fraction

The runoff and albedo parameterizations are also provided and can be changed via keyword argument;
additional sources may be required in your model if the soil model will be composed with other
component models.

Roughness lengths and soil emissivity are currently treated as constants; these can be passed in as Floats
by kwarg; otherwise the default values are used.
"""
function EnergyHydrology{FT}(
    domain,
    forcing,
    toml_dict::CP.ParamDict;
    prognostic_land_components = (:soil,),
    albedo::AbstractSoilAlbedoParameterization = CLMTwoBandSoilAlbedo{FT}(;
        clm_soil_albedo_parameters(domain.space.surface)...,
    ),
    runoff::Runoff.AbstractRunoffModel = Runoff.TOPMODELRunoff(
        toml_dict,
        f_max = topmodel_fmax(domain.space.surface, FT),
    ),
    retention_parameters = soil_vangenuchten_parameters(
        domain.space.subsurface,
        FT,
    ),
    composition_parameters = soil_composition_parameters(
        domain.space.subsurface,
        FT,
    ),
    S_s = ClimaCore.Fields.zeros(domain.space.subsurface) .+ FT(1e-3),
    z_0m = toml_dict["soil_momentum_roughness_length"],
    z_0b = toml_dict["soil_scalar_roughness_length"],
    emissivity = toml_dict["emissivity_bare_soil"],
    additional_sources = (),
) where {FT <: AbstractFloat}
    # TODO: Move runoff scalar parameters to ClimaParams, possibly use types in retention, composition,
    #  roughness, and emissivity.
    top_bc = AtmosDrivenFluxBC(
        forcing.atmos,
        forcing.radiation,
        runoff;
        prognostic_land_components,
    )
    bottom_bc = EnergyWaterFreeDrainage()
    boundary_conditions = (; top = top_bc, bottom = bottom_bc)
    # sublimation and subsurface runoff are added automatically
    sources = (additional_sources..., PhaseChange{FT}())
    parameters = EnergyHydrologyParameters(
        toml_dict;
        retention_parameters...,
        composition_parameters...,
        albedo,
        S_s,
        z_0m,
        z_0b,
        emissivity,
    )
    return EnergyHydrology{FT}(;
        parameters,
        domain,
        boundary_conditions,
        sources,
    )
end


"""
    RichardsModel{FT}(domain, forcing;
                         runoff =  ClimaLand.Soil.Runoff.TOPMODELRunoff{FT}(f_over = FT(3.28),
                                                                            R_sb = FT(1.484e-4 / 1000),
                                                                            f_max = topmodel_fmax(domain.space.surface,FT),
                                                                            ),
                         retention_parameters = soil_vangenuchten_parameters(domain.space.subsurface, FT),
                         S_s = ClimaCore.Fields.zeros(domain.space.subsurface) .+ 1e-3,
                         )

Creates a RichardsModel model with the given float type `FT`, `domain` and `forcing`.
Here, `forcing` should be a `NamedTuple` containing a field `atmos` with the atmospheric forcing.

Default spatially varying parameters (for retention curve parameters and specific storativity) are provided but can be
changed with keyword arguments.

The runoff parameterization is also provided and can be changed via keyword argument.
"""
function RichardsModel{FT}(
    domain,
    forcing;
    runoff::Runoff.AbstractRunoffModel = Runoff.TOPMODELRunoff{FT}(
        f_over = FT(3.28), # extract from EPS
        R_sb = FT(1.484e-4 / 1000),# extract from EPS
        f_max = topmodel_fmax(domain.space.surface, FT),
    ),
    retention_parameters = soil_vangenuchten_parameters(
        domain.space.subsurface,
        FT,
    ),
    S_s = ClimaCore.Fields.zeros(domain.space.subsurface) .+ 1e-3,
) where {FT}
    # TODO: Move scalar parameters to ClimaParams and obtain from earth_param_set, possibly use types in retention argument.
    top_bc = RichardsAtmosDrivenFluxBC(forcing.atmos, runoff)
    bottom_bc = WaterFluxBC((p, t) -> 0.0)
    boundary_conditions = (; top = top_bc, bottom = bottom_bc)
    parameters = RichardsParameters(retention_parameters..., S_s)
    return RichardsModel{FT}(;
        parameters,
        domain,
        boundary_conditions,
        sources = (),
        lateral_flow = false,
    )
end

end
