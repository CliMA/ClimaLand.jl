module Canopy
using DocStringExtensions
using Thermodynamics
using ClimaLand
using LazyBroadcast: lazy
using ClimaCore
using ClimaCore.MatrixFields
import ClimaCore.MatrixFields: @name, ⋅
import ClimaUtilities.TimeManager: ITime, date
import LinearAlgebra: I
using ClimaLand: AbstractRadiativeDrivers, AbstractAtmosphericDrivers
import ..Parameters as LP
import LinearAlgebra: dot

import ClimaLand:
    name,
    prognostic_vars,
    prognostic_types,
    auxiliary_vars,
    auxiliary_types,
    auxiliary_domain_names,
    prognostic_domain_names,
    initialize_prognostic,
    initialize_auxiliary,
    make_update_boundary_fluxes,
    make_update_aux,
    make_compute_exp_tendency,
    make_compute_imp_tendency,
    make_compute_jacobian,
    get_drivers,
    total_liq_water_vol_per_area!,
    total_energy_per_area!

using ClimaLand.Domains: Point, Plane, SphericalSurface
export SharedCanopyParameters, CanopyModel, set_canopy_prescribed_field!
include("./component_models.jl")
include("./ground_drivers.jl")
include("./PlantHydraulics.jl")
using .PlantHydraulics
include("./stomatalconductance.jl")
include("./photosynthesis.jl")
include("./pmodel.jl")
include("./radiation.jl")
include("./solar_induced_fluorescence.jl")
include("./pfts.jl")
include("./canopy_energy.jl")
include("./canopy_parameterizations.jl")
using Dates
include("./autotrophic_respiration.jl")

"""
    SharedCanopyParameters{FT <: AbstractFloat, PSE}

A place to store shared parameters that are required by multiple canopy components.
$(DocStringExtensions.FIELDS)
"""
struct SharedCanopyParameters{FT <: AbstractFloat, PSE}
    "Roughness length for momentum (m)"
    z_0m::FT
    "Roughness length for scalars (m)"
    z_0b::FT
    "Earth param set"
    earth_param_set::PSE
end

"""
     CanopyModel{FT, AR, RM, PM, SM, PHM, EM, SM, A, R, S, PS, D} <: ClimaLand.AbstractImExModel{FT}

The model struct for the canopy, which contains
- the canopy model domain (a point for site-level simulations, or
an extended surface (plane/spherical surface) for regional or global simulations.
- subcomponent model type for radiative transfer. This is of type
`AbstractRadiationModel`.
- subcomponent model type for photosynthesis. This is of type
`AbstractPhotosynthesisModel`, and currently only the `FarquharModel`
is supported.
- subcomponent model type for stomatal conductance. This is of type
 `AbstractStomatalConductanceModel` and currently only the `MedlynModel`
is supported
- subcomponent model type for plant hydraulics. This is of type
 `AbstractPlantHydraulicsModel` and currently only a version which
prognostically solves Richards equation in the plant is available.
- subcomponent model type for canopy energy. This is of type
 `AbstractCanopyEnergyModel` and currently we support a version where
  the canopy temperature is prescribed, and one where it is solved for
  prognostically.
- subcomponent model type for canopy SIF.
  prognostically.
- canopy model parameters, which include parameters that are shared
between canopy model components or those needed to compute boundary
fluxes.
- The boundary conditions, which contain:
    - The atmospheric conditions, which are either prescribed
      (of type `PrescribedAtmosphere`) or computed via a coupled simulation
      (of type `CoupledAtmosphere`).
    - The radiative flux conditions, which are either prescribed
      (of type `PrescribedRadiativeFluxes`) or computed via a coupled simulation
      (of type `CoupledRadiativeFluxes`).
    - The ground conditions, which are either prescribed or prognostic

Note that the canopy height is specified as part of the
PlantHydraulicsModel, along with the area indices of the leaves, roots, and
stems. Eventually, when plant biomass becomes a prognostic variable (by
integrating with a carbon model), some parameters specified here will be
treated differently.

$(DocStringExtensions.FIELDS)
"""
struct CanopyModel{FT, AR, RM, PM, SM, PHM, EM, SIFM, B, PS, D} <:
       ClimaLand.AbstractImExModel{FT}
    "Autotrophic respiration model, a canopy component model"
    autotrophic_respiration::AR
    "Radiative transfer model, a canopy component model"
    radiative_transfer::RM
    "Photosynthesis model, a canopy component model"
    photosynthesis::PM
    "Stomatal conductance model, a canopy component model"
    conductance::SM
    "Plant hydraulics model, a canopy component model"
    hydraulics::PHM
    "Energy balance model, a canopy component model"
    energy::EM
    "SIF model, a canopy component model"
    sif::SIFM
    "Boundary Conditions"
    boundary_conditions::B
    "Shared canopy parameters between component models"
    parameters::PS
    "Canopy model domain"
    domain::D
end

"""
    CanopyModel{FT}(;
        autotrophic_respiration::AbstractAutotrophicRespirationModel{FT},
        radiative_transfer::AbstractRadiationModel{FT},
        photosynthesis::AbstractPhotosynthesisModel{FT},
        conductance::AbstractStomatalConductanceModel{FT},
        hydraulics::AbstractPlantHydraulicsModel{FT},
        energy::AbstractCanopyEnergyModel{FT},
        sif::AbstractSIFModel{FT},
        boundary_conditions::B,
        parameters::SharedCanopyParameters{FT, PSE},
        domain::Union{
            ClimaLand.Domains.Point,
            ClimaLand.Domains.Plane,
            ClimaLand.Domains.SphericalSurface,
        },
        energy = PrescribedCanopyTempModel{FT}(),
    ) where {FT, PSE}

An outer constructor for the `CanopyModel`. The primary
constraints this applies are (1) ensuring that the domain is 1d or 2d
(a ``surface" domain of a column, box, or sphere) and (2) ensuring
consistency between the PlantHydraulics model and the general canopy model,
since these are built separately.
"""
function CanopyModel{FT}(;
    autotrophic_respiration::AbstractAutotrophicRespirationModel{FT},
    radiative_transfer::AbstractRadiationModel{FT},
    photosynthesis::AbstractPhotosynthesisModel{FT},
    conductance::AbstractStomatalConductanceModel{FT},
    hydraulics::AbstractPlantHydraulicsModel{FT},
    energy = PrescribedCanopyTempModel{FT}(),
    sif = Lee2015SIFModel{FT}(),
    boundary_conditions::B,
    parameters::SharedCanopyParameters{FT, PSE},
    domain::Union{
        ClimaLand.Domains.Point,
        ClimaLand.Domains.Plane,
        ClimaLand.Domains.SphericalSurface,
    },
) where {FT, B, PSE}
    if typeof(energy) <: PrescribedCanopyTempModel{FT}
        @info "Using the PrescribedAtmosphere air temperature as the canopy temperature"
        @assert typeof(boundary_conditions.atmos) <: PrescribedAtmosphere{FT}
    end

    args = (
        autotrophic_respiration,
        radiative_transfer,
        photosynthesis,
        conductance,
        hydraulics,
        energy,
        sif,
        boundary_conditions,
        parameters,
        domain,
    )
    return CanopyModel{FT, typeof.(args)...}(args...)
end

ClimaLand.name(::CanopyModel) = :canopy

"""
    canopy_components(::CanopyModel)

Returns the names of the components of the CanopyModel.

These names are used for storing prognostic and auxiliary variables
in a hierarchical manner within the state vectors.

These names must match the field names of the CanopyModel struct.
"""
canopy_components(::CanopyModel) = (
    :hydraulics,
    :conductance,
    :photosynthesis,
    :radiative_transfer,
    :autotrophic_respiration,
    :energy,
    :sif,
)

"""
    prognostic_vars(canopy::CanopyModel)

Returns the prognostic variables for the canopy model by
looping over each sub-component name in `canopy_components`.

This relies on the propertynames of `CanopyModel` being the same
as those returned by `canopy_components`.
"""
function prognostic_vars(canopy::CanopyModel)
    components = canopy_components(canopy)
    prognostic_list = map(components) do model
        prognostic_vars(getproperty(canopy, model))
    end
    return NamedTuple{components}(prognostic_list)
end

"""
    prognostic_types(canopy::CanopyModel)

Returns the prognostic types for the canopy model by
looping over each sub-component name in `canopy_components`.

This relies on the propertynames of `CanopyModel` being the same
as those returned by `canopy_components`.
"""
function prognostic_types(canopy::CanopyModel)
    components = canopy_components(canopy)
    prognostic_list = map(components) do model
        prognostic_types(getproperty(canopy, model))
    end
    return NamedTuple{components}(prognostic_list)
end

"""
    auxiliary_vars(canopy::CanopyModel)

Returns the auxiliary variables for the canopy model by
looping over each sub-component name in `canopy_components`.

This relies on the propertynames of `CanopyModel` being the same
as those returned by `canopy_components`.
"""
function auxiliary_vars(canopy::CanopyModel)
    components = canopy_components(canopy)
    auxiliary_list = map(components) do model
        auxiliary_vars(getproperty(canopy, model))
    end
    return NamedTuple{components}(auxiliary_list)
end

"""
    auxiliary_types(canopy::CanopyModel)

Returns the auxiliary types for the canopy model by
looping over each sub-component name in `canopy_components`.

This relies on the propertynames of `CanopyModel` being the same
as those returned by `canopy_components`.
"""
function auxiliary_types(canopy::CanopyModel)
    components = canopy_components(canopy)
    auxiliary_list = map(components) do model
        auxiliary_types(getproperty(canopy, model))
    end
    return NamedTuple{components}(auxiliary_list)
end

"""
    filter_nt(nt::NamedTuple)

Removes all key/value pairs of a NamedTuple where the value is `nothing`.
Note that NamedTuples are immutable, so rather than updating the input
in-place, this creates a new NamedTuple with the filtered key/value pairs.

This results in unnecessary allocations because a new object is being
created, and we may want to implement a better solution in the future.
"""
function filter_nt(nt::NamedTuple)
    pairs = []
    for (k, v) in (zip(keys(nt), values(nt)))
        ~(isnothing(v)) ? push!(pairs, k => (filter_nt(v))) : (;)
    end
    return NamedTuple(pairs)
end

"""
    filter_nt(nt)

Base case for `filter_nt` recursion, used when this function is called on
a NamedTuple with no nested NamedTuples.
"""
filter_nt(nt) = nt

"""
    initialize_prognostic(
        model::CanopyModel{FT},
        coords,
    ) where {FT}

Creates the prognostic state vector of the `CanopyModel` and returns
it as a ClimaCore.Fields.FieldVector.

The input `state` is usually a ClimaCore Field object.

This function loops over the components of the `CanopyModel` and appends
each component models prognostic state vector into a single state vector,
structured by component name.
"""
function initialize_prognostic(model::CanopyModel{FT}, coords) where {FT}
    components = canopy_components(model)
    Y_state_list = map(components) do (component)
        submodel = getproperty(model, component)
        getproperty(initialize_prognostic(submodel, coords), component)
    end
    # `Y_state_list` contains `nothing` for components with no prognostic
    #  variables, which we need to filter out before constructing `Y`
    Y = ClimaCore.Fields.FieldVector(;
        name(model) => filter_nt(NamedTuple{components}(Y_state_list)),
    )
    return Y
end

"""
    initialize_auxiliary(
        model::CanopyModel{FT},
        coords,
    ) where {FT}

Creates the auxiliary state vector of the `CanopyModel` and returns
 it as a ClimaCore.Fields.FieldVector.

The input `coords` is usually a ClimaCore Field object.

This function loops over the components of the `CanopyModel` and appends
each component models auxiliary state vector into a single state vector,
structured by component name.
"""
function initialize_auxiliary(model::CanopyModel{FT}, coords) where {FT}
    components = canopy_components(model)
    p_state_list = map(components) do (component)
        submodel = getproperty(model, component)
        getproperty(initialize_auxiliary(submodel, coords), component)
    end
    # `p_state_list` contains `nothing` for components with no auxiliary
    #  variables, which we need to filter out before constructing `p`
    # We also add in boundary variables here.
    p = (;
        name(model) => (;
            filter_nt(NamedTuple{components}(p_state_list))...,
            initialize_boundary_vars(model, coords)...,
        )
    )
    p = ClimaLand.add_dss_buffer_to_aux(p, model.domain)
    return p
end

"""
    initialize_boundary_vars(model::CanopyModel{FT}, coords)

Add boundary condition-related variables to the cache.
This calls functions defined in canopy_boundary_fluxes.jl
which dispatch on the boundary condition type to add the correct variables.
"""
function initialize_boundary_vars(model::CanopyModel{FT}, coords) where {FT}
    vars = boundary_vars(model.boundary_conditions, ClimaLand.TopBoundary())
    types = boundary_var_types(
        model,
        model.boundary_conditions,
        ClimaLand.TopBoundary(),
    )
    domains = boundary_var_domain_names(
        model.boundary_conditions,
        ClimaLand.TopBoundary(),
    )
    additional_aux = map(zip(types, domains)) do (T, domain)
        zero_instance = ClimaCore.RecursiveApply.rzero(T)
        f = map(_ -> zero_instance, getproperty(coords, domain))
        fill!(ClimaCore.Fields.field_values(f), zero_instance)
        f
    end
    return NamedTuple{vars}(additional_aux)
end

"""
     ClimaLand.make_update_aux(canopy::CanopyModel{FT,
                                                  <:AutotrophicRespirationModel,
                                                  <:Union{BeerLambertModel, TwoStreamModel},
                                                  <:FarquharModel,
                                                  <:MedlynConductanceModel,
                                                  <:PlantHydraulicsModel,},
                              ) where {FT}

Creates the `update_aux!` function for the `CanopyModel`; a specific
method for `update_aux!` for the case where the canopy model components
are of the type in the parametric type signature: `AutotrophicRespirationModel`, `AbstractRadiationModel`,
`FarquharModel`, `MedlynConductanceModel`, and `PlantHydraulicsModel`.

Please note that the plant hydraulics model has auxiliary variables
that are updated in its prognostic `compute_exp_tendency!` function.
While confusing, this is better for performance as it saves looping
over the state vector multiple times.

The other sub-components rely heavily on each other,
so the version of the `CanopyModel` with these subcomponents
has a single update_aux! function, given here.
"""
function ClimaLand.make_update_aux(
    canopy::CanopyModel{
        FT,
        <:AutotrophicRespirationModel,
        <:Union{BeerLambertModel, TwoStreamModel},
        <:Union{FarquharModel, OptimalityFarquharModel},
        <:MedlynConductanceModel,
        <:PlantHydraulicsModel,
        <:AbstractCanopyEnergyModel,
    },
) where {FT}
    function update_aux!(p, Y, t)

        # Extend to other fields when necessary
        # Update the prescribed fields to the current time `t`,
        # prior to updating the rest of the auxiliary state to
        # the current time, as they depend on prescribed fields.
        set_canopy_prescribed_field!(canopy.hydraulics, p, t)

        # Shortcut names
        An = p.canopy.photosynthesis.An
        ψ = p.canopy.hydraulics.ψ
        ϑ_l = Y.canopy.hydraulics.ϑ_l
        fa = p.canopy.hydraulics.fa
        par_d = p.canopy.radiative_transfer.par_d
        nir_d = p.canopy.radiative_transfer.nir_d
        cosθs = p.drivers.cosθs
        area_index = p.canopy.hydraulics.area_index
        LAI = area_index.leaf
        SAI = area_index.stem

        bc = canopy.boundary_conditions

        # update radiative transfer
        (; G_Function, Ω, λ_γ_PAR) = canopy.radiative_transfer.parameters
        @. p.canopy.radiative_transfer.ϵ =
            canopy.radiative_transfer.parameters.ϵ_canopy *
            (1 - exp(-(LAI + SAI))) #from CLM 5.0, Tech note 4.20
        RT = canopy.radiative_transfer
        compute_PAR!(par_d, RT, bc.radiation, p, t)
        compute_NIR!(nir_d, RT, bc.radiation, p, t)

        compute_fractional_absorbances!(
            p,
            RT,
            LAI,
            ground_albedo_PAR(
                Val(bc.prognostic_land_components),
                bc.ground,
                Y,
                p,
                t,
            ),
            ground_albedo_NIR(
                Val(bc.prognostic_land_components),
                bc.ground,
                Y,
                p,
                t,
            ),
        )

        # update plant hydraulics aux
        hydraulics = canopy.hydraulics
        n_stem = hydraulics.n_stem
        n_leaf = hydraulics.n_leaf
        PlantHydraulics.lai_consistency_check.(n_stem, n_leaf, area_index)
        (; retention_model, conductivity_model, S_s, ν) = hydraulics.parameters
        # We can index into a field of Tuple{FT} to extract a field of FT
        # using the following notation: field.:index
        @inbounds @. ψ.:1 = PlantHydraulics.water_retention_curve(
            retention_model,
            PlantHydraulics.effective_saturation(ν, ϑ_l.:1),
            ν,
            S_s,
        )
        # Inside of a loop, we need to use a single dollar sign
        # for indexing into Fields of Tuples in non broadcasted
        # expressions, and two dollar signs for
        # for broadcasted expressions using the macro @.
        # field.:($index) .= value # works
        # @ field.:($$index) = value # works
        @inbounds for i in 1:(n_stem + n_leaf - 1)
            ip1 = i + 1
            @. ψ.:($$ip1) = PlantHydraulics.water_retention_curve(
                retention_model,
                PlantHydraulics.effective_saturation(ν, ϑ_l.:($$ip1)),
                ν,
                S_s,
            )

            areai = getproperty(area_index, hydraulics.compartment_labels[i])
            areaip1 =
                getproperty(area_index, hydraulics.compartment_labels[ip1])

            # Compute the flux*area between the current compartment `i`
            # and the compartment above.
            @. fa.:($$i) =
                PlantHydraulics.water_flux(
                    hydraulics.compartment_midpoints[i],
                    hydraulics.compartment_midpoints[ip1],
                    ψ.:($$i),
                    ψ.:($$ip1),
                    PlantHydraulics.hydraulic_conductivity(
                        conductivity_model,
                        ψ.:($$i),
                    ),
                    PlantHydraulics.hydraulic_conductivity(
                        conductivity_model,
                        ψ.:($$ip1),
                    ),
                ) * PlantHydraulics.harmonic_mean(areaip1, areai)
        end
        # We update the fa[n_stem+n_leaf] element once we have computed transpiration

        # Update Rd, An, Vcmax25 (if applicable to model) in place, GPP
        update_photosynthesis!(p, Y, canopy.photosynthesis, canopy)

        # update SIF
        update_SIF!(p, Y, canopy.sif, canopy)

        # update stomatal conductance
        update_canopy_conductance!(p, Y, canopy.conductance, canopy)

        # update autotrophic respiration
        update_autotrophic_respiration!(
            p,
            Y,
            canopy.autotrophic_respiration,
            canopy,
        )
    end
    return update_aux!
end

"""
    make_compute_exp_tendency(canopy::CanopyModel)

Creates and returns the compute_exp_tendency! for the `CanopyModel`.
"""
function make_compute_exp_tendency(
    canopy::CanopyModel{
        FT,
        <:AutotrophicRespirationModel,
        <:Union{BeerLambertModel, TwoStreamModel},
        <:Union{FarquharModel, OptimalityFarquharModel},
        <:MedlynConductanceModel,
        <:PlantHydraulicsModel,
        <:Union{PrescribedCanopyTempModel, BigLeafEnergyModel},
    },
) where {FT}
    components = canopy_components(canopy)
    compute_exp_tendency_list = map(
        x -> make_compute_exp_tendency(getproperty(canopy, x), canopy),
        components,
    )
    function compute_exp_tendency!(dY, Y, p, t)
        for f! in compute_exp_tendency_list
            f!(dY, Y, p, t)
        end

    end
    return compute_exp_tendency!
end

"""
    make_compute_imp_tendency(canopy::CanopyModel)

Creates and returns the compute_imp_tendency! for the `CanopyModel`.
"""
function make_compute_imp_tendency(
    canopy::CanopyModel{
        FT,
        <:AutotrophicRespirationModel,
        <:Union{BeerLambertModel, TwoStreamModel},
        <:Union{FarquharModel, OptimalityFarquharModel},
        <:MedlynConductanceModel,
        <:PlantHydraulicsModel,
        <:Union{PrescribedCanopyTempModel, BigLeafEnergyModel},
    },
) where {FT}
    components = canopy_components(canopy)
    compute_imp_tendency_list = map(
        x -> make_compute_imp_tendency(getproperty(canopy, x), canopy),
        components,
    )
    function compute_imp_tendency!(dY, Y, p, t)
        for f! in compute_imp_tendency_list
            f!(dY, Y, p, t)
        end

    end
    return compute_imp_tendency!
end

"""
    ClimaLand.make_compute_jacobian(canopy::CanopyModel)

Creates and returns the compute_jacobian! for the `CanopyModel`.
"""
function ClimaLand.make_compute_jacobian(
    canopy::CanopyModel{
        FT,
        <:AutotrophicRespirationModel,
        <:Union{BeerLambertModel, TwoStreamModel},
        <:Union{FarquharModel, OptimalityFarquharModel},
        <:MedlynConductanceModel,
        <:PlantHydraulicsModel,
        <:Union{PrescribedCanopyTempModel, BigLeafEnergyModel},
    },
) where {FT}
    components = canopy_components(canopy)
    update_jacobian_list = map(
        x -> make_compute_jacobian(getproperty(canopy, x), canopy),
        components,
    )
    function compute_jacobian!(W, Y, p, dtγ, t)
        for f! in update_jacobian_list
            f!(W, Y, p, dtγ, t)
        end

    end
    return compute_jacobian!
end


function ClimaLand.get_drivers(model::CanopyModel)
    return (
        model.boundary_conditions.atmos,
        model.boundary_conditions.radiation,
    )
end
include("./canopy_boundary_fluxes.jl")
#Make the canopy model broadcastable
Base.broadcastable(C::CanopyModel) = tuple(C)

"""
    ClimaLand.total_energy_per_area!(
        surface_field,
        model::CanopyModel,
        Y,
        p,
        t,
)

A function which updates `surface_field` in place with the value for
the total energy per unit ground area for the `CanopyModel`.

This acts by calling the method for the energy component of
the canopy model.
"""
function ClimaLand.total_energy_per_area!(
    surface_field,
    model::CanopyModel,
    Y,
    p,
    t,
)
    ClimaLand.total_energy_per_area!(surface_field, model.energy, Y, p, t)
end

"""
    ClimaLand.total_liq_water_vol_per_area!(
        surface_field,
        model::CanopyModel,
        Y,
        p,
        t,
)

A function which updates `surface_field` in place with the value for
the total liquid water volume per unit ground area for the `CanopyModel`.

This acts by calling the method for the PlantHydraulics component of
the canopy model.
"""
function ClimaLand.total_liq_water_vol_per_area!(
    surface_field,
    model::CanopyModel,
    Y,
    p,
    t,
)
    ClimaLand.total_liq_water_vol_per_area!(
        surface_field,
        model.hydraulics,
        Y,
        p,
        t,
    )
end

end
