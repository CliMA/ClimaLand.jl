module Canopy
using DocStringExtensions
using Thermodynamics
using ClimaLSM
using ClimaCore
using ClimaLSM: AbstractRadiativeDrivers, AbstractAtmosphericDrivers
import ..Parameters as LSMP

import ClimaLSM:
    AbstractExpModel,
    name,
    domain,
    prognostic_vars,
    prognostic_types,
    auxiliary_vars,
    auxiliary_types,
    initialize_prognostic,
    initialize_auxiliary,
    make_update_aux,
    make_compute_exp_tendency,
    surface_temperature,
    surface_specific_humidity,
    surface_air_density,
    surface_evaporative_scaling,
    surface_height

using ClimaLSM.Domains: Point, Plane, SphericalSurface
export SharedCanopyParameters, CanopyModel
include("./component_models.jl")
include("./PlantHydraulics.jl")
using .PlantHydraulics
include("./stomatalconductance.jl")
include("./photosynthesis.jl")
include("./radiation.jl")
include("./canopy_parameterizations.jl")

"""
    SharedCanopyParameters{FT <: AbstractFloat, PSE}

A place to store shared parameters that are required by all canopy components.
$(DocStringExtensions.FIELDS)
"""
struct SharedCanopyParameters{FT <: AbstractFloat, PSE}
    "Leaf Area Index (m2/m2)"
    LAI::FT
    "Canopy height (m)"
    h_c::FT
    "Roughness length for momentum (m)"
    z_0m::FT
    "Roughness length for scalars (m)"
    z_0b::FT
    "Earth param set"
    earth_param_set::PSE
end

"""
     CanopyModel{FT, RM, PM, SM, PHM, A, R, PS, D} <: AbstractExpModel{FT}

The model struct for the canopy, which contains
- the canopy model domain (a point for site-level simulations, or
an extended surface (plane/spherical surface) for regional or global simulations.
- subcomponent model type for radiative transfer. This is of type
`AbstractRadiationModel` and currently only the `BeerLambertModel` is
supported.
- subcomponent model type for photosynthesis. This is of type
`AbstractPhotosynthesisModel`, and currently only the `FarquharModel`
is supported.
- subcomponent model type for stomatal conductance. This is of type
 `AbstractStomatalConductanceModel` and currently only the `MedlynModel`
is supported
- subcomponent model type for plant hydraulics. This is of type
 `AbstractPlantHydraulicsModel` and currently only a version which
prognostically solves Richards equation in the plant is available.
- canopy model parameters, which include parameters that are shared
between canopy model components or those needed to compute boundary
fluxes.
- The atmospheric conditions, which are either prescribed
(of type `PrescribedAtmosphere`) or computed via a coupled simulation
(of type `CoupledAtmosphere`).
- The radiative flux conditions, which are either prescribed
(of type `PrescribedRadiativeFluxes`) or computed via a coupled simulation
(of type `CoupledRadiativeFluxes`).

$(DocStringExtensions.FIELDS)
"""
struct CanopyModel{FT, RM, PM, SM, PHM, A, R, PS, D} <: AbstractExpModel{FT}
    "Radiative transfer model, a canopy component model"
    radiative_transfer::RM
    "Photosynthesis model, a canopy component model"
    photosynthesis::PM
    "Stomatal conductance model, a canopy component model"
    conductance::SM
    "Plant hydraulics model, a canopy component model"
    hydraulics::PHM
    "Atmospheric forcing: prescribed or coupled"
    atmos::A
    "Radiative forcing: prescribed or coupled"
    radiation::R
    "Shared canopy parameters between component models"
    parameters::PS
    "Canopy model domain"
    domain::D
end

"""
    CanopyModel{FT}(;
        radiative_transfer::AbstractRadiationModel{FT},
        photosynthesis::AbstractPhotosynthesisModel{FT},
        conductance::AbstractStomatalConductanceModel{FT},
        hydraulics::AbstractPlantHydraulicsModel{FT},
        atmos::AbstractAtmosphericDrivers{FT},
        radiation::AbstractRadiativeDrivers{FT},
        parameters::SharedCanopyParameters{FT, PSE},
        domain::Union{
            ClimaLSM.Domains.Point,
            ClimaLSM.Domains.Plane,
            ClimaLSM.Domains.SphericalSurface,
        },
    ) where {FT, PSE}

An outer constructor for the `CanopyModel`. The primary
constraints this applies are (1) ensuring that the domain is 1d or 2d
(a ``surface" domain of a column, box, or sphere) and (2) ensuring
consistency between the PlantHydraulics model and the general canopy model,
since these are built separately.
"""
function CanopyModel{FT}(;
    radiative_transfer::AbstractRadiationModel{FT},
    photosynthesis::AbstractPhotosynthesisModel{FT},
    conductance::AbstractStomatalConductanceModel{FT},
    hydraulics::AbstractPlantHydraulicsModel{FT},
    atmos::AbstractAtmosphericDrivers{FT},
    radiation::AbstractRadiativeDrivers{FT},
    parameters::SharedCanopyParameters{FT, PSE},
    domain::Union{
        ClimaLSM.Domains.Point,
        ClimaLSM.Domains.Plane,
        ClimaLSM.Domains.SphericalSurface,
    },
) where {FT, PSE}
    if parameters.LAI != hydraulics.parameters.area_index[:leaf]
        throw(
            AssertionError(
                "The leaf area index must be the same between the plant hydraulics and shared canopy parameters.",
            ),
        )
    end
    args = (
        radiative_transfer,
        photosynthesis,
        conductance,
        hydraulics,
        atmos,
        radiation,
        parameters,
        domain,
    )
    return CanopyModel{FT, typeof.(args)...}(args...)
end

ClimaLSM.name(::CanopyModel) = :canopy
ClimaLSM.domain(::CanopyModel) = :surface

"""
    canopy_components(::CanopyModel)

Returns the names of the components of the CanopyModel.

These names are used for storing prognostic and auxiliary variables
in a hierarchical manner within the state vectors.

These names must match the field names of the CanopyModel struct.
"""
canopy_components(::CanopyModel) =
    (:hydraulics, :conductance, :photosynthesis, :radiative_transfer)

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
    initialize_prognostic(
        model::CanopyModel{FT},
        coords::ClimaCore.Fields.Field,
    ) where {FT}

Creates the prognostic state vector of the `CanopyModel` and returns
it as a ClimaCore.Fields.FieldVector.

This function loops over the components of the `CanopyModel` and appends
each component models prognostic state vector into a single state vector,
structured by component name.
"""
function initialize_prognostic(
    model::CanopyModel{FT},
    coords::ClimaCore.Fields.Field,
) where {FT}
    components = canopy_components(model)
    Y_state_list = map(components) do (component)
        submodel = getproperty(model, component)
        zero_state = map(_ -> zero(FT), coords)
        getproperty(initialize_prognostic(submodel, zero_state), component)
    end
    Y = ClimaCore.Fields.FieldVector(;
        name(model) => NamedTuple{components}(Y_state_list),
    )
    return Y
end

"""
    initialize_auxiliary(
        model::CanopyModel{FT},
        coords::ClimaCore.Fields.Field,
    ) where {FT}

Creates the auxiliary state vector of the `CanopyModel` and returns
 it as a ClimaCore.Fields.FieldVector.

This function loops over the components of the `CanopyModel` and appends
each component models auxiliary state vector into a single state vector,
structured by component name.
"""
function initialize_auxiliary(
    model::CanopyModel{FT},
    coords::ClimaCore.Fields.Field,
) where {FT}
    components = canopy_components(model)
    p_state_list = map(components) do (component)
        submodel = getproperty(model, component)
        zero_state = map(_ -> zero(FT), coords)
        getproperty(initialize_auxiliary(submodel, zero_state), component)
    end
    p = ClimaCore.Fields.FieldVector(;
        name(model) => NamedTuple{components}(p_state_list),
    )
    return p
end

"""
     ClimaLSM.make_update_aux(canopy::CanopyModel{FT, <:BeerLambertModel,
                                                  <:FarquharModel,
                                                  <:MedlynConductanceModel,
                                                  <:PlantHydraulicsModel,},
                              ) where {FT}

Creates the `update_aux!` function for the `CanopyModel`; a specific
method for `update_aux!` for the case where the canopy model components
are of the type in the parametric type signature: `BeerLambertModel`,
`FarquharModel`, `MedlynConductanceModel`, and `PlantHydraulicsModel`.

Please note that the plant hydraulics model has auxiliary variables
that are updated in its prognostic `compute_exp_tendency!` function.
While confusing, this is better for performance as it saves looping
over the state vector multiple times.

The other sub-components rely heavily on each other,
so the version of the `CanopyModel` with these subcomponents
has a single update_aux! function, given here.
"""
function ClimaLSM.make_update_aux(
    canopy::CanopyModel{
        FT,
        <:BeerLambertModel,
        <:FarquharModel,
        <:MedlynConductanceModel,
        <:PlantHydraulicsModel,
    },
) where {FT}
    function update_aux!(p, Y, t)
        # variables being updated
        APAR = p.canopy.radiative_transfer.apar
        PAR = p.canopy.radiative_transfer.par
        β = p.canopy.hydraulics.β
        medlyn_factor = p.canopy.conductance.medlyn_term
        gs = p.canopy.conductance.gs
        transpiration = p.canopy.conductance.transpiration
        An = p.canopy.photosynthesis.An
        GPP = p.canopy.photosynthesis.GPP
        ψ = p.canopy.hydraulics.ψ
        ϑ_l = Y.canopy.hydraulics.ϑ_l
        fa = p.canopy.hydraulics.fa

        #unpack parameters
        earth_param_set = canopy.parameters.earth_param_set
        c = FT(LSMP.light_speed(earth_param_set))
        h = FT(LSMP.planck_constant(earth_param_set))
        N_a = FT(LSMP.avogadro_constant(earth_param_set))
        (; ld, Ω, ρ_leaf, λ_γ) = canopy.radiative_transfer.parameters
        (; LAI) = canopy.parameters
        energy_per_photon = h * c / λ_γ
        R = FT(LSMP.gas_constant(earth_param_set))
        thermo_params = canopy.parameters.earth_param_set.thermo_params
        (; g1, g0, Drel) = canopy.conductance.parameters
        (; sc, ψc) = canopy.photosynthesis.parameters

        # Current atmospheric conditions
        θs::FT = canopy.radiation.θs(t, canopy.radiation.orbital_data)
        c_co2::FT = canopy.atmos.c_co2(t)
        P::FT = canopy.atmos.P(t)
        u::FT = canopy.atmos.u(t)
        T::FT = canopy.atmos.T(t)
        h::FT = canopy.atmos.h
        q::FT = canopy.atmos.q(t)

        # update radiative transfer
        PAR .= compute_PAR(canopy.radiative_transfer, canopy.radiation, t)
        K = extinction_coeff(ld, θs)
        @. APAR = plant_absorbed_ppfd(
            PAR / (energy_per_photon * N_a),
            ρ_leaf,
            K,
            LAI,
            Ω,
        )

        # update plant hydraulics aux
        hydraulics = canopy.hydraulics
        n_stem = hydraulics.n_stem
        n_leaf = hydraulics.n_leaf
        (; vg_α, vg_n, vg_m, S_s, ν, K_sat, area_index) = hydraulics.parameters
        @inbounds @. ψ[1] = PlantHydraulics.water_retention_curve(
            vg_α,
            vg_n,
            vg_m,
            PlantHydraulics.effective_saturation(ν, ϑ_l[1]),
            ν,
            S_s,
        )

        @inbounds for i in 1:(n_stem + n_leaf - 1)

            @. ψ[i + 1] = PlantHydraulics.water_retention_curve(
                vg_α,
                vg_n,
                vg_m,
                PlantHydraulics.effective_saturation(ν, ϑ_l[i + 1]),
                ν,
                S_s,
            )
            # Compute the flux*area between the current compartment `i`
            # and the compartment above.
            @. fa[i] =
                PlantHydraulics.flux(
                    hydraulics.compartment_midpoints[i],
                    hydraulics.compartment_midpoints[i + 1],
                    ψ[i],
                    ψ[i + 1],
                    vg_α,
                    vg_n,
                    vg_m,
                    ν,
                    S_s,
                    K_sat[hydraulics.compartment_labels[i]],
                    K_sat[hydraulics.compartment_labels[i + 1]],
                ) * (
                    area_index[hydraulics.compartment_labels[i]] +
                    area_index[hydraulics.compartment_labels[i + 1]]
                ) / 2
        end
        @. β = moisture_stress(ψ[n_stem + n_leaf], sc, ψc)
        # We update the fa[n_stem+n_leaf] element once we have computed transpiration, below

        # update photosynthesis and conductance terms
        medlyn_factor .= medlyn_term.(g1, T, P, q, Ref(thermo_params))
        An .=
            compute_photosynthesis.(
                Ref(canopy.photosynthesis),
                T,
                medlyn_factor,
                APAR,
                c_co2,
                β,
                R,
            )
        @. GPP = compute_GPP(An, K, LAI, Ω)
        @. gs = medlyn_conductance(g0, Drel, medlyn_factor, An, c_co2)

        # Compute transpiration
        (canopy_transpiration, shf, lhf) =
            canopy_surface_fluxes(canopy.atmos, canopy, Y, p, t)
        transpiration .= canopy_transpiration
        # Transpiration is per unit ground area, not leaf area.
        fa[n_stem + n_leaf] .= PlantHydraulics.transpiration_per_ground_area(
            hydraulics.transpiration,
            Y,
            p,
            t,
        )

    end
    return update_aux!
end

"""
    make_compute_exp_tendency(canopy::CanopyModel)

Creates and returns the compute_exp_tendency! for the `CanopyModel`.

This allows for prognostic variables in each canopy component, and
specifies that they will be stepped explicitly.
"""
function make_compute_exp_tendency(canopy::CanopyModel)
    components = canopy_components(canopy)
    compute_exp_tendency_list = map(
        x -> make_compute_exp_tendency(getproperty(canopy, x), canopy),
        components,
    )
    function compute_exp_tendency!(dY, Y, p, t)
        # aux vars are updated in `compute_exp_tendency` functions
        for f! in compute_exp_tendency_list
            f!(dY, Y, p, t)
        end
    end
    return compute_exp_tendency!
end

"""
    canopy_surface_fluxes(atmos::PrescribedAtmosphere{FT},
                          model::CanopyModel,
                          Y,
                          p,
                          t::FT) where {FT}

Computes canopy transpiration using Monin-Obukhov Surface Theory,
the prescribed atmospheric conditions, and the canopy conductance.

Please note that in the future the SurfaceFluxes.jl code will compute
fluxes taking into account the canopy conductance, so that
what is returned by `surface_fluxes` is correct. At present, it does not,
so we are adjusting for it after the fact here in both ET and LHF.
"""
function canopy_surface_fluxes(
    atmos::PrescribedAtmosphere{FT},
    model::CanopyModel,
    Y,
    p,
    t::FT,
) where {FT}
    conditions = surface_fluxes(atmos, model, Y, p, t) # per unit m^2 of leaf
    earth_param_set = model.parameters.earth_param_set
    R = FT(LSMP.gas_constant(earth_param_set))
    ρ_liq = FT(LSMP.ρ_cloud_liq(earth_param_set))
    P::FT = model.atmos.P(t)
    T::FT = model.atmos.T(t)

    # here is where we adjust evaporation for stomatal conductance = 1/r_sfc
    r_ae = 1 / (conditions.Ch * abs(atmos.u(t))) # [s/m]

    leaf_conductance = p.canopy.conductance.gs
    canopy_conductance =
        upscale_leaf_conductance.(
            leaf_conductance,
            model.parameters.LAI,
            T,
            R,
            P,
        )
    r_sfc = @. 1 / (canopy_conductance) # [s/m]
    r_eff = r_ae .+ r_sfc
    canopy_transpiration = @. conditions.vapor_flux * r_ae / r_eff

    # we also need to correct the LHF
    canopy_lhf = @. conditions.lhf * r_ae / r_eff
    return canopy_transpiration, conditions.shf, canopy_lhf
end


"""
    ClimaLSM.surface_temperature(model::CanopyModel, Y, p, t)

A helper function which returns the surface temperature for the canopy
model, which is stored in the aux state.
"""
function ClimaLSM.surface_temperature(model::CanopyModel, Y, p, t)
    return model.atmos.T(t)
end

"""
    ClimaLSM.surface_height(model::CanopyModel, Y, _...)

A helper function which returns the surface height for the canopy
model, which is stored in the parameter struct.
"""
function ClimaLSM.surface_height(model::CanopyModel, _...)
    return model.parameters.h_c
end

"""
    ClimaLSM.surface_specific_humidity(model::CanopyModel, Y, p)

A helper function which returns the surface specific humidity for the canopy
model, which is stored in the aux state.
"""
function ClimaLSM.surface_specific_humidity(
    model::CanopyModel,
    Y,
    p,
    T_sfc,
    ρ_sfc,
)
    thermo_params =
        LSMP.thermodynamic_parameters(model.parameters.earth_param_set)
    return Thermodynamics.q_vap_saturation_generic.(
        Ref(thermo_params),
        T_sfc,
        ρ_sfc,
        Ref(Thermodynamics.Liquid()),
    )
end

"""
    ClimaLSM.surface_air_density(model::CanopyModel, Y, p)

A helper function which computes and returns the surface air density for the canopy
model.
"""
function ClimaLSM.surface_air_density(
    atmos::PrescribedAtmosphere,
    model::CanopyModel,
    Y,
    p,
    t,
    T_sfc,
)
    thermo_params =
        LSMP.thermodynamic_parameters(model.parameters.earth_param_set)
    ts_in = construct_atmos_ts(atmos, t, thermo_params)
    return compute_ρ_sfc.(Ref(thermo_params), Ref(ts_in), T_sfc)
end

"""
    ClimaLSM.surface_evaporative_scaling(model::CanopyModel, Y, p)

A helper function which computes and returns the surface evaporative scaling
 factor for the canopy model.
"""
function ClimaLSM.surface_evaporative_scaling(
    model::CanopyModel{FT},
    Y,
    p,
) where {FT}
    beta = FT(1.0)
    return beta
end

end
