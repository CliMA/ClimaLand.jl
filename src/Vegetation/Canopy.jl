module Canopy
using Thermodynamics
using ClimaLSM
using ClimaCore
using ClimaLSM:
    AbstractModel, AbstractRadiativeDrivers, AbstractAtmosphericDrivers
import ..Parameters as LSMP

import ClimaLSM:
    name,
    domain,
    prognostic_vars,
    prognostic_types,
    auxiliary_vars,
    auxiliary_types,
    initialize_prognostic,
    initialize_auxiliary,
    make_update_aux,
    make_rhs,
    surface_temperature,
    surface_specific_humidity,
    surface_air_density,
    surface_evaporative_scaling,
    surface_height

using ClimaLSM.Domains: Point, Plane, SphericalSurface
using DocStringExtensions
export SharedCanopyParameters, CanopyModel, DiagnosticTranspiration
include("./PlantHydraulics.jl")
using .PlantHydraulics
import .PlantHydraulics: transpiration, AbstractTranspiration
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
     CanopyModel{FT, RM, PM, SM, PHM, A, R, PS, D} <: AbstractModel{FT}

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
struct CanopyModel{FT, RM, PM, SM, PHM, A, R, PS, D} <: AbstractModel{FT}
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
    if hydraulics.domain != domain
        throw(
            AssertionError(
                "The provided canopy model domain must be the same as the hydraulics model domain.",
            ),
        )
    end
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

"""
    CanopyModel{FT}(hydraulics::PlantHydraulicsModel{FT}) where {FT}

An outer constructor for the CanopyModel which reduces to the standalone
PlantHydraulicsModel.

For testing purposes.
"""
function CanopyModel{FT}(hydraulics::PlantHydraulicsModel{FT}) where {FT}
    args = (
        nothing,
        nothing,
        nothing,
        hydraulics,
        nothing,
        nothing,
        nothing,
        hydraulics.domain,
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
        if submodel isa Nothing
            FT[]
        else
            zero_state = map(_ -> zero(FT), coords)
            getproperty(initialize_prognostic(submodel, zero_state), component)
        end
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
        if submodel isa Nothing
            FT[]
        else
            zero_state = map(_ -> zero(FT), coords)
            getproperty(initialize_auxiliary(submodel, zero_state), component)
        end
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
that are updated in its prognostic `rhs!` function. While confusing, this
is better for performance as it saves looping over the state vector
multiple times. As such, there is no plant hydraulics update_aux! function
ever defined or used here.

The other sub-components rely heavily on each other and will never be run in
standalone fashion, so they do not have separate `update_aux!` functions,
but rather, the version of the `CanopyModel` with these subcomponents 
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
        earth_param_set = canopy.parameters.earth_param_set
        R = FT(LSMP.gas_constant(earth_param_set))
        c = FT(LSMP.light_speed(earth_param_set))
        h = FT(LSMP.planck_constant(earth_param_set))
        N_a = FT(LSMP.avogadro_constant(earth_param_set))
        thermo_params = canopy.parameters.earth_param_set.thermo_params

        (;
            Vcmax25,
            Γstar25,
            ΔHJmax,
            ΔHVcmax,
            ΔHΓstar,
            f,
            ΔHRd,
            To,
            θj,
            ϕ,
            mechanism,
            sc,
            ψc,
            oi,
            Kc25,
            Ko25,
            ΔHkc,
            ΔHko,
        ) = canopy.photosynthesis.parameters
        (; g1, g0, Drel) = canopy.conductance.parameters
        (; ld, Ω, ρ_leaf, λ_γ) = canopy.radiative_transfer.parameters
        (; LAI) = canopy.parameters
        energy_per_photon = h * c / λ_γ

        top_index = canopy.hydraulics.n_stem + canopy.hydraulics.n_leaf
        c_co2::FT = canopy.atmos.c_co2(t)
        P::FT = canopy.atmos.P(t)
        u::FT = canopy.atmos.u(t)
        T::FT = canopy.atmos.T(t)
        h::FT = canopy.atmos.h
        q::FT = canopy.atmos.q(t)
        SW_d::FT = canopy.radiation.SW_d(t)
        LW_d::FT = canopy.radiation.LW_d(t)
        θs::FT = canopy.radiation.θs(t)

        # compute VPD
        es =
            Thermodynamics.saturation_vapor_pressure.(
                Ref(thermo_params),
                T,
                Ref(Thermodynamics.Liquid()),
            )
        ea =
            Thermodynamics.partial_pressure_vapor.(
                Ref(thermo_params),
                P,
                Thermodynamics.PhasePartition.(q),
            )
        VPD = es .- ea
        PAR = @.(SW_d / (energy_per_photon * N_a) / 2)

        K = extinction_coeff.(ld, θs)
        APAR = plant_absorbed_ppfd.(PAR, ρ_leaf, K, LAI, Ω)
        β = moisture_stress.(p.canopy.hydraulics.ψ[top_index], sc, ψc)
        Jmax = max_electron_transport.(Vcmax25, ΔHJmax, T, To, R)
        J = electron_transport.(APAR, Jmax, θj, ϕ)
        Vcmax = compute_Vcmax.(Vcmax25, T, To, R, ΔHVcmax)
        Γstar = co2_compensation.(Γstar25, ΔHΓstar, T, To, R)
        m = medlyn_term.(g1, VPD)
        ci = intercellular_co2.(c_co2, Γstar, m)
        Aj = light_assimilation.(Ref(mechanism), J, ci, Γstar)
        Kc = MM_Kc.(Kc25, ΔHkc, T, To, R)
        Ko = MM_Ko.(Ko25, ΔHko, T, To, R)
        Ac = rubisco_assimilation.(Ref(mechanism), Vcmax, ci, Γstar, Kc, Ko, oi)
        Rd = dark_respiration.(Vcmax25, β, f, ΔHRd, T, To, R)

        p.canopy.radiative_transfer.apar .= APAR
        p.canopy.photosynthesis.An .= net_photosynthesis.(Ac, Aj, Rd, β)
        p.canopy.photosynthesis.GPP .=
            compute_GPP.(p.canopy.photosynthesis.An, K, LAI, Ω)
        p.canopy.conductance.gs .=
            medlyn_conductance.(g0, Drel, m, p.canopy.photosynthesis.An, c_co2)
        p.canopy.conductance.medlyn_term .= m
        (evapotranspiration, shf, lhf) =
            canopy_surface_fluxes(canopy.atmos, canopy, Y, p, t)
        p.canopy.hydraulics.fa[top_index] .= evapotranspiration ./ LAI # this should be per leaf area
        # note, somewhat confusingly, that the other aux variables for plant hydraulics are updated in the RHS
    end
    return update_aux!
end

"""
    make_rhs(canopy::CanopyModel)

Creates and returns the rhs! for the `CanopyModel`. 

This allows for prognostic variables in each canopy component.
"""
function make_rhs(canopy::CanopyModel)
    components = canopy_components(canopy)
    rhs_function_list = map(components) do (component)
        submodel = getproperty(canopy, component)
        if submodel isa Nothing
            (dY, Y, p, t) -> nothing
        else
            make_rhs(submodel)
        end
    end
    function rhs!(dY, Y, p, t)
        for f! in rhs_function_list
            f!(dY, Y, p, t)
        end
    end
    return rhs!
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



"""
    DiagnosticTranspiration{FT} <: AbstractTranspiration{FT}

A concrete type used for dispatch when computing the transpiration
from the leaves, in the case where transpiration is computed
diagnostically, as a function of prognostic variables and parameters,
and stored in `p` during the `update_aux!` step.
"""
struct DiagnosticTranspiration{FT} <: AbstractTranspiration{FT} end

"""
    transpiration(transpiration::DiagnosticTranspiration, Y, p, t)

Returns the transpiration computed diagnostically using local conditions.
In this case, it just returns the value which was computed and stored in 
the `aux` state during the update_aux! step.
"""
function transpiration(transpiration::DiagnosticTranspiration, Y, p, t)

    @inbounds return p.canopy.hydraulics.fa[length(
        propertynames(p.canopy.hydraulics.fa),
    )]
end


end
