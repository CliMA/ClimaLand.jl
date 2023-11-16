module Snow

using UnPack
using DocStringExtensions
import ...Parameters as LSMP
using ClimaCore
using Thermodynamics
using ClimaLSM
using ClimaLSM:
    AbstractAtmosphericDrivers,
    AbstractRadiativeDrivers,
    liquid_precipitation,
    snow_precipitation,
    surface_fluxes,
    net_radiation,
    construct_atmos_ts,
    compute_ρ_sfc,
    AbstractModel,
    heaviside,
    PrescribedAtmosphere

import ..Parameters as LSMP
import ClimaLSM.Domains
import ClimaLSM:
    make_update_aux,
    make_compute_exp_tendency,
    prognostic_vars,
    auxiliary_vars,
    name,
    prognostic_types,
    prognostic_domain_names,
    auxiliary_types,
    auxiliary_domain_names,
    initialize_vars,
    initialize,
    initialize_auxiliary,
    make_set_initial_aux_state,
    surface_temperature,
    surface_air_density,
    surface_specific_humidity,
    surface_evaporative_scaling,
    surface_height,
    surface_albedo,
    surface_emissivity
export SnowParameters, SnowModel

"""
    SnowParameters{FT <: AbstractFloat, PSE}

A struct for storing parameters of the `SnowModel`.
$(DocStringExtensions.FIELDS)
"""
struct SnowParameters{FT <: AbstractFloat, PSE}
    "Density of snow (kg/m^3)"
    ρ_snow::FT
    "Roughness length over snow for momentum (m)"
    z_0m::FT
    "Roughness length over snow for scalars (m)"
    z_0b::FT
    "Albedo of snow (unitless)"
    α_snow::FT
    "Emissivity of snow (unitless)"
    ϵ_snow::FT
    "Volumetric holding capacity of water in snow (unitless)"
    θ_r::FT
    "Hydraulic conductivity of wet snow (m/s)"
    Ksat::FT
    "Thermal conductivity of ice (W/m/K)"
    κ_ice::FT
    "Timestep of the model"
    Δt::FT
    "Critical depth (m) under which the snowpack temperature is the ground surface temperature"
    z_crit::FT
    "Clima-wide parameters"
    earth_param_set::PSE
end

"""
   SnowParameters{FT}(Δt;
                      ρ_snow = FT(200),
                      z_0m = FT(0.0024),
                      z_0b = FT(0.00024),
                      α_snow = FT(0.8),
                      ϵ_snow = FT(0.99),
                      θ_r = FT(0.08),
                      Ksat = FT(1e-3),
                      κ_ice = FT(2.21),
                      z_crit = FT(0.0508),
                      earth_param_set::PSE) where {FT, PSE}

An outer constructor for `SnowParameters` which supplies defaults for
all arguments but `earth_param_set`.
"""
function SnowParameters{FT}(
    Δt;
    ρ_snow = FT(200),
    z_0m = FT(0.0024),
    z_0b = FT(0.00024),
    α_snow = FT(0.8),
    ϵ_snow = FT(0.99),
    θ_r = FT(0.08),
    Ksat = FT(1e-3),
    κ_ice = FT(2.21),
    z_crit = FT(0.0508),
    earth_param_set::PSE,
) where {FT <: AbstractFloat, PSE}
    return SnowParameters{FT, PSE}(
        ρ_snow,
        z_0m,
        z_0b,
        α_snow,
        ϵ_snow,
        θ_r,
        Ksat,
        κ_ice,
        Δt,
        z_crit,
        earth_param_set,
    )
end

struct SnowModel{
    FT,
    PS <: SnowParameters{FT},
    ATM <: AbstractAtmosphericDrivers{FT},
    RAD <: AbstractRadiativeDrivers{FT},
    D,
} <: AbstractModel{FT}
    "Parameters required by the snow model"
    parameters::PS
    "The atmospheric drivers: Prescribed or Coupled"
    atmos::ATM
    "The radiation drivers: Prescribed or Coupled"
    radiation::RAD
    "The domain of the model"
    domain::D
end
function SnowModel(;
    parameters::SnowParameters{FT, PSE},
    domain::ClimaLSM.Domains.AbstractDomain,
    atmosphere::ATM,
    radiation::RAD,
) where {FT, PSE, ATM, RAD}
    args = (parameters, atmosphere, radiation, domain)
    SnowModel{FT, typeof.(args)...}(args...)
end


prognostic_types(::SnowModel{FT}) where {FT} = (FT, FT)
prognostic_vars(::SnowModel) = (:S, :U)

prognostic_domain_names(::SnowModel) = (:surface, :surface)

auxiliary_types(::SnowModel{FT}) where {FT} =
    (FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT)
auxiliary_domain_names(::SnowModel) = 
    (:surface,:surface,:surface,:surface,:surface,:surface,:surface,:surface,:surface,:surface,:surface,:surface,:surface,:surface,)
auxiliary_vars(::SnowModel) = (
    :q_l,
    :T,
    :T_sfc,
    :ρ_sfc,
    :q_sfc,
    :evaporation,
    :turbulent_energy_flux,
    :R_n,
    :energy_runoff,
    :water_runoff,
    :total_energy_flux,
    :total_water_flux,
    :snow_energy_flux,
    :snow_water_flux,
)


ClimaLSM.name(::SnowModel) = :snow

function ClimaLSM.make_update_aux(model::SnowModel{FT}) where {FT}
    function update_aux!(p, Y, t)
        parameters = model.parameters

        _T_freeze = FT(LSMP.T_freeze(parameters.earth_param_set))
        _LH_f0 = FT(LSMP.LH_f0(parameters.earth_param_set))
        _ρ_liq = FT(LSMP.ρ_cloud_liq(parameters.earth_param_set))

        p.snow.q_l .=
            snow_liquid_mass_fraction.(Y.snow.U, Y.snow.S, Ref(parameters))
        p.snow.T .=
            snow_bulk_temperature.(
                Y.snow.U,
                Y.snow.S,
                p.snow.q_l,
                Ref(parameters),
            )
        p.snow.T_sfc .= snow_surface_temperature.(p.snow.T)
        p.snow.ρ_sfc .=
            surface_air_density(model.atmos, model, Y, p, t, p.snow.T_sfc)        # We are currently computing evaporation as if it is entirely over ice.
        # In the future, it'll be a mass weighted average over ice and liquid water.
        # This requires two calls to surface fluxes.
        thermo_params =
            LSMP.thermodynamic_parameters(parameters.earth_param_set)

        p.snow.q_sfc .=
            Thermodynamics.q_vap_saturation_generic.(
                Ref(thermo_params),
                p.snow.T_sfc,
                p.snow.ρ_sfc,
                Ref(Thermodynamics.Ice()),
            )

        turbulent_fluxes = surface_fluxes(model.atmos, model, Y, p, t)


        @. p.snow.turbulent_energy_flux =
            turbulent_fluxes.lhf + turbulent_fluxes.shf
        @. p.snow.evaporation = turbulent_fluxes.vapor_flux
        p.snow.R_n .= net_radiation(model.radiation, model, Y, p, t)


        p.snow.water_runoff .=
            compute_water_runoff.(
                Y.snow.S,
                p.snow.q_l,
                p.snow.T,
                Ref(parameters),
            )

        I_liq =
            volumetric_internal_energy_liq(eltype(Y.snow.S), model.parameters)
        p.snow.energy_runoff .= p.snow.water_runoff .* I_liq

        P_liq = liquid_precipitation(model.atmos, p, t) # volume flux of liquid water
        P_snow = snow_precipitation(model.atmos, p, t) # volume flux of liquid water

        # Heaviside function because we cannot change the water content of the snow by
        # sublimation, evap, melting, or rain
        # unless there is already snow on the ground. 
        @. p.snow.total_water_flux =
            -P_snow +
            (- p.snow.evaporation + p.snow.water_runoff) *
            heaviside(Y.snow.S - eps(FT))
        @. p.snow.snow_water_flux =
            clip_dSdt(Y.snow.S, p.snow.total_water_flux, model.parameters.Δt)

        # I think we want dU/dt to include energy of falling snow.
        # otherwise snow can fall but energy wont change
        ρe_falling_snow = -_LH_f0 * _ρ_liq # per unit vol of liquid water
        # We are assuming that the internal energy of rain and sensible heat portion of ice is negligible.

        @. p.snow.total_energy_flux =
            -P_snow * ρe_falling_snow +
            (
                -p.snow.turbulent_energy_flux - p.snow.R_n +
                p.snow.energy_runoff
            ) * heaviside(Y.snow.S - eps(FT))
        @. p.snow.snow_energy_flux = clip_dUdt(
            Y.snow.U,
            Y.snow.S,
            p.snow.total_energy_flux,
            p.snow.total_water_flux,
            model.parameters.Δt,
        )

    end
end

function ClimaLSM.make_compute_exp_tendency(model::SnowModel)
    function compute_exp_tendency!(dY, Y, p, t)
        @. dY.snow.S = p.snow.snow_water_flux
        @. dY.snow.U = p.snow.snow_energy_flux
    end
    return compute_exp_tendency!
end

function clip_dSdt(S, dSdt, Δt)
    if S + dSdt * Δt < 0
        return -S / Δt
    else
        return dSdt
    end
end

function clip_dUdt(U, S, dUdt, dSdt, Δt)
    if (U + dUdt * Δt) > 0
        return -U / Δt
    elseif S + dSdt * Δt < 0
        return -U / Δt
    else
        return dUdt
    end
end

include("./snow_parameterizations.jl")
end
