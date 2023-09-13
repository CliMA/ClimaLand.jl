using Thermodynamics
using ClimaCore
using DocStringExtensions
using SurfaceFluxes
import SurfaceFluxes.Parameters as SFP
using StaticArrays
import ..Parameters as LSMP
export AbstractAtmosphericDrivers,
    AbstractRadiativeDrivers,
    PrescribedAtmosphere,
    PrescribedRadiativeFluxes,
    compute_ρ_sfc,
    construct_atmos_ts,
    surface_fluxes,
    net_radiation,
    surface_fluxes_at_a_point,
    liquid_precipitation,
    snow_precipitation,
    vapor_pressure_deficit,
    displacement_height

"""
     AbstractAtmosphericDrivers{FT <: AbstractFloat}

An abstract type of atmospheric drivers of land models.
"""
abstract type AbstractAtmosphericDrivers{FT <: AbstractFloat} end

"""
     AbstractRadiativeDrivers{FT <: AbstractFloat}

An abstract type of radiative drivers of land models.
"""
abstract type AbstractRadiativeDrivers{FT <: AbstractFloat} end

"""
    PrescribedAtmosphere{FT, LP, SP, TA, UA, QA, RA, CA, DT} <: AbstractAtmosphericDrivers{FT}

Container for holding prescribed atmospheric drivers and other
information needed for computing turbulent surface fluxes when
driving land models in standalone mode.

Since not all models require co2 concentration, the default for that
is `nothing`.
$(DocStringExtensions.FIELDS)
"""
struct PrescribedAtmosphere{FT, LP, SP, TA, UA, QA, RA, CA, DT} <:
       AbstractAtmosphericDrivers{FT}
    "Precipitation (m/s) function of time: positive by definition"
    liquid_precip::LP
    "Snow precipitation (m/s) function of time: positive by definition"
    snow_precip::SP
    "Prescribed atmospheric temperature (function of time)  at the reference height (K)"
    T::TA
    "Prescribed wind speed (function of time)  at the reference height (m/s)"
    u::UA
    "Prescribed specific humidity (function of time)  at the reference height (_)"
    q::QA
    "Prescribed air pressure (function of time)  at the reference height (Pa)"
    P::RA
    "CO2 concentration in atmosphere (mol/mol)"
    c_co2::CA
    "Reference time - the datetime corresponding to t=0 for the simulation"
    ref_time::DT
    "Reference height (m), relative to surface elevation"
    h::FT
    "Minimum wind speed (gustiness; m/s)"
    gustiness::FT
    function PrescribedAtmosphere(
        liquid_precip,
        snow_precip,
        T,
        u,
        q,
        P,
        ref_time,
        h::FT;
        gustiness = FT(1),
        c_co2 = nothing,
    ) where {FT}
        args = (liquid_precip, snow_precip, T, u, q, P, c_co2, ref_time)
        return new{typeof(h), typeof.(args)...}(args..., h, gustiness)
    end

end

"""
    compute_ρ_sfc(thermo_params, ts_in, T_sfc)

Computes the density of air at the surface, given the temperature
at the surface T_sfc, the thermodynamic state of the atmosphere,
ts_in, and a set of Clima.Thermodynamics parameters thermo_params.

This assumes the ideal gas law and hydrostatic balance to
extrapolate to the surface.

"""
function compute_ρ_sfc(thermo_params, ts_in, T_sfc)
    T_int = Thermodynamics.air_temperature(thermo_params, ts_in)
    Rm_int = Thermodynamics.gas_constant_air(thermo_params, ts_in)
    ρ_air = Thermodynamics.air_density(thermo_params, ts_in)
    ρ_sfc =
        ρ_air *
        (T_sfc / T_int)^(Thermodynamics.cv_m(thermo_params, ts_in) / Rm_int)
    return ρ_sfc
end

"""
    construct_atmos_ts(
        atmos::PrescribedAtmosphere,
        t::FT,
        thermo_params,
    ) where {FT}

A helper function which constructs a Clima.Thermodynamics
thermodynamic state given a PrescribedAtmosphere, a time
at which the state is needed, and a set of Clima.Thermodynamics
parameters thermo_params.
"""
function construct_atmos_ts(
    atmos::PrescribedAtmosphere,
    t::FT,
    thermo_params,
) where {FT}
    P::FT = atmos.P(t)
    T::FT = atmos.T(t)
    q::FT = atmos.q(t)
    ts_in = Thermodynamics.PhaseEquil_pTq(thermo_params, P, T, q)
    return ts_in
end


"""
    surface_fluxes(atmos::PrescribedAtmosphere{FT},
                   model::AbstractModel{FT},
                   Y::ClimaCore.Fields.FieldVector,
                   p::NamedTuple,
                   t::FT
                   ) where {FT <: AbstractFloat}

Computes the turbulent surface flux terms at the ground for a standalone simulation,
including turbulent energy fluxes as well as the water vapor flux
(in units of m^3/m^2/s of water).
Positive fluxes indicate flow from the ground to the atmosphere.

It solves for these given atmospheric conditions, stored in `atmos`,
model parameters, and the surface conditions.
"""
function surface_fluxes(
    atmos::PrescribedAtmosphere{FT},
    model::AbstractModel{FT},
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    t::FT,
) where {FT <: AbstractFloat}
    T_sfc = surface_temperature(model, Y, p, t)
    ρ_sfc = surface_air_density(atmos, model, Y, p, t, T_sfc)
    q_sfc = surface_specific_humidity(model, Y, p, T_sfc, ρ_sfc)
    β_sfc = surface_evaporative_scaling(model, Y, p)
    h_sfc = surface_height(model, Y, p)
    r_sfc = surface_resistance(model, Y, p, t)
    d_sfc = displacement_height(model, Y, p)
    return surface_fluxes_at_a_point.(
        T_sfc,
        q_sfc,
        ρ_sfc,
        β_sfc,
        h_sfc,
        r_sfc,
        d_sfc,
        t,
        Ref(model.parameters),
        Ref(atmos),
    )
end

"""
    surface_fluxes_at_a_point(T_sfc::FT,
                              q_sfc::FT,
                              ρ_sfc::FT,
                              β_sfc::FT,
                              h_sfc::FT,
                              r_sfc::FT,
                              d_sfc::FT,
                              t::FT,
                              parameters,
                              atmos::PA,
                              ) where {FT <: AbstractFloat, PA <: PrescribedAtmosphere{FT}}

Computes turbulent surface fluxes at a point on a surface given
(1) the surface temperature, specific humidity, and air density,
(2) the time at which the fluxes are needed,
(3) a factor β_sfc  which scales the evaporation from the potential rate
    (used in bucket models), and/or the surface resistance r_sfc (used
    in more complex land models),
(4) the parameter set for the model, which must have fields `earth_param_set`,
and roughness lengths `z_0m, z_0b`.
(5) the prescribed atmospheric state, stored in `atmos`,
(6) the displacement height for the model.

This returns an energy flux and a liquid water volume flux, stored in
a tuple with self explanatory keys.
"""
function surface_fluxes_at_a_point(
    T_sfc::FT,
    q_sfc::FT,
    ρ_sfc::FT,
    β_sfc::FT,
    h_sfc::FT,
    r_sfc::FT,
    d_sfc::FT,
    t::FT,
    parameters::P,
    atmos::PA,
) where {FT <: AbstractFloat, P, PA <: PrescribedAtmosphere{FT}}
    (; z_0m, z_0b, earth_param_set) = parameters
    thermo_params = LSMP.thermodynamic_parameters(earth_param_set)

    u::FT = atmos.u(t)
    h::FT = atmos.h

    ts_in = construct_atmos_ts(atmos, t, thermo_params)
    ts_sfc = Thermodynamics.PhaseEquil_ρTq(thermo_params, ρ_sfc, T_sfc, q_sfc)

    # SurfaceFluxes.jl expects a relative difference between where u = 0
    # and the atmosphere height. Here, we assume h and h_sfc are measured
    # relative to a common reference. Then d_sfc + h_sfc + z_0m is the apparent
    # source of momentum, and
    # Δh ≈ h - d_sfc - h_sfc is the relative height difference between the
    # apparent source of momentum and the atmosphere height.

    # In this we have neglected z_0m and z_0b (i.e. assumed they are small
    # compared to Δh).
    state_sfc = SurfaceFluxes.StateValues(FT(0), SVector{2, FT}(0, 0), ts_sfc)
    state_in = SurfaceFluxes.StateValues(
        h - d_sfc - h_sfc,
        SVector{2, FT}(u, 0),
        ts_in,
    )

    # State containers
    sc = SurfaceFluxes.ValuesOnly(
        state_in,
        state_sfc,
        z_0m,
        z_0b,
        beta = β_sfc,
        gustiness = atmos.gustiness,
    )
    surface_flux_params = LSMP.surface_fluxes_parameters(earth_param_set)
    conditions = SurfaceFluxes.surface_conditions(
        surface_flux_params,
        sc;
        tol_neutral = SFP.cp_d(surface_flux_params) / 100000,
    )
    _LH_v0::FT = LSMP.LH_v0(earth_param_set)
    _ρ_liq = LSMP.ρ_cloud_liq(earth_param_set)

    cp_m = Thermodynamics.cp_m(thermo_params, ts_in)
    T_in = Thermodynamics.air_temperature(thermo_params, ts_in)
    ΔT = T_in - T_sfc
    r_ae = 1 / (conditions.Ch * SurfaceFluxes.windspeed(sc))
    ρ_air = Thermodynamics.air_density(thermo_params, ts_in)

    E0 = SurfaceFluxes.evaporation(surface_flux_params, sc, conditions.Ch)
    E = E0 * r_ae / (r_sfc + r_ae)
    Ẽ = E / _ρ_liq
    H = -ρ_air * cp_m * ΔT / r_ae
    LH = _LH_v0 * E
    return (lhf = LH, shf = H, vapor_flux = Ẽ, r_ae = r_ae)
end

"""
    PrescribedRadiativeFluxes{FT, SW, LW, DT, T, OD} <: AbstractRadiativeDrivers{FT}

Container for the prescribed radiation functions needed to drive land models in standalone mode.
$(DocStringExtensions.FIELDS)
"""
struct PrescribedRadiativeFluxes{FT, SW, LW, DT, T, OD} <:
       AbstractRadiativeDrivers{FT}
    "Downward shortwave radiation function of time (W/m^2): positive indicates towards surface"
    SW_d::SW
    "Downward longwave radiation function of time (W/m^2): positive indicates towards surface"
    LW_d::LW
    "Reference time - the datetime corresponding to t=0 for the simulation"
    ref_time::DT
    "Sun zenith angle, in radians"
    θs::T
    "Orbital Data for Insolation.jl"
    orbital_data::OD
    function PrescribedRadiativeFluxes(
        FT,
        SW_d,
        LW_d,
        ref_time;
        θs = nothing,
        orbital_data,
    )
        args = (SW_d, LW_d, ref_time, θs, orbital_data)
        @assert !isnothing(orbital_data)
        return new{FT, typeof.(args)...}(args...)
    end
end

"""
    net_radiation(radiation::PrescribedRadiativeFluxes{FT},
                  model::AbstractModel{FT},
                  Y::ClimaCore.Fields.FieldVector,
                  p::NamedTuple,
                  t::FT
                  ) where {FT <: AbstractFloat}

Computes net radiative fluxes for a prescribed incoming
longwave and shortwave radiation.

This returns an energy flux.
"""
function net_radiation(
    radiation::PrescribedRadiativeFluxes{FT},
    model::AbstractModel{FT},
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    t::FT,
) where {FT <: AbstractFloat}
    LW_d::FT = radiation.LW_d(t)
    SW_d::FT = radiation.SW_d(t)
    earth_param_set = model.parameters.earth_param_set
    _σ = LSMP.Stefan(earth_param_set)
    T_sfc = surface_temperature(model, Y, p, t)
    α_sfc = surface_albedo(model, Y, p)
    ϵ_sfc = surface_emissivity(model, Y, p)
    # Recall that the user passed the LW and SW downwelling radiation,
    # where positive values indicate toward surface, so we need a negative sign out front
    # in order to inidicate positive R_n  = towards atmos.
    R_n = @.(-(1 - α_sfc) * SW_d - ϵ_sfc * (LW_d - _σ * T_sfc^4))
    return R_n
end

"""
    liquid_precipitation(atmos::PrescribedAtmosphere, p, t)

Returns the liquid precipitation (m/s) at the surface.
"""
function liquid_precipitation(atmos::PrescribedAtmosphere, p, t)
    return atmos.liquid_precip(t)
end

"""
    snow_precipitation(atmos::PrescribedAtmosphere, p, t)

Returns the precipitation in snow (m of liquid water/s) at the surface.
"""
function snow_precipitation(atmos::PrescribedAtmosphere, p, t)
    return atmos.snow_precip(t)
end

"""
    surface_temperature(model::AbstractModel, Y, p, t)

A helper function which returns the surface temperature for a given
model, needed because different models compute and store surface temperature in
different ways and places.

Extending this function for your model is only necessary if you need to
compute surface fluxes and radiative fluxes at the surface using
the functions in this file.
"""
function surface_temperature(model::AbstractModel, Y, p, t) end

"""
    surface_resistance(model::AbstractModel, Y, p, t)

A helper function which returns the surface resistance for a given
model, needed because different models compute and store surface resistance in
different ways and places.

Extending this function for your model is only necessary if you need to
compute surface fluxes and radiative fluxes at the surface using
the functions in this file.

The default is 0, which is no additional resistance aside from the usual
aerodynamic resistance from MOST.
"""
function surface_resistance(model::AbstractModel{FT}, Y, p, t) where {FT}
    return FT(0)
end


"""
    surface_air_density(
                        atmos::AbstractAtmosphericDrivers,
                        model::AbstractModel,
                        Y,
                        p,
                        t,
                        T_sfc,
                        )

A helper function which returns the surface air density for a given
model, needed because different models compute and store surface air density
 in different ways and places.

We additionally include the `atmos` type as an argument because
the surface air density computation may change between a coupled simulation
and a prescibed atmos simulation.

Extending this function for your model is only necessary if you need to
compute surface fluxes and radiative fluxes at the surface using
the functions in this file.
"""
function surface_air_density(
    atmos::AbstractAtmosphericDrivers,
    model::AbstractModel,
    Y,
    p,
    t,
    T_sfc,
) end

"""
    surface_specific_humidity(model::AbstractModel, Y, p, T_sfc, ρ_sfc)

A helper function which returns the surface specific humidity for a given
model, needed because different models compute and store q_sfc in
different ways and places.

Extending this function for your model is only necessary if you need to
compute surface fluxes and radiative fluxes at the surface using
the functions in this file.
"""
function surface_specific_humidity(model::AbstractModel, Y, p, T_sfc, ρ_sfc) end

"""
    surface_evaporative_scaling(model::AbstractModel{FT}, Y, p) where {FT}

A helper function which returns the surface evaporative scaling factor
 for a given model, needed because different models compute and store β_sfc in
different ways and places. Currently, this factor is 1 for all models
besides the bucket model, so we have chosen a default of 1.

Extending this function for your model is only necessary if you need to
compute surface fluxes and radiative fluxes at the surface using
the functions in this file.
"""
function surface_evaporative_scaling(model::AbstractModel{FT}, Y, p) where {FT}
    return FT(1)
end


"""
    surface_albedo(model::AbstractModel, Y, p)

A helper function which returns the surface albedo
 for a given model, needed because different models compute and store α_sfc in
different ways and places.

Extending this function for your model is only necessary if you need to
compute surface fluxes and radiative fluxes at the surface using
the functions in this file.
"""
function surface_albedo(model::AbstractModel, Y, p) end

"""
    surface_emissivity(model::AbstractModel, Y, p)

A helper function which returns the surface emissivity
 for a given model, needed because different models compute and store ϵ_sfc in
different ways and places.

Extending this function for your model is only necessary if you need to
compute surface fluxes and radiative fluxes at the surface using
the functions in this file.
"""
function surface_emissivity(model::AbstractModel, Y, p) end


"""
    surface_height(model::AbstractModel, Y, p)

A helper function which returns the surface height (canopy height+elevation)
 for a given model, needed because different models compute and store h_sfc in
different ways and places.

Extending this function for your model is only necessary if you need to
compute surface fluxes and radiative fluxes at the surface using
the functions in this file.
"""
function surface_height(model::AbstractModel, Y, p) end

"""
    displacement_height(model::AbstractModel, Y, p)

A helper function which returns the displacement height 
 for a given model; the default is zero.

Extending this function for your model is only necessary if you need to
compute surface fluxes and radiative fluxes at the surface using
the functions in this file.
"""
function displacement_height(model::AbstractModel{FT}, Y, p) where {FT}
    return FT(0)
end

"""
    vapor_pressure_deficit(T_air, P_air, q_air, thermo_params)

Computes the vapor pressure deficit for air with temperature T_air,
pressure P_air, and specific humidity q_air, using thermo_params,
a Thermodynamics.jl param set.
"""
function vapor_pressure_deficit(T_air, P_air, q_air, thermo_params)
    es = Thermodynamics.saturation_vapor_pressure(
        thermo_params,
        T_air,
        Thermodynamics.Liquid(),
    )
    ea = Thermodynamics.partial_pressure_vapor(
        thermo_params,
        P_air,
        Thermodynamics.PhasePartition(q_air),
    )
    return es - ea
end
