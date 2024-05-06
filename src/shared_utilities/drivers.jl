import ClimaUtilities.TimeVaryingInputs:
    TimeVaryingInput, AbstractTimeVaryingInput
using Thermodynamics
using ClimaCore
using DocStringExtensions
using SurfaceFluxes
import SurfaceFluxes.Parameters as SFP
using StaticArrays
import ..Parameters as LP
export AbstractAtmosphericDrivers,
    AbstractRadiativeDrivers,
    PrescribedAtmosphere,
    PrescribedPrecipitation,
    CoupledAtmosphere,
    PrescribedRadiativeFluxes,
    CoupledRadiativeFluxes,
    compute_ρ_sfc,
    set_atmos_ts!,
    turbulent_fluxes,
    net_radiation,
    turbulent_fluxes_at_a_point,
    liquid_precipitation,
    snow_precipitation,
    vapor_pressure_deficit,
    displacement_height,
    make_update_drivers

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
    PrescribedAtmosphere{FT, CA, DT} <: AbstractAtmosphericDrivers{FT}

Container for holding prescribed atmospheric drivers and other
information needed for computing turbulent surface fluxes when
driving land models in standalone mode.

The default CO2 concentration is a constant as a function of time, equal to
4.2e-4 mol/mol.

Since not all models require co2 concentration, the default for that
is `nothing`.
$(DocStringExtensions.FIELDS)
"""
struct PrescribedAtmosphere{
    FT,
    LP <: Union{Nothing, AbstractTimeVaryingInput},
    SP <: Union{Nothing, AbstractTimeVaryingInput},
    TA <: Union{Nothing, AbstractTimeVaryingInput},
    UA <: Union{Nothing, AbstractTimeVaryingInput},
    QA <: Union{Nothing, AbstractTimeVaryingInput},
    RA <: Union{Nothing, AbstractTimeVaryingInput},
    CA <: Union{Nothing, AbstractTimeVaryingInput},
    DT,
    TP,
} <: AbstractAtmosphericDrivers{FT}
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
    "Thermodynamic parameters"
    thermo_params::TP
    function PrescribedAtmosphere(
        liquid_precip,
        snow_precip,
        T,
        u,
        q,
        P,
        ref_time,
        h::FT,
        earth_param_set;
        gustiness = FT(1),
        c_co2 = TimeVaryingInput((t) -> 4.2e-4),
    ) where {FT}
        thermo_params = LP.thermodynamic_parameters(earth_param_set)
        args = (liquid_precip, snow_precip, T, u, q, P, c_co2, ref_time)
        return new{typeof(h), typeof.(args)..., typeof(thermo_params)}(
            args...,
            h,
            gustiness,
            thermo_params,
        )
    end
end

"""
    PrescribedPrecipitation{FT, LP} <: AbstractAtmosphericDrivers{FT}
Container for holding prescribed precipitation driver
for models which only require precipitation (RichardsModel).
$(DocStringExtensions.FIELDS)
"""
struct PrescribedPrecipitation{FT, LP <: AbstractTimeVaryingInput} <:
       AbstractAtmosphericDrivers{FT}
    "Precipitation (m/s) function of time: positive by definition"
    liquid_precip::LP
end

PrescribedPrecipitation{FT}(liquid_precip) where {FT} =
    PrescribedPrecipitation{FT, typeof(liquid_precip)}(liquid_precip)

"""
    CoupledRadiativeFluxes{FT} <: AbstractRadiativeDrivers{FT}

To be used when coupling to an atmosphere model.
"""
struct CoupledRadiativeFluxes{FT} <: AbstractRadiativeDrivers{FT} end

"""
    CoupledAtmosphere{FT} <: AbstractAtmosphericDrivers{FT}

To be used when coupling to an atmosphere model.
"""
struct CoupledAtmosphere{FT} <: AbstractAtmosphericDrivers{FT} end

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
    set_atmos_ts!(ts_in, atmos::PrescribedAtmosphere{FT}, p)

Fill the pre-allocated ts_in `Field` with a thermodynamic state computed from the
atmosphere.
"""
function set_atmos_ts!(ts_in, atmos::PrescribedAtmosphere{FT}, p) where {FT}
    P = p.drivers.P
    T = p.drivers.T
    q = p.drivers.q
    ts_in .= Thermodynamics.PhaseEquil_pTq.(atmos.thermo_params, P, T, q)
    return nothing
end

"""
    turbulent_fluxes(atmos::PrescribedAtmosphere,
                   model::AbstractModel,
                   Y::ClimaCore.Fields.FieldVector,
                   p::NamedTuple,
                   t
                   )

Computes the turbulent surface flux terms at the ground for a standalone simulation,
including turbulent energy fluxes as well as the water vapor flux
(in units of m^3/m^2/s of water).
Positive fluxes indicate flow from the ground to the atmosphere.

It solves for these given atmospheric conditions, stored in `atmos`,
model parameters, and the surface conditions.
"""
function turbulent_fluxes(
    atmos::PrescribedAtmosphere,
    model::AbstractModel,
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    t,
)
    T_sfc = surface_temperature(model, Y, p, t)
    ρ_sfc = surface_air_density(atmos, model, Y, p, t, T_sfc)
    q_sfc = surface_specific_humidity(model, Y, p, T_sfc, ρ_sfc)
    β_sfc = surface_evaporative_scaling(model, Y, p)
    h_sfc = surface_height(model, Y, p)
    r_sfc = surface_resistance(model, Y, p, t)
    d_sfc = displacement_height(model, Y, p)
    u_air = p.drivers.u
    h_air = atmos.h

    return turbulent_fluxes_at_a_point.(
        T_sfc,
        q_sfc,
        ρ_sfc,
        β_sfc,
        h_sfc,
        r_sfc,
        d_sfc,
        p.drivers.thermal_state,
        u_air,
        h_air,
        atmos.gustiness,
        model.parameters.z_0m,
        model.parameters.z_0b,
        Ref(model.parameters.earth_param_set),
    )
end

"""
    turbulent_fluxes_at_a_point(T_sfc::FT,
                                q_sfc::FT,
                                ρ_sfc::FT,
                                β_sfc::FT,
                                h_sfc::FT,
                                r_sfc::FT,
                                d_sfc::FT,
                                ts_in,
                                u::FT,
                                h::FT,
                                gustiness::FT,
                                z_0m::FT,
                                z_0b::FT,
                                earth_param_set::EP,
                               ) where {FT <: AbstractFloat, P}

Computes turbulent surface fluxes at a point on a surface given
(1) the surface temperature (T_sfc), specific humidity (q_sfc),
    and air density (ρ_sfc),
(2) Other surface properties, such as the factor
    β_sfc  which scales the evaporation from the potential rate
    (used in bucket models), and the surface resistance r_sfc (used
    in more complex land models), and the topographical height of the surface (h_sfc)
(3) the roughness lengths `z_0m, z_0b`, and the Earth parameter set for the model
    `earth_params`.
(4) the prescribed atmospheric state, `ts_in`, u, h the height
    at which these measurements are made, and the gustiness parameter (m/s).
(5) the displacement height for the model d_sfc

This returns an energy flux and a liquid water volume flux, stored in
a tuple with self explanatory keys.
"""
function turbulent_fluxes_at_a_point(
    T_sfc::FT,
    q_sfc::FT,
    ρ_sfc::FT,
    β_sfc::FT,
    h_sfc::FT,
    r_sfc::FT,
    d_sfc::FT,
    ts_in,
    u::FT,
    h::FT,
    gustiness::FT,
    z_0m::FT,
    z_0b::FT,
    earth_param_set::EP,
) where {FT <: AbstractFloat, EP}
    thermo_params = LP.thermodynamic_parameters(earth_param_set)
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
        gustiness = gustiness,
    )
    surface_flux_params = LP.surface_fluxes_parameters(earth_param_set)
    conditions = SurfaceFluxes.surface_conditions(
        surface_flux_params,
        sc;
        tol_neutral = SFP.cp_d(surface_flux_params) / 100000,
    )
    _LH_v0::FT = LP.LH_v0(earth_param_set)
    _ρ_liq::FT = LP.ρ_cloud_liq(earth_param_set)

    cp_m::FT = Thermodynamics.cp_m(thermo_params, ts_in)
    T_in::FT = Thermodynamics.air_temperature(thermo_params, ts_in)
    ΔT = T_in - T_sfc
    r_ae::FT = 1 / (conditions.Ch * SurfaceFluxes.windspeed(sc))
    ρ_air::FT = Thermodynamics.air_density(thermo_params, ts_in)

    E0::FT = SurfaceFluxes.evaporation(surface_flux_params, sc, conditions.Ch)
    E = E0 * r_ae / (r_sfc + r_ae)
    Ẽ = E / _ρ_liq
    H = -ρ_air * cp_m * ΔT / r_ae
    LH = _LH_v0 * E
    return (lhf = LH, shf = H, vapor_flux = Ẽ, r_ae = r_ae)
end

"""
    turbulent_fluxes(atmos::CoupledAtmosphere,
                    model::AbstractModel,
                    Y,
                    p,
                    t)

Computes the turbulent surface fluxes terms at the ground for a coupled simulation.
"""
function ClimaLand.turbulent_fluxes(
    atmos::CoupledAtmosphere,
    model::AbstractModel,
    Y,
    p,
    t,
)
    # coupler has done its thing behind the scenes already
    model_name = ClimaLand.name(model)
    model_cache = getproperty(p, model_name)
    return model_cache.turbulent_fluxes
end


"""
    PrescribedRadiativeFluxes{FT, SW, LW, DT, T} <: AbstractRadiativeDrivers{FT}

Container for the prescribed radiation functions needed to drive land models in standalone mode.
$(DocStringExtensions.FIELDS)
"""
struct PrescribedRadiativeFluxes{
    FT,
    SW <: Union{Nothing, AbstractTimeVaryingInput},
    LW <: Union{Nothing, AbstractTimeVaryingInput},
    DT,
    T,
} <: AbstractRadiativeDrivers{FT}
    "Downward shortwave radiation function of time (W/m^2): positive indicates towards surface"
    SW_d::SW
    "Downward longwave radiation function of time (W/m^2): positive indicates towards surface"
    LW_d::LW
    "Reference time - the datetime corresponding to t=0 for the simulation"
    ref_time::DT
    "Sun zenith angle, in radians"
    θs::T
    function PrescribedRadiativeFluxes(FT, SW_d, LW_d, ref_time; θs = nothing)
        args = (SW_d, LW_d, ref_time, θs)
        return new{FT, typeof.(args)...}(args...)
    end
end

"""
    net_radiation(radiation::PrescribedRadiativeFluxes{FT},
                  model::AbstractModel{FT},
                  Y::ClimaCore.Fields.FieldVector,
                  p::NamedTuple,
                  t,
                  ) where {FT}

Computes net radiative fluxes for a prescribed incoming
longwave and shortwave radiation.

This returns an energy flux.
"""
function net_radiation(
    radiation::PrescribedRadiativeFluxes,
    model::AbstractModel,
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    t,
)
    LW_d = p.drivers.LW_d
    SW_d = p.drivers.SW_d
    earth_param_set = model.parameters.earth_param_set
    _σ = LP.Stefan(earth_param_set)
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
    net_radiation(radiation::CoupledRadiativeFluxes,
                  model::AbstractModel,
                  Y,
                  p,
                  t)

Computes the net radiative flux at the ground for a coupled simulation.
Your model cache must contain the field `R_n`.
"""
function ClimaLand.net_radiation(
    radiation::CoupledRadiativeFluxes,
    model::AbstractModel,
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    t,
)
    # coupler has done its thing behind the scenes already
    model_name = ClimaLand.name(model)
    model_cache = getproperty(p, model_name)
    return model_cache.R_n
end


"""
    liquid_precipitation(atmos::AbstractAtmosphericDrivers, p, t)

Returns the liquid precipitation (m/s) at the surface.
"""
function liquid_precipitation(atmos::AbstractAtmosphericDrivers, p, t)
    return p.drivers.P_liq
end

"""
    snow_precipitation(atmos::AbstractAtmosphericDrivers, p, t)

Returns the precipitation in snow (m of liquid water/s) at the surface.
"""
function snow_precipitation(atmos::AbstractAtmosphericDrivers, p, t)
    return p.drivers.P_snow
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
                        atmos::PrescribedAtmosphere,
                        model::AbstractModel,
                        Y,
                        p,
                        t,
                        T_sfc,
                        )

A helper function which returns the surface air density; this assumes that
the `model` has a property called `parameters` containing `earth_param_set`.

We additionally include the `atmos` type as an argument because
the surface air density computation will change between a coupled simulation
and a prescibed atmos simulation.

Extending this function for your model is only necessary if you need to
compute the air density in a different way.
"""
function surface_air_density(
    atmos::PrescribedAtmosphere,
    model::AbstractModel,
    Y,
    p,
    t,
    T_sfc,
)
    thermo_params =
        LP.thermodynamic_parameters(model.parameters.earth_param_set)
    return compute_ρ_sfc.(thermo_params, p.drivers.thermal_state, T_sfc)
end


"""
    ClimaLand.surface_air_density(
                    atmos::CoupledAtmosphere,
                    model::AbstractModel,
                    Y,
                    p,
                    _...,
                )
Returns the air density at the surface in the case of a coupled simulation.

This requires the field `ρ_sfc` to be present in the cache `p` under the name
of the model.
"""
function surface_air_density(
    atmos::CoupledAtmosphere,
    model::AbstractModel,
    Y,
    p,
    _...,
)
    model_name = ClimaLand.name(model)
    model_cache = getproperty(p, model_name)
    return model_cache.ρ_sfc
end

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

"""
    initialize_drivers(a::PrescribedAtmosphere{FT}, coords) where {FT}

Creates and returns a NamedTuple for the `PrescribedAtmosphere` driver,
with variables `P_liq`, `P_snow`, and air temperature `T`, pressure `P`,
horizontal wind speed `u`, specific humidity `q`, and CO2 concentration
`c_co2`.
"""
function initialize_drivers(a::PrescribedAtmosphere{FT}, coords) where {FT}
    keys = (:P_liq, :P_snow, :T, :P, :u, :q, :c_co2, :thermal_state)
    # The thermal state is a different type
    types = ([FT for k in keys[1:(end - 1)]]..., Thermodynamics.PhaseEquil{FT})
    domain_names = ([:surface for k in keys]...,)
    model_name = :drivers
    # intialize_vars packages the variables as a named tuple,
    # as part of a named tuple with `model_name` as the key.
    # Here we just want the variable named tuple itself
    vars =
        ClimaLand.initialize_vars(keys, types, domain_names, coords, model_name)
    return vars.drivers
end

"""
    initialize_drivers(a::PrescribedPrecipitation{FT}, coords) where {FT}
Creates and returns a NamedTuple for the `PrescribedPrecipitation` driver,
with variable `P_liq`.
"""
function initialize_drivers(a::PrescribedPrecipitation{FT}, coords) where {FT}
    keys = (:P_liq,)
    types = (FT,)
    domain_names = (:surface,)
    model_name = :drivers
    # intialize_vars packages the variables as a named tuple,
    # as part of a named tuple with `model_name` as the key.
    # Here we just want the variable named tuple itself
    vars =
        ClimaLand.initialize_vars(keys, types, domain_names, coords, model_name)
    return vars.drivers
end

"""
    initialize_drivers(r::CoupledAtmosphere{FT}, coords) where {FT}

Creates and returns a NamedTuple for the `CoupledAtmosphere` driver,
with variables `P_liq`, and `P_snow`. This is intended to be used in coupled
simulations with ClimaCoupler.jl
"""
function ClimaLand.initialize_drivers(
    a::CoupledAtmosphere{FT},
    coords,
) where {FT}
    keys = (:P_liq, :P_snow)
    types = ([FT for k in keys]...,)
    domain_names = ([:surface for k in keys]...,)
    model_name = :drivers
    # intialize_vars packages the variables as a named tuple,
    # as part of a named tuple with `model_name` as the key.
    # Here we just want the variable named tuple itself
    vars =
        ClimaLand.initialize_vars(keys, types, domain_names, coords, model_name)
    return vars.drivers
end

"""
    initialize_drivers(r::PrescribedRadiativeFluxes{FT}, coords) where {FT}

Creates and returns a NamedTuple for the `PrescribedRadiativeFluxes` driver,
 with variables `SW_d`, `LW_d`, and zenith angle `θ_s`.
"""
function initialize_drivers(r::PrescribedRadiativeFluxes{FT}, coords) where {FT}
    keys = (:SW_d, :LW_d, :θs)
    types = ([FT for k in keys]...,)
    domain_names = ([:surface for k in keys]...,)
    model_name = :drivers
    # intialize_vars packages the variables as a named tuple,
    # as part of a named tuple with `model_name` as the key.
    # Here we just want the variable named tuple itself
    vars =
        ClimaLand.initialize_vars(keys, types, domain_names, coords, model_name)
    return vars.drivers
end

"""
    initialize_drivers(d::Union{AbstractAtmosphericDrivers, AbstractRadiativeDrivers, Nothing}, coords)

Creates and returns a NamedTuple with `nothing` when no driver cache variables are needed.
"""
function initialize_drivers(
    d::Union{AbstractAtmosphericDrivers, AbstractRadiativeDrivers, Nothing},
    coords,
)
    return (;)
end

"""
    initialize_drivers(a::Union{AbstractAtmosphericDrivers, Nothing},
                       r::Union{AbstractRadiativeDrivers, Nothing},
                       coords)

Creates and returns a NamedTuple with the cache variables required by the
atmospheric and radiative drivers.

If no forcing is required, `a` and `r` are type `Nothing` and an
empty NamedTuple is returned.
"""
function initialize_drivers(
    a::Union{AbstractAtmosphericDrivers, Nothing},
    r::Union{AbstractRadiativeDrivers, Nothing},
    coords,
)
    atmos_drivers = initialize_drivers(a, coords)
    radiation_drivers = initialize_drivers(r, coords)
    merge(atmos_drivers, radiation_drivers)
end

"""
    make_update_drivers(a::Union{AbstractAtmosphericDrivers, Nothing},
                          r::Union{AbstractRadiativeDrivers, Nothing},
                         )

Creates and returns a function which updates the atmospheric
and radiative forcing variables ("drivers").
"""
function make_update_drivers(
    a::Union{AbstractAtmosphericDrivers, Nothing},
    r::Union{AbstractRadiativeDrivers, Nothing},
)
    update_atmos! = make_update_drivers(a)
    update_radiation! = make_update_drivers(r)
    function update_drivers!(p, t)
        update_atmos!(p, t)
        update_radiation!(p, t)
    end
    return update_drivers!
end

"""
    make_update_drivers(d::Union{AbstractAtmosphericDrivers, AbstractRadiativeDrivers, Nothing})

Creates and returns a function which updates the driver variables
in the case of no driver variables. This is also the default.
"""
make_update_drivers(
    d::Union{AbstractAtmosphericDrivers, AbstractRadiativeDrivers, Nothing},
) = (p, t) -> nothing

"""
    make_update_drivers(a::PrescribedAtmosphere{FT}) where {FT}

Creates and returns a function which updates the driver variables
in the case of a PrescribedAtmosphere.
"""
function make_update_drivers(a::PrescribedAtmosphere{FT}) where {FT}
    function update_drivers!(p, t)
        evaluate!(p.drivers.P_liq, a.liquid_precip, t)
        evaluate!(p.drivers.P_snow, a.snow_precip, t)
        evaluate!(p.drivers.T, a.T, t)
        evaluate!(p.drivers.P, a.P, t)
        evaluate!(p.drivers.u, a.u, t)
        evaluate!(p.drivers.q, a.q, t)
        evaluate!(p.drivers.c_co2, a.c_co2, t)
        set_atmos_ts!(p.drivers.thermal_state, a, p)
    end
    return update_drivers!
end

"""
    make_update_drivers(a::PrescribedPrecipitation{FT}) where {FT}
Creates and returns a function which updates the driver variables
in the case of a PrescribedPrecipitation.
"""
function make_update_drivers(a::PrescribedPrecipitation{FT}) where {FT}
    function update_drivers!(p, t)
        evaluate!(p.drivers.P_liq, a.liquid_precip, t)
    end
    return update_drivers!
end

"""
    make_update_drivers(r::PrescribedRadiativeFluxes{FT}) where {FT}

Creates and returns a function which updates the driver variables
in the case of a PrescribedRadiativeFluxes.
"""
function make_update_drivers(r::PrescribedRadiativeFluxes{FT}) where {FT}
    function update_drivers!(p, t)
        evaluate!(p.drivers.SW_d, r.SW_d, t)
        evaluate!(p.drivers.LW_d, r.LW_d, t)
        if !isnothing(r.θs)
            p.drivers.θs .= FT.(r.θs(t, r.ref_time))
        else
            p.drivers.θs .= FT(0)
        end
    end
    return update_drivers!
end
