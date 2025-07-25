import ClimaUtilities.TimeVaryingInputs:
    TimeVaryingInput,
    AbstractTimeVaryingInput,
    LinearInterpolation,
    PeriodicCalendar
import ClimaUtilities.Regridders: InterpolationsRegridder
import ClimaUtilities.FileReaders: NCFileReader, read
import ClimaUtilities.TimeManager: ITime, date
using Thermodynamics
using ClimaCore
using Dates
using DocStringExtensions
using Insolation
using SurfaceFluxes
import SurfaceFluxes.Parameters as SFP
using StaticArrays
import ..Parameters as LP
export AbstractAtmosphericDrivers,
    AbstractRadiativeDrivers,
    PrescribedAtmosphere,
    PrescribedPrecipitation,
    PrescribedSoilOrganicCarbon,
    PrescribedGroundConditions,
    CoupledAtmosphere,
    PrescribedRadiativeFluxes,
    CoupledRadiativeFluxes,
    compute_ρ_sfc,
    set_atmos_ts!,
    turbulent_fluxes!,
    net_radiation!,
    turbulent_fluxes_at_a_point,
    vapor_pressure_deficit,
    displacement_height,
    relative_humidity,
    specific_humidity_from_dewpoint,
    make_update_drivers,
    prescribed_lai_era5,
    prescribed_forcing_era5,
    prescribed_analytic_forcing,
    default_zenith_angle

"""
     AbstractClimaLandDrivers{FT <: AbstractFloat}

An abstract type of drivers (forcing) for land models.
"""
abstract type AbstractClimaLandDrivers{FT <: AbstractFloat} end

"""
    initialize_drivers(::AbstractClimaLandDrivers, coords)

Creates and returns a default empty NamedTuple for AbstractClimaLandDrivers.
More generally this should return a named tuple of the driver fields, which will
then be stored in the cache under `p.drivers`.
"""
initialize_drivers(::AbstractClimaLandDrivers, coords) = (;)


"""
    make_update_drivers(::AbstractClimaLandDrivers)

Creates and returns a function which updates the driver variables
in the default case of no drivers. More generally, this should return
a function which updates the driver fields stored in `p.drivers`.
"""
function make_update_drivers(::AbstractClimaLandDrivers)
    update_drivers!(p, t) = nothing
    return update_drivers!
end


"""
     AbstractAtmosphericDrivers{FT}

An abstract type of atmospheric drivers of land models.
"""
abstract type AbstractAtmosphericDrivers{FT} <: AbstractClimaLandDrivers{FT} end

"""
     AbstractRadiativeDrivers{FT}

An abstract type of radiative drivers of land models.
"""
abstract type AbstractRadiativeDrivers{FT} <: AbstractClimaLandDrivers{FT} end


"""
     PrescribedSoilOrganicCarbon{FT}

A type for prescribing soil organic carbon.
$(DocStringExtensions.FIELDS)
"""
struct PrescribedSoilOrganicCarbon{FT, SOC <: AbstractTimeVaryingInput} <:
       AbstractClimaLandDrivers{FT}
    "Soil organic carbon, function of time and space: kg C/m^3"
    soc::SOC
end

PrescribedSoilOrganicCarbon{FT}(soc) where {FT} =
    PrescribedSoilOrganicCarbon{FT, typeof(soc)}(soc)

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
    "Start date - the datetime corresponding to t=0 for the simulation"
    start_date::DT
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
        start_date,
        h::FT,
        earth_param_set;
        gustiness = FT(1),
        c_co2 = TimeVaryingInput((t) -> 4.2e-4),
    ) where {FT}
        thermo_params = LP.thermodynamic_parameters(earth_param_set)
        args = (liquid_precip, snow_precip, T, u, q, P, c_co2, start_date)
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
    CoupledRadiativeFluxes{
        FT,
        F <: Union{Function, Nothing},
        T,
    } <: AbstractRadiativeDrivers{FT}

To be used when coupling to an atmosphere model. Either both `θs` and `start_date`
must be `nothing`, or both must not be `nothing``.

During the driver update, cosθs is unchanged if `θs` is `nothing`. This behavior differs from
the `PrescribedRadiativeFluxes` where `cosθs` set to `NaN` if `θs` is `nothing`.
Otherwise, `θs` recieves the following arguments:
(time_from_start, `start_date`), and is expected to return zenith angle at the given time.
$(DocStringExtensions.FIELDS)
"""
struct CoupledRadiativeFluxes{FT, F <: Union{Function, Nothing}, T} <:
       AbstractRadiativeDrivers{FT}
    """Function that fills a climacore field with the zenith angle given the following arguments:
    (time_from_start, `start_date`)"""
    θs::F
    "Start date - the datetime corresponding to t=0 for the simulation"
    start_date::T
    function CoupledRadiativeFluxes{FT, F, T}(
        θs::F,
        start_date::T,
    ) where {FT, F, T}
        (
            (isnothing(θs) && isnothing(start_date)) ||
            ((!isnothing(θs) && !isnothing(start_date)))
        ) || error(
            "CoupledRadiativeFluxes: `θs` and start_date` must both be `nothing` or both not `nothing`.",
        )
        new{FT, F, T}(θs, start_date)
    end
end

# This constructor is here to maintain backwards compatability with ClimaCoupler
# If this constructor is used, the coupler should calulate and update the zenith angle
CoupledRadiativeFluxes{FT}() where {FT} =
    CoupledRadiativeFluxes{FT, Nothing, Nothing}(nothing, nothing)

CoupledRadiativeFluxes(::Type{FT}, args...) where {FT} =
    CoupledRadiativeFluxes{FT}(args...)

"""
    CoupledRadiativeFluxes{FT}(
        start_date::Dates.DateTime;
        latitude,
        longitude,
        insol_params=LP.LandParameters(FT).insol_params,
    )
Creates a `CoupledRadiativeFluxes` object with a default zenith angle function that uses Insolation.jl
to compute the zenith angle at a given time and location.
"""
function CoupledRadiativeFluxes{FT}(
    start_date::DT;
    latitude::LT,
    longitude::LT,
    insol_params::IP = LP.LandParameters(FT).insol_params,
) where {FT, DT, LT, IP}
    zenith_angle =
        (t, s) -> default_zenith_angle(
            t,
            s;
            latitude = latitude,
            longitude = longitude,
            insol_params = insol_params,
        )
    return CoupledRadiativeFluxes{FT, typeof(zenith_angle), DT}(
        zenith_angle,
        start_date,
    )
end

"""
    default_zenith_angle(
        t::T,
        start_date::Dates.DateTime;
        latitude::LT,
        longitude::LT,
        insol_params::Insolation.Parameters.InsolationParameters{FT},
    )

Calculate zenith angle with Insolation for the given start date, insolation parameters, latitude,
and longitude.

`latitude` and `longitude` can be a collections or a Number.
"""
function default_zenith_angle(
    t::T,
    start_date::Dates.DateTime;
    latitude::LT,
    longitude::LT,
    insol_params,
) where {T, LT}
    FT = eltype(latitude)
    current_datetime =
        T <: ITime ? date(t) : start_date + Dates.Second(round(t))
    d, δ, η_UTC =
        FT.(
            Insolation.helper_instantaneous_zenith_angle(
                current_datetime,
                start_date,
                insol_params,
            )
        )
    # Reduces allocations by throwing away unwanted values
    zenith_only = (args...) -> Insolation.instantaneous_zenith_angle(args...)[1]
    return zenith_only.(d, δ, η_UTC, longitude, latitude)
end

"""
    CoupledAtmosphere{FT} <: AbstractAtmosphericDrivers{FT}

To be used when coupling to an atmosphere model. Contains fields that
are used to compute surface fluxes in the coupled setup.

When constructed without a space, the struct doesn't contain anything,
but it still acts as a flag that fluxes have been updated by the coupler
and don't need to be recomputed.
When constructed with a space, the struct contains the fields needed to compute
surface fluxes in the coupled setup, which are accessed by ClimaCoupler.
"""
struct CoupledAtmosphere{FT} <: AbstractAtmosphericDrivers{FT}
    "Atmospheric horizontal wind velocity vector at the reference height (m/s)"
    u::Union{Nothing, Fields.Field}
    "Atmospheric reference height (m), relative to surface elevation"
    h::Union{Nothing, Fields.Field}
    "Minimum wind speed (gustiness; m/s), which is always a spatial constant"
    gustiness::Union{Nothing, FT}
    "Atmospheric thermodynamic state at the reference height"
    thermal_state::Union{Nothing, Fields.Field}
    "Create a `CoupledAtmosphere` with default values"
    function CoupledAtmosphere{FT}() where {FT}
        return new{FT}(nothing, nothing, nothing, nothing)
    end
    function CoupledAtmosphere{FT}(space) where {FT}
        return new{FT}(
            Fields.zeros(SVector{2, FT}, space),
            Fields.zeros(space),
            FT(1), # gustiness is always a spatial constant, for now
            Fields.zeros(Thermodynamics.PhaseEquil{FT}, space),
        )
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
    turbulent_fluxes!(dest,
                      atmos::PrescribedAtmosphere,
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
function turbulent_fluxes!(
    dest,
    atmos::PrescribedAtmosphere,
    model::AbstractModel,
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    t,
)
    T_sfc = surface_temperature(model, Y, p, t)
    ρ_sfc = surface_air_density(atmos, model, Y, p, t, T_sfc)
    q_sfc = surface_specific_humidity(atmos, model, Y, p, T_sfc, ρ_sfc)
    β_sfc = surface_evaporative_scaling(model, Y, p)
    h_sfc = surface_height(model, Y, p)
    r_sfc = surface_resistance(model, Y, p, t)
    d_sfc = displacement_height(model, Y, p)
    u_air = p.drivers.u
    h_air = atmos.h

    dest .=
        turbulent_fluxes_at_a_point.(
            Val(false), # return_extra_fluxes
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
            model.parameters.earth_param_set,
        )
    return nothing
end

"""
    coupler_compute_turbulent_fluxes!(dest, atmos::CoupledAtmosphere, model::AbstractModel, Y::ClimaCore.Fields.FieldVector, p::NamedTuple, t)

This function computes the turbulent surface fluxes for a coupled simulation.
This function is very similar to the default method of `turbulent_fluxes!`,
but it is used with a `CoupledAtmosphere` which contains all the necessary
atmosphere fields to compute the surface fluxes, rather than some being stored in `p`.

This function is intended to be called by ClimaCoupler.jl when computing
fluxes for a coupled simulation with the integrated land model.
"""
function coupler_compute_turbulent_fluxes!(
    dest,
    atmos::CoupledAtmosphere,
    model::AbstractModel,
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    t,
)
    T_sfc = surface_temperature(model, Y, p, t)
    ρ_sfc = surface_air_density(atmos, model, Y, p, t, T_sfc)
    q_sfc = surface_specific_humidity(atmos, model, Y, p, T_sfc, ρ_sfc)
    β_sfc = surface_evaporative_scaling(model, Y, p)
    h_sfc = surface_height(model, Y, p)
    r_sfc = surface_resistance(model, Y, p, t)
    d_sfc = displacement_height(model, Y, p)

    dest .=
        turbulent_fluxes_at_a_point.(
            Val(true), # return_extra_fluxes
            T_sfc,
            q_sfc,
            ρ_sfc,
            β_sfc,
            h_sfc,
            r_sfc,
            d_sfc,
            atmos.thermal_state,
            atmos.u,
            atmos.h,
            atmos.gustiness,
            model.parameters.z_0m,
            model.parameters.z_0b,
            model.parameters.earth_param_set,
        )
    return nothing
end

"""
    turbulent_fluxes_at_a_point(return_extra_fluxes, args...)

This is a wrapper function that allows us to dispatch on the type of `return_extra_fluxes`
as we compute the turbulent fluxes pointwise. This is needed because space for the
extra fluxes is only allocated in the cache when running with a `CoupledAtmosphere`.
The function `compute_turbulent_fluxes_at_a_point` does the actual flux computation.

The `return_extra_fluxes` argument indicates whether to return the following:
- momentum fluxes (`ρτxz`, `ρτyz`)
- buoyancy flux (`buoy_flux`)
"""
function turbulent_fluxes_at_a_point(return_extra_fluxes::Val{false}, args...)
    (LH, SH, Ẽ, r_ae, _, _, _) = compute_turbulent_fluxes_at_a_point(args...)
    return (lhf = LH, shf = SH, vapor_flux = Ẽ, r_ae = r_ae)
end
function turbulent_fluxes_at_a_point(return_extra_fluxes::Val{true}, args...)
    (LH, SH, Ẽ, r_ae, ρτxz, ρτyz, buoy_flux) =
        compute_turbulent_fluxes_at_a_point(args...)
    return (
        lhf = LH,
        shf = SH,
        vapor_flux = Ẽ,
        r_ae = r_ae,
        ρτxz = ρτxz,
        ρτyz = ρτyz,
        buoy_flux = buoy_flux,
    )
end

"""
    compute_turbulent_fluxes_at_a_point(T_sfc::FT,
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
                                earth_param_set::EP;
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

Please note that this function, if r_sfc is set to zero, makes no alteration
to the output of the `SurfaceFluxes.surface_conditions` call, aside
from unit conversion of evaporation from a mass to a liquid water volume
flux. When r_sfc is nonzero, an additional resistance is applied in series
to the vapor flux (and hence also the latent heat flux).
"""
function compute_turbulent_fluxes_at_a_point(
    T_sfc::FT,
    q_sfc::FT,
    ρ_sfc::FT,
    β_sfc::FT,
    h_sfc::FT,
    r_sfc::FT,
    d_sfc::FT,
    ts_in,
    u::Union{FT, SVector{2, FT}},
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
    # u is already a vector when we get it from a coupled atmosphere, otherwise we need to make it one
    if u isa FT
        u = SVector{2, FT}(u, 0)
    end
    state_in = SurfaceFluxes.StateValues(h - d_sfc - h_sfc, u, ts_in)
    # The following line wont work on GPU
    #    h - d_sfc - h_sfc < 0 &&
    #        @error("Surface height is larger than atmos height in surface fluxes")
    # State containers
    states = SurfaceFluxes.ValuesOnly(
        state_in,
        state_sfc,
        z_0m,
        z_0b,
        beta = β_sfc,
        gustiness = gustiness,
    )
    surface_flux_params = LP.surface_fluxes_parameters(earth_param_set)
    scheme = SurfaceFluxes.PointValueScheme()
    conditions =
        SurfaceFluxes.surface_conditions(surface_flux_params, states, scheme)
    _LH_v0::FT = LP.LH_v0(earth_param_set)
    _ρ_liq::FT = LP.ρ_cloud_liq(earth_param_set)

    # aerodynamic resistance
    r_ae::FT = 1 / (conditions.Ch * SurfaceFluxes.windspeed(states))

    # latent heat flux
    E0::FT =
        SurfaceFluxes.evaporation(surface_flux_params, states, conditions.Ch) # mass flux at potential evaporation rate
    E = E0 * r_ae / (r_sfc + r_ae) # mass flux accounting for additional surface resistance
    LH = _LH_v0 * E # Latent heat flux

    # sensible heat flux
    SH = SurfaceFluxes.sensible_heat_flux(
        surface_flux_params,
        conditions.Ch,
        states,
        scheme,
    )

    # vapor flux in volume of liquid water with density 1000kg/m^3
    Ẽ = E / _ρ_liq

    return (
        LH,
        SH,
        Ẽ,
        r_ae,
        conditions.ρτxz,
        conditions.ρτyz,
        conditions.buoy_flux,
    )
end

"""
    turbulent_fluxes!(dest,
                    atmos::CoupledAtmosphere,
                    model::AbstractModel,
                    Y,
                    p,
                    t)

Computes the turbulent surface fluxes terms at the ground for a coupled simulation.
In this case, the coupler has already computed turbulent fluxes and updated
them in each of the component models, so this function does nothing.
"""
function ClimaLand.turbulent_fluxes!(
    dest,
    atmos::CoupledAtmosphere,
    model::AbstractModel,
    Y,
    p,
    t,
)
    # coupler has done its thing behind the scenes already
    return nothing
end


"""
    PrescribedRadiativeFluxes{FT, SW, LW, DT, T, TP} <: AbstractRadiativeDrivers{FT}

Container for the prescribed radiation functions needed to drive land models in standalone mode.

Note that some models require the zenith angle AND diffuse fraction. The diffuse fraction
may be provided directly (of type TimeVaryingInput), or it may be computed empirically.
In the latter case, it requires the thermodynamic parameters as well to compute.

Therefore, the allowed combinations are:
1. Zenith angle and diffuse fraction not needed: zenith angle=nothing, thermo_params=nothing, diffuse fraction=nothing
2. Zenith angle provided and diffuse fraction computed empirically:  thermo params used, diffuse fraction=nothing
3. Zenith angle provided and diffuse fraction provided: thermo_params not used, diffuse fraction TimeVaryingInput
$(DocStringExtensions.FIELDS)
"""
struct PrescribedRadiativeFluxes{
    FT,
    SW <: Union{Nothing, AbstractTimeVaryingInput},
    FD <: Union{Nothing, AbstractTimeVaryingInput},
    LW <: Union{Nothing, AbstractTimeVaryingInput},
    DT,
    T,
    TP,
} <: AbstractRadiativeDrivers{FT}
    "Downward shortwave radiation function of time (W/m^2): positive indicates towards surface"
    SW_d::SW
    "Diffuse Fraction of shortwave radiation (unitless, [0,1])"
    frac_diff::FD
    "Downward longwave radiation function of time (W/m^2): positive indicates towards surface"
    LW_d::LW
    "Start date - the datetime corresponding to t=0 for the simulation"
    start_date::DT
    "Sun zenith angle, in radians"
    θs::T
    "Thermodynamic parameters"
    thermo_params::TP
    function PrescribedRadiativeFluxes(
        FT,
        SW_d,
        LW_d,
        start_date;
        θs = nothing,
        earth_param_set = nothing,
        frac_diff = nothing,
    )
        if isnothing(earth_param_set)
            thermo_params = nothing
        else
            thermo_params = LP.thermodynamic_parameters(earth_param_set)
        end
        args = (SW_d, frac_diff, LW_d, start_date, θs, thermo_params)
        return new{FT, typeof.(args)...}(args...)
    end
end

"""
    net_radiation!(dest::ClimaCore.Fields.Field,
                   radiation::PrescribedRadiativeFluxes{FT},
                   model::AbstractModel{FT},
                   Y::ClimaCore.Fields.FieldVector,
                   p::NamedTuple,
                   t,
                   ) where {FT}


Computes net radiative  fluxes for a prescribed incoming  longwave and shortwave
radiation.
"""
function net_radiation!(
    dest::ClimaCore.Fields.Field,
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
    @. dest = -(1 - α_sfc) * SW_d - ϵ_sfc * (LW_d - _σ * T_sfc^4)
    return nothing
end


"""
    net_radiation!(dest,
                  radiation::CoupledRadiativeFluxes,
                  model::AbstractModel,
                  Y,
                  p,
                  t)

Computes the net radiative flux at the ground for a coupled simulation.
Your model cache must contain the field `R_n`.
"""
function ClimaLand.net_radiation!(
    dest,
    radiation::CoupledRadiativeFluxes,
    model::AbstractModel,
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    t,
)
    # coupler has done its thing behind the scenes already
    model_name = ClimaLand.name(model)
    model_cache = getproperty(p, model_name)
    dest .= model_cache.R_n
    return nothing
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
    surface_air_density(
                        atmos::CoupledAtmosphere,
                        model::AbstractModel,
                        Y,
                        p,
                        t,
                        T_sfc,
                        )

A helper function which returns the surface air density; this assumes that
the `model` has a property called `parameters` containing `earth_param_set`.

This method is similar to the general method above, except in this case
we get the thermodynamic parameters from the `atmos` object. This is used
when running with a coupled atmosphere.
"""
function surface_air_density(
    atmos::CoupledAtmosphere,
    model::AbstractModel,
    Y,
    p,
    t,
    T_sfc,
)
    thermo_params =
        LP.thermodynamic_parameters(model.parameters.earth_param_set)
    return compute_ρ_sfc.(thermo_params, atmos.thermal_state, T_sfc)
end

"""
    surface_specific_humidity(atmos, model::AbstractModel, Y, p, T_sfc, ρ_sfc)

A helper function which returns the surface specific humidity for a given
model, needed because different models compute and store q_sfc in
different ways and places.

Extending this function for your model is only necessary if you need to
compute surface fluxes and radiative fluxes at the surface using
the functions in this file.
"""
function surface_specific_humidity(
    atmos,
    model::AbstractModel,
    Y,
    p,
    T_sfc,
    ρ_sfc,
) end

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
    relative_humidity(T_air, P_air, q_air, thermo_params)

Computes the vapor pressure deficit for air with temperature T_air,
pressure P_air, and specific humidity q_air, using thermo_params,
a Thermodynamics.jl param set.
"""
function relative_humidity(
    T_air::FT,
    P_air::FT,
    q_air::FT,
    thermo_params,
) where {FT}
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
    return es / ea
end

"""
    specific_humidity_from_dewpoint(dewpoint_temperature, temperature, air_pressure, earth_param_set)

Estimates the specific humidity given the dewpoint temperature, temperature of the air
in Kelvin, and air pressure in Pa, along with the ClimaLand earth_param_set. This is useful
for creating the PrescribedAtmosphere - which needs specific humidity - from ERA5 reanalysis data.

We first compute the relative humidity using the Magnus formula, then the saturated vapor pressure, and then
we compute q from vapor pressure and saturated vapor pressure.

For more information on the Magnus formula, see e.g.
Lawrence, Mark G. (1 February 2005).
"The Relationship between Relative Humidity and the Dewpoint Temperature in Moist Air:
A Simple Conversion and Applications".
Bulletin of the American Meteorological Society. 86 (2): 225–234.
"""
function specific_humidity_from_dewpoint(
    T_dew_air::data_FT,
    T_air::data_FT,
    P_air::data_FT,
    earth_param_set,
) where {data_FT <: Real}
    thermo_params = LP.thermodynamic_parameters(earth_param_set)
    _T_freeze = LP.T_freeze(earth_param_set)
    sim_FT = typeof(_T_freeze)
    # Obtain the relative humidity. This function requires temperatures in Celsius
    rh::sim_FT = rh_from_dewpoint(
        sim_FT(T_dew_air) - _T_freeze,
        sim_FT(T_air) - _T_freeze,
    )

    # In the future, we will use Thermodynamics functions for all of the below. Note that the type of T_air must match the type of the parameters in thermo_params.
    esat::sim_FT = Thermodynamics.saturation_vapor_pressure(
        thermo_params,
        sim_FT(T_air),
        Thermodynamics.Liquid(),
    )

    e = rh * esat
    q = data_FT(0.622 * e / (P_air - 0.378 * e))
    return q
end

"""
    rh_from_dewpoint(Td_C, T_C)

Returns the relative humidity given the dewpoint temperature in Celsius and the
air temperature in Celsius, using the Magnus formula.

For more information on the Magnus formula, see e.g.
Lawrence, Mark G. (1 February 2005).
"The Relationship between Relative Humidity and the Dewpoint Temperature in Moist Air:
A Simple Conversion and Applications".
Bulletin of the American Meteorological Society. 86 (2): 225–234.
"""
function rh_from_dewpoint(Td_C::FT, T_C::FT) where {FT <: Real}
    c = FT(243.04) # C
    b = FT(17.625) # unitless
    γ = Td_C * b / (c + Td_C)
    rh = exp(γ - b * T_C / (c + T_C))
    return rh
end

"""
An abstract type of ground conditions for the canopy model;
we support prescribed (canopy run in standalone mode)
or prognostic ground conditions (integrated land models).
"""
abstract type AbstractGroundConditions{FT} <: AbstractClimaLandDrivers{FT} end

"""
     PrescribedGroundConditions <: AbstractGroundConditions

A container for holding prescribed ground conditions needed by the canopy model
when running the canopy in standalone mode, including the soil pressure, surface
temperature, albedo, and emissivity.
$(DocStringExtensions.FIELDS)
"""
struct PrescribedGroundConditions{
    FT,
    F1 <: AbstractTimeVaryingInput,
    F2 <: AbstractTimeVaryingInput,
} <: AbstractGroundConditions{FT}
    "Prescribed soil potential (m) in the root zone as a function of time"
    ψ::F1
    "Prescribed ground surface temperature (K) as a function of time"
    T::F2
    "Ground albedo for PAR"
    α_PAR::FT
    "Ground albedo for NIR"
    α_NIR::FT
    "Ground emissivity"
    ϵ::FT
end

"""
     function PrescribedGroundConditions{FT}(;
         ψ::TimeVaryingInput,
         T::TimeVaryingInput,
         α_PAR::FT,
         α_NIR::FT,
         ϵ::FT
     ) where {FT}

An outer constructor for the PrescribedGroundConditions allowing the user to
specify the ground parameters by keyword arguments.
"""
function PrescribedGroundConditions{FT}(;
    ψ = TimeVaryingInput((t) -> 0.0),
    T = TimeVaryingInput((t) -> 298.0),
    α_PAR = FT(0.2),
    α_NIR = FT(0.4),
    ϵ = FT(0.99),
) where {FT <: AbstractFloat}
    return PrescribedGroundConditions{FT, typeof(ψ), typeof(T)}(
        ψ,
        T,
        α_PAR,
        α_NIR,
        ϵ,
    )
end

"""
     PrognosticSoilConditions <: AbstractGroundConditions

 A type of AbstractGroundConditions to use when the soil model is prognostic and
of type `EnergyHydrology`. `PrognosticSoilConditions` functions as a flag and is used for dispatch
with the canopy model.
"""
struct PrognosticSoilConditions{FT} <: AbstractGroundConditions{FT} end

"""
     PrognosticGroundConditions <: Canopy.AbstractGroundConditions

A type of AbstractGroundConditions to use when the soil model is prognostic and
of type `EnergyHydrology`, and the snow model is prognostic and included.

Note that this struct is linked with the EnergyHydrology/SnowModel model. If we ever had a different
soil model, we might need to construct a different `PrognosticGroundConditions` because
the fields may be stored in different places.
"""
struct PrognosticGroundConditions{FT} <: AbstractGroundConditions{FT} end

"""
    initialize_drivers(a::PrescribedGroundConditions{FT}, coords) where {FT}

Creates and returns a NamedTuple for the `PrescribedGroundConditions` driver,
with variables `ψ` (matric potential in the soil in the root zone), `T` (temperature
of the surface of the ground).
"""
function initialize_drivers(
    a::PrescribedGroundConditions{FT},
    coords,
) where {FT}
    keys = (:ψ, :T_ground)
    types = (FT, FT)
    domain_names = (:surface, :surface)
    model_name = :drivers
    # intialize_vars packages the variables as a named tuple,
    # as part of a named tuple with `model_name` as the key.
    # Here we just want the variable named tuple itself
    vars =
        ClimaLand.initialize_vars(keys, types, domain_names, coords, model_name)
    return vars.drivers
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
    initialize_drivers(r::AbstractRadiativeDrivers{FT}, coords) where {FT}

Creates and returns a NamedTuple for the radiative fluxes driver,
with variables `SW_d`, `LW_d`, cosine of the zenith angle `cosθs`, and the
diffuse fraction of shortwave radiation `frac_diff`.

We require the same variables for both prescribed and coupled radiative drivers,
so we can use the same function for both of them.
"""
function initialize_drivers(r::AbstractRadiativeDrivers{FT}, coords) where {FT}
    keys = (:SW_d, :LW_d, :cosθs, :frac_diff)
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
    initialize_drivers(r::PrescribedSoilOrganicCarbon{FT}, coords) where {FT}

Creates and returns a NamedTuple for the `PrescribedSoilOrganicCarbon` driver,
 with variable `soc`.
"""
function initialize_drivers(
    r::PrescribedSoilOrganicCarbon{FT},
    coords,
) where {FT}
    keys = (:soc,)
    types = (FT,)
    domain_names = (:subsurface,)
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
with variables `P_liq`, `P_snow`, `c_co2`, `T`, `P`, and `q`.

This is intended to be used in coupled simulations with ClimaCoupler.jl
"""
function ClimaLand.initialize_drivers(
    a::CoupledAtmosphere{FT},
    coords,
) where {FT}
    keys = (:P_liq, :P_snow, :c_co2, :T, :P, :q)
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
    initialize_drivers(driver_tuple::Tuple,
                       coords)

Creates and returns a NamedTuple with the cache variables required by the
model drivers.

If no forcing is required,  driver_tuple is an empty tuple, and an
empty NamedTuple is returned.
"""
function initialize_drivers(driver_tuple::Tuple, coords)
    if isempty(driver_tuple)
        return (;)
    else
        tmp = map(driver_tuple) do (driver)
            nt = initialize_drivers(driver, coords)
        end
        merge(tmp...)
    end

end

"""
    make_update_drivers(driver_tuple)

Creates and returns a function which updates the forcing variables ("drivers"). If no drivers are being used, driver_tuple is empty, and the update
function does nothing.
"""
function make_update_drivers(driver_tuple::Tuple)
    if isempty(driver_tuple)
        return (p, t) -> nothing
    else
        update_driver_list = map(driver_tuple) do (driver)
            make_update_drivers(driver)
        end
        function update_drivers!(p, t)
            for ud! in update_driver_list
                ud!(p, t)
            end
        end
        return update_drivers!
    end
end

"""
    make_update_drivers(a::PrescribedGroundConditions{FT}) where {FT}

Creates and returns a function which updates the driver variables
in the case of a PrescribedGroundConditions.
"""
function make_update_drivers(a::PrescribedGroundConditions{FT}) where {FT}
    function update_drivers!(p, t)
        evaluate!(p.drivers.ψ, a.ψ, t)
        evaluate!(p.drivers.T_ground, a.T, t)
    end
    return update_drivers!
end

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
    make_update_drivers(r::CoupledRadiativeFluxes{FT}) where {FT}

Creates and returns a function which updates the driver variables
in the case of a CoupledRadiativeFluxes.

When `r.θs` is `nothing`, the cosine zenith angle
not changed, and should be updated by the coupler. This differs from the behavior of
`PrescribedRadiativeFluxes`, where the cosine zenith angle is set to `NaN` if `θs` is `nothing`.

Otherwise, the cosine zenith angle is computed using `cos.(r.θs(t, r.start_date))`.
"""
make_update_drivers(r::CoupledRadiativeFluxes{FT, Nothing}) where {FT} =
    (p, t) -> nothing
function make_update_drivers(r::CoupledRadiativeFluxes{FT}) where {FT}
    function update_drivers!(p, t)
        p.drivers.cosθs .= cos.(r.θs(t, r.start_date))
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
        # Next we update the zenith angle and diffuse fraction of light. These are either
        # both required, or neither required.
        if !isnothing(r.θs)
            p.drivers.cosθs .= FT.(cos.(r.θs(t, r.start_date)))
            if !isnothing(r.frac_diff)
                # If the diffuse fraction is a TimeVaryingInput, we use that directly.
                evaluate!(p.drivers.frac_diff, r.frac_diff, t)
            else
                # Otherwise, we use an empirical function that requires atmospheric
                # quantities, thermo params, and zenith angle as input.
                @assert !isnothing(r.thermo_params)
                # We assume that temperature, pressure, and specific humidity have already updated
                # and that PrescribedRadiation is only used with PrescribedAtmos.
                # If this is not the case, this function would use the vaues of T, P, q
                # from the previous step.
                # Since this function is empirical and our steps are O(10 minutes), this is
                # not a big issue.
                @. p.drivers.frac_diff = empirical_diffuse_fraction(
                    t,
                    r.start_date,
                    p.drivers.T,
                    p.drivers.P,
                    p.drivers.q,
                    p.drivers.SW_d,
                    p.drivers.cosθs,
                    r.thermo_params,
                )
            end
        elseif isnothing(r.thermo_params) & isnothing(r.frac_diff)
            # These should not be accessed or used. Set to NaN.
            p.drivers.cosθs .= FT(NaN)
            p.drivers.frac_diff .= FT(NaN)
        else
            @error(
                "Zenith angle is not provided, but either diffuse fraction or thermo params are. If zenith angle is nothing, both of these must be nothing. See doc string for PrescribedRadiativeFluxes"
            )
        end
    end

    return update_drivers!
end

"""
    make_update_drivers(d::PrescribedSoilOrganicCarbon{FT}) where {FT}
Creates and returns a function which updates the driver variables
in the case of a PrescribedSoilOrganicCarbon.
"""
function make_update_drivers(d::PrescribedSoilOrganicCarbon{FT}) where {FT}
    function update_drivers!(p, t)
        evaluate!(p.drivers.soc, d.soc, t)
    end
    return update_drivers!
end

"""
     prescribed_forcing_era5(era5_ncdata_path,
                             surface_space,
                             start_date,
                             earth_param_set,
                             FT;
                             gustiness=1,
                             max_wind_speed = nothing,
                             c_co2 = TimeVaryingInput((t) -> 4.2e-4),
                             time_interpolation_method = LinearInterpolation(PeriodicCalendar()),
                             regridder_type = :InterpolationsRegridder)

A helper function which constructs the `PrescribedAtmosphere` and `PrescribedRadiativeFluxes`
from a file path pointing to the ERA5 data in a netcdf file, the surface_space, the start date,
and the earth_param_set.

The argument `era5_ncdata_path` is either a list of nc files, each with all of the variables required, but with different time intervals in the different files, or else it is a single file with all the variables.

########## WARNING ##########

High wind speed anomalies (10-100x increase and decrease over a period of a several hours) appear in the ERA5
reanalysis data. These generate very large surface fluxes (due to wind speeds up to 300 m/s), which lead to instability. The kwarg max_wind_speed,
with a value give in m/s,
is used to clip these if it is not `nothing`.
See: https://confluence.ecmwf.int/display/CKB/ERA5%3A+large+10m+winds
"""
function prescribed_forcing_era5(
    era5_ncdata_path,
    surface_space,
    start_date,
    earth_param_set,
    FT;
    gustiness = 1,
    max_wind_speed = nothing,
    c_co2 = TimeVaryingInput((t) -> 4.2e-4),
    time_interpolation_method = LinearInterpolation(PeriodicCalendar()),
    regridder_type = :InterpolationsRegridder,
)
    # Pass a list of files in all cases
    era5_ncdata_path isa String && (era5_ncdata_path = [era5_ncdata_path])

    # Precip is provide as a mass flux; convert to volume flux of liquid water with ρ =1000 kg/m^3
    precip = TimeVaryingInput(
        [era5_ncdata_path, era5_ncdata_path],
        ["mtpr", "msr"],
        surface_space;
        start_date,
        regridder_type,
        file_reader_kwargs = (; preprocess_func = (data) -> -data / 1000,),
        method = time_interpolation_method,
        compose_function = (mtpr, msr) -> min.(mtpr .- msr, Float32(0)),
    )
    # Precip is provide as a mass flux; convert to volume flux of liquid water with ρ =1000 kg/m^3
    snow_precip = TimeVaryingInput(
        era5_ncdata_path,
        "msr",
        surface_space;
        start_date,
        regridder_type,
        file_reader_kwargs = (; preprocess_func = (data) -> -data / 1000,),
        method = time_interpolation_method,
    )
    if max_wind_speed isa Nothing
        wind_function = (u, v) -> sqrt.(u .^ 2 .+ v .^ 2)
    else
        wind_function = (u, v) -> min.(sqrt.(u .^ 2 .+ v .^ 2), max_wind_speed)
    end

    u_atmos = TimeVaryingInput(
        [era5_ncdata_path, era5_ncdata_path],
        ["u10", "v10"],
        surface_space;
        start_date,
        regridder_type,
        compose_function = wind_function,
        method = time_interpolation_method,
    )
    specific_humidity(Td, T, P; params = earth_param_set) =
        ClimaLand.specific_humidity_from_dewpoint.(Td, T, P, params)
    q_atmos = TimeVaryingInput(
        [era5_ncdata_path, era5_ncdata_path, era5_ncdata_path],
        ["d2m", "t2m", "sp"],
        surface_space;
        start_date,
        regridder_type,
        compose_function = specific_humidity,
        method = time_interpolation_method,
    )
    P_atmos = TimeVaryingInput(
        era5_ncdata_path,
        "sp",
        surface_space;
        start_date,
        regridder_type,
        method = time_interpolation_method,
    )

    T_atmos = TimeVaryingInput(
        era5_ncdata_path,
        "t2m",
        surface_space;
        start_date,
        regridder_type,
        method = time_interpolation_method,
    )
    h_atmos = FT(10)

    atmos = PrescribedAtmosphere(
        precip,
        snow_precip,
        T_atmos,
        u_atmos,
        q_atmos,
        P_atmos,
        start_date,
        h_atmos,
        earth_param_set;
        gustiness = FT(gustiness),
        c_co2 = c_co2,
    )

    SW_d = TimeVaryingInput(
        era5_ncdata_path,
        "msdwswrf",
        surface_space;
        start_date,
        regridder_type,
        method = time_interpolation_method,
    )
    function compute_diffuse_fraction(total, direct)
        diff = max(total - direct, Float32(0))
        return min(diff / (total + eps(Float32)), Float32(1))
    end
    function compute_diffuse_fraction_broadcasted(total, direct)
        return @. compute_diffuse_fraction(total, direct)
    end

    frac_diff = TimeVaryingInput(
        [era5_ncdata_path, era5_ncdata_path],
        ["msdwswrf", "msdrswrf"],
        surface_space;
        start_date,
        regridder_type,
        method = time_interpolation_method,
        compose_function = compute_diffuse_fraction_broadcasted,
    )
    LW_d = TimeVaryingInput(
        era5_ncdata_path,
        "msdwlwrf",
        surface_space;
        start_date,
        regridder_type,
        method = time_interpolation_method,
    )
    zenith_angle =
        (t, s) -> default_zenith_angle(
            t,
            s;
            latitude = ClimaCore.Fields.coordinate_field(surface_space).lat,
            longitude = ClimaCore.Fields.coordinate_field(surface_space).long,
            insol_params = earth_param_set.insol_params,
        )

    radiation = PrescribedRadiativeFluxes(
        FT,
        SW_d,
        LW_d,
        start_date;
        θs = zenith_angle,
        earth_param_set = earth_param_set,
        frac_diff = frac_diff,
    )
    return (; atmos, radiation)
end


"""
     prescribed_lai_era5(era5_lai_ncdata_path,
                         era5_lai_cover_ncdata_path,
                         surface_space,
                         start_date,
                         earth_param_set;
                         time_interpolation_method = LinearInterpolation(PeriodicCalendar()),
                         regridder_type = :InterpolationsRegridder)

A helper function which constructs the TimeVaryingInput object for Leaf Area Index, from a
file path pointing to the ERA5 LAI data in a netcdf file, a file path pointing to the ERA5
LAI cover data in a netcdf file, the surface_space, the start date, and the earth_param_set.

This currently one works when a single file is passed for both the era5 lai and era5 lai cover data.
"""
function prescribed_lai_era5(
    era5_lai_ncdata_path,
    era5_lai_cover_ncdata_path,
    surface_space,
    start_date;
    time_interpolation_method = LinearInterpolation(PeriodicCalendar()),
    regridder_type = :InterpolationsRegridder,
)
    hvc_ds = NCFileReader(era5_lai_cover_ncdata_path, "cvh")
    lvc_ds = NCFileReader(era5_lai_cover_ncdata_path, "cvl")
    hv_cover = read(hvc_ds)
    lv_cover = read(lvc_ds)
    close(hvc_ds)
    close(lvc_ds)
    compose_function = let hv_cover = hv_cover, lv_cover = lv_cover
        (lai_hv, lai_lv) -> lai_hv .* hv_cover .+ lai_lv .* lv_cover
    end
    return TimeVaryingInput(
        era5_lai_ncdata_path,
        ["lai_hv", "lai_lv"],
        surface_space;
        start_date,
        regridder_type,
        method = time_interpolation_method,
        compose_function = compose_function,
    )
end

"""
     prescribed_lai_modis(modis_lai_ncdata_path,
                          surface_space,
                          start_date,
                          earth_param_set;
                          time_interpolation_method =
                                        LinearInterpolation(PeriodicCalendar()))
                          regridder_type = :InterpolationsRegridder)

A helper function which constructure the TimeVaryingInput object for Lead Area
Index from a file path pointint to the MODIS LAI data in a netcdf file, the
surface_space, the start date, and the earth_param_set.
"""
function prescribed_lai_modis(
    modis_lai_ncdata_path,
    surface_space,
    start_date;
    time_interpolation_method = LinearInterpolation(PeriodicCalendar()),
    regridder_type = :InterpolationsRegridder,
)
    return TimeVaryingInput(
        modis_lai_ncdata_path,
        ["lai"],
        surface_space;
        start_date,
        regridder_type,
        method = time_interpolation_method,
    )
end

"""
     prescribed_analytic_forcing(FT = Float32;
                                 earth_param_set = LP.LandParameters(FT),
                                 start_date = DateTime(2005),
                                 SW_d = (t) -> 0,
                                 LW_d = (t) -> 5.67e-8 * 280.0^4.0,
                                 precip = (t) -> 0, # no precipitation
                                 T_atmos = (t) -> 280.0,
                                 u_atmos = (t) -> 1.0,
                                 q_atmos = (t) -> 0.0, # no atmos water
                                 h_atmos = FT(1e-8),
                                 P_atmos = (t) -> 101325,
                                 atmos = PrescribedAtmosphere(
                                     TimeVaryingInput(precip),
                                     TimeVaryingInput(precip),
                                     TimeVaryingInput(T_atmos),
                                     TimeVaryingInput(u_atmos),
                                     TimeVaryingInput(q_atmos),
                                     TimeVaryingInput(P_atmos),
                                     start_date,
                                     h_atmos,
                                     earth_param_set,
                                 ),
                                 radiation = PrescribedRadiativeFluxes(
                                     FT,
                                     TimeVaryingInput(SW_d),
                                     TimeVaryingInput(LW_d),
                                     start_date,
                                 ),
                             )

A helper function which constructs the `PrescribedAtmosphere` and `PrescribedRadiativeFluxes`
for a simple analytic case.
"""
function prescribed_analytic_forcing(
    FT = Float32;
    earth_param_set = LP.LandParameters(FT),
    start_date = DateTime(2005),
    SW_d = (t) -> 0,
    LW_d = (t) -> 5.67e-8 * 280.0^4.0,
    precip = (t) -> 0, # no precipitation
    T_atmos = (t) -> 280.0,
    u_atmos = (t) -> 1.0,
    q_atmos = (t) -> 0.0, # no atmos water
    h_atmos = FT(1e-8),
    P_atmos = (t) -> 101325,
    atmos = PrescribedAtmosphere(
        TimeVaryingInput(precip),
        TimeVaryingInput(precip),
        TimeVaryingInput(T_atmos),
        TimeVaryingInput(u_atmos),
        TimeVaryingInput(q_atmos),
        TimeVaryingInput(P_atmos),
        start_date,
        h_atmos,
        earth_param_set,
    ),
    radiation = PrescribedRadiativeFluxes(
        FT,
        TimeVaryingInput(SW_d),
        TimeVaryingInput(LW_d),
        start_date,
    ),
)
    return atmos, radiation
end


"""
    empirical_diffuse_fraction(td::FT, T::FT, P, q, SW_d::FT, cosθs::FT, thermo_params) where {FT}

Computes the fraction of diffuse radiation (`diff_frac`) as a function
of the solar zenith angle (`θs`), the total surface downwelling shortwave radiation (`SW_d`),
the air temperature (`T`), air pressure (`P`), specific humidity (`q`), and the day of the year
(`td`).

See Appendix A of Braghiere, "Evaluation of turbulent fluxes of CO2, sensible heat,
and latent heat as a function of aerosol optical depth over the course of deforestation
in the Brazilian Amazon" 2013.

Note that cos(θs) is equal to zero when θs = π/2, and this is a coefficient
of k₀, which we divide by in this expression. This can amplify small errors
when θs is near π/2.

This formula is empirical and can yied negative numbers depending on the
input, which, when dividing by something very near zero,
can become large negative numbers.

Because of that, we cap the returned value to lie within [0,1].
"""
function empirical_diffuse_fraction(
    t,
    start_date,
    T::FT,
    P::FT,
    q::FT,
    SW_d::FT,
    cosθs::FT,
    thermo_params,
) where {FT}
    if t isa ITime
        DOY = Dates.dayofyear(date(t))
    else
        DOY = Dates.dayofyear(start_date + Dates.Second(floor(Int64, t)))
    end
    RH = ClimaLand.relative_humidity(T, P, q, thermo_params)
    k₀ = FT(1370 * (1 + 0.033 * cos(2π * DOY / 365))) * cosθs
    kₜ = SW_d / k₀
    if kₜ ≤ 0.3
        diff_frac = FT(
            kₜ * (1 - 0.232 * kₜ + 0.0239 * cosθs - 6.82e-4 * T + 0.0195 * RH),
        )
    elseif kₜ ≤ 0.78
        diff_frac = FT(
            kₜ *
            (1.329 - 1.716 * kₜ + 0.267 * cosθs - 3.57e-3 * T + 0.106 * RH),
        )
    else
        diff_frac =
            FT(kₜ * (0.426 * kₜ + 0.256 * cosθs - 3.49e-3 * T + 0.0734 * RH))
    end
    return max(min(diff_frac, FT(1)), FT(0))
end
