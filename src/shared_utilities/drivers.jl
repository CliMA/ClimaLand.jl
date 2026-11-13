import ClimaUtilities.TimeVaryingInputs:
    TimeVaryingInput,
    AbstractTimeVaryingInput,
    LinearInterpolation,
    PeriodicCalendar,
    Flat
import Interpolations: Constant
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
import SurfaceFluxes.UniversalFunctions as UF
using StaticArrays
using NVTX
import ..Parameters as LP
export AbstractAtmosphericDrivers,
    AbstractRadiativeDrivers,
    AbstractAtmosphericDrivers,
    PrescribedAtmosphere,
    PrescribedPrecipitation,
    PrescribedGroundConditions,
    CoupledAtmosphere,
    PrescribedRadiativeFluxes,
    CoupledRadiativeFluxes,
    turbulent_fluxes!,
    net_radiation!,
    turbulent_fluxes_at_a_point,
    make_update_drivers,
    prescribed_forcing_era5,
    prescribed_perturbed_temperature_era5,
    prescribed_perturbed_rh_era5,
    prescribed_analytic_forcing,
    default_cos_zenith_angle,
    prescribed_forcing_crujra,
    CycledYearsInput,
    prescribed_co2_annual,
    read_annual_co2

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
    PrescribedAtmosphere{FT, CA, DT} <: AbstractAtmosphericDrivers{FT}

Container for holding prescribed atmospheric drivers and other
information needed for computing turbulent surface fluxes when
driving land models in standalone mode.

The default CO2 concentration is a constant as a function of time, equal to
4.2e-4 mol/mol. Since not all models require CO2 concentration, the 
default for that is `nothing`.

The atmospheric height passed is, by default, the height of the forcing above
the ground surface: `h_atmos = z_atmos - z_ground`. In this case, the reference
level is `ReferenceAboveSurface`. However, in some cases, it is more natural
to provide the atmospheric forcing at the height relative to the apparent sink for
momentum, which lies at a height `d+z_0m` above the ground surface, where 
`d` is the displacement height 
and `z_0m` is the roughness length for momentum. If a canopy is present,
these refer to the canopy, but if only bare ground is present, the
displacement height and roughness are that of the bare ground. In that case,
one must change the reference to be `ReferenceAboveApparentSink`.

The keyword `gustiness` is a SurfaceFluxes gustiness model, or a number,
which is stored as a `SurfaceFluxes.ConstantGustinessSpec` with that 
constant minimum wind speed [m/s]; the default is a
`SurfaceFluxes.FlooredDeardorffGustinessSpec` with a floor of 1 m/s, the
larger of that floor and the Deardorff convective gustiness.
$(DocStringExtensions.FIELDS)
"""
struct PrescribedAtmosphere{
    FT,
    LP <: AbstractTimeVaryingInput,
    SP <: AbstractTimeVaryingInput,
    TA <: AbstractTimeVaryingInput,
    UA <: AbstractTimeVaryingInput,
    QA <: AbstractTimeVaryingInput,
    RA <: AbstractTimeVaryingInput,
    CA <: AbstractTimeVaryingInput,
    DT,
    G,
    RL,
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
    "Height of the forcing (m), above the ground or above the apparent sink for momentum, according to `reference`"
    h::FT
    "Gustiness model of the effective wind speed (a `SurfaceFluxes.AbstractGustinessSpec`)"
    gustiness::G
    "Convention for the height `h`: `SurfaceFluxes.ReferenceAboveSurface()` (the default) or `SurfaceFluxes.ReferenceAboveApparentSink()`"
    reference::RL
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
        toml_dict::CP.ParamDict;
        gustiness = SurfaceFluxes.FlooredDeardorffGustinessSpec(FT(1)),
        reference = SurfaceFluxes.ReferenceAboveSurface(),
        c_co2 = TimeVaryingInput((t) -> 4.2e-4),
    ) where {FT}
        earth_param_set = LP.LandParameters(toml_dict)
        thermo_params = LP.thermodynamic_parameters(earth_param_set)
        args = (liquid_precip, snow_precip, T, u, q, P, c_co2, start_date)
        gustiness isa Number &&
            (gustiness = SurfaceFluxes.ConstantGustinessSpec(FT(gustiness)))
        return new{
            typeof(h),
            typeof.(args)...,
            typeof(gustiness),
            typeof(reference),
            typeof(thermo_params),
        }(
            args...,
            h,
            gustiness,
            reference,
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

To be used when coupling to an atmosphere model. Either both `cosθs` and `start_date`
must be `nothing`, or both must not be `nothing``.

During the driver update, cosθs is unchanged if `cosθs` is `nothing`. This behavior differs from
the `PrescribedRadiativeFluxes` where `cosθs` set to `NaN` if `cosθs` is `nothing`.
Otherwise, `cosθs` recieves the following arguments:
(`time_from_start`, `start_date`), and is expected to return zenith angle at the given time.
$(DocStringExtensions.FIELDS)
"""
struct CoupledRadiativeFluxes{FT, F <: Union{Function, Nothing}, T} <:
       AbstractRadiativeDrivers{FT}
    """Function that fills a climacore field with the zenith angle given the following arguments:
    (time_from_start, `start_date`)"""
    cosθs::F
    "Start date - the datetime corresponding to t=0 for the simulation"
    start_date::T
    function CoupledRadiativeFluxes{FT, F, T}(
        cosθs::F,
        start_date::T,
    ) where {FT, F, T}
        (
            (isnothing(cosθs) && isnothing(start_date)) ||
            ((!isnothing(cosθs) && !isnothing(start_date)))
        ) || error(
            "CoupledRadiativeFluxes: `cosθs` and start_date` must both be `nothing` or both not `nothing`.",
        )
        new{FT, F, T}(cosθs, start_date)
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
        toml_dict::CP.ParamDict,
    )

Creates a `CoupledRadiativeFluxes` object with a default zenith angle function that uses Insolation.jl
to compute the zenith angle at a given time and location.
"""
function CoupledRadiativeFluxes{FT}(
    start_date::DT;
    latitude::LT,
    longitude::LT,
    toml_dict::CP.ParamDict,
) where {FT, DT, LT}
    insol_params = LP.LandParameters(toml_dict).insol_params
    cos_zenith_angle =
        (t, s) -> default_cos_zenith_angle(
            t,
            s;
            latitude = latitude,
            longitude = longitude,
            insol_params = insol_params,
        )
    return CoupledRadiativeFluxes{FT, typeof(cos_zenith_angle), DT}(
        cos_zenith_angle,
        start_date,
    )
end

"""
    default_cos_zenith_angle(
        t::T,
        start_date::Dates.DateTime;
        latitude::LT,
        longitude::LT,
        insol_params::Insolation.Parameters.InsolationParameters{FT},
    )

Calculate cosine of the zenith angle with Insolation for the given start date, insolation parameters, latitude,
and longitude. Note that Insolation.jl returns 0 for the cosine when the
sun is below the horizon.

`latitude` and `longitude` can be a collections or a Number.
"""
function default_cos_zenith_angle(
    t::T,
    start_date::Dates.DateTime;
    latitude::LT,
    longitude::LT,
    insol_params,
) where {T, LT}
    current_datetime = if T <: ITime
        # ITime may not have an epoch, so use start_date as fallback
        isnothing(t.epoch) ? start_date + t.counter * t.period : date(t)
    else
        start_date + Dates.Second(round(t))
    end

    # Reduces allocations by throwing away unwanted values
    zenith_only = (args...) -> Insolation.insolation(args...).μ
    return @. lazy(
        zenith_only(current_datetime, latitude, longitude, insol_params),
    )
end

"""
    get_solar_noon_zenith(day_utc_noon::DateTime, lat::FT, lon::FT, params)::FT where {FT<:AbstractFloat}

Utility function for using Insolation.jl to get the solar zenith angle at a particular coordinate's
solar noon, deriving its solar noon time from UTC noon of that day using the longitude.
"""
function get_solar_noon_zenith(
    day_utc_noon::Dates.DateTime,
    lat::FT,
    lon::FT,
    params,
)::FT where {FT <: AbstractFloat}
    # Get the UTC time of solar noon: 1 degree of longitude is 4 minutes offset
    # (There is an EquationOfTime component to this, but it will affect the output negligbly)
    coord_noon = day_utc_noon - Minute(round(Int, FT(4) * lon))
    return FT(Insolation.insolation(coord_noon, lat, lon, params).μ)
end

"""
    solar_noon_zenith_cosine(
        t::T,
        start_date::Dates.DateTime;
        latitude::LT,
        longitude::LT,
        insol_params::Insolation.Parameters.InsolationParameters{FT},
    )

Calculate cosine of the zenith angle with Insolation for solar noon of the given start date, insolation parameters, latitude,
and longitude. Note that Insolation.jl returns 0 for the cosine when the
sun is below the horizon.

`lat` (latitude) and `lon` (longitude) can be a collections or a Number.
"""
function solar_noon_zenith_cosine(
    t::T,
    start_date::Dates.DateTime;
    lat::LT,
    lon::LT,
    insol_params,
) where {T, LT}
    current_datetime = if T <: ITime
        # ITime may not have an epoch, so use start_date as fallback
        isnothing(t.epoch) ? start_date + t.counter * t.period : date(t)
    else
        start_date + Dates.Second(round(t))
    end
    day_utc_noon = DateTime(
        year(current_datetime),
        month(current_datetime),
        day(current_datetime),
        12,
    )
    return @. lazy(get_solar_noon_zenith(day_utc_noon, lat, lon, insol_params))
end

"""
    CoupledAtmosphere{FT, T <: Union{FT, Fields.Field}, RL} <: AbstractAtmosphericDrivers{FT}
    CoupledAtmosphere{FT, T}(h, gustiness; reference = SurfaceFluxes.ReferenceAboveApparentSink())

To be used when coupling to an atmosphere model. Contains fields that
are used to compute surface fluxes in the coupled setup.

The height `h` of the atmospheric state is measured according to `reference`:
by default from the apparent sink for momentum, `d + z_0m` above the ground,
the convention of an atmosphere model that does not resolve the canopy; pass
`reference = SurfaceFluxes.ReferenceAboveSurface()` for a height above the
ground.

When constructed without a space, the struct doesn't contain anything,
but it still acts as a flag that fluxes have been updated by the coupler
and don't need to be recomputed.
When constructed with a space, the struct contains the fields needed to compute
surface fluxes in the coupled setup, which are accessed by ClimaCoupler.

The gustiness is a number, the minimum wind speed [m/s] that the coupler uses
in its flux computation; land models driven by a `CoupledAtmosphere` inside
ClimaLand use it as the floor of a
`SurfaceFluxes.FlooredDeardorffGustinessSpec` (see `gustiness_spec`).
"""
struct CoupledAtmosphere{FT, T <: Union{FT, Fields.Field}, RL} <:
       AbstractAtmosphericDrivers{FT}
    "Height of the atmospheric state (m), according to `reference`"
    h::T
    "Convention for the height `h` (a `SurfaceFluxes.AbstractReferenceLevel`)"
    reference::RL
    "Minimum wind speed (gustiness; m/s), which is always a spatial constant"
    gustiness::FT
end

function CoupledAtmosphere{FT, T}(
    h::T,
    gustiness::FT;
    reference = SurfaceFluxes.ReferenceAboveApparentSink(),
) where {FT, T}
    return CoupledAtmosphere{FT, T, typeof(reference)}(h, reference, gustiness)
end

"""
    gustiness_spec(atmos::AbstractAtmosphericDrivers)

Return the gustiness model of an atmospheric driver as a
`SurfaceFluxes.AbstractGustinessSpec`: the `gustiness` field of a
`PrescribedAtmosphere`, and a `SurfaceFluxes.FlooredDeardorffGustinessSpec` with
the minimum wind speed that a `CoupledAtmosphere` holds for the coupler.
"""
gustiness_spec(atmos::PrescribedAtmosphere) = atmos.gustiness
gustiness_spec(atmos::CoupledAtmosphere) =
    SurfaceFluxes.FlooredDeardorffGustinessSpec(atmos.gustiness)

"""
    compute_ρ_sfc(surface_flux_params, T_air, P_air, q_air, Δh, T_sfc, q_vap_sfc = nothing)

Computes the density of air at the displacement height, where SurfaceFluxes
places the surface state, given the temperature at the surface T_sfc, the
atmospheric temperature T_air, pressure P_air, specific humidity q_air, and the
height `Δh` of the atmospheric state above the displacement height (see
`compute_Δz_eff`), using surface_flux_params. The hydrostatic extrapolation
uses the surface specific humidity `q_vap_sfc` where given and `q_air`
otherwise. A negative `T_sfc` gives `NaN`, which the NaN check reports, where
the extrapolation would give a finite, unphysical density.
"""
function compute_ρ_sfc(
    surface_flux_params,
    T_air,
    P_air,
    q_air,
    Δh,
    T_sfc,
    q_vap_sfc = nothing,
)
    T_sfc = T_sfc < 0 ? oftype(T_sfc, NaN) : T_sfc
    thermo_params =
        SurfaceFluxes.Parameters.thermodynamics_params(surface_flux_params)
    ρ_air = Thermodynamics.air_density(thermo_params, T_air, P_air, q_air)
    q_liq_air = zero(q_air)
    q_ice_air = zero(q_air)
    return SurfaceFluxes.surface_density(
        surface_flux_params,
        T_air,
        ρ_air,
        T_sfc,
        Δh,
        q_air,
        q_liq_air,
        q_ice_air,
        q_vap_sfc,
    )
end

return_momentum_fluxes(atmos::PrescribedAtmosphere) = false
return_momentum_fluxes(atmos::CoupledAtmosphere) = true

"""
    turbulent_fluxes!(dest,
                      atmos::AbstractAtmosphericDrivers,
                      model::AbstractModel,
                      Y,
                      p,
                      t;
                      h_atmos = atmos.h,
                      u_atmos = p.drivers.u,
                      T_atmos = p.drivers.T,
                      q_atmos = p.drivers.q,
                      gustiness = gustiness_spec(atmos),
                      reference = atmos.reference
                      )

Computes the turbulent surface flux terms at the ground,
including turbulent energy fluxes as well as the water vapor flux
(in units of m^3/m^2/s of water).
Positive fluxes indicate flow from the ground to the atmosphere.

It solves for these given atmospheric conditions,
model parameters, and the surface conditions. The elements of `dest` select,
by name, which of the quantities of `turbulent_fluxes_at_a_point` are stored
(see `select_fluxes`).

The reference height `h_atmos`, the wind `u_atmos`, temperature `T_atmos`, and
specific humidity `q_atmos` at it, and the gustiness model (a
`SurfaceFluxes.AbstractGustinessSpec`, or a field of them) default to those of
the atmospheric forcing (`gustiness_spec`). The `reference` is the convention for
the atmosphere height (relative to surface, or relative to the apparent sink for momentum).

Integrated models pass the
sub-canopy reference height, attenuated wind, canopy-air temperature and
humidity, and gustiness model for the surfaces beneath a canopy (see
`Canopy.subcanopy_forcing`).
"""
function turbulent_fluxes!(
    dest,
    atmos::AbstractAtmosphericDrivers,
    model::AbstractModel,
    Y,
    p,
    t;
    h_atmos = atmos.h,
    u_atmos = p.drivers.u,
    T_atmos = p.drivers.T,
    q_atmos = p.drivers.q,
    gustiness = gustiness_spec(atmos),
    reference = atmos.reference,
)

    T_sfc = component_temperature(model, Y, p) # guess
    # Guess
    q_sfc = component_specific_humidity(
        model,
        Y,
        p,
        (; h_atmos, T_atmos, q_atmos, reference),
    )
    roughness_model = surface_roughness_model(model, Y, p)
    update_T_sfc = get_update_surface_temperature_function(model, Y, p)
    update_q_sfc = get_update_surface_humidity_function(model, Y, p)
    z_sfc = surface_elevation(model, Y, p)
    displ = surface_displacement_height(model, Y, p)
    update_∂T_sfc∂T = get_∂T_sfc∂T_function(model, Y, p)
    update_∂q_sfc∂T = get_∂q_sfc∂T_function(model, Y, p)
    earth_param_set = get_earth_param_set(model)
    stored = Val(fieldnames(eltype(dest)))
    dest .= turbulent_fluxes_at_a_point.(
        stored,
        p.drivers.P,
        T_atmos,
        q_atmos, # q_tot
        u_atmos,
        h_atmos,
        T_sfc,
        q_sfc,
        roughness_model,
        update_T_sfc,
        update_q_sfc,
        z_sfc,
        displ,
        update_∂T_sfc∂T,
        update_∂q_sfc∂T,
        gustiness,
        reference,
        earth_param_set,
    )
    return nothing
end

"""
    select_fluxes(::Val{names}, fluxes)

Return the NamedTuple of the fields `names` of the NamedTuple `fluxes` returned
by `turbulent_fluxes_from_output`. Each model stores the subset it needs in its
cache: all store the energy and vapor fluxes and their temperature derivatives,
and the soil, snow, and canopy also store the surface state and similarity
scales of the flux solve for the screen-level diagnostics.
"""
@inline select_fluxes(::Val{names}, fluxes) where {names} =
    NamedTuple{names}(fluxes)

"""
    flux_names(return_extra_fluxes::Val)

Return, as a `Val`, the names of the quantities that
`turbulent_fluxes_at_a_point` returns for `return_extra_fluxes`: the energy
fluxes `lhf` and `shf`, the vapor flux `vapor_flux`, and the temperature
derivatives `∂lhf∂T` and `∂shf∂T`, followed for `Val(true)` by the
momentum fluxes `ρτxz` and `ρτyz` and the buoyancy flux `buoyancy_flux`.
"""
flux_names(::Val{false}) = Val((:lhf, :shf, :vapor_flux, :∂lhf∂T, :∂shf∂T))
flux_names(::Val{true}) = Val((
    :lhf,
    :shf,
    :vapor_flux,
    :∂lhf∂T,
    :∂shf∂T,
    :ρτxz,
    :ρτyz,
    :buoyancy_flux,
))

"""
    flux_tuple_type(::Val{names}, ::Type{FT})

Return the `NamedTuple` type with fields `names`, all of type `FT`; the
element type of the cache variables holding the turbulent fluxes.
"""
flux_tuple_type(::Val{names}, ::Type{FT}) where {names, FT} =
    NamedTuple{names, NTuple{length(names), FT}}
"""
    turbulent_fluxes_at_a_point(return_extra_fluxes, P_atmos, T_atmos, q_tot_atmos,
                                u_atmos, h_atmos, T_sfc_guess, q_vap_sfc_guess,
                                roughness_model, update_T_sfc, update_q_vap_sfc, z_sfc,
                                displ, update_∂T_sfc∂T, update_∂q_sfc∂T, gustiness, reference,
                                earth_param_set)
    turbulent_fluxes_at_a_point(stored::Val{names}, args...)

Computes turbulent surface fluxes at a point on a surface given
(1) the prescribed atmospheric conditions, `P_atmos`, `T_atmos`, `q_tot_atmos`,
    `u_atmos`, and `h_atmos` the relative height at which these measurements are made, according
    to `reference`.
    `u_atmos` is a can be a scalar if only the horizontal wind speed is known,
    otherwise a tuple can be passed with the (u,v) components of the wind.
(2) the surface temperature (`T_sfc_guess`), specific humidity (`q_vap_sfc_guess`).
    If no `update_T_sfc` and `update_q_val_sfc` functions are passed, these ``guesses"
    are the values used in the solve. Otherwise, the update functions are used to
    modify the ``guess" values to the actual surface values.
(3) Other surface properties, such the roughness model (see SurfaceFluxes.jl
    for possible types), the ground elevation,
    the displacement height.
(4) functions for computing the partial derivatives of `T_sfc` and `q_vap_sfc` with respect
    to the actual temperature (of the leaves, soil, or snow). These functions are component
    specific.
(5) the parameter set.

With `return_extra_fluxes = Val(false)`, this returns the NamedTuple
`(; lhf, shf, vapor_flux, ∂lhf∂T, ∂shf∂T)` of the energy fluxes, the
liquid water volume flux, and the derivatives of the energy fluxes with respect to
the component temperature; with `Val(true)`, it also returns the momentum flux
components in the horizontal directions, `ρτxz` and `ρτyz`, and the buoyancy
flux `buoyancy_flux` (see `flux_names`). ClimaCoupler evaluates this method
directly into its flux fields, with the `reference` of its `CoupledAtmosphere`.
The form without `reference` is deprecated; it measures `h_atmos` from the
ground (`SurfaceFluxes.ReferenceAboveSurface`).

With a `Val` of a tuple of names, it returns the fields `names` of the full
NamedTuple of `turbulent_fluxes_from_output`, which also holds the surface
temperature `T_sfc` and specific humidity `q_sfc` at which the fluxes were
evaluated and the friction velocity `ustar`, stability parameter `ζ`, and
effective height `Δz_eff` of the Monin-Obukhov profiles. `turbulent_fluxes!`
passes the field names of its destination, so that each model stores the
subset it needs; space for the extra fluxes is only allocated in the cache
when running with a `CoupledAtmosphere`.
"""
@inline turbulent_fluxes_at_a_point(
    return_extra_fluxes::Union{Val{true}, Val{false}},
    args...,
) = turbulent_fluxes_at_a_point(flux_names(return_extra_fluxes), args...)

@inline function turbulent_fluxes_at_a_point(
    stored::Val{names},
    P_atmos,
    T_atmos,
    q_tot_atmos,
    u_atmos,
    h_atmos,
    T_sfc_guess,
    q_vap_sfc_guess,
    roughness_model,
    update_T_sfc::UT,
    update_q_vap_sfc::UQ,
    z_sfc,
    displ,
    update_∂T_sfc∂T::UDT,
    update_∂q_sfc∂T::UDQ,
    gustiness,
    reference,
    earth_param_set,
) where {names, UT, UQ, UDT, UDQ}
    output = surface_fluxes_at_a_point(
        T_sfc_guess,
        q_vap_sfc_guess,
        update_T_sfc,
        update_q_vap_sfc,
        P_atmos,
        T_atmos,
        q_tot_atmos,
        u_atmos,
        h_atmos,
        z_sfc,
        displ,
        roughness_model,
        gustiness,
        reference,
        earth_param_set,
    )
    # The momentum and buoyancy fluxes are cheap, so they are evaluated for
    # every selection and kept only where named
    fluxes = turbulent_fluxes_from_output(
        Val(true),
        output,
        T_sfc_guess,
        q_vap_sfc_guess,
        update_∂T_sfc∂T,
        update_∂q_sfc∂T,
        P_atmos,
        T_atmos,
        q_tot_atmos,
        h_atmos,
        displ,
        roughness_model,
        reference,
        earth_param_set,
    )
    return select_fluxes(stored, fluxes)
end

# Deprecated form without `reference`, which predates the reference conventions
# and is called by the ClimaCoupler bucket: `h_atmos` is the height above the
# ground, as it was before `reference` was added
@inline turbulent_fluxes_at_a_point(
    return_extra_fluxes::Union{Val{true}, Val{false}},
    P_atmos,
    T_atmos,
    q_tot_atmos,
    u_atmos,
    h_atmos,
    T_sfc_guess,
    q_vap_sfc_guess,
    roughness_model,
    update_T_sfc::UT,
    update_q_vap_sfc::UQ,
    z_sfc,
    displ,
    update_∂T_sfc∂T::UDT,
    update_∂q_sfc∂T::UDQ,
    gustiness,
    earth_param_set,
) where {UT, UQ, UDT, UDQ} = turbulent_fluxes_at_a_point(
    return_extra_fluxes,
    P_atmos,
    T_atmos,
    q_tot_atmos,
    u_atmos,
    h_atmos,
    T_sfc_guess,
    q_vap_sfc_guess,
    roughness_model,
    update_T_sfc,
    update_q_vap_sfc,
    z_sfc,
    displ,
    update_∂T_sfc∂T,
    update_∂q_sfc∂T,
    gustiness,
    SurfaceFluxes.ReferenceAboveSurface(),
    earth_param_set,
)

"""
    surface_flux_config(roughness_model, gustiness, reference)

Return the SurfaceFluxes configuration used for the turbulent fluxes of all
land surfaces: the given roughness and gustiness models, moist thermodynamics,
no roughness sublayer correction, the `MaxHeatFluxStabilityCap` in stable
conditions, and a type indicating the reference convention for atmospheric height.

The stability cap holds the exchange coefficients at their values at the
stability `ζ_p` at which the Monin-Obukhov sensible heat flux at fixed wind
speed is maximal. `ζ_p` depends only on the ratio of the effective forcing
height to the momentum roughness length: about 0.15–0.4 over tall forests,
0.4–0.8 over grass, 1–1.2 over bare soil, and 1.4–1.6 over snow.
Without a cap, the Monin-Obukhov exchange collapses once the bulk Richardson
number becomes supercritical, which decouples the surface from the
atmosphere at night (runaway cooling). The cap has no free parameters. See the
SurfaceFluxes documentation.
"""
@inline surface_flux_config(roughness_model, gustiness, reference) =
    SurfaceFluxes.SurfaceFluxConfig(
        roughness_model,
        gustiness,
        SurfaceFluxes.MoistModel(),
        SurfaceFluxes.NoRoughnessSubLayer(),
        SurfaceFluxes.MaxHeatFluxStabilityCap(),
        reference,
    )

"""
    surface_fluxes_at_a_point(T_sfc_guess, q_vap_sfc_guess, update_T_sfc, update_q_vap_sfc,
                              P_atmos, T_atmos, q_tot_atmos, u_atmos, h_atmos, z_sfc, displ,
                              roughness_model, gustiness, reference, earth_param_set)

Solve the Monin-Obukhov similarity equations with SurfaceFluxes.jl at a point
and return its output. The surface temperature and specific humidity start
from `T_sfc_guess` and `q_vap_sfc_guess` and are updated within the iterations
by the callbacks `update_T_sfc` and `update_q_vap_sfc` (or held fixed if these
are `nothing`). The atmospheric pressure, temperature, specific humidity, and
wind (a speed or a horizontal vector) are given at the height according to the `reference`,
`h_atmos`, and the ground is at elevation `z_sfc` with displacement height
`displ`.

Called from `turbulent_fluxes_at_a_point` and the surface temperature solves of
the soil and snow models, which use its output with
`turbulent_fluxes_from_output`.
"""
@inline function surface_fluxes_at_a_point(
    T_sfc_guess,
    q_vap_sfc_guess,
    update_T_sfc::UT,
    update_q_vap_sfc::UQ,
    P_atmos,
    T_atmos,
    q_tot_atmos,
    u_atmos,
    h_atmos,
    z_sfc,
    displ,
    roughness_model,
    gustiness,
    reference,
    earth_param_set,
) where {UT, UQ}
    thermo_params = LP.thermodynamic_parameters(earth_param_set)
    surface_flux_params = LP.surface_fluxes_parameters(earth_param_set)
    _grav = LP.grav(earth_param_set)
    config = surface_flux_config(roughness_model, gustiness, reference)
    positional_default_args = (
        scheme = SurfaceFluxes.PointValueScheme(),
        solver_opts = nothing,
        flux_specs = nothing,
    )
    # A coupled atmosphere provides a wind vector; a prescribed one a speed
    u = u_atmos isa Real ? (u_atmos, zero(u_atmos)) : u_atmos
    u_sfc = (zero(first(u)), zero(first(u)))
    q_liq_atmos = zero(q_tot_atmos)
    q_ice_atmos = zero(q_tot_atmos)
    ρ_atmos =
        Thermodynamics.air_density(thermo_params, T_atmos, P_atmos, q_tot_atmos)
    return SurfaceFluxes.surface_fluxes(
        surface_flux_params,
        T_atmos,
        q_tot_atmos,
        q_liq_atmos,
        q_ice_atmos,
        ρ_atmos,
        T_sfc_guess,
        q_vap_sfc_guess,
        _grav * z_sfc,
        h_atmos, # above the ground or apparent sink, according to `reference`
        displ,
        u,
        u_sfc,
        nothing, # roughness inputs
        config,
        positional_default_args...,
        update_T_sfc,
        update_q_vap_sfc,
    )
end

"""
     compute_Δz_eff(h_atmos, z_0m, displ, reference::SurfaceFluxes.ReferenceAboveSurface)

The surface density and geopotential difference in SurfaceFluxes.jl is computed using z_atmos - z_ground - d.
Under the reference `ReferenceAboveSurface`, `h_atmos = z_atmos - z_ground`.
"""
function compute_Δz_eff(
    h_atmos,
    z_0m,
    displ,
    reference::SurfaceFluxes.ReferenceAboveSurface,
)
    return h_atmos - displ
end

"""
     compute_Δz_eff(h_atmos, z_0m, displ, reference::SurfaceFluxes.ReferenceAboveApparentSink)

The surface density and geopotential difference in SurfaceFluxes.jl is computed using z_atmos - z_ground - d.
Under the reference `ReferenceAboveApparentSink`, `h_atmos = z_atmos - (z_ground + d + z_0m)`.
"""
function compute_Δz_eff(
    h_atmos,
    z_0m,
    displ,
    reference::SurfaceFluxes.ReferenceAboveApparentSink,
)
    return h_atmos + z_0m
end

"""
    turbulent_fluxes_from_output(return_extra_fluxes, output, T_sfc_guess, q_vap_sfc_guess,
                                 update_∂T_sfc∂T, update_∂q_sfc∂T, P_atmos, T_atmos,
                                 q_tot_atmos, h_atmos, displ, roughness_model, reference, earth_param_set)

Return the NamedTuple of `turbulent_fluxes_at_a_point` from the SurfaceFluxes.jl
`output` of `surface_fluxes_at_a_point`: the latent and sensible heat fluxes,
the vapor flux in volume of liquid water, the approximate derivatives of the
heat fluxes with respect to the component temperature (evaluated with
`update_∂T_sfc∂T` and `update_∂q_sfc∂T` at the surface temperature and
humidity `T_sfc_guess` and `q_vap_sfc_guess`, or omitted when both are
`nothing`), the surface temperature `T_sfc` and specific humidity `q_sfc` at
which the fluxes were evaluated, the friction velocity `ustar`, the stability
parameter `ζ`, and the effective height `Δz_eff`, followed for
`return_extra_fluxes = Val(true)` by the momentum fluxes and the buoyancy
flux. The atmospheric state is given at height `h_atmos` according to the
`reference` convention, the surface has displacement height `displ`, and the
roughness model is as given. Models that solve for their surface temperature
within the Monin-Obukhov iterations use it to obtain the fluxes from that
solve.
"""
@inline function turbulent_fluxes_from_output(
    return_extra_fluxes::Val,
    output,
    T_sfc_guess,
    q_vap_sfc_guess,
    update_∂T_sfc∂T::UDT,
    update_∂q_sfc∂T::UDQ,
    P_atmos,
    T_atmos,
    q_tot_atmos,
    h_atmos,
    displ,
    roughness_model,
    reference,
    earth_param_set,
) where {UDT, UDQ}
    surface_flux_params = LP.surface_fluxes_parameters(earth_param_set)
    _ρ_liq = LP.ρ_cloud_liq(earth_param_set)
    # The surface state applies at the displacement height, so the density is
    # extrapolated over the effective height, as in the flux solve
    z_0m = SurfaceFluxes.momentum_roughness(
        roughness_model,
        nothing,
        surface_flux_params,
        nothing,
    )
    Δz_eff = compute_Δz_eff(h_atmos, z_0m, displ, reference)
    needs_ρ_sfc = !(
        isnothing(update_∂T_sfc∂T) &&
        isnothing(update_∂q_sfc∂T) &&
        return_extra_fluxes isa Val{false}
    )
    ρ_sfc =
        needs_ρ_sfc ?
        compute_ρ_sfc(
            surface_flux_params,
            T_atmos,
            P_atmos,
            q_tot_atmos,
            Δz_eff,
            output.T_sfc,
            output.q_vap_sfc,
        ) : zero(output.T_sfc)
    derivatives = flux_temperature_derivatives(
        update_∂T_sfc∂T,
        update_∂q_sfc∂T,
        ρ_sfc,
        output.ustar,
        output.g_h,
        q_vap_sfc_guess,
        T_sfc_guess,
        earth_param_set,
    )
    fluxes = (;
        lhf = output.lhf,
        shf = output.shf,
        vapor_flux = output.evaporation / _ρ_liq, # volume of liquid water
        derivatives...,
        T_sfc = output.T_sfc,
        q_sfc = output.q_vap_sfc,
        ustar = output.ustar,
        # Stability parameter at which the exchange coefficients were
        # evaluated (capped in stable conditions)
        ζ = output.ζ_eff,
        Δz_eff,
    )
    return with_extra_fluxes(
        return_extra_fluxes,
        fluxes,
        output,
        surface_flux_params,
        ρ_sfc,
    )
end

# Approximate derivatives of the turbulent heat fluxes with respect to T_sfc;
# omitted when neither derivative callback is supplied
@inline flux_temperature_derivatives(::Nothing, ::Nothing, args...) = (;)
@inline function flux_temperature_derivatives(
    update_∂T_sfc∂T::UDT,
    update_∂q_sfc∂T::UDQ,
    ρ_sfc,
    u_star,
    g_h,
    q_vap_sfc_guess,
    T_sfc_guess,
    earth_param_set,
) where {UDT, UDQ}
    thermo_params = LP.thermodynamic_parameters(earth_param_set)
    _LH_v0 = LP.LH_v0(earth_param_set)
    cp_d = Thermodynamics.Parameters.cp_d(thermo_params)
    # This assumes that q_vap_sfc_guess is the saturated value
    ∂lhf∂T =
        ρ_sfc *
        g_h *
        _LH_v0 *
        update_∂q_sfc∂T(
            u_star,
            g_h,
            q_vap_sfc_guess,
            T_sfc_guess,
            earth_param_set,
        )
    ∂shf∂T = ρ_sfc * g_h * cp_d * update_∂T_sfc∂T(u_star, g_h, earth_param_set)
    return (; ∂lhf∂T, ∂shf∂T)
end

"""
    with_extra_fluxes(return_extra_fluxes, fluxes, output, surface_flux_params, ρ_sfc)

Return the NamedTuple `fluxes`, followed for `Val(true)` by the momentum fluxes
`ρτxz`, `ρτyz` of the SurfaceFluxes.jl `output` and the buoyancy flux at the
surface air density `ρ_sfc`.
"""
@inline with_extra_fluxes(::Val{false}, fluxes, args...) = fluxes
@inline function with_extra_fluxes(
    ::Val{true},
    fluxes,
    output,
    surface_flux_params,
    ρ_sfc,
)
    buoyancy_flux = SurfaceFluxes.buoyancy_flux(
        surface_flux_params,
        output.shf,
        output.lhf,
        output.T_sfc,
        ρ_sfc,
        output.q_vap_sfc,
    )
    return (; fluxes..., ρτxz = output.ρτxz, ρτyz = output.ρτyz, buoyancy_flux)
end

"""
    surface_temperature_newton_update(inputs, param_set, thermo_params, g_h, ρ_sfc, q_sfc,
                                      ∂q_sfc∂T, T_below, r, T_max, ϵ, σ, SW_n, LW_d)

Return the Newton update `clamp(T + ΔT, T / 2, T_max)` [K] of the temperature
`T = inputs.T_sfc_guess` of a surface with zero heat capacity, connected by the
thermal resistance `r` [m² K/W] to a layer at temperature `T_below` [K], with
`ΔT = -f(T)/f'(T)` and

    r f(T) = r (SW_n + LW_n(T) + L(T) + H(T)) + (T - T_below),

where all fluxes are positive upward and `LW_n = -ϵ (LW_d - σ T⁴)`. The latent
and sensible heat fluxes `L` and `H` are evaluated with the SurfaceFluxes.jl
`inputs`, the heat conductance `g_h` [m/s], the surface air density `ρ_sfc`
[kg/m³], and the surface specific humidity `q_sfc` [kg/kg]; `f'(T)` uses its
derivative `∂q_sfc∂T` [1/K] and neglects those of `g_h` and `ρ_sfc`. The cap
`T_max` [K] applies within the Monin-Obukhov iterations, so the fluxes of the
solve are those at the capped temperature; the energy that would warm the
surface above `T_max` goes into phase change. The floor `T / 2` limits a
downward step to half the current temperature, so that every iterate stays
positive and the saturation humidity and surface density remain defined.

Called from the SurfaceFluxes.jl `update_T` callbacks of the soil skin and the
snow surface temperature.
"""
function surface_temperature_newton_update(
    inputs,
    param_set,
    thermo_params,
    g_h,
    ρ_sfc,
    q_sfc,
    ∂q_sfc∂T,
    T_below,
    r,
    T_max,
    ϵ,
    σ,
    SW_n,
    LW_d,
)
    T_sfc = inputs.T_sfc_guess
    E = SurfaceFluxes.evaporation(
        param_set,
        inputs,
        g_h,
        SurfaceFluxes.interior_vapor_specific_humidity(inputs),
        q_sfc,
        ρ_sfc,
        inputs.moisture_model,
    )
    L = SurfaceFluxes.latent_heat_flux(
        param_set,
        inputs,
        E,
        inputs.moisture_model,
    )
    H = SurfaceFluxes.sensible_heat_flux(
        param_set,
        inputs,
        g_h,
        inputs.T_int,
        T_sfc,
        ρ_sfc,
        E,
    )
    _LH_v0 = Thermodynamics.Parameters.LH_v0(thermo_params)
    cp_d = Thermodynamics.Parameters.cp_d(thermo_params)
    LW_n = -ϵ * (LW_d - σ * T_sfc^4)
    F = SW_n + LW_n + L + H
    ∂F∂T = 4 * ϵ * σ * T_sfc^3 + ρ_sfc * g_h * (_LH_v0 * ∂q_sfc∂T + cp_d)
    ΔT = -(r * F + (T_sfc - T_below)) / (r * ∂F∂T + 1)
    return clamp(T_sfc + ΔT, T_sfc / 2, T_max)
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
    SW <: AbstractTimeVaryingInput,
    FD <: Union{Nothing, AbstractTimeVaryingInput},
    LW <: AbstractTimeVaryingInput,
    DT,
    T <: Union{Nothing, Function},
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
    cosθs::T
    "Thermodynamic parameters"
    thermo_params::TP
    function PrescribedRadiativeFluxes(
        FT,
        SW_d,
        LW_d,
        start_date;
        cosθs = nothing,
        toml_dict::Union{CP.ParamDict, Nothing} = nothing,
        frac_diff = nothing,
    )
        if isnothing(toml_dict)
            thermo_params = nothing
        else
            earth_param_set = LP.LandParameters(toml_dict)
            thermo_params = LP.thermodynamic_parameters(earth_param_set)
        end
        args = (SW_d, frac_diff, LW_d, start_date, cosθs, thermo_params)
        return new{FT, typeof.(args)...}(args...)
    end
end

"""
    net_radiation!(dest::ClimaCore.Fields.Field,
                   radiation::AbstractRadiativeDrivers,
                   model::AbstractModel{FT},
                   Y::ClimaCore.Fields.FieldVector,
                   p::NamedTuple,
                   t,
                   ) where {FT}


Computes net radiative flux assuming that the incoming shortwave
radiation is stored in `p.drivers.SW_d`, and the incoming longwave
radiation is stored in `p.drivers.LW_d`.  This can be because
- the prescribed functions are called to update these values
- the coupler has updated these values.
"""
function net_radiation!(
    dest::ClimaCore.Fields.Field,
    radiation::AbstractRadiativeDrivers,
    model::AbstractModel,
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    t,
)
    LW_d = p.drivers.LW_d
    SW_d = p.drivers.SW_d
    earth_param_set = model.parameters.earth_param_set
    _σ = LP.Stefan(earth_param_set)
    T_sfc = component_temperature(model, Y, p)
    α_sfc = surface_albedo(model, Y, p)
    ϵ_sfc = surface_emissivity(model, Y, p)
    @. dest = -(1 - α_sfc) * SW_d - ϵ_sfc * (LW_d - _σ * T_sfc^4)
    return nothing
end

"""
    component_temperature(model::AbstractModel, Y, p)

A helper function which returns the component temperature for a given
model. This is used in the turbulent_fluxes! function as the value of
surface temperature and used to compute long wave radiative fluxes.

Extending this function for your model is only necessary if you need to
compute surface fluxes and radiative fluxes at the surface using
the functions in this file. If your model uses different temperatures to
compute turbulent fluxes vs long wave radiative fluxes, please see
the documentation for an update function in the `turbulent_fluxes!` call.
"""
function component_temperature(model::AbstractModel, Y, p) end

"""
    component_specific_humidity(model::AbstractModel, Y, p)
    component_specific_humidity(model::AbstractModel, Y, p, forcing)

A helper function which returns the surface specific humidity for a given
model, needed because different models compute and store q_sfc in
different ways and places.

The second form takes the forcing `(; h_atmos, T_atmos, q_atmos)` of the flux
solve, the reference height and the air temperature and specific humidity at
it, which differ from the atmospheric forcing below a canopy. Models whose
surface humidity depends on the forcing extend it; for other models, it
returns the first form.

Extending this function for your model is only necessary if you need to
compute surface fluxes using the functions in this file. Please see
the documentation for an update function in the `turbulent_fluxes!` call.
"""
function component_specific_humidity(model::AbstractModel, Y, p) end
component_specific_humidity(model::AbstractModel, Y, p, forcing) =
    component_specific_humidity(model, Y, p)

"""
    surface_roughness_model(model::AbstractModel, Y, p)

Returns the SurfaceFluxes roughness model for `model`.

Extending this function for your model is only necessary if you need to
compute surface fluxes using the functions in this file.
"""
function surface_roughness_model(model::AbstractModel, Y, p) end

"""
    surface_displacement_height(model::AbstractModel, Y, p)

Returns the displacement height for `model`.

Extending this function for your model is only necessary if you need to
compute surface fluxes using the functions in this file.
"""
function surface_displacement_height(model::AbstractModel, Y, p)
    FT = FTfromY(Y)
    return FT(0)
end

"""
    get_update_surface_temperature_function(model::AbstractModel, Y, p)

Returns the SurfaceFluxes `update_T_sfc` function for `model`.

This is only required if the output of `component_temperature` does not coincide
with the temperature that should be used to compute turbulent fluxes.

Extending this function for your model is only necessary if you need to
compute surface fluxes using the functions in this file.
"""
function get_update_surface_temperature_function(model::AbstractModel, Y, p) end

"""
    get_update_surface_humidity_function(model::AbstractModel, Y, p)

Returns the SurfaceFluxes `update_q_vap_sfc` function for `model`.

This is only required if the output of `component_specific_humidity` does not
coincide with the specific humidity that should be used to compute turbulent
fluxes.

Extending this function for your model is only necessary if you need to
compute surface fluxes using the functions in this file.
"""
function get_update_surface_humidity_function(model::AbstractModel, Y, p) end

"""
    get_∂T_sfc∂T_function(model::AbstractModel, Y, p)

Returns a function which computes the partial derivative of the surface temperature
used in `turbulent_fluxes!` with respect to the component temperature. The expected
arguments are (`u_star`, `g_h`, `earth_param_set`).

This is only required if the output of `component_temperature` does not coincide
with the temperature that should be used to compute turbulent fluxes, and if your
model requires derivatives of the fluxes with respect to the component temperature.
A model whose cache stores no flux derivatives returns `nothing` from this
function and from `get_∂q_sfc∂T_function`, and `turbulent_fluxes_from_output`
then omits the derivatives (see `Soil.EnergyHydrology`).

Extending this function for your model is only necessary if you need to
compute surface fluxes using the functions in this file.
"""
function get_∂T_sfc∂T_function(model::AbstractModel, Y, p)
    update_∂T_sfc∂T_at_a_point(u_star, g_h, earth_param_set) = one(g_h)
    return update_∂T_sfc∂T_at_a_point
end

"""
    get_∂q_sfc∂T_function(model::AbstractModel, Y, p)

Returns a function which computes the partial derivative of the surface humidity
used in `turbulent_fluxes!` with respect to the component temperature. The expected
arguments are (`u_star`, `g_h`, `q_sat`, `T_sfc`, `earth_param_set`).

This is only required if the output of `component_specific_humidity` does not coincide
with the humidity that should be used to compute turbulent fluxes, and if your
model requires derivatives of the fluxes with respect to the component temperature.
A model whose cache stores no flux derivatives returns `nothing` from this
function and from `get_∂T_sfc∂T_function` (see `Soil.EnergyHydrology`).

Extending this function for your model is only necessary if you need to
compute surface fluxes using the functions in this file.
"""
function get_∂q_sfc∂T_function(model::AbstractModel, Y, p)

    function update_∂q_sfc∂T_at_a_point(
        u_star,
        g_h,
        q_sat,
        T_sfc,
        earth_param_set,
    )
        _T_freeze = LP.T_freeze(earth_param_set)
        phase = ifelse(
            T_sfc > _T_freeze,
            Thermodynamics.Liquid(),
            Thermodynamics.Ice(),
        )
        return ClimaLand.partial_q_sat_partial_T(
            q_sat,
            T_sfc,
            phase,
            earth_param_set,
        )
    end
    return update_∂q_sfc∂T_at_a_point
end
"""
    surface_elevation(model::AbstractModel, Y, p)

Elevation of the ground surface [m], the level from which 
SurfaceFluxes measures the displacement height. 

Zero for all components until topography is represented.
"""
function surface_elevation(model::AbstractModel, Y, p)
    FT = FTfromY(Y)
    return FT(0)
end

# Deprecated name of `surface_elevation`, which the ClimaCoupler bucket calls
Base.@deprecate(
    surface_height(model::AbstractModel, Y, p),
    surface_elevation(model, Y, p),
    false,
)

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
    specific_humidity_from_dewpoint(dewpoint_temperature, temperature, air_pressure, earth_param_set)

Estimates the specific humidity given the dewpoint temperature, temperature of the air
in Kelvin, and air pressure in Pa, along with the ClimaLand earth_param_set. This is useful
for creating the PrescribedAtmosphere - which needs specific humidity - from ERA5 reanalysis data.
"""
function specific_humidity_from_dewpoint(
    T_dew_air::data_FT,
    T_air::data_FT,
    P_air::data_FT,
    earth_param_set,
) where {data_FT <: Real}
    thermo_params = LP.thermodynamic_parameters(earth_param_set)
    sim_FT = typeof(LP.T_freeze(earth_param_set))

    # Obtain true rh using vapor pressure
    e_sat_dew = Thermodynamics.saturation_vapor_pressure(
        thermo_params,
        sim_FT(T_dew_air),
        Thermodynamics.Liquid(),
    )
    e_sat_T = Thermodynamics.saturation_vapor_pressure(
        thermo_params,
        sim_FT(T_air),
        Thermodynamics.Liquid(),
    )
    rh::sim_FT = e_sat_dew / e_sat_T

    q = Thermodynamics.q_vap_from_RH(
        thermo_params,
        sim_FT(P_air),
        sim_FT(T_air),
        rh,
        Thermodynamics.Liquid(),
    )

    return q
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
when running the canopy in standalone mode, including the surface temperature,
albedo, and emissivity, and soil water content, porosity, residual water fraction,
and hydrology closure model.

Note that internally we enforce that the soil water content `θ` must be within the
range `(θ\\_r, ν]`.

$(DocStringExtensions.FIELDS)
"""
struct PrescribedGroundConditions{
    FT,
    F1 <: AbstractTimeVaryingInput,
    F2 <: AbstractTimeVaryingInput,
    C,
} <: AbstractGroundConditions{FT}
    "Prescribed soil water content (m^3/m^3) in the root zone as a function of time"
    θ::F1
    "Prescribed ground surface temperature (K) as a function of time"
    T::F2
    "Ground albedo for PAR"
    α_PAR::FT
    "Ground albedo for NIR"
    α_NIR::FT
    "Ground emissivity"
    ϵ::FT
    "Soil porosity (m^3/m^3"
    ν::FT
    "The soil residual water fraction (m^3/m^3)"
    θ_r::FT
    "The soil hydrology closure model: van Genuchten or Brooks and Corey"
    hydrology_cm::C
end

"""
    function PrescribedGroundConditions{FT}(;
        ν = FT(0.4),
        θ = TimeVaryingInput((t) -> ν),
        T = TimeVaryingInput((t) -> 298.0),
        α_PAR = FT(0.2),
        α_NIR = FT(0.4),
        ϵ = FT(0.99),
        θ_r = FT(0),
        hydrology_cm = vanGenuchten{FT}(; α = FT(2), n = FT(1.5)),
    ) where {FT <: AbstractFloat}

An outer constructor for the PrescribedGroundConditions allowing the user to
specify the ground parameters by keyword arguments.

The defaults should be overridden for physically meaningful results.
"""
function PrescribedGroundConditions{FT}(;
    ν = FT(0.4),
    θ = TimeVaryingInput((t) -> ν),
    T = TimeVaryingInput((t) -> 298.0),
    α_PAR = FT(0.2),
    α_NIR = FT(0.4),
    ϵ = FT(0.99),
    θ_r = FT(0),
    hydrology_cm = vanGenuchten{FT}(; α = FT(2), n = FT(1.5)),
) where {FT <: AbstractFloat}
    return PrescribedGroundConditions{
        FT,
        typeof(θ),
        typeof(T),
        typeof(hydrology_cm),
    }(
        θ,
        T,
        α_PAR,
        α_NIR,
        ϵ,
        ν,
        θ_r,
        hydrology_cm,
    )
end

"""
     PrognosticGroundConditions <: Canopy.AbstractGroundConditions

A type of AbstractGroundConditions to use when running the CanopyModel
as part of an integrated model, i.e. with a prognostic `EnergyHydrology`
soil model and either with or without a snow model.

Note that this struct is linked with the `EnergyHydrology`/`SnowModel` models. If we ever had a different
soil model, we might need to construct a different `PrognosticGroundConditions` because
the fields may be stored in different places.
"""
struct PrognosticGroundConditions{FT} <: AbstractGroundConditions{FT} end

"""
    initialize_drivers(a::PrescribedGroundConditions{FT}, coords) where {FT}

Creates and returns a NamedTuple for the `PrescribedGroundConditions` driver,
with variables `θ` (water content in the soil in the root zone), `T` (temperature
of the surface of the ground).
"""
function initialize_drivers(
    a::PrescribedGroundConditions{FT},
    coords,
) where {FT}
    keys = (:θ, :T_ground)
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
    keys = (:P_liq, :P_snow, :T, :P, :u, :q, :c_co2)
    types = ([FT for _ in keys]...,)
    domain_names = ([:surface for _ in keys]...,)
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
    types = ([FT for _ in keys]...,)
    domain_names = ([:surface for _ in keys]...,)
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
with variables `P_liq`, `P_snow`, `c_co2`, `T`, `P`, `q`, and `u`.

This is intended to be used in coupled simulations with ClimaCoupler.jl.
"""
function ClimaLand.initialize_drivers(
    a::CoupledAtmosphere{FT},
    coords,
) where {FT}
    keys = (:P_liq, :P_snow, :c_co2, :T, :P, :q, :u)
    types = ([FT for _ in 1:(length(keys) - 1)]..., SVector{2, FT})
    domain_names = ([:surface for _ in keys]...,)
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
            initialize_drivers(driver, coords)
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
        NVTX.@annotate function update_drivers!(p, t)
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
        evaluate!(p.drivers.θ, a.θ, t)
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

When `r.cosθs` is `nothing`, the cosine zenith angle
not changed, and should be updated by the coupler. This differs from the behavior of
`PrescribedRadiativeFluxes`, where the cosine zenith angle is set to `NaN` if `cosθs` is `nothing`.

Otherwise, the cosine zenith angle is computed using `r.cosθs(t, r.start_date)`.
"""
make_update_drivers(r::CoupledRadiativeFluxes{FT, Nothing}) where {FT} =
    (p, t) -> nothing
function make_update_drivers(r::CoupledRadiativeFluxes{FT}) where {FT}
    function update_drivers!(p, t)
        p.drivers.cosθs .= r.cosθs(t, r.start_date)
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
        if !isnothing(r.cosθs)
            p.drivers.cosθs .= FT.(r.cosθs(t, r.start_date))
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
    prescribed_forcing_era5(start_date,
                            stop_date,
                            surface_space,
                            toml_dict::CP.ParamDict,
                            FT;
                            use_lowres_forcing = false,
                            gustiness = SurfaceFluxes.FlooredDeardorffGustinessSpec(FT(1)),
                            reference = SurfaceFluxes.ReferenceAboveApparentSink(),
                            max_wind_speed = nothing,
                            c_co2 = TimeVaryingInput((t) -> 4.2e-4),
                            time_interpolation_method = LinearInterpolation(PeriodicCalendar()),
                            regridder_type = :InterpolationsRegridder,
                            context = nothing,
                            )

A helper function which constructs the `PrescribedAtmosphere` and
`PrescribedRadiativeFluxes` from ERA5 data stored in NetCDF files.

There are two versions of the ERA5 data available through ClimaLand:
- a high resolution version (1° x 1°) available for the years 1979-2024
- a low resolution version (8° x 8°) available only for the year 2008

The high resolution version will be used if `use_lowres_forcing = false` (the default).
This artifact is used for global runs on compute clusters, but is too large to be used for local testing,
and requires you to have acquired the data in advance.
The low resolution version will be used if `use_lowres_forcing = true`.
If the simulation dates are outside of 2008, the 2008 data will be reused for each year of simulation.
This artifact is recommended for local testing or quick runs where accuracy is less critical.

The method for temporal interpolation is controlled via the `time_interpolation_method` kwarg.
We suggest `LinearInterpolation(PeriodicCalendar())`, which implies linear interpolation in time;
the inner argument implies how extrapolation outside the bounds of the data is handled. For example,
the ERA5 forcing data we use is hourly, which implies Dec 31 of the last year of the data, at midnight,
is not in the data. With the PeriodicCalendar()
option, the interpolated value in the data at Jan 1 at timestamp 00 of the first year of the simulation will
be used. For the low-resolution forcing data, which only exists for 2008, multi-year runs will repeat
the forcing. If this behavior is not what you want, you can change the extrapolation argument to
error `LinearInterpolation(Throw())` or extrapolate by using the last value ``LinearInterpolation(Flat())`.
More information is available in the ClimaUtilities documentation:
https://clima.github.io/ClimaUtilities.jl/dev/inputs/#Extrapolation-boundary-conditions.

The ClimaLand default is to use nearest neighbor spatial interpolation for low resolution forcing,
and linear spatial interpolation for high resolution forcing.

!!! warning "Clipped values"
    High wind speed anomalies (10-100x increase and decrease over a period of a
    several hours) appear in the ERA5 reanalysis data. These generate very large
    surface fluxes (due to wind speeds up to 300 m/s), which lead to
    instability. The kwarg `max_wind_speed`, with a value give in m/s, is used to
    clip these if it is not `nothing`. See
    [here](https://confluence.ecmwf.int/display/CKB/ERA5%3A+large+10m+winds).

!!! note "Full high resolution dataset available on clima cluster only"
    The full 40 year dataset of high resolution ERA5 data is only available on the
    clima cluster.

!!! note "Reference heights"
    The wind is the ERA5 10 m wind, `sqrt(u10² + v10²)`, and the reference
    height of the forcing is set to 10 m. The temperature and humidity are the
    ERA5 2 m fields (`t2m`, `d2m`), which are themselves products of ERA5's own
    surface-layer scheme, not 10 m values. The Monin-Obukhov solve therefore
    evaluates the temperature and humidity differences to the surface over
    2 m but the wind over 10 m, which weakens the diagnosed stability and
    instability. The hourly-mean components also give a scalar wind speed
    somewhat below the mean of the instantaneous speed.
"""
function prescribed_forcing_era5(
    start_date,
    stop_date,
    surface_space,
    toml_dict::CP.ParamDict,
    FT;
    use_lowres_forcing = false,
    gustiness = SurfaceFluxes.FlooredDeardorffGustinessSpec(FT(1)),
    reference = SurfaceFluxes.ReferenceAboveApparentSink(),
    max_wind_speed = nothing,
    c_co2 = TimeVaryingInput((t) -> 4.2e-4),
    time_interpolation_method = LinearInterpolation(PeriodicCalendar()),
    regridder_type = :InterpolationsRegridder,
    context = nothing,
)
    # Determine the ERA5 dataset and interpolation method based on the simulation years and resolution
    if use_lowres_forcing
        era5_ncdata_path =
            ClimaLand.Artifacts.era5_land_forcing_data2008_lowres_path(;
                context,
            )
        interpolation_method = Interpolations.Constant()
        # We can use the 2008 forcing outside of that year, but we need periodic boundaries for the time
        if year(start_date) != 2008 ||
           (year(stop_date) != 2008 && stop_date != DateTime(2009, 1, 1))
            time_interpolation_method !=
            LinearInterpolation(PeriodicCalendar()) &&
                @warn "Using low resolution ERA5 forcing outside of 2008 requires periodic time interpolation. Overriding `time_interpolation_method`."
            time_interpolation_method = LinearInterpolation(PeriodicCalendar())
        end
    else
        era5_ncdata_path = ClimaLand.Artifacts.find_era5_year_paths(
            start_date,
            stop_date;
            context,
        )
        interpolation_method = Interpolations.Constant()
    end

    # Pass a list of files in all cases
    era5_ncdata_path isa String && (era5_ncdata_path = [era5_ncdata_path])

    earth_param_set = LP.LandParameters(toml_dict)
    _ρ_liq = LP.ρ_cloud_liq(earth_param_set)
    # Precip is provided as a mass flux; convert to volume flux of liquid water with ρ = 1000 kg/m^3
    precip = TimeVaryingInput(
        [era5_ncdata_path, era5_ncdata_path],
        ["mtpr", "msr"],
        surface_space;
        start_date,
        regridder_type,
        regridder_kwargs = (; interpolation_method),
        file_reader_kwargs = (; preprocess_func = (data) -> -data / _ρ_liq,),
        method = time_interpolation_method,
        compose_function = (mtpr, msr) -> min.(mtpr .- msr, Float32(0)),
    )
    # Precip is provided as a mass flux; convert to volume flux of liquid water with ρ = 1000 kg/m^3
    snow_precip = TimeVaryingInput(
        era5_ncdata_path,
        "msr",
        surface_space;
        start_date,
        regridder_type,
        regridder_kwargs = (; interpolation_method),
        file_reader_kwargs = (; preprocess_func = (data) -> -data / _ρ_liq,),
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
        regridder_kwargs = (; interpolation_method),
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
        regridder_kwargs = (; interpolation_method),
        compose_function = specific_humidity,
        method = time_interpolation_method,
    )
    P_atmos = TimeVaryingInput(
        era5_ncdata_path,
        "sp",
        surface_space;
        start_date,
        regridder_type,
        regridder_kwargs = (; interpolation_method),
        method = time_interpolation_method,
    )

    T_atmos = TimeVaryingInput(
        era5_ncdata_path,
        "t2m",
        surface_space;
        start_date,
        regridder_type,
        regridder_kwargs = (; interpolation_method),
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
        toml_dict;
        gustiness,
        reference,
        c_co2 = c_co2,
    )

    SW_d = TimeVaryingInput(
        era5_ncdata_path,
        "msdwswrf",
        surface_space;
        start_date,
        regridder_type,
        regridder_kwargs = (; interpolation_method),
        file_reader_kwargs = (;
            preprocess_func = (data) -> max.(data, Float32(0)),
        ),
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
        regridder_kwargs = (; interpolation_method),
        method = time_interpolation_method,
        compose_function = compute_diffuse_fraction_broadcasted,
    )
    LW_d = TimeVaryingInput(
        era5_ncdata_path,
        "msdwlwrf",
        surface_space;
        start_date,
        regridder_type,
        regridder_kwargs = (; interpolation_method),
        file_reader_kwargs = (;
            preprocess_func = (data) -> max.(data, Float32(0)),
        ),
        method = time_interpolation_method,
    )
    cos_zenith_angle =
        (t, s) -> default_cos_zenith_angle(
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
        cosθs = cos_zenith_angle,
        toml_dict = toml_dict,
        frac_diff = frac_diff,
    )
    return (; atmos, radiation)
end

"""
     prescribed_analytic_forcing(FT = Float32;
                                 toml_dict::CP.ParamDict,
                                 start_date = DateTime(2005),
                                 SW_d = (t) -> 0,
                                 LW_d = (t) -> 5.67e-8 * 280.0^4.0,
                                 precip = (t) -> 0, # no precipitation
                                 T_atmos = (t) -> 280.0,
                                 u_atmos = (t) -> 1.0,
                                 q_atmos = (t) -> 0.0, # no atmos water
                                 h_atmos = FT(1),
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
                                     LP.LandParameters(toml_dict),
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

Note that this uses the reference (SurfaceFluxes.ReferenceAboveSurface),
where `h_atmos` is provided relative to the ground.
"""
function prescribed_analytic_forcing(
    FT = Float32;
    toml_dict::CP.ParamDict,
    start_date = DateTime(2005),
    SW_d = (t) -> 0,
    LW_d = (t) -> 5.67e-8 * 280.0^4.0,
    precip = (t) -> 0, # no precipitation
    T_atmos = (t) -> 280.0,
    u_atmos = (t) -> 1.0,
    q_atmos = (t) -> 0.0, # no atmos water
    h_atmos = FT(1), # m
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
        toml_dict,
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
of the solar zenith angle (`cosθs`), the total surface downwelling shortwave radiation (`SW_d`),
the air temperature (`T`), air pressure (`P`), specific humidity (`q`), and the day of the year
(`td`).

See Appendix A of Braghiere, "Evaluation of turbulent fluxes of CO2, sensible heat,
and latent heat as a function of aerosol optical depth over the course of deforestation
in the Brazilian Amazon" 2013.

Note that cosθs is a coefficient of k₀, which we divide by in this expression.
This can amplify small errors when cosθs is near 0.

This formula is empirical and can yield negative numbers depending on the
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
        # If t is an ITime, and it has an epoch, use that instead of start_date to compute
        # the current simulation date
        DOY = Dates.dayofyear(
            isnothing(t.epoch) ?
            start_date + Dates.Second(floor(Int64, float(t))) : date(t),
        )
    else
        DOY = Dates.dayofyear(start_date + Dates.Second(floor(Int64, t)))
    end
    # Pass explicit zero liquid/ice fractions to avoid Thermodynamics.jl MethodError with default Integer arguments
    RH = Thermodynamics.relative_humidity(thermo_params, T, P, q) # Note: This used RH over liquid before; not it's RH over ice when below freezing
    # TODO Replace magic numbers with appropriate constants; for example 1370 probably is TSI and should come from ClimaParams for consistency with rest of model
    # TODO Nondimensionalize equations appropriately; here we have dimensional constants that are not clearly identified as such, easily leading to errors
    cosθs = max(cosθs, eps(FT))
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

"""
    prescribed_forcing_crujra(start_date,
                              stop_date,
                              surface_space,
                              toml_dict::CP.ParamDict,
                              FT;
                              gustiness = SurfaceFluxes.FlooredDeardorffGustinessSpec(FT(1)),
                              reference = SurfaceFluxes.ReferenceAboveApparentSink(),
                              c_co2 = TimeVaryingInput((t) -> 4.2e-4),
                              version = "2.5",
                              crujra_path = nothing,
                              cycle_years = nothing,
                              time_interpolation_method = isnothing(cycle_years) ?
                                  LinearInterpolation() : LinearInterpolation(Flat()),
                              sw_time_interpolation = :coszen,
                              state_time_interpolation = :harmonic,
                              flux_time_interpolation = :interval,
                              updateat = Dates.Hour(3),
                              sw_quadrature_step = Dates.Minute(15),
                              regridder_type = :InterpolationsRegridder,
                              context = nothing,
                             )

Construct `PrescribedAtmosphere` and `PrescribedRadiativeFluxes` drivers from
the CRU-JRA reanalysis NetCDF data.

Data version
- `version` is the CRU-JRA version, `"2.5"` (the `crujra_forcing_data`
    artifact, 1901–2023) or `"4.10"` (1901–2025, the forcing of the TRENDY v15 /
    Global Carbon Budget 2026 protocol). See `ClimaLand.Artifacts.crujra_years`.
- `crujra_path`, if given, is the directory of the yearly files
    `crujra_forcing_data_YYYY_0.5x0.5.nc`, used instead of the artifact. v4.10
    is not an artifact yet, so it needs `crujra_path`.

Time handling
- CRU-JRA gives temperature, humidity, pressure and wind as instantaneous values
    at 00, 06, 12 and 18 UTC, and shortwave and longwave radiation and
    precipitation as means over the 6 hours that start at the time stamp. The
    drivers are updated every `updateat` and held in between; `updateat` must
    be the driver update interval of the simulation (`updateat` of
    `LandSimulation`, 3 hours by default; 1 hour is recommended with CRU-JRA).
- Shortwave: with `sw_time_interpolation = :coszen` (the default), each
    6-hour mean is distributed within its interval in proportion to the cosine
    of the solar zenith angle, as the `coszen` time interpolation of the CESM
    data atmosphere (CDEPS) that CLM uses for CRU-JRA. The shortwave given at an
    update is the mean of this distribution over the time until the next update,
    and the cosine of the zenith angle is its shortwave-weighted mean over the
    same window (see [`ClimaLand.CosZenithWeightedInput`](@ref)). The 6-hour
    energy is conserved exactly when the updates fall on multiples of
    `sw_quadrature_step` after 00 UTC. The diffuse fraction is constant over
    each 6-hour interval.
- Longwave radiation and precipitation: with `flux_time_interpolation =
    :interval` (the default), each 6-hour mean is held constant over its
    interval, and the value given at an update is the mean over the time until
    the next update (see [`ClimaLand.IntervalMeanInput`](@ref)), so the energy
    and water of each interval stay in that interval.
- Temperature, humidity, pressure and wind: with `state_time_interpolation =
    :harmonic` (the default), the value given at an update is taken at the
    middle of the hold window from a constant plus diurnal and semidiurnal
    harmonics through the four 6-hourly samples around it (see
    [`ClimaLand.HoldWindowStateInput`](@ref)). Linear interpolation between
    6-hourly samples cuts the afternoon temperature peak, which lowers the
    daytime vapor pressure deficit. `state_time_interpolation = :cubic` uses
    the Catmull-Rom cubic through the same samples, and `:linear_midpoint`
    interpolates linearly, also at the middle of the window.
- `:linear` for any of the three keywords recovers the earlier behaviour for
    those variables: the records are interpolated linearly in time with
    `time_interpolation_method` (defaults to `LinearInterpolation()`, no
    periodic extrapolation), as if they were instantaneous values, and taken
    at the update time (for the shortwave, the zenith angle is the one at the
    update time).
- With `cycle_years`, e.g. `1901:1920`, the forcing of a year `y` is read from
    year `first(cycle_years) + mod(y - first(cycle_years), length(cycle_years))`,
    on the same day and time, as for the pre-industrial climate of the TRENDY
    protocol (1901–1920 cycled); a February 29 that maps to a non-leap year uses
    February 28. The default `time_interpolation_method` is then
    `LinearInterpolation(Flat())`, so that the last 6 hours of the last cycled
    year use its last record. The cycling is applied to the 6-hourly records,
    before the time handling above, so the samples that the interpolation takes
    after the last cycled year come from the first one, as the model sees them.

How it differs from `prescribed_forcing_era5`
- CRU-JRA provides scalar 10 m wind speed (`"wind"`), so no u/v combination is
    needed (ERA5 combines `u10`/`v10`).
- Only a single resolution is supported here; ERA5 exposes both high- and
    low-resolution artifacts and switches interpolation strategies accordingly.
- The default time interpolation is non-periodic; ERA5 defaults to
    `LinearInterpolation(PeriodicCalendar())` to enable seamless multi-year
    reuse of limited data (especially for low-res forcing).

!!! note "Full dataset availability"
        The full CRU-JRA dataset used by ClimaLand is only available on the clima
        cluster.
"""
function prescribed_forcing_crujra(
    start_date,
    stop_date,
    surface_space,
    toml_dict::CP.ParamDict,
    FT;
    gustiness = SurfaceFluxes.FlooredDeardorffGustinessSpec(FT(1)),
    reference = SurfaceFluxes.ReferenceAboveApparentSink(),
    c_co2 = TimeVaryingInput((t) -> 4.2e-4),
    version = "2.5",
    crujra_path = nothing,
    cycle_years = nothing,
    time_interpolation_method = isnothing(cycle_years) ? LinearInterpolation() :
                                LinearInterpolation(Flat()),
    sw_time_interpolation = :coszen,
    state_time_interpolation = :harmonic,
    flux_time_interpolation = :interval,
    updateat = Dates.Hour(3),
    sw_quadrature_step = Dates.Minute(15),
    regridder_type = :InterpolationsRegridder,
    context = nothing,
)
    sw_time_interpolation in (:coszen, :linear) ||
        error("`sw_time_interpolation` must be `:coszen` or `:linear`")
    state_time_interpolation in
    (:cubic, :harmonic, :linear_midpoint, :linear) || error(
        "`state_time_interpolation` must be `:cubic`, `:harmonic`, `:linear_midpoint` or `:linear`",
    )
    flux_time_interpolation in (:interval, :linear) ||
        error("`flux_time_interpolation` must be `:interval` or `:linear`")

    data_start_date, data_stop_date =
        isnothing(cycle_years) ? (start_date, stop_date) :
        (DateTime(first(cycle_years)), DateTime(last(cycle_years), 12, 31))
    crujra_ncdata_path = ClimaLand.Artifacts.find_crujra_year_paths(
        data_start_date,
        data_stop_date;
        context,
        version,
        crujra_path,
    )
    cycle(input) =
        isnothing(cycle_years) ? input :
        CycledYearsInput(input, start_date, cycle_years)

    interpolation_method = Interpolations.Constant()

    crujra_ncdata_path isa String && (crujra_ncdata_path = [crujra_ncdata_path])

    earth_param_set = LP.LandParameters(toml_dict)
    _ρ_liq = LP.ρ_cloud_liq(earth_param_set)

    # Replaces missing values with NaN.
    clean_missing = data -> coalesce.(data, Float32(NaN))

    # Precip is provided as a mass flux; convert to volume flux of liquid water with ρ = 1000 kg/m^3
    precip = TimeVaryingInput(
        [crujra_ncdata_path, crujra_ncdata_path],
        ["mtpr", "msr"],
        surface_space;
        start_date,
        regridder_type,
        regridder_kwargs = (; interpolation_method),
        file_reader_kwargs = (;
            preprocess_func = (data) -> -clean_missing(data) / _ρ_liq
        ),
        method = time_interpolation_method,
        compose_function = (mtpr, msr) -> min.(mtpr .- msr, Float32(0)),
    )

    snow_precip = TimeVaryingInput(
        crujra_ncdata_path,
        "msr",
        surface_space;
        start_date,
        regridder_type,
        regridder_kwargs = (; interpolation_method),
        file_reader_kwargs = (;
            preprocess_func = (data) -> -clean_missing(data) / _ρ_liq
        ),
        method = time_interpolation_method,
    )

    u_atmos = TimeVaryingInput(
        crujra_ncdata_path,
        "wind",
        surface_space;
        start_date,
        regridder_type,
        regridder_kwargs = (; interpolation_method),
        file_reader_kwargs = (; preprocess_func = clean_missing),
        method = time_interpolation_method,
    )
    specific_humidity(Td, T, P; params = earth_param_set) =
        ClimaLand.specific_humidity_from_dewpoint.(Td, T, P, params)
    q_atmos = TimeVaryingInput(
        [crujra_ncdata_path, crujra_ncdata_path, crujra_ncdata_path],
        ["d2m", "t2m", "sp"],
        surface_space;
        start_date,
        regridder_type,
        regridder_kwargs = (; interpolation_method),
        compose_function = specific_humidity,
        file_reader_kwargs = (; preprocess_func = clean_missing),
        method = time_interpolation_method,
    )
    P_atmos = TimeVaryingInput(
        crujra_ncdata_path,
        "sp",
        surface_space;
        start_date,
        regridder_type,
        regridder_kwargs = (; interpolation_method),
        file_reader_kwargs = (; preprocess_func = clean_missing),
        method = time_interpolation_method,
    )

    T_atmos = TimeVaryingInput(
        crujra_ncdata_path,
        "t2m",
        surface_space;
        start_date,
        regridder_type,
        regridder_kwargs = (; interpolation_method),
        file_reader_kwargs = (; preprocess_func = clean_missing),
        method = time_interpolation_method,
    )
    h_atmos = FT(10)

    SW_d = TimeVaryingInput(
        crujra_ncdata_path,
        "msdwswrf",
        surface_space;
        start_date,
        regridder_type,
        regridder_kwargs = (; interpolation_method),
        file_reader_kwargs = (; preprocess_func = clean_missing),
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
        [crujra_ncdata_path, crujra_ncdata_path],
        ["msdwswrf", "msdrswrf"],
        surface_space;
        start_date,
        regridder_type,
        regridder_kwargs = (; interpolation_method),
        file_reader_kwargs = (; preprocess_func = clean_missing),
        method = time_interpolation_method,
        compose_function = compute_diffuse_fraction_broadcasted,
    )
    LW_d = TimeVaryingInput(
        crujra_ncdata_path,
        "msdwlwrf",
        surface_space;
        start_date,
        regridder_type,
        regridder_kwargs = (; interpolation_method),
        file_reader_kwargs = (; preprocess_func = clean_missing),
        method = time_interpolation_method,
    )
    # Cycle the 6-hourly records first, so that the time handling below sees
    # the cycled sequence of records (after the last cycled year, the first one)
    precip,
    snow_precip,
    u_atmos,
    q_atmos,
    P_atmos,
    T_atmos,
    SW_d,
    frac_diff,
    LW_d = cycle.((
        precip,
        snow_precip,
        u_atmos,
        q_atmos,
        P_atmos,
        T_atmos,
        SW_d,
        frac_diff,
        LW_d,
    ))
    latitude = ClimaCore.Fields.coordinate_field(surface_space).lat
    longitude = ClimaCore.Fields.coordinate_field(surface_space).long
    insol_params = earth_param_set.insol_params
    if state_time_interpolation != :linear
        method =
            state_time_interpolation == :linear_midpoint ? :linear :
            state_time_interpolation
        state(input; nonnegative = false) = HoldWindowStateInput(
            input,
            start_date,
            zero(latitude);
            updateat,
            method,
            nonnegative,
        )
        T_atmos = state(T_atmos)
        q_atmos = state(q_atmos; nonnegative = true)
        P_atmos = state(P_atmos)
        u_atmos = state(u_atmos; nonnegative = true)
    end
    if flux_time_interpolation == :interval
        interval_mean(input) =
            IntervalMeanInput(input, start_date, zero(latitude); updateat)
        precip = interval_mean(precip)
        snow_precip = interval_mean(snow_precip)
        LW_d = interval_mean(LW_d)
    end
    atmos = PrescribedAtmosphere(
        precip,
        snow_precip,
        T_atmos,
        u_atmos,
        q_atmos,
        P_atmos,
        start_date,
        h_atmos,
        toml_dict;
        gustiness,
        reference,
        c_co2 = c_co2,
    )
    if sw_time_interpolation == :coszen
        SW_d = CosZenithWeightedInput(
            SW_d,
            start_date,
            latitude,
            longitude,
            insol_params;
            updateat,
            quadrature_step = sw_quadrature_step,
        )
        frac_diff = IntervalStartInput(frac_diff, start_date)
        cos_zenith_angle = WindowCosZenith(SW_d)
    else
        cos_zenith_angle =
            (t, s) -> default_cos_zenith_angle(
                t,
                s;
                latitude,
                longitude,
                insol_params,
            )
    end

    radiation = PrescribedRadiativeFluxes(
        FT,
        SW_d,
        LW_d,
        start_date;
        cosθs = cos_zenith_angle,
        toml_dict = toml_dict,
        frac_diff = frac_diff,
    )
    return (; atmos, radiation)
end

"""
    CycledYearsInput{I, D} <: AbstractTimeVaryingInput

A time-varying input that evaluates `input` in a cycle of years: at a date of
year `y`, `input` is evaluated on the same day and time of year
`first(years) + mod(y - first(years), length(years))`. February 29 uses
February 28 when the cycled year is not a leap year.

Used for the pre-industrial climate of the TRENDY protocol, which cycles the
forcing of 1901–1920 (see `prescribed_forcing_crujra`).
$(DocStringExtensions.FIELDS)
"""
struct CycledYearsInput{I <: AbstractTimeVaryingInput, D} <:
       AbstractTimeVaryingInput
    "Time-varying input that covers `years`"
    input::I
    "Date of `t = 0` of the simulation"
    start_date::D
    "Years that are cycled"
    years::UnitRange{Int}
end

"""
    cycled_date(date, years)

Return the date in the cycle of `years` that corresponds to `date`, as defined
in [`CycledYearsInput`](@ref).
"""
function cycled_date(date::DateTime, years)
    y = Dates.year(date)
    cycled_year = first(years) + mod(y - first(years), length(years))
    return date + Dates.Year(cycled_year - y)
end

function TimeVaryingInputs.evaluate!(
    dest,
    input::CycledYearsInput,
    date::DateTime,
    args...;
    kwargs...,
)
    TimeVaryingInputs.evaluate!(
        dest,
        input.input,
        cycled_date(date, input.years),
        args...;
        kwargs...,
    )
    return nothing
end

TimeVaryingInputs.evaluate!(
    dest,
    input::CycledYearsInput,
    t::Union{Number, ITime},
    args...;
    kwargs...,
) = TimeVaryingInputs.evaluate!(
    dest,
    input,
    input.start_date + Dates.Millisecond(round(Int64, 1000 * float(t))),
    args...;
    kwargs...,
)


"""
    simulation_date(start_date, t)

Return the date of the simulation time `t`, given as seconds since
`start_date`, as an `ITime`, or as a `DateTime`.
"""
simulation_date(start_date, t::DateTime) = t
simulation_date(start_date, t::ITime) =
    isnothing(t.epoch) ? start_date + t.counter * t.period : date(t)
simulation_date(start_date, t::Number) =
    start_date + Dates.Millisecond(round(Int64, 1000 * float(t)))

"""
    hold_length(updateat)

Return the driver update interval `updateat` (a `Dates.Period`, seconds, or an
`ITime`) as a number of milliseconds.
"""
hold_length(updateat::Dates.Period) = Dates.Millisecond(updateat)
hold_length(updateat::ITime) =
    Dates.Millisecond(updateat.counter * updateat.period)
hold_length(updateat::Real) = Dates.Millisecond(round(Int64, 1000 * updateat))
hold_length(updateat) = error(
    "`updateat` must be a single update interval (a `Dates.Period`, seconds, or an `ITime`)",
)

"""
    IntervalStartInput{I, D} <: AbstractTimeVaryingInput

A time-varying input that is constant over each interval of length `period`
(counted from 00 UTC) and equal to `input` at the start of the interval. Used
for quantities given as means over intervals stamped at their start, such as
the diffuse fraction of the 6-hourly CRU-JRA shortwave radiation.

Construct it with `IntervalStartInput(input, start_date; period = Dates.Hour(6))`.
$(DocStringExtensions.FIELDS)
"""
struct IntervalStartInput{I <: AbstractTimeVaryingInput, D} <:
       AbstractTimeVaryingInput
    "Time-varying input stamped at the start of each interval"
    input::I
    "Date of `t = 0` of the simulation"
    start_date::D
    "Length of the intervals"
    period::Dates.Millisecond
end

IntervalStartInput(input, start_date; period = Dates.Hour(6)) =
    IntervalStartInput(input, start_date, Dates.Millisecond(period))

function TimeVaryingInputs.evaluate!(
    dest,
    input::IntervalStartInput,
    t,
    args...;
    kwargs...,
)
    interval_start = floor(simulation_date(input.start_date, t), input.period)
    TimeVaryingInputs.evaluate!(
        dest,
        input.input,
        interval_start,
        args...;
        kwargs...,
    )
    return nothing
end

"""
    CosZenithWeightedInput{I, D, LA, LO, P, C} <: AbstractTimeVaryingInput

A time-varying input for a solar flux given as means over intervals of length
`period` (counted from 00 UTC) stamped at the start of each interval, such as
the 6-hourly CRU-JRA shortwave radiation, for drivers that are updated every
`updateat` and held in between.

Within an interval, the mean `F̄` is distributed in proportion to the cosine of
the solar zenith angle `μ(s)` (zero when the sun is below the horizon), as in
the `coszen` time interpolation of the CESM data models (CDEPS):

    F(s) = F̄ μ(s) / μ̄,

where `μ̄` is the mean of `μ` over the interval. Evaluated at time `t`, the
input returns the mean of `F` over the hold window `[t, t + updateat)`, which
may span several intervals. Means are taken with the midpoint rule on a grid of
spacing `gcd(updateat, quadrature_step)` (15 minutes by default). When the
updates fall on this grid (counted from 00 UTC), the windows and the interval
means use the same points, and the energy of each interval is conserved to
rounding error. Otherwise it is conserved up to the quadrature error.

If the sun is below the horizon at all points of an interval (`μ̄ = 0`), for
example in polar night, `F = F̄` over that interval, so that a nonzero mean is
not lost.

The shortwave-weighted mean of `μ` over the last window, `Σ F μ / Σ F` (or the
mean of `μ` if `Σ F = 0`), is stored in `window_cosθs`; see
[`ClimaLand.WindowCosZenith`](@ref).

Construct it with
`CosZenithWeightedInput(input, start_date, latitude, longitude, insol_params;
period = Dates.Hour(6), updateat = Dates.Hour(3), quadrature_step =
Dates.Minute(15))`, where `input` evaluated at the start of an interval gives
its mean, `latitude` and `longitude` (degrees) are fields or arrays shaped like
the destination, and `updateat` is the driver update interval of the simulation.
A warning is shown once if the input is evaluated at a spacing other than
`updateat`.
$(DocStringExtensions.FIELDS)
"""
struct CosZenithWeightedInput{I <: AbstractTimeVaryingInput, D, LA, LO, P, C} <:
       AbstractTimeVaryingInput
    "Time-varying input of the interval means, stamped at the start of each interval"
    input::I
    "Date of `t = 0` of the simulation"
    start_date::D
    "Latitude (degrees)"
    latitude::LA
    "Longitude (degrees)"
    longitude::LO
    "Insolation parameters"
    insol_params::P
    "Length of the intervals"
    period::Dates.Millisecond
    "Hold window: time between driver updates"
    updateat::Dates.Millisecond
    "Spacing of the quadrature points"
    quadrature_step::Dates.Millisecond
    "Cache: mean flux of the current interval"
    interval_mean::C
    "Cache: mean cosine of the zenith angle over the current interval"
    interval_mean_cosθs::C
    "Scratch: cosine of the zenith angle at a quadrature point"
    cosθs::C
    "Scratch: sum of the cosine of the zenith angle over the window"
    window_sum_cosθs::C
    "Shortwave-weighted mean cosine of the zenith angle over the last window"
    window_cosθs::C
    "Scratch destination used by `WindowCosZenith`"
    scratch::C
    "Start of the interval of the cached values"
    interval_start::Base.RefValue{DateTime}
    "Start of the last evaluated window"
    window_start::Base.RefValue{DateTime}
end

function CosZenithWeightedInput(
    input,
    start_date,
    latitude,
    longitude,
    insol_params;
    period = Dates.Hour(6),
    updateat = Dates.Hour(3),
    quadrature_step = Dates.Minute(15),
)
    period = Dates.Millisecond(period)
    updateat = hold_length(updateat)
    quadrature_step = Dates.Millisecond(quadrature_step)
    Dates.value(updateat) > 0 || error("`updateat` must be positive")
    Dates.value(quadrature_step) > 0 &&
    Dates.value(period) % Dates.value(quadrature_step) == 0 ||
        error("`quadrature_step` must divide the interval length $period")
    step = Dates.Millisecond(
        gcd(Dates.value(updateat), Dates.value(quadrature_step)),
    )
    caches = ntuple(_ -> zero(latitude), 6)
    return CosZenithWeightedInput(
        input,
        start_date,
        latitude,
        longitude,
        insol_params,
        period,
        updateat,
        step,
        caches...,
        Ref(typemin(DateTime)),
        Ref(typemin(DateTime)),
    )
end

cos_zenith(date, latitude, longitude, insol_params) =
    Insolation.insolation(date, latitude, longitude, insol_params).μ

cos_zenith_weight(cosθs, mean_cosθs) =
    ifelse(mean_cosθs > 0, cosθs / max(mean_cosθs, eps(mean_cosθs)), one(cosθs))

quadrature_offset(input::CosZenithWeightedInput) =
    Dates.Millisecond(div(Dates.value(input.quadrature_step), 2))

function update_interval!(
    input::CosZenithWeightedInput,
    interval_start,
    args...;
    kwargs...,
)
    (; latitude, longitude, insol_params, interval_mean_cosθs) = input
    TimeVaryingInputs.evaluate!(
        input.interval_mean,
        input.input,
        interval_start,
        args...;
        kwargs...,
    )
    n = div(input.period, input.quadrature_step)
    offset = quadrature_offset(input)
    interval_mean_cosθs .= 0
    for i in 0:(n - 1)
        date_i = interval_start + i * input.quadrature_step + offset
        @. interval_mean_cosθs +=
            cos_zenith(date_i, latitude, longitude, insol_params)
    end
    interval_mean_cosθs ./= n
    input.interval_start[] = interval_start
    return nothing
end

function check_update_spacing!(input::CosZenithWeightedInput, window_start)
    last_start = input.window_start[]
    if last_start != typemin(DateTime) &&
       window_start > last_start &&
       window_start - last_start != input.updateat
        @warn "CosZenithWeightedInput: drivers updated $(window_start - last_start) apart, but `updateat` is $(input.updateat); the shortwave energy is not conserved" maxlog =
            1
    end
    if floor(window_start, input.quadrature_step) != window_start
        @warn "CosZenithWeightedInput: update at $window_start is not on the $(input.quadrature_step) quadrature grid; the shortwave energy is conserved only approximately" maxlog =
            1
    end
    input.window_start[] = window_start
    return nothing
end

function TimeVaryingInputs.evaluate!(
    dest,
    input::CosZenithWeightedInput,
    t,
    args...;
    kwargs...,
)
    (; latitude, longitude, insol_params, cosθs) = input
    (; interval_mean, interval_mean_cosθs) = input
    (; window_sum_cosθs, window_cosθs) = input
    window_start = simulation_date(input.start_date, t)
    check_update_spacing!(input, window_start)
    n = div(input.updateat, input.quadrature_step)
    offset = quadrature_offset(input)
    dest .= 0
    window_cosθs .= 0
    window_sum_cosθs .= 0
    for j in 0:(n - 1)
        date_j = window_start + j * input.quadrature_step + offset
        interval_start = floor(date_j, input.period)
        interval_start == input.interval_start[] ||
            update_interval!(input, interval_start, args...; kwargs...)
        @. cosθs = cos_zenith(date_j, latitude, longitude, insol_params)
        @. dest += interval_mean * cos_zenith_weight(cosθs, interval_mean_cosθs)
        @. window_cosθs +=
            interval_mean *
            cos_zenith_weight(cosθs, interval_mean_cosθs) *
            cosθs
        @. window_sum_cosθs += cosθs
    end
    @. window_cosθs = ifelse(
        dest > 0,
        window_cosθs / max(dest, eps(dest)),
        window_sum_cosθs / n,
    )
    dest ./= n
    return nothing
end

"""
    WindowCosZenith{I <: CosZenithWeightedInput} <: Function

The cosine of the solar zenith angle held with the shortwave radiation of a
[`ClimaLand.CosZenithWeightedInput`](@ref): called as `f(t, start_date)`, it
returns the shortwave-weighted mean cosine of the zenith angle over the hold
window that starts at `t`. Pass it as `cosθs` to `PrescribedRadiativeFluxes`
so that the canopy sees the zenith angle of the light it receives.
"""
struct WindowCosZenith{I <: CosZenithWeightedInput} <: Function
    input::I
end

function (f::WindowCosZenith)(t, start_date)
    input = f.input
    window_start = simulation_date(start_date, t)
    window_start == input.window_start[] ||
        TimeVaryingInputs.evaluate!(input.scratch, input, window_start)
    return input.window_cosθs
end

"""
    IntervalMeanInput{I, D, C} <: AbstractTimeVaryingInput

A time-varying input for a flux given as means over intervals of length
`period` (counted from 00 UTC) stamped at the start of each interval, such as
the 6-hourly CRU-JRA longwave radiation and precipitation, for drivers that are
updated every `updateat` and held in between.

The flux is taken as constant within each interval, equal to the interval mean.
Evaluated at time `t`, the input returns the mean of this piecewise-constant
flux over the hold window `[t, t + updateat)`, weighting each interval by its
overlap with the window. Each interval keeps its own amount (energy or water)
exactly, for any `updateat`, as long as the windows tile time.

Construct it with `IntervalMeanInput(input, start_date, scratch; period =
Dates.Hour(6), updateat = Dates.Hour(3))`, where `input` evaluated at the start
of an interval gives its mean and `scratch` is a field or array shaped like the
destination.
$(DocStringExtensions.FIELDS)
"""
struct IntervalMeanInput{I <: AbstractTimeVaryingInput, D, C} <:
       AbstractTimeVaryingInput
    "Time-varying input of the interval means, stamped at the start of each interval"
    input::I
    "Date of `t = 0` of the simulation"
    start_date::D
    "Length of the intervals"
    period::Dates.Millisecond
    "Hold window: time between driver updates"
    updateat::Dates.Millisecond
    "Scratch: mean of one interval"
    scratch::C
end

function IntervalMeanInput(
    input,
    start_date,
    scratch;
    period = Dates.Hour(6),
    updateat = Dates.Hour(3),
)
    updateat = hold_length(updateat)
    Dates.value(updateat) > 0 || error("`updateat` must be positive")
    return IntervalMeanInput(
        input,
        start_date,
        Dates.Millisecond(period),
        updateat,
        scratch,
    )
end

function TimeVaryingInputs.evaluate!(
    dest,
    input::IntervalMeanInput,
    t,
    args...;
    kwargs...,
)
    (; period, updateat, scratch) = input
    window_start = simulation_date(input.start_date, t)
    window_end = window_start + updateat
    interval_start = floor(window_start, period)
    FT = eltype(dest)
    dest .= 0
    while interval_start < window_end
        overlap =
            min(window_end, interval_start + period) -
            max(window_start, interval_start)
        weight = FT(Dates.value(overlap) / Dates.value(updateat))
        TimeVaryingInputs.evaluate!(
            scratch,
            input.input,
            interval_start,
            args...;
            kwargs...,
        )
        @. dest += weight * scratch
        interval_start += period
    end
    return nothing
end

"""
    HoldWindowStateInput{I, D, C} <: AbstractTimeVaryingInput

A time-varying input for a state variable sampled at the start of intervals of
length `period` (counted from 00 UTC), such as the 6-hourly CRU-JRA
temperature, humidity, pressure and wind, which are instantaneous values at
00, 06, 12 and 18 UTC, for drivers that are updated every `updateat` and held
in between.

Evaluated at time `t`, the input returns the value at the middle of the hold
window, `t + updateat / 2`, which approximates the mean over the window. With
`method = :cubic` (the default), the value is the Catmull-Rom cubic through
the four samples around it (at `k - period`, `k`, `k + period` and `k + 2
period`, with `k` the start of the interval), which passes through the samples
and, unlike linear interpolation, does not flatten the diurnal cycle between
them. With `method = :harmonic`, it is the sum of a constant, a diurnal and a
semidiurnal harmonic through the same four samples, `a + b cos θ + c sin θ +
d cos 2θ`, with `θ` the phase in a cycle of four periods (one day for 6-hourly
samples); this follows the diurnal cycle of temperature more closely, but does
not reproduce linear trends between the samples. Samples outside the dates of the data are replaced by the nearest sample
of the interval. With `method = :linear`, the value is interpolated linearly
between the two samples of the interval. If `nonnegative` is true, negative
values (cubic overshoots of humidity or wind speed) are set to zero.

Construct it with `HoldWindowStateInput(input, start_date, scratch; period =
Dates.Hour(6), updateat = Dates.Hour(3), method = :cubic, nonnegative =
false)`, where `input` evaluated at a sample time gives the sample and
`scratch` is a field or array shaped like the destination.
$(DocStringExtensions.FIELDS)
"""
struct HoldWindowStateInput{I <: AbstractTimeVaryingInput, D, C} <:
       AbstractTimeVaryingInput
    "Time-varying input of the samples"
    input::I
    "Date of `t = 0` of the simulation"
    start_date::D
    "Spacing of the samples"
    period::Dates.Millisecond
    "Hold window: time between driver updates"
    updateat::Dates.Millisecond
    "Interpolation between the samples: `:cubic` or `:linear`"
    method::Symbol
    "Whether negative values are set to zero"
    nonnegative::Bool
    "Cache: samples at `k - period`, `k`, `k + period` and `k + 2 period`"
    samples::NTuple{4, C}
    "Start `k` of the interval of the cached samples"
    interval_start::Base.RefValue{DateTime}
end

function HoldWindowStateInput(
    input,
    start_date,
    scratch;
    period = Dates.Hour(6),
    updateat = Dates.Hour(3),
    method = :cubic,
    nonnegative = false,
)
    method in (:cubic, :harmonic, :linear) ||
        error("`method` must be `:cubic`, `:harmonic` or `:linear`")
    updateat = hold_length(updateat)
    Dates.value(updateat) > 0 || error("`updateat` must be positive")
    return HoldWindowStateInput(
        input,
        start_date,
        Dates.Millisecond(period),
        updateat,
        method,
        nonnegative,
        ntuple(_ -> zero(scratch), 4),
        Ref(typemin(DateTime)),
    )
end

in_data_range(input, date) =
    hasproperty(input, :data_handler) ? date in input : true
in_data_range(input::CycledYearsInput, date) =
    in_data_range(input.input, cycled_date(date, input.years))

function update_samples!(
    input::HoldWindowStateInput,
    interval_start,
    args...;
    kwargs...,
)
    (; period, samples) = input
    dates = interval_start .+ period .* (-1, 0, 1, 2)
    # Outside the data, fall back on the nearest sample of the interval
    fallback = (dates[2], dates[2], dates[3], dates[3])
    for i in 1:4
        date = in_data_range(input.input, dates[i]) ? dates[i] : fallback[i]
        TimeVaryingInputs.evaluate!(
            samples[i],
            input.input,
            date,
            args...;
            kwargs...,
        )
    end
    input.interval_start[] = interval_start
    return nothing
end

function TimeVaryingInputs.evaluate!(
    dest,
    input::HoldWindowStateInput,
    t,
    args...;
    kwargs...,
)
    (; period, samples) = input
    midpoint =
        simulation_date(input.start_date, t) +
        Dates.Millisecond(div(Dates.value(input.updateat), 2))
    interval_start = floor(midpoint, period)
    interval_start == input.interval_start[] ||
        update_samples!(input, interval_start, args...; kwargs...)
    FT = eltype(dest)
    s = FT(Dates.value(midpoint - interval_start) / Dates.value(period))
    p0, p1, p2, p3 = samples
    if input.method == :cubic
        # Catmull-Rom weights
        w0 = (-s^3 + 2s^2 - s) / 2
        w1 = (3s^3 - 5s^2 + 2) / 2
        w2 = (-3s^3 + 4s^2 + s) / 2
        w3 = (s^3 - s^2) / 2
        @. dest = w0 * p0 + w1 * p1 + w2 * p2 + w3 * p3
    elseif input.method == :harmonic
        # a + b cos θ + c sin θ + d cos 2θ through the samples at θ = -π/2, 0,
        # π/2 and π (a 6-hour period is a quarter of the diurnal cycle)
        sinθ, cosθ = sincospi(s / 2)
        cos2θ = cospi(s)
        w0 = (1 - 2sinθ - cos2θ) / 4
        w1 = (1 + 2cosθ + cos2θ) / 4
        w2 = (1 + 2sinθ - cos2θ) / 4
        w3 = (1 - 2cosθ + cos2θ) / 4
        @. dest = w0 * p0 + w1 * p1 + w2 * p2 + w3 * p3
    else
        @. dest = (1 - s) * p1 + s * p2
    end
    input.nonnegative && (@. dest = max(dest, zero(dest)))
    return nothing
end

"""
    read_annual_co2(filepath)

Read a text file of annual global mean atmospheric CO₂ and return the years and
the CO₂ mole fractions (mol/mol).

Each line holds a year and a CO₂ value in ppm separated by white space, as in
the `global_co2_ann_1700_2025.txt` file of the TRENDY/Global Carbon Budget
protocol. Empty lines and lines starting with `#` are skipped. The years must be
consecutive.
"""
function read_annual_co2(filepath)
    years = Int[]
    co2 = Float64[]
    for line in eachline(filepath)
        stripped = strip(line)
        (isempty(stripped) || startswith(stripped, '#')) && continue
        parts = split(stripped)
        push!(years, parse(Int, parts[1]))
        push!(co2, parse(Float64, parts[2]) * 1e-6)
    end
    isempty(years) && error("No CO₂ values found in $filepath")
    years == first(years):last(years) ||
        error("The years in $filepath are not consecutive")
    return years, co2
end

"""
    prescribed_co2_annual(filepath, start_date; interpolation = :step, fixed_year = nothing)

Return a `TimeVaryingInput` of the atmospheric CO₂ mole fraction (mol/mol) from a
file of annual global means (see [`read_annual_co2`](@ref)), for the `c_co2`
keyword of `PrescribedAtmosphere`, `prescribed_forcing_era5` and
`prescribed_forcing_crujra`. `start_date` is the date of `t = 0` of the
simulation.

With `interpolation = :step` (the default), the CO₂ of a year is used from
January 1 to December 31 of that year, so the mean over each year equals the
value in the file, as in the TRENDY protocol, which gives one value per year.
With `interpolation = :linear`, the values are placed at the middle of each year
and interpolated linearly in between. In both cases, dates before the first year
or after the last year of the file use the first or the last value.

With `fixed_year`, the CO₂ of that year is used at all times, for example
`fixed_year = 1700` for the pre-industrial CO₂ of the TRENDY simulations S0, S4
and S6.
"""
function prescribed_co2_annual(
    filepath,
    start_date::Dates.DateTime;
    interpolation::Symbol = :step,
    fixed_year::Union{Nothing, Integer} = nothing,
)
    years, co2 = read_annual_co2(filepath)
    first_year, last_year = first(years), last(years)
    if !isnothing(fixed_year)
        first_year <= fixed_year <= last_year ||
            error("Year $fixed_year is not in $filepath")
        fixed_co2 = co2[fixed_year - first_year + 1]
        return TimeVaryingInput((t) -> fixed_co2)
    end
    if interpolation == :step
        function co2_of_year(t)
            date = start_date + Dates.Millisecond(round(Int64, 1000 * float(t)))
            year = clamp(Dates.year(date), first_year, last_year)
            return co2[year - first_year + 1]
        end
        return TimeVaryingInput(co2_of_year)
    elseif interpolation == :linear
        mid_year(y) =
            Dates.DateTime(y) + (Dates.DateTime(y + 1) - Dates.DateTime(y)) ÷ 2
        times = [
            Dates.value(Dates.Millisecond(mid_year(y) - start_date)) / 1000
            for y in years
        ]
        return TimeVaryingInput(
            times,
            co2;
            method = LinearInterpolation(Flat()),
        )
    else
        error("interpolation must be :step or :linear, got :$interpolation")
    end
end
