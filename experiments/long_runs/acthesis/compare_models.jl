using ClimaLand
using ClimaLand.Snow:
    AbstractDensityModel,
    AbstractAlbedoModel,
    SnowModel,
    snow_bulk_density
import ClimaComms
import ClimaLand.Snow:
    extra_prog_vars,
    extra_prog_types,
    extra_prog_domain_names,
    update_density_and_depth!,
    update_snow_albedo!,
    compute_extra_prog_tendency!

"""
    Anderson1976{FT <: AbstractFloat} <: AbstractDensityModel{FT}

Establishes the density parameterization where snow density
compacts according to the seminal works of Anderson in the 1970s (and specifically the
numerical implementation documented in the Snow17 model).
This uses the default constants defined in the Snow17 documentation, but allows for redefinition of the model constants.
For consistent translation from Snow17 runs and literature, constants are preserved in their original units (the associated `Anderson1976` functions
handle the subsequent unit converstions). As Snow17 evolves unitless relative densities ρ = ρ_snow/ρ_water, some of these constants or rates
are also unitless. Documentation can be found at https://www.weather.gov/media/owp/oh/hrl/docs/22snow17.pdf.
"""
Base.@kwdef struct Anderson1976{FT <: AbstractFloat} <: AbstractDensityModel{FT}
    "rate of (unitless) density increases from overhead load at 0ᵒC (0.026 hr⁻¹ cm⁻¹; cm of water-equivalent load)"
    c1::FT = FT(0.026)
    "compression constant estimated by Kojima in 1967 (21 cm³/g)"
    c2::FT = FT(21.0)
    "the \"fractional settling\" rate (the change in unitless density, per hour, representing destructive metamorphism) at 0ᵒC for densities less than the critical density `ρ_d` (default value is 0.005 hr⁻¹)"
    c3::FT = FT(0.005)
    "the coefficient of temperature-dependent adjustments to the (unitless) fractional settling (0.1 ᵒC⁻¹)"
    c4::FT = FT(0.10)
    "the multiplicative scaling of the settling rate when no water is present in the snowpack (sources vary on whether this factor should be 0 or 1, default is set to 0)"
    c5::FT = FT(0)
    "the multiplicative scaling of the settling rate when water is present in the snowpack (default is 2)"
    c5_wet::FT = FT(2)
    "the coefficient of unitless-density-dependent adjustments to the (unitless) fractional settling (default value is 23)"
    cx::FT = FT(23)
    "the critical (unitless) density `ρ_d` (default is 0.15)"
    ρ_d::FT = FT(0.15)
    "The scaling factor of how the settling ratios depend on temperature (default is 0.08 ᵒC⁻¹)"
    T_factor::FT = FT(0.08)
    "The maximum possible unitless density (default is 0.45)"
    max_density::FT = FT(0.45)
    "The minimum possible unitless density (default is 0.05)"
    min_density::FT = FT(0.05)
    "The minimum allowable water (in ice form) in a snowpack (default is 1mm; 0.001 m)"
    min_ice_lim::FT = FT(0.001)
    "Rate of how density of new snowfall varies with the temperature (default is 0.0017)"
    new_snowfall_rate::FT = FT(0.0017)
    "The threshold at which the minimum density of new snowfall occurs (default is -15 C)"
    new_snowfall_temp::FT = FT(-15)
    "The power scaling on how the density of new snowfall scales with temperature (default is 1.5)"
    new_snowfall_power::FT = FT(1.5)
end

#Define the additional prognostic variables needed for using this parameterization:
ClimaLand.Snow.extra_prog_vars(m::Anderson1976) = (:Z,)
ClimaLand.Snow.extra_prog_types(m::Anderson1976{FT}) where {FT} = (FT,)
ClimaLand.Snow.extra_prog_domain_names(m::Anderson1976) = (:surface,)

"""
    compact_density(
        density::Anderson1976,
        W_i::FT,
        sliq::FT,
        ρ::FT,
        Δt_t::FT,
        Tsnow::FT
    ) where {FT}

Returns the new compressed density ρ_new from the input ρ under the compaction model defined
by Snow17 for a given `Anderson1976` parameterization configuration.
"""
function compact_density(
    density::Anderson1976,
    W_i::FT,
    sliq::FT,
    ρ::FT,
    Δt_t::FT,
    Tsnow::FT,
)::FT where {FT}
    c5 = (sliq > 0) ? density.c5_wet : density.c5
    cx = (ρ > density.ρ_d) ? density.cx : FT(0)

    #Snow17 includes 0.1 * W_i with W_i in millimeters, but we have 100 * W_i as we provide W_i in meters.
    B = 100 * W_i * density.c1 * Δt_t * exp(density.T_factor * Tsnow - density.c2 * ρ)
    factor = (W_i > density.min_ice_lim) ? (exp(B) - 1) / B : FT(1)
    A = exp(
        density.c3 *
        c5 *
        Δt_t *
        exp(density.c4 * Tsnow - cx * (ρ - density.ρ_d)),
    )
    factor = factor * A

    ρ_new = clamp(ρ*factor, density.min_density, density.max_density)
    return ρ_new
end

"""
    newsnow_density(air_temp::FT)::FT where {FT}
Estimates the density of newly fallen snow as a function of air temperature (in celsius), in relative units (i.e. as a fraction of the density
of water), according to the Snow17 model.
Used in computing the time derivative of snow depth under the Anderson1976 density parameterization.
"""
function newsnow_density(air_temp::FT, density::Anderson1976{FT})::FT where {FT}
    if air_temp < density.new_snowfall_temp
        return density.min_density
    else
        return density.min_density + density.new_snowfall_rate * (air_temp - density.new_snowfall_temp)^density.new_snowfall_rate
    end
end

"""
    newsnow_temp(air_temp::FT)::FT where {FT}
Estimates the temperature of newly fallen snow as a function of air temperature, according to the Snow17 model.
Used in computing the time derivative of snow depth under the Anderson1976 density parameterization.
"""
function newsnow_temp(air_temp::FT)::FT where {FT}
    return (air_temp > 0) ? FT(0) : air_temp
end

"""
    get_dzdt(
        density::Anderson1976{FT},
        airT::FT,
        snowT::FT,
        Psnow::FT,
        S::FT,
        Z::FT,
        q_l::FT,
        dt::FT,
    ) where {FT}

Returns the change in snow depth (rate) given the current model state and the `Anderson1976`
density paramterization, which estimates contributions from new snowfall, as well as compaction
effects from the existing and newly fallen snow.
"""
function get_dzdt(
    density::Anderson1976,
    airT::FT,
    snowT::FT,
    Psnow::FT,
    S::FT,
    Z::FT,
    q_l::FT,
    dt::FT,
)::FT where {FT}
    #change in z =  (contribution from new snowfall; liquid equivalent divided by its density)
    #             + (decrease from how the old snowpack compresses from itself and new snowfall on top)

    #The contribution from new snowfall:
    ρ_newsnow =
        compact_density(
            density,
            Psnow * dt,
            FT(0),
            newsnow_density(airT, density),
            dt/3600, #in hours
            newsnow_temp(airT),
        )
    Δz = Psnow * dt / ρ_newsnow

    #The change in height from the compaction of the older snow:
    W_ice = (1 - q_l) * S
    ρ_est = 
        compact_density(
            density,
            Psnow * dt + W_ice,
            q_l * S,
            W_ice / Z,
            dt/3600, #in hours
            snowT,
        )
    dz_compaction = 
        W_ice < density.min_ice_lim ? -Z : W_ice / ρ_est - Z

    return (dz_compaction + Δz) / dt
end

"""
    anderson_clip_dZdt(S::FT, Z::FT, dSdt::FT, dZdt::FT, Δt::FT, ρ_min_frac::FT)::FT

A helper function which clips the tendency of Z such that
its behavior is consistent with that of S: if all snow melts
within a timestep, we clip the tendency of S so that it does
not become negative, and here we also clip the tendency of Z
so that depth does not become negative. Additionally, if the
tendencies of Z and S are such that we would encounter Z < S
(rho_snow > rho_liq), we also clip the tendency.
"""
function anderson_clip_dZdt(
    S::FT,
    Z::FT,
    dSdt::FT,
    dZdt::FT,
    Δt::FT,
    ρ_min_frac::FT,
)::FT where {FT}
    #This function is simple for now since we let Y.snow.Z also be the per-ground-area value, just like Y.snow.S is
    new_S_ground_area = dSdt * Δt + S
    new_Z = dZdt * Δt + Z #also ground-area presently
    #new_scf = (something), if we pivot to let Y.snow.Z be the per-snow-area value

    new_safe_z = clamp(new_Z, new_S_ground_area, new_S_ground_area / ρ_min_frac) #new_S_ground_area -> (new_S_ground_area / new_scf) if if Y.snow.Z becomes per-snow-area, which can be unstable
    return (new_safe_z - Z) / Δt
end

"""
    update_density_and_depth!(ρ_snow, z_snow, density::Anderson1976, Y, p, earth_param_set)

Updates the snow density and depth in place given the current model state. Extends the update_density_and_depth!
function for the Anderson1976 type.
"""
function update_density_and_depth!(
    ρ_snow,
    z_snow,
    density::Anderson1976,
    Y,
    p,
    earth_param_set,
)
    @. z_snow = clamp(Y.snow.Z, Y.snow.S, Y.snow.S / density.min_density) #assume we are doing z as ground-area
    @. ρ_snow = snow_bulk_density(Y.snow.S, z_snow, earth_param_set) #make sure they are consistent types (ground-area vs snow-area)
end

"""
    compute_extra_prog_tendency!(
        density::Anderson1976,
        model::SnowModel,
        dY,
        Y,
        p,
        t
    )

Updates all prognostic variables associated with density/depth given the current model state and the `Anderson1976`
density paramterization.
"""
function ClimaLand.Snow.compute_extra_prog_tendency!(
    density::Anderson1976,
    model::SnowModel,
    dY,
    Y,
    p,
    t,
)
    dY.snow.Z .= get_dzdt.(
        Ref(density),
        p.drivers.T .- model.parameters.earth_param_set.T_freeze,
        p.snow.T .- model.parameters.earth_param_set.T_freeze,
        abs.(p.drivers.P_snow),
        Y.snow.S,
        p.snow.z_snow,
        p.snow.q_l,
        model.parameters.Δt
    )

    # Now we clip the tendency so that Z stays within approximately physical bounds.
    @. dY.snow.Z = anderson_clip_dZdt(
        Y.snow.S,
        Y.snow.Z,
        dY.snow.S, #assumes dY.snow.S is updated (and clipped) before dY.snow.Z
        dY.snow.Z,
        model.parameters.Δt,
        density.min_density,
    )
end

"""
    HTESSELAlbedoModel{FT <: AbstractFloat} <: AbstractAlbedoModel{FT}
    
Establishes the albedo parameterization where snow albedo evolves according to the HTESSEL model.
This uses the default constants defined in the documentation, but allows for redefinition of the model constants.
Temporal rates have been scaled to work with the expression of time in ClimaLand in seconds.
Documentation can be found at https://www.ecmwf.int/en/elibrary/74297-new-snow-scheme-htessel-description-and-offline-validation.
"""
Base.@kwdef struct HTESSELAlbedoModel{FT <: AbstractFloat} <: AbstractAlbedoModel{FT}
    "The minimum albedo (0-1)"
    α_min::FT = FT(0.5)
    "The maximum albedo (0-1)"
    α_max::FT = FT(0.85)
    "The rate of decay of albedo for the exponential decay (dimensionless, multiplying Δt / 1 day)"
    k_exp::FT = FT(0.24)
    "The rate of decay of albedo for the linear decay (dimensionless, multiplying Δt / 1 day)"
    k_lin::FT = FT(0.008)
    "The critical amount of snowfall (m) to reset the snow albedo"
    P_thresh::FT = FT(0.01)
    "The temperature threshold to distinguish between exponential and linear decay curves"
    T_thresh::FT = FT(-2)
end

#Define the additional prognostic variables needed for using this parameterization:
ClimaLand.Snow.extra_prog_vars(m::HTESSELAlbedoModel) = (:A,)
ClimaLand.Snow.extra_prog_types(m::HTESSELAlbedoModel{FT}) where {FT} = (FT,)
ClimaLand.Snow.extra_prog_domain_names(m::HTESSELAlbedoModel) = (:surface,)

"""
    get_dαdt(
        albedo::HTESSELAlbedoModel,
        α::FT,
        P_snow::FT,
        dSWEdt::FT,
        airT::FT,
        dt::FT,
        earth_param_set,
    )

Establihes the numerical behavior for the albedo growth and decay for an HTESSELAlbedoModel.
The presence of melt (dSWEdt being negative) or an air temperature greater than the model's threshold will
dictate evolution according to an exponential decay curve, and otherwise using a linear decay curve. Any
precipitation will reset the albedo in linear proportion to the amount of snowfall compared to the dictated
precipitation reset threshold.
"""
function get_dαdt(
    albedo::HTESSELAlbedoModel,
    α::FT,
    P_snow::FT,
    dSWEdt::FT,
    snowT::FT,
    dt::FT,
    earth_param_set,
)
    _DT_ = FT(earth_param_set.insol_params.day) #86400 seconds, from earth_param_set
    T_freeze = FT(earth_param_set.T_freeze)
    p_scaled = abs(P_snow) * dt / albedo.P_thresh
    return ifelse(
        p_scaled > 0,
        min(1, p_scaled)*(albedo.α_max - α)/dt,
        ifelse(
            (dSWEdt < 0) || (snowT - T_freeze > albedo.T_thresh),
            (albedo.α_min - α)*(1 - exp(-albedo.k_exp * dt / _DT_)) / dt,
            max(-albedo.k_lin / _DT_, (albedo.α_min - α)/dt)
        )
    )
end

"""
    update_snow_albedo!(α, m::HTESSELAlbedoModel, Y, p, t, earth_param_set,)

Updates the snow albedo in place given the current model state. Dispatched form
specifically for the HTESSELAlbedoModel type.
"""
function ClimaLand.Snow.update_snow_albedo!(
    α,
    m::HTESSELAlbedoModel,
    Y,
    p,
    t,
    earth_param_set,
)
    @. α = clamp(Y.snow.A, m.α_min, m.α_max)
end

"""
    compute_extra_prog_tendency!(albedo::HTESSELAlbedoModel, model::SnowModel, dY, Y, p, t)

Updates all prognostic variables associated with snow albedo given the current model state and the `HTESSELAlbedoModel`
albedo paramterization.
"""
function ClimaLand.Snow.compute_extra_prog_tendency!(
    albedo::HTESSELAlbedoModel,
    model::SnowModel,
    dY,
    Y,
    p,
    t,
)
    dY.snow.A .= get_dαdt.(
        Ref(albedo),
        Y.snow.A,
        p.drivers.P_snow,
        dY.snow.S, #assumes it has been set/clipped by now, which it has been
        p.snow.T,
        model.parameters.Δt,
        model.parameters.earth_param_set
    )
end