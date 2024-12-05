export snow_surface_temperature,
    snow_depth,
    specific_heat_capacity,
    snow_thermal_conductivity,
    snow_bulk_temperature,
    snow_liquid_mass_fraction,
    maximum_liquid_mass_fraction,
    runoff_timescale,
    compute_water_runoff,
    energy_from_q_l_and_swe,
    energy_from_T_and_swe,
    snow_cover_fraction

"""
    snow_cover_fraction(x::FT; α = FT(1e-3))::FT where {FT}

Returns the snow cover fraction, assuming it is a heaviside
function at 1e-3 meters.

In the future we can play around with other forms.
"""
function snow_cover_fraction(x::FT; α = FT(1e-3))::FT where {FT}
    return heaviside(x - α)
end

"""
    snow_depth(model::AbstractDensityModel, Y, p, params)

A default function which returns the snow depth through the scaling of SWE by the snow density.
This can be extended for additional types of parameterizations.
"""
function snow_depth(density::AbstractDensityModel{FT}, Y, p, params) where {FT}
    ρ_l = FT(LP.ρ_cloud_liq(params.earth_param_set))
    return @. ρ_l * Y.snow.S / p.snow.ρ_snow
end


"""
    ClimaLand.surface_height(
        model::SnowModel{FT},
        Y,
        p,
    ) where {FT}

Returns the surface height of the `Snow` model; surface (ground) elevation
and ignores snow depth (CLM).

Once topography or land ice is incorporated, this will need to change to
z_sfc + land_ice_depth. Note that land ice can 
be ~1-3 km thick on Greenland/

In order to compute surface fluxes, this cannot be larger than the
height of the atmosphere measurement location (z_atmos > z_land always).

This assumes that the surface elevation is zero. 
"""
function ClimaLand.surface_height(model::SnowModel{FT}, Y, p) where {FT}
    return FT(0)
end


"""
    surface_albedo(model::SnowModel, Y, p)

A helper function which computes and returns the snow albedo.
"""
function ClimaLand.surface_albedo(model::SnowModel, _...)
    return model.parameters.α_snow
end

"""
    surface_emissivity(model::SnowModel, Y, p)

A helper function which computes and returns the snow emissivity.
"""
function ClimaLand.surface_emissivity(model::SnowModel, _...)
    return model.parameters.ϵ_snow
end

"""
    ClimaLand.surface_temperature(model::SnowModel, Y, p)

a helper function which returns the surface temperature for the snow
model, which is stored in the aux state.
"""
function ClimaLand.surface_temperature(model::SnowModel, Y, p, t)
    return p.snow.T_sfc
end


"""
    ClimaLand.surface_specific_humidity(model::SnowModel, Y, p)

Computes and returns the specific humidity over snow as a weighted
fraction of the saturated specific humidity over liquid and frozen
water.
"""
function ClimaLand.surface_specific_humidity(
    model::SnowModel,
    Y,
    p,
    T_sfc,
    ρ_sfc,
)
    thermo_params =
        LP.thermodynamic_parameters(model.parameters.earth_param_set)
    qsat_over_ice =
        Thermodynamics.q_vap_saturation_generic.(
            Ref(thermo_params),
            T_sfc,
            ρ_sfc,
            Ref(Thermodynamics.Ice()),
        )
    qsat_over_liq =
        Thermodynamics.q_vap_saturation_generic.(
            Ref(thermo_params),
            T_sfc,
            ρ_sfc,
            Ref(Thermodynamics.Liquid()),
        )
    q_l = p.snow.q_l
    return @. qsat_over_ice * (1 - q_l) + q_l * (qsat_over_liq)
end



"""
    snow_surface_temperature(T::FT) where {FT}

Returns the snow surface temperature assuming it is the same
as the bulk temperature T.
"""
snow_surface_temperature(T::FT) where {FT} = T


"""
    specific_heat_capacity(q_l::FT,
                           parameters::SnowParameters{FT}
                           ) where {FT}

Computes the specific heat capacity of the snow, neglecting
any contribution from air in the pore spaces, given
the liquid water mass fraction q_l and other parameters.
"""
function specific_heat_capacity(
    q_l::FT,
    parameters::SnowParameters{FT},
) where {FT}
    q_i = 1 - q_l
    _cp_i = FT(LP.cp_i(parameters.earth_param_set))

    _cp_l = FT(LP.cp_l(parameters.earth_param_set))

    cp_s = q_l * _cp_l + q_i * _cp_i
    return cp_s
end

"""
    snow_thermal_conductivity(ρ_snow::FT,
                         parameters::SnowParameters{FT},
                         ) where {FT}

Computes the thermal conductivity, given the density
of the snow, according to Equation
5.33 from Bonan's textbook, which in turn is taken from
Jordan (1991).

We have adjusted the original equation to make the coefficients
non-dimensional by multiplying by the first by x = ρ_ice/ρ_ice
and the second by x², with ρ_ice in kg/m³.

When ρ_snow = ρ_ice, we recover κ_snow = κ_ice.
"""
function snow_thermal_conductivity(
    ρ_snow::FT,
    parameters::SnowParameters{FT},
) where {FT}
    _κ_air = FT(LP.K_therm(parameters.earth_param_set))
    _ρ_ice = FT(LP.ρ_cloud_ice(parameters.earth_param_set))
    κ_ice = parameters.κ_ice
    return _κ_air +
           (FT(0.07) * (ρ_snow / _ρ_ice) + FT(0.93) * (ρ_snow / _ρ_ice)^2) *
           (κ_ice - _κ_air)
end

"""
    snow_bulk_temperature(U::FT,
                          SWE::FT,
                          q_l::FT,
                          parameters::SnowParameters{FT}) where {FT}

Computes the bulk snow temperature from the snow water equivalent SWE,
energy per unit area U, liquid water mass fraction q_l, and specific heat
capacity c_s, along with other needed parameters.

If there is no snow (U = SWE = 0), the bulk temperature is the reference temperature,
which is 273.16K.
"""
function snow_bulk_temperature(
    U::FT,
    SWE::FT,
    q_l::FT,
    parameters::SnowParameters{FT},
) where {FT}
    _ρ_l = FT(LP.ρ_cloud_liq(parameters.earth_param_set))
    _T_ref = FT(LP.T_0(parameters.earth_param_set))
    _LH_f0 = FT(LP.LH_f0(parameters.earth_param_set))
    cp_s = specific_heat_capacity(q_l, parameters)
    _ρcD_g = parameters.ρcD_g
    return _T_ref +
           (U + (1 - q_l) * _LH_f0 * _ρ_l * SWE) / (_ρ_l * SWE * cp_s + _ρcD_g)

end

"""
    snow_bulk_density(SWE::FT, z::FT, parameters::SnowParameters{FT}) where {FT}
Returns the snow density given the current model state when depth and SWE are available.
"""
function snow_bulk_density(
    SWE::FT,
    z::FT,
    parameters::SnowParameters{FT},
)::FT where {FT}
    ρ_l = FT(LP.ρ_cloud_liq(parameters.earth_param_set))
    if SWE <= eps(FT) #if there is no snowpack, aka avoid NaN
        return FT(ρ_l)
    end
    ρ_new = SWE / z * ρ_l
    return ρ_new
end

"""
    snow_liquid_mass_fraction(U::FT, SWE::FT, parameters::SnowParameters{FT}) where {FT}

Computes the snow liquid water mass fraction, given the snow water equivalent SWE,
snow energy per unit area U, and other needed parameters.
"""
function snow_liquid_mass_fraction(
    U::FT,
    SWE::FT,
    parameters::SnowParameters{FT},
) where {FT}
    _ρ_l = FT(LP.ρ_cloud_liq(parameters.earth_param_set))
    _T_ref = FT(LP.T_0(parameters.earth_param_set))
    _T_freeze = FT(LP.T_freeze(parameters.earth_param_set))
    _LH_f0 = FT(LP.LH_f0(parameters.earth_param_set))
    _cp_i = FT(LP.cp_i(parameters.earth_param_set))
    _cp_l = FT(LP.cp_l(parameters.earth_param_set))

    _ρcD_g = parameters.ρcD_g
    Uminus =
        (_ρ_l * SWE * _cp_i + _ρcD_g) * (_T_freeze - _T_ref) -
        _ρ_l * SWE * _LH_f0
    Uplus = (_ρ_l * SWE * _cp_l + _ρcD_g) * (_T_freeze - _T_ref)
    if U < Uminus
        return FT(0)
    elseif U > Uplus
        return FT(1)
    else
        return (U - Uminus) / max((Uplus - Uminus), eps(FT))
    end
end


"""
    maximum_liquid_mass_fraction(T::FT, ρ_snow::FT, parameters::SnowParameters{FT}) where {FT}

Computes the maximum liquid water mass fraction, given the bulk temperature of the snow T,
the density of the snow ρ_snow, and parameters.
"""
function maximum_liquid_mass_fraction(
    T::FT,
    ρ_snow::FT,
    parameters::SnowParameters{FT},
) where {FT}
    _ρ_l = FT(LP.ρ_cloud_liq(parameters.earth_param_set))
    _T_freeze = FT(LP.T_freeze(parameters.earth_param_set))
    if T > _T_freeze
        return FT(0)
    else
        return parameters.θ_r * _ρ_l / ρ_snow
    end
end


"""
    runoff_timescale(z::FT, Ksat::FT, Δt::FT) where {FT}

Computes the timescale for liquid water to percolate and leave the snowpack,
given the depth of the snowpack z and the hydraulic conductivity Ksat.
"""
function runoff_timescale(z::FT, Ksat::FT, Δt::FT) where {FT}
    τ = max(Δt, z / Ksat)
    return τ
end

"""
    volumetric_internal_energy_liq(FT, parameters)

Computes the volumetric internal energy of the liquid water
in the snowpack.

Since liquid water can only exist in the snowpack at
the freezing point, this is a constant.
"""
function volumetric_internal_energy_liq(FT, parameters)
    _ρ_l = FT(LP.ρ_cloud_liq(parameters.earth_param_set))
    _T_freeze = FT(LP.T_freeze(parameters.earth_param_set))
    _cp_l = FT(LP.cp_l(parameters.earth_param_set))
    _T_ref = FT(LP.T_0(parameters.earth_param_set))

    I_liq = _ρ_l * _cp_l * (_T_freeze .- _T_ref)
    return I_liq
end


"""
    compute_energy_runoff(S::FT, q_l::FT, T::FT, parameters) where {FT}

Computes the rate of change in the snow water equivalent S due to loss of
liquid water (runoff) from the snowpack.

Runoff occurs as the snow melts and exceeds the water holding capacity.
"""
function compute_water_runoff(
    S::FT,
    q_l::FT,
    T::FT,
    ρ_snow::FT,
    z::FT,
    parameters,
) where {FT}
    τ = runoff_timescale(z, parameters.Ksat, parameters.Δt)
    q_l_max::FT = maximum_liquid_mass_fraction(T, ρ_snow, parameters)
    return -(q_l - q_l_max) * heaviside(q_l - q_l_max) / τ * S /
           max(1 - q_l, FT(0.01))
end


"""
    energy_from_q_l_and_swe(S::FT, q_l::FT, parameters) where {FT}

A helper function for compute the snow energy per unit area, given snow
water equivalent S, liquid fraction q_l,  and snow model parameters.

Note that liquid water can only exist at the freezing point in this model,
so temperature is not required as an input.
"""
function energy_from_q_l_and_swe(S::FT, q_l::FT, parameters) where {FT}
    _T_freeze = FT(LP.T_freeze(parameters.earth_param_set))
    _ρ_l = FT(LP.ρ_cloud_liq(parameters.earth_param_set))
    _T_ref = FT(LP.T_0(parameters.earth_param_set))
    _LH_f0 = FT(LP.LH_f0(parameters.earth_param_set))

    c_snow = specific_heat_capacity(q_l, parameters)
    _ρcD_g = parameters.ρcD_g
    return _ρ_l * S * (c_snow * (_T_freeze - _T_ref) - (1 - q_l) * _LH_f0) +
           _ρcD_g * (_T_freeze - _T_ref)
end


"""
    energy_from_T_and_swe(S::FT, T::FT, parameters) where {FT}

A helper function for compute the snow energy per unit area, given snow
water equivalent S, bulk temperature T, and snow model parameters.

If T = T_freeze, we return the energy as if q_l = 0.

"""
function energy_from_T_and_swe(S::FT, T::FT, parameters) where {FT}
    _ρ_l = FT(LP.ρ_cloud_liq(parameters.earth_param_set))
    _T_ref = FT(LP.T_0(parameters.earth_param_set))
    _T_freeze = FT(LP.T_freeze(parameters.earth_param_set))
    _LH_f0 = FT(LP.LH_f0(parameters.earth_param_set))
    _cp_i = FT(LP.cp_i(parameters.earth_param_set))
    _cp_l = FT(LP.cp_l(parameters.earth_param_set))
    _ρcD_g = parameters.ρcD_g
    if T <= _T_freeze
        return (_ρ_l * S * _cp_i + _ρcD_g) * (T - _T_ref) - _ρ_l * S * _LH_f0
    else
        T > _T_freeze
        return (_ρ_l * S * _cp_l + _ρcD_g) * (T - _T_ref)
    end

end

"""
    update_density!(density::AbstractDensityModel, params::SnowParameters, Y, p)
Updates the snow density given the current model state. Default for all model types,
can be extended for alternative density paramterizations.
"""
function update_density!(
    density::AbstractDensityModel,
    params::SnowParameters,
    Y,
    p,
)
    p.snow.ρ_snow .=
        snow_bulk_density.(Y.snow.S, snow_depth(density, Y, p, params), params)
end

"""
    update_density!(density::ConstantDensityModel, params::SnowParameters, Y, p)
Extends the update_density! function for the ConstantDensityModel type.
"""
function update_density!(
    density::ConstantDensityModel,
    params::SnowParameters,
    Y,
    p,
)
    p.snow.ρ_snow .= density.ρ_snow
end

"""
    update_density_prog!(density::ConstantDensityModel{FT}, model::SnowModel{FT}, Y, p) where {FT}
Updates all prognostic variables associated with density/depth given the current model state.
This is the default method for all density model types, which can be extended for alternative paramterizations.
"""
function update_density_prog!(
    density::AbstractDensityModel,
    model::SnowModel,
    dY,
    Y,
    p,
)
    return nothing
end

"""
    snow_depth(density::Anderson1976, Y, p, params)

An extension of the `snow_depth` function to the Anderson1976 density parameterization, which includes the prognostic depth variable
and thus does not need to derive snow depth from SWE and density.
This is sufficient to enable dynamics of the auxillary variable ρ_snow without extension of update_density!, and avoids
redundant computations in the computation of runoff.
"""
snow_depth(density::Anderson1976, Y, p, params) = Y.snow.Z

"""
    compact_density(W_i::FT, sliq::FT, ρ::FT, Δt_t::FT, Tsnow::FT, density::Anderson1976) where {FT}
Returns the new compressed density ρ_new from the input ρ under the compaction model defined
by Snow17 for a given `Anderson1976` parameterization configuration.
"""
function compact_density(
    W_i::FT,
    sliq::FT,
    ρ::FT,
    Δt_t::FT,
    Tsnow::FT,
    density::Anderson1976,
)::FT where {FT}
    c5 = (sliq > 0) ? FT(2) : density.c5
    cx = (ρ > density.ρ_d) ? density.cx : FT(0)

    #Snow17 includes 0.1 * W_i with W_i in millimeters, but we have 100 * W_i as we provide W_i in meters.
    B = 100 * W_i * density.c1 * Δt_t * exp(FT(0.08) * Tsnow - density.c2 * ρ)
    factor = (W_i > 0.001) ? (exp(B) - 1) / B : FT(1)
    A = exp(
        density.c3 *
        c5 *
        Δt_t *
        exp(density.c4 * Tsnow - cx * (ρ - density.ρ_d)),
    )
    factor = factor * A

    ρ_new = minimum([ρ * factor, FT(0.45)])
    ρ_new = (ρ_new < FT(0.05)) ? ρ : ρ_new
    return ρ_new
end

"""
    newsnow_density(air_temp::FT)::FT where {FT}
Estimates the density of newly fallen snow as a function of air temperature (in celsius), in relative units (i.e. as a fraction of the density
of water), according to the Snow17 model.
Used in computing the time derivative of snow depth under the Anderson1976 density parameterization.
"""
function newsnow_density(air_temp::FT)::FT where {FT}
    if air_temp < -15
        return FT(0.05)
    else
        return FT(0.05) + FT(0.0017) * FT((air_temp + FT(15.0))^1.5)
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
    dzdt(density::Anderson1976{FT}, model::SnowModel{FT}, Y, p, t) where {FT}
Returns the change in snow depth (rate) given the current model state and the `Anderson1976`
density paramterization, which estimates contributions from new snowfall, as well as compaction
effects from the existing and newly fallen snow.
"""
function dzdt(density::Anderson1976{FT}, model::SnowModel{FT}, Y, p) where {FT}
    Δt_t = model.parameters.Δt / 3600 #in hours (when Δt in SnowParameters is ::FT instead of ::Period)

    #Contribution from new snowfall: (uses parameterization of newsnow_temp, newsnow_density)
    air_temp = p.drivers.T .- model.parameters.earth_param_set.T_freeze
    snowfall = abs.(p.drivers.P_snow) .* model.parameters.Δt
    T_newsnow = newsnow_temp.(p.drivers.T)
    ρ_precip = newsnow_density.(air_temp) #unitless
    ρ_newsnow =
        compact_density.(
            snowfall,
            FT(0),
            ρ_precip,
            Δt_t,
            T_newsnow,
            Ref(density),
        )
    Δz = snowfall ./ ρ_newsnow

    #Estimate resulting change from compaction effects of old snow (and the incoming snow on top):
    W_ice = (1 .- p.snow.q_l) .* Y.snow.S
    ρ_ice = W_ice ./ Y.snow.Z #unitless
    sliq = p.snow.q_l .* Y.snow.S
    Tsnow = p.snow.T .- model.parameters.earth_param_set.T_freeze
    W_use = W_ice .+ snowfall
    parent(ρ_ice)[parent(W_ice) .<= eps(FT)] .= FT(0.1) #any nonzero value for no-snowpack to avoid NaN and return 0.0
    ρ_est = compact_density.(W_use, sliq, ρ_ice, Δt_t, Tsnow, Ref(density))
    dz_compaction = (W_ice ./ ρ_est) .- Y.snow.Z

    return (dz_compaction .+ Δz) ./ model.parameters.Δt
end

"""
    clip_dZdt(S::FT, Z::FT, dSdt::FT, dZdt::FT, Δt::FT)::FT
A helper function which clips the tendency of Z such that
its behavior is consistent with that of S: if all snow melts
within a timestep, we clip the tendency of S so that it does
not become negative, and here we also clip the tendency of Z
so that depth does not become negative. Additionally, if the
tendencies of Z and S are such that we would encounter Z < S
(rho_snow > rho_liq), we also clip the tendency.
"""
function clip_dZdt(S::FT, Z::FT, dSdt::FT, dZdt::FT, Δt::FT)::FT where {FT}
    if (S + dSdt * Δt) <= eps(FT)
        return -Z / Δt
    elseif (Z + dZdt * Δt) < (S + dSdt * Δt)
        return (S - Z) / Δt + dSdt
    else
        return dZdt
    end
end

"""
    update_density_prog!(density::Anderson1976, model::SnowModel, Y, p)
Extends update_density_prog! for the Anderson1976 density paramteterization type.
"""
function update_density_prog!(density::Anderson1976, model::SnowModel, dY, Y, p)
    dY.snow.Z .=
        clip_dZdt.(
            Y.snow.S,
            Y.snow.Z,
            dY.snow.S, #this assumes dY.snow.S is updated (and clipped) before dY.snow.Z
            dzdt(density, model, Y, p),
            model.parameters.Δt,
        )
end
