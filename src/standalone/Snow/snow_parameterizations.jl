export snow_surface_temperature,
    snow_depth,
    specific_heat_capacity,
    snow_thermal_conductivity,
    snow_bulk_temperature,
    snow_liquid_mass_fraction,
    maximum_liquid_mass_fraction,
    runoff_timescale
"""
    snow_surface_temperature(T::FT) where {FT}

Returns the snow surface temperature assuming it is the same
as the bulk temperature T.
"""
snow_surface_temperature(T::FT) where {FT} = T


"""
    snow_depth(SWE::FT, ρ_snow::FT, ρ_l::FT) where {FT}

Returns the snow depth given SWE, snow density ρ_snow, and
the density of liquid water ρ_l.

"""
snow_depth(SWE::FT, ρ_snow::FT, ρ_l::FT) where {FT} = SWE * ρ_l / ρ_snow


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
    _cp_i = FT(LSMP.cp_i(parameters.earth_param_set))

    _cp_l = FT(LSMP.cp_l(parameters.earth_param_set))

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
"""
function snow_thermal_conductivity(
    ρ_snow::FT,
    parameters::SnowParameters{FT},
) where {FT}
    _κ_air = FT(LSMP.K_therm(parameters.earth_param_set))
    κ_ice = parameters.κ_ice
    return _κ_air +
           (FT(7.75e-5) * ρ_snow + FT(1.105e-6) * ρ_snow^2) * (κ_ice - _κ_air)
end

"""
    snow_bulk_temperature(U::FT,
                          SWE::FT,
                          q_l::FT,
                          c_s::FT,
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
    _ρ_i = FT(LSMP.ρ_cloud_ice(parameters.earth_param_set))
    _ρ_l = FT(LSMP.ρ_cloud_liq(parameters.earth_param_set))
    _T_ref = FT(LSMP.T_0(parameters.earth_param_set))
    _LH_f0 = FT(LSMP.LH_f0(parameters.earth_param_set))
    cp_s = specific_heat_capacity(q_l, parameters)
    _T_freeze = FT(LSMP.T_freeze(parameters.earth_param_set))
    # to prevent dividing by zero
    SWE = max(SWE, eps(FT))
    T = _T_ref + (U / (_ρ_l * SWE) + (1 - q_l) * _LH_f0) / cp_s

end

"""
    snow_liquid_mass_fraction(U::FT, SWE::FT, parameters::SnowParameters{FT}) where {FT}

Computes the snow liquid water mass fraction, given the snow water equivalent SWE,
snow energy per unit area U, and other needed parameters.

If there is no snow (U = SWE = 0), the liquid mass fraction is 1.
"""
function snow_liquid_mass_fraction(
    U::FT,
    SWE::FT,
    parameters::SnowParameters{FT},
) where {FT}
    _ρ_i = FT(LSMP.ρ_cloud_ice(parameters.earth_param_set))
    _ρ_l = FT(LSMP.ρ_cloud_liq(parameters.earth_param_set))
    _T_ref = FT(LSMP.T_0(parameters.earth_param_set))
    _T_freeze = FT(LSMP.T_freeze(parameters.earth_param_set))
    _LH_f0 = FT(LSMP.LH_f0(parameters.earth_param_set))
    _cp_i = FT(LSMP.cp_i(parameters.earth_param_set))
    _cp_l = FT(LSMP.cp_l(parameters.earth_param_set))

    ΔT = _T_freeze - _T_ref
    # to prevent dividing by zero
    SWE = max(SWE, eps(FT))
    energy_per_mass = U / (_ρ_l * SWE)
    if energy_per_mass < _cp_i * ΔT - _LH_f0
        return FT(0)
    elseif energy_per_mass > _cp_l * ΔT
        return FT(1)
    else
        return (energy_per_mass + _LH_f0 - _cp_i * ΔT) /
               ((_cp_l - _cp_i) * ΔT + _LH_f0)
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
    _ρ_l = FT(LSMP.ρ_cloud_liq(parameters.earth_param_set))
    _T_freeze = FT(LSMP.T_freeze(parameters.earth_param_set))

    if T > _T_freeze
        return FT(0)
    else
        return parameters.θ_r * _ρ_l / ρ_snow
    end
end

"""
    runoff_timescale(z::FT, Ksat::FT, Δt) where {FT}

Computes the timescale for liquid water to percolate and leave the snowpack,
given the depth of the snowpack z and the hydraulic conductivity Ksat.
"""
function runoff_timescale(z::FT, Ksat::FT, Δt) where {FT}
    τ = max(Δt, z / Ksat)
    return τ
end
