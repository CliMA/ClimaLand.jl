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
    energy_from_temperature_and_swe

"""
    ClimaLSM.surface_height(
        model::SnowModel{FT},
        Y,
        p,
    ) where {FT}

Returns the surface height of the `Snow` model.
"""
function ClimaLSM.surface_height(model::SnowModel{FT}, Y, p) where {FT}
    z_sfc = ClimaCore.Fields.coordinate_field(model.domain.space.surface).z
    return z_sfc
end


"""
    surface_albedo(model::SnowModel, Y, p)

A helper function which computes and returns the snow albedo.
"""
function ClimaLSM.surface_albedo(model::SnowModel, _...)
    return model.parameters.α_snow
end

"""
    surface_emissivity(model::SnowModel, Y, p)

A helper function which computes and returns the snow emissivity.
"""
function ClimaLSM.surface_emissivity(model::SnowModel, _...)
    return model.parameters.ϵ_snow
end

"""
    ClimaLSM.surface_temperature(model::SnowModel, Y, p)

a helper function which returns the surface temperature for the snow 
model, which is stored in the aux state.
"""
function ClimaLSM.surface_temperature(model::SnowModel, Y, p, t)
    return p.snow.T_sfc
end

"""
    ClimaLSM.surface_specific_humidity(model::BucketModel, Y, p)

a helper function which returns the surface specific
humidity for the snow
model, which is stored in the aux state.
"""
function ClimaLSM.surface_specific_humidity(model::SnowModel, Y, p, _...)
    return p.snow.q_sfc
end

"""
    ClimaLSM.surface_air_density(
        atmos::PrescribedAtmosphere,
        model::SnowModel,
        Y,
        p,
        t,
        T_sfc,
    )
a helper function which computes and returns the surface air density for the
snow model.
"""
function ClimaLSM.surface_air_density(
    atmos::PrescribedAtmosphere,
    model::SnowModel,
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
    ClimaLSM.surface_temperature(model::SnowModel, Y, p)

a helper function which computes and returns the surface 
evaporative scaling factor for the snow model.
"""
function ClimaLSM.surface_evaporative_scaling(
    model::SnowModel{FT},
    _...,
) where {FT}
    return FT(1.0)
end


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
function snow_depth(SWE::FT, ρ_snow::FT, ρ_l::FT)::FT where {FT}
    return SWE * ρ_l / ρ_snow
end



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
    z = snow_depth(SWE, parameters.ρ_snow, _ρ_l)
    if z < parameters.z_crit
        return _T_freeze # Replace with T_ground: then this is like the bucket model
    else
        return _T_ref + (U / (_ρ_l * SWE) + (1 - q_l) * _LH_f0) / cp_s
    end

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
    _ρ_i = FT(LSMP.ρ_cloud_ice(parameters.earth_param_set))
    _ρ_l = FT(LSMP.ρ_cloud_liq(parameters.earth_param_set))
    _T_ref = FT(LSMP.T_0(parameters.earth_param_set))
    _T_freeze = FT(LSMP.T_freeze(parameters.earth_param_set))
    _LH_f0 = FT(LSMP.LH_f0(parameters.earth_param_set))
    _cp_i = FT(LSMP.cp_i(parameters.earth_param_set))
    _cp_l = FT(LSMP.cp_l(parameters.earth_param_set))

    ΔT = _T_freeze - _T_ref
    z = snow_depth(SWE, parameters.ρ_snow, _ρ_l)
    if z < parameters.z_crit
        return FT(0) # Not sure what to put here, but this is equivalent to the bucket model.
    else
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
)::FT where {FT}
    _ρ_l = FT(LSMP.ρ_cloud_liq(parameters.earth_param_set))
    _T_freeze = FT(LSMP.T_freeze(parameters.earth_param_set))

    if T > _T_freeze
        return FT(0)
    else
        return parameters.θ_r * _ρ_l / ρ_snow
    end
end

"""
    energy_at_maximum_q_l(S::FT, q_l_max::FT, parameters)::FT where {FT}

Helper function which computes the energy per unit area, U, given the current 
snow water equivalent S, at the maximum liquid mass fraction permitted by the
snowpack.
"""
function energy_at_maximum_q_l(S::FT, q_l_max::FT, parameters)::FT where {FT}
    _ρ_l = FT(LSMP.ρ_cloud_liq(parameters.earth_param_set))
    _T_freeze = FT(LSMP.T_freeze(parameters.earth_param_set))
    _T_ref = FT(LSMP.T_0(parameters.earth_param_set))
    _LH_f0 = FT(LSMP.LH_f0(parameters.earth_param_set))

    c_snow = specific_heat_capacity(q_l_max, parameters)
    return _ρ_l * S * (c_snow * (_T_freeze - _T_ref) - (1 - q_l_max) * _LH_f0)
end


"""
    runoff_timescale(z::FT, Ksat::FT, Δt::FT) where {FT}
>>>>>>> 6b1c3289 (snow rhs)

Computes the timescale for liquid water to percolate and leave the snowpack,
given the depth of the snowpack z and the hydraulic conductivity Ksat.
"""
function runoff_timescale(z::FT, Ksat::FT, Δt) where {FT}
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
    _ρ_l = FT(LSMP.ρ_cloud_liq(parameters.earth_param_set))
    _T_freeze = FT(LSMP.T_freeze(parameters.earth_param_set))
    _cp_l = FT(LSMP.cp_l(parameters.earth_param_set))
    _T_ref = FT(LSMP.T_0(parameters.earth_param_set))

    I_liq = _ρ_l * _cp_l * (_T_freeze .- _T_ref)
    return I_liq
end


"""
    compute_energy_runoff(S::FT, q_l::FT, T::FT, parameters) where {FT}

Computes the rate of change in the snow water equivalent S due to loss of
liquid water (runoff) from the snowpack.

Runoff occurs as the snow melts and exceeds the water holding capacity.
"""
function compute_water_runoff(S::FT, q_l::FT, T::FT, parameters) where {FT}
    _ρ_l = FT(LSMP.ρ_cloud_liq(parameters.earth_param_set))
    ρ_snow::FT = parameters.ρ_snow
    Ksat::FT = parameters.Ksat
    Δt::FT = parameters.Δt
    depth = snow_depth(S, ρ_snow, _ρ_l)
    τ = runoff_timescale(depth, Ksat, Δt)
    q_l_max::FT = maximum_liquid_mass_fraction(T, ρ_snow, parameters)
    return -(q_l - q_l_max) * heaviside(q_l - q_l_max) / τ * S /
           (1 - q_l + eps(FT))
end


"""
    energy_from_temperature_and_swe(S::FT, T::FT, parameters) where {FT}

A helper function for computing the snow energy per unit area, given snow 
water equivalent S, bulk temperature T, and snow model parameters.
"""
function energy_from_temperature_and_swe(S::FT, T::FT, parameters) where {FT}
    _T_freeze = FT(LSMP.T_freeze(parameters.earth_param_set))
    if T > _T_freeze
        throw(
            AssertionError(
                "Snow temperature must be below the freezing point if completely frozen",
            ),
        )
    end

    _ρ_l = FT(LSMP.ρ_cloud_liq(parameters.earth_param_set))
    _T_ref = FT(LSMP.T_0(parameters.earth_param_set))
    _LH_f0 = FT(LSMP.LH_f0(parameters.earth_param_set))

    c_snow = specific_heat_capacity(FT(0), parameters)

    return _ρ_l * S * (c_snow * (T - _T_ref) - _LH_f0)
end

"""
    energy_from_q_l_and_swe(S::FT, q_l::FT, parameters) where {FT}

A helper function for compute the snow energy per unit area, given snow 
water equivalent S, liquid fraction q_l,  and snow model parameters.

Note that liquid water can only exist at the freezing point in this model,
so temperature is not required as an input.
"""
function energy_from_q_l_and_swe(S::FT, q_l::FT, parameters) where {FT}
    _T_freeze = FT(LSMP.T_freeze(parameters.earth_param_set))
    _ρ_l = FT(LSMP.ρ_cloud_liq(parameters.earth_param_set))
    _T_ref = FT(LSMP.T_0(parameters.earth_param_set))
    _LH_f0 = FT(LSMP.LH_f0(parameters.earth_param_set))

    c_snow = specific_heat_capacity(q_l, parameters)

    return _ρ_l * S * (c_snow * (_T_freeze - _T_ref) - (1 - q_l) * _LH_f0)
end
