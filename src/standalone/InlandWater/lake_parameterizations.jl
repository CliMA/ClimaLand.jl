
# ─── Lake physics parameterizations ──────────────────────────────────────────
# Pure physics functions that don't depend on the model struct.

"""
    lake_energy_at_freezing(q_l::FT, depth::FT, params::SlabLakeParameters{FT}) where {FT}

Returns the slab energy per unit area (J m⁻²) at the freezing point for a
given liquid fraction `q_l`. When `q_l = 0` the slab is fully frozen; when
`q_l = 1` it is fully liquid at the freezing temperature.
"""
function lake_energy_at_freezing(
    q_l::FT,
    depth::FT,
    params::SlabLakeParameters{FT},
) where {FT}
    earth_param_set = params.earth_param_set
    _ρ_l = LP.ρ_cloud_liq(earth_param_set)
    _T_ref = LP.T_0(earth_param_set)
    _T_freeze = LP.T_freeze(earth_param_set)
    _LH_f0 = LP.LH_f0(earth_param_set)
    _cp_i = LP.cp_i(earth_param_set)
    _cp_l = LP.cp_l(earth_param_set)
    cp = q_l * _cp_l + (1 - q_l) * _cp_i
    return _ρ_l * depth * cp * (_T_freeze - _T_ref) -
           _ρ_l * depth * (1 - q_l) * _LH_f0
end

"""
    lake_liquid_fraction(U::FT, depth::FT, params::SlabLakeParameters{FT}) where {FT}

Diagnoses the liquid fraction `q_l ∈ [0, 1]` from the prognostic slab
energy `U` (J m⁻²). Returns 0 when the lake is fully frozen, 1 when fully
liquid, and interpolates linearly during the phase-change interval.
"""
function lake_liquid_fraction(
    U::FT,
    depth::FT,
    params::SlabLakeParameters{FT},
) where {FT}
    U_ice = lake_energy_at_freezing(FT(0), depth, params)
    U_liq = lake_energy_at_freezing(FT(1), depth, params)
    if U <= U_ice
        return FT(0)
    elseif U >= U_liq
        return FT(1)
    else
        return (U - U_ice) / (U_liq - U_ice)
    end
end

"""
    lake_energy_from_temperature(T::FT, depth::FT, params::SlabLakeParameters{FT}) where {FT}

Returns the slab energy per unit area (J m⁻²) for a given temperature `T` (K).
Below the freezing point the slab is treated as pure ice; above it, as pure
liquid.
"""
function lake_energy_from_temperature(
    T::FT,
    depth::FT,
    params::SlabLakeParameters{FT},
) where {FT}
    earth_param_set = params.earth_param_set
    _ρ_l = LP.ρ_cloud_liq(earth_param_set)
    _T_ref = LP.T_0(earth_param_set)
    _T_freeze = LP.T_freeze(earth_param_set)
    _LH_f0 = LP.LH_f0(earth_param_set)
    _cp_i = LP.cp_i(earth_param_set)
    _cp_l = LP.cp_l(earth_param_set)
    if T <= _T_freeze
        return _ρ_l * depth * _cp_i * (T - _T_ref) -
               _ρ_l * depth * _LH_f0
    else
        return _ρ_l * depth * _cp_l * (T - _T_ref)
    end
end

"""
    lake_temperature(U::FT, q_l::FT, depth::FT, params::SlabLakeParameters{FT}) where {FT}

Diagnoses the lake temperature (K) from the prognostic slab energy `U`
(J m⁻²) and liquid fraction `q_l`. During the phase-change interval
(`0 < q_l < 1`) the temperature is clamped to the freezing point.
"""
function lake_temperature(
    U::FT,
    q_l::FT,
    depth::FT,
    params::SlabLakeParameters{FT},
) where {FT}
    earth_param_set = params.earth_param_set
    _ρ_l = LP.ρ_cloud_liq(earth_param_set)
    _T_ref = LP.T_0(earth_param_set)
    _T_freeze = LP.T_freeze(earth_param_set)
    _LH_f0 = LP.LH_f0(earth_param_set)
    _cp_i = LP.cp_i(earth_param_set)
    _cp_l = LP.cp_l(earth_param_set)
    if q_l < eps(FT)
        return _T_ref + (U + _ρ_l * depth * _LH_f0) / (_ρ_l * depth * _cp_i)
    elseif q_l > FT(1) - eps(FT)
        return _T_ref + U / (_ρ_l * depth * _cp_l)
    else
        return _T_freeze
    end
end

"""
    lake_volumetric_internal_energy(T::FT, q_l::FT, earth_param_set) where {FT}

Returns the volumetric internal energy (J m⁻³) of the lake water at
temperature `T` (K) and liquid fraction `q_l`, referenced to T₀. Used to
compute the energy carried away by runoff.
"""
function lake_volumetric_internal_energy(
    T::FT,
    q_l::FT,
    earth_param_set,
) where {FT}
    _ρ_l = LP.ρ_cloud_liq(earth_param_set)
    _T_ref = LP.T_0(earth_param_set)
    _LH_f0 = LP.LH_f0(earth_param_set)
    _cp_i = LP.cp_i(earth_param_set)
    _cp_l = LP.cp_l(earth_param_set)
    cp = q_l * _cp_l + (1 - q_l) * _cp_i
    return _ρ_l * (cp * (T - _T_ref) - (1 - q_l) * _LH_f0)
end

"""
    lake_runoff_energy_flux(runoff::FT, T::FT, q_l::FT, earth_param_set) where {FT}

Returns the energy flux (W m⁻²) carried by lake runoff. The runoff
volume flux (m s⁻¹) is multiplied by the volumetric internal energy
of the lake water at its current temperature and liquid fraction.
"""
function lake_runoff_energy_flux(
    runoff::FT,
    T::FT,
    q_l::FT,
    earth_param_set,
) where {FT}
    return runoff * lake_volumetric_internal_energy(T, q_l, earth_param_set)
end

"""
    lake_surface_albedo(q_l::FT, params::SlabLakeParameters{FT}) where {FT}

Returns the lake surface albedo as a liquid-fraction-weighted blend of
the open-water (`liquid_albedo`) and ice (`ice_albedo`) albedos.
"""
lake_surface_albedo(q_l::FT, params::SlabLakeParameters{FT}) where {FT} =
    q_l * params.liquid_albedo + (FT(1) - q_l) * params.ice_albedo

"""
    lake_sediment_heat_flux(T_lake::FT, T_soil::FT, κ_soil::FT, Δz_soil::FT, params::SlabLakeParameters{FT}) where {FT}

Compute sediment heat flux (W m⁻²) using a series-conductance model:
    Q_sed = -G_eff * (T_lake - T_soil)
where G_eff = 1 / (1/G + Δz_soil/(2*κ_soil)) combines the lake-side
conductance `G` (from `params`) and the half-cell soil conductance.

In standalone runs, this is not computed and set to zero in simulations
"""
function lake_sediment_heat_flux(
    T_lake::FT,
    T_soil::FT,
    κ_soil::FT,
    Δz_soil::FT,
    params::SlabLakeParameters{FT},
) where {FT}
    G = params.conductance
    G_eff = FT(1) / (FT(1) / G + Δz_soil / (FT(2) * κ_soil))
    return -G_eff * (T_lake - T_soil)
end

"""
    lake_specific_humidity(T_sfc::FT, q_l::FT, ρ_sfc::FT, earth_param_set) where {FT}

Returns the surface specific humidity (kg kg⁻¹) over the lake as a
liquid-fraction-weighted blend of the saturation specific humidities
over ice and liquid water at temperature `T_sfc` and surface air
density `ρ_sfc`.
"""
function lake_specific_humidity(
    T_sfc::FT,
    q_l::FT,
    ρ_sfc::FT,
    earth_param_set,
) where {FT}
    thermo_params = LP.thermodynamic_parameters(earth_param_set)
    qsat_ice = Thermodynamics.q_vap_saturation(
        thermo_params,
        T_sfc,
        ρ_sfc,
        Thermodynamics.Ice(),
    )
    qsat_liq = Thermodynamics.q_vap_saturation(
        thermo_params,
        T_sfc,
        ρ_sfc,
        Thermodynamics.Liquid(),
    )
    return (FT(1) - q_l) * qsat_ice + q_l * qsat_liq
end
