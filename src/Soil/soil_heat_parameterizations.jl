## To do: use clima parameters; remove hardwired parameters

export volumetric_internal_energy,
    temperature_from_ρe_int, volumetric_internal_energy_liq
"""
    temperature_from_ρe_int(ρe_int::FT, θ_i::FT, ρc_s::FT) where {FT}

A pointwise function for computing the temperature from the volumetric
internal energy, volumetric ice content, and volumetric heat capacity of
the soil.
"""
function temperature_from_ρe_int(ρe_int::FT, θ_i::FT, ρc_s::FT) where {FT}
    ρ_i = FT(916.7)
    T_ref = FT(273.16)
    LH_f0 = FT(333600.0)
    T = T_ref + (ρe_int + θ_i * ρ_i * LH_f0) / ρc_s
    return T
end

"""
    volumetric_internal_energy(θ_i::FT, ρc_s::FT, T::FT) where {FT}

A pointwise function for computing the volumetric internal energy of the soil,
given the volumetric ice content, volumetric heat capacity, and temperature.
"""
function volumetric_internal_energy(θ_i::FT, ρc_s::FT, T::FT) where {FT}
    ρ_i = FT(916.7)
    LH_f0 = FT(333600.0)
    T_ref = FT(273.16)
    ρe_int = ρc_s * (T - T_ref) - θ_i * ρ_i * LH_f0
    return ρe_int
end

"""
    volumetric_internal_energy_liq(T::FT) where {FT}

A pointwise function for computing the volumetric internal energy
of the liquid water in the soil, given the temperature T.
"""
function volumetric_internal_energy_liq(T::FT) where {FT}
    T_ref = FT(273.16)
    ρ_l = FT(1000.0)
    ρcp_l = FT(4181.0 * ρ_l)
    ρe_int_l = ρcp_l * (T - T_ref)
    return ρe_int_l
end
