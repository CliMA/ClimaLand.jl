export volumetric_liquid_fraction,
    effective_saturation,
    matric_potential,
    inverse_matric_potential,
    pressure_head,
    hydraulic_conductivity,
    impedance_factor,
    viscosity_factor,
    dψdϑ,
    dry_soil_layer_thickness,
    soil_resistance,
    soil_tortuosity,
    is_saturated
"""
    volumetric_liquid_fraction(ϑ_l::FT, ν_eff::FT, θ_r::FT) where {FT}

A pointwise function returning the volumetric liquid fraction
given the augmented liquid fraction and the effective porosity.
"""
function volumetric_liquid_fraction(ϑ_l::FT, ν_eff::FT, θ_r::FT) where {FT}
    ϑ_l_safe = max(ϑ_l, θ_r + eps(FT))
    if ϑ_l_safe < ν_eff
        θ_l = ϑ_l_safe
    else
        θ_l = ν_eff
    end
    return θ_l
end

"""
    effective_saturation(porosity::FT, ϑ_l::FT, θr::FT) where {FT}

A point-wise function computing the effective saturation.
"""
function effective_saturation(porosity::FT, ϑ_l::FT, θr::FT) where {FT}
    ϑ_l_safe = max(ϑ_l, θr + eps(FT))
    S_l = (ϑ_l_safe - θr) / (porosity - θr)
    return S_l
end


"""
     matric_potential(cm::vanGenuchten{FT}, S::FT) where {FT}

A point-wise function returning the matric potential, using the
van Genuchten formulation.
"""
function matric_potential(cm::vanGenuchten{FT}, S::FT) where {FT}
    (; α, m, n) = cm
    ψ = -((S^(-FT(1) / m) - FT(1)) * α^(-n))^(FT(1) / n)
    return ψ
end

"""
     inverse_matric_potential(cm::vanGenuchten{FT}, ψ::FT) where {FT}

A point-wise function returning the effective saturation, given
the matric potential, using the
van Genuchten formulation.
"""
function inverse_matric_potential(cm::vanGenuchten{FT}, ψ::FT) where {FT}
    ψ > 0 && error("Matric potential is positive")
    (; α, m, n) = cm
    S = (FT(1) + (α * abs(ψ))^n)^(-m)
    return S
end


"""
     approximate_ψ_S_slope(cm::vanGenuchten)

An estimate of the slope of the absolute value of the logψ-logS curve.
Following Lehmann, Assouline, and Or (2008), we linearize the ψ(S) curve about the inflection point (where d²ψ/dS² = 0, at S = (1+m)^(-m)). 
"""
function approximate_ψ_S_slope(cm::vanGenuchten)
    m = cm.m
    n = cm.n
    return (1 + m) / (n * m * m)
end



"""
    pressure_head(
        cm::vanGenuchten{FT},
        θ_r::FT,
        ϑ_l::FT,
        ν_eff::FT,
        S_s::FT,
    ) where {FT}

A point-wise function returning the pressure head in
variably saturated soil, using the van Genuchten matric potential 
if the soil is not saturated, and
an approximation of the positive pressure in the soil if the
soil is saturated.
"""
function pressure_head(
    cm::vanGenuchten{FT},
    θ_r::FT,
    ϑ_l::FT,
    ν_eff::FT,
    S_s::FT,
) where {FT}
    ϑ_l_safe = max(ϑ_l, θ_r + eps(FT))
    S_l_eff = effective_saturation(ν_eff, ϑ_l_safe, θ_r)
    if S_l_eff <= FT(1.0)
        ψ = matric_potential(cm, S_l_eff)
    else
        ψ = (ϑ_l - ν_eff) / S_s
    end
    return ψ
end

"""
   dψdϑ(cm::vanGenuchten{FT}, ϑ_l, ν, θ_r, S_s)

Computes and returns the derivative of the pressure head 
with respect to ϑ for the van Genuchten formulation.
"""
function dψdϑ(cm::vanGenuchten{FT}, ϑ_l, ν, θ_r, S_s) where {FT}
    ϑ_l_safe = max(ϑ_l, θ_r + eps(FT))
    S = effective_saturation(ν, ϑ_l_safe, θ_r)
    (; α, m, n) = cm
    if S < 1.0
        return FT(
            1.0 / (α * m * n) / (ν - θ_r) *
            (S^(-1 / m) - 1)^(1 / n - 1) *
            S^(-1 / m - 1),
        )
    else
        return 1 / S_s
    end
end


"""
     hydraulic_conductivity(cm::vanGenuchten{FT}, K_sat::FT, S::FT) where {FT}

A point-wise function returning the hydraulic conductivity, using the
van Genuchten formulation.
"""
function hydraulic_conductivity(
    cm::vanGenuchten{FT},
    K_sat::FT,
    S::FT,
) where {FT}
    (; m) = cm
    if S < FT(1)
        K = sqrt(S) * (FT(1) - (FT(1) - S^(FT(1) / m))^m)^FT(2)
    else
        K = FT(1)
    end
    return K * K_sat
end


"""
     matric_potential(cm::BrooksCorey{FT}, S::FT) where {FT}

A point-wise function returning the matric potential, using the
Brooks and Corey formulation.
"""
function matric_potential(cm::BrooksCorey{FT}, S::FT) where {FT}
    (; c, ψb) = cm
    ψ = ψb * S^(-1 / c)
    return ψ
end

"""
     inverse_matric_potential(cm::BrooksCorey{FT}, ψ::FT) where {FT}

A point-wise function returning the effective saturation, given
the matric potential, using the
Brooks and Corey formulation.
"""
function inverse_matric_potential(cm::BrooksCorey{FT}, ψ::FT) where {FT}
    ψ > 0 && error("Matric potential is positive")
    (; c, ψb) = cm
    S = (ψ / ψb)^(-c)
    return S
end


"""
     approximate_ψ_S_slope(cm::BrooksCorey)

The slope of the logψ-logS curve for the Brooks and Corey
model.
"""
function approximate_ψ_S_slope(cm::BrooksCorey)
    return 1 / cm.c
end


"""
   dψdϑ(cm::BrooksCorey{FT}, ϑ, ν, θ_r, S_s)

Computes and returns the derivative of the pressure head 
with respect to ϑ for the Brooks and Corey formulation.
"""
function dψdϑ(cm::BrooksCorey{FT}, ϑ, ν, θ_r, S_s) where {FT}
    S = effective_saturation(ν, ϑ, θ_r)
    (; ψb, c) = cm
    if S < 1.0
        return -ψb / (c * (ν - θ_r)) * S^(-(1 + 1 / c))
    else
        return 1 / S_s
    end
end

"""
     hydraulic_conductivity(cm::BrooksCorey{FT}, K_sat::FT, S::FT) where {FT}

A point-wise function returning the hydraulic conductivity, using the
Brooks and Corey formulation.
"""
function hydraulic_conductivity(
    cm::BrooksCorey{FT},
    K_sat::FT,
    S::FT,
) where {FT}
    (; c) = cm
    if S < FT(1)
        K = S^(2 / c + 3)
    else
        K = FT(1)
    end
    return K * K_sat
end



"""
    pressure_head(
        cm::BrooksCorey{FT},
        θ_r::FT,
        ϑ_l::FT,
        ν_eff::FT,
        S_s::FT,
    ) where {FT}

A point-wise function returning the pressure head in
variably saturated soil, using the Brooks and Corey matric potential 
if the soil is not saturated, and
an approximation of the positive pressure in the soil if the
soil is saturated.
"""
function pressure_head(
    cm::BrooksCorey{FT},
    θ_r::FT,
    ϑ_l::FT,
    ν_eff::FT,
    S_s::FT,
) where {FT}
    S_l_eff = effective_saturation(ν_eff, ϑ_l, θ_r)
    if S_l_eff <= FT(1.0)
        ψ = matric_potential(cm, S_l_eff)
    else
        # For Brooks and Corey, S = 1 does not correspond to ψ = 0
        ψ = (ϑ_l - ν_eff) / S_s + cm.ψb
    end
    return ψ
end


"""
    impedance_factor(
        f_i::FT,
        Ω::FT
    ) where {FT}
Returns the multiplicative factor reducing conductivity when 
a fraction of ice `f_i` is present.

Only for use with the `EnergyHydrology` model.
"""
function impedance_factor(f_i::FT, Ω::FT) where {FT}
    gamma = FT(10.0^(-Ω * f_i))
    return gamma
end


"""
    viscosity_factor(
        T::FT,
        γ::FT,
        γT_ref::FT,
    ) where {FT}

Returns the multiplicative factor which accounts for
the temperature dependence of the conductivity.

Only for use with the `EnergyHydrology` model.
"""
function viscosity_factor(T::FT, γ::FT, γT_ref::FT) where {FT}
    factor = FT(γ * (T - γT_ref))
    Theta = FT(exp(factor))
    return Theta
end

"""
    ice_fraction(θ_l::FT, θ_i::FT, ν::FT, θ_r::FT)::FT where {FT}

Computes and returns the ice fraction, which is the
fraction  of the vapor flux that is due to sublimation, and 
the fraction of the humidity in the air due to ice, as

f = S_i/(S_i+S_l)

This same fraction is used to estimate the specific humidity, i.e.
q = q_over_ice * f + q_over_water * (1-f).
"""
function ice_fraction(θ_l::FT, θ_i::FT, ν::FT, θ_r::FT)::FT where {FT}
    S_l = effective_saturation(ν, θ_l, θ_r)
    S_i = effective_saturation(ν, θ_i, θ_r)
    f = S_i / (S_i + S_l)
    return f
end

"""
    soil_tortuosity(θ_l::FT, θ_i::FT, ν::FT) where {FT}

Computes the tortuosity of water vapor in a porous medium,
as a function of porosity `ν` and the volumetric liquid water
and ice contents, `θ_l` and `θ_i`.

See Equation (1) of : Shokri, N., P. Lehmann, and
D. Or (2008), Effects of hydrophobic layers on evaporation from
porous media, Geophys. Res. Lett., 35, L19407, doi:10.1029/
2008GL035230.
"""
function soil_tortuosity(θ_l::FT, θ_i::FT, ν::FT) where {FT}
    safe_θ_a = max(ν - θ_l - θ_i, eps(FT))
    return safe_θ_a^(FT(2.5)) / ν
end

"""
    soil_resistance(θ_l::FT,
                    θ_i::FT,
                    hydrology_cm::C,
                    ν::FT,
                    θ_r::FT,
                    d_ds::FT,
                    earth_param_set::EP,
                   ) where {FT, EP, C}

Computes the resistance of the top of the soil column to
water vapor diffusion, as a function of the surface 
volumetric liquid water fraction `θ_l`, the augmented
liquid water fraction `ϑ_l`,  the volumetric ice water
fraction `θ_i`, and other soil parameters.
"""
function soil_resistance(
    θ_l::FT,
    θ_i::FT,
    hydrology_cm::C,
    ν::FT,
    θ_r::FT,
    d_ds::FT,
    earth_param_set::EP,
) where {FT, EP, C}
    (; S_c) = hydrology_cm
    _D_vapor = FT(LP.D_vapor(earth_param_set))
    S_w = effective_saturation(ν, θ_l + θ_i, θ_r)
    τ_a = soil_tortuosity(θ_l, θ_i, ν)
    dsl::FT = dry_soil_layer_thickness(S_w, S_c, d_ds)
    r_soil = dsl / (_D_vapor * τ_a) # [s\m]
    return r_soil
end

"""
    dry_soil_layer_thickness(S_w::FT, S_c::FT, d_ds::FT)::FT where {FT}

Returns the maximum dry soil layer thickness that can develop under vapor flux; 
this is used when computing the soil resistance to vapor flux according to
Swenson et al (2012)/Sakaguchi and Zeng (2009).
"""
function dry_soil_layer_thickness(S_w::FT, S_c::FT, d_ds::FT)::FT where {FT}
    return S_w < S_c ? d_ds * (S_c - S_w) / S_c : FT(0)
end

"""
    is_saturated(twc::FT, ν::FT) where {FT}

A helper function which can be used to indicate whether a layer of soil is 
saturated based on if the total volumetric water content, `twc` is greater
than porosity `ν`.
"""
function is_saturated(twc::FT, ν::FT) where {FT}
    return ClimaLand.heaviside(twc - ν)
end
