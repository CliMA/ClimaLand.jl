export volumetric_liquid_fraction,
    effective_saturation,
    matric_potential,
    inverse_matric_potential,
    pressure_head,
    hydraulic_conductivity,
    impedance_factor,
    viscosity_factor,
    dψdϑ,
    dry_soil_layer_thickness
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
    S_l_eff = effective_saturation(ν_eff, ϑ_l, θ_r)
    if S_l_eff <= FT(1.0)
        ψ = matric_potential(cm, S_l_eff)
    else
        ψ = (ϑ_l - ν_eff) / S_s
    end
    return ψ
end

"""
   dψdϑ(cm::vanGenuchten{FT}, ϑ, ν, θ_r, S_s)

Computes and returns the derivative of the pressure head 
with respect to ϑ for the van Genuchten formulation.
"""
function dψdϑ(cm::vanGenuchten{FT}, ϑ, ν, θ_r, S_s) where {FT}
    S = effective_saturation(ν, ϑ, θ_r)
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
    dry_soil_layer_thickness(S_l_sfc::FT, S_c::FT, d_ds::FT) where {FT}

Returns the maximum dry soil layer thickness that can develop under evaporation; 
this is used when computing the soil resistance to evaporation according to
Swenson et al (2012).
"""
function dry_soil_layer_thickness(S_l_sfc::FT, S_c::FT, d_ds::FT) where {FT}
    return S_l_sfc < S_c ? d_ds * (S_c - S_l_sfc) / S_c : FT(0)
end
