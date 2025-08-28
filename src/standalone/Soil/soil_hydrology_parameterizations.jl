export volumetric_liquid_fraction,
    effective_saturation,
    matric_potential,
    inverse_matric_potential,
    pressure_head,
    hydraulic_conductivity,
    impedance_factor,
    viscosity_factor,
    dψdϑ

"""
    volumetric_liquid_fraction(ϑ_l::FT, ν_eff::FT, θ_r::FT) where {FT}

A pointwise function returning the volumetric liquid fraction
given the augmented liquid fraction and the effective porosity.
The output is guaranteed to be in (θ_r, ν_eff].

For Richards model, ν_eff = ν; and the clipping below is not required,
(and will do nothing; ν_eff_safe = ν_eff), but we leave it in for a 
simpler interface.
"""
function volumetric_liquid_fraction(ϑ_l::FT, ν_eff::FT, θ_r::FT) where {FT}
    ϑ_l_safe = max(ϑ_l, θ_r + sqrt(eps(FT)))
    ν_eff_safe = max(ν_eff, θ_r + sqrt(eps(FT)))

    if ϑ_l_safe < ν_eff_safe
        θ_l = ϑ_l_safe
    else
        θ_l = ν_eff_safe
    end
    return θ_l
end

"""
    effective_saturation(ν_eff::FT, ϑ_l::FT, θr::FT) where {FT}

A point-wise function computing the effective saturation given the
effective porosity, augmented liquid fraction, and residual water
fraction as input.

For Richards model, or any other parameterization where ice is not
relevant, ν_eff = ν; and the clipping below is not required,
(and will do nothing; ν_eff_safe = ν_eff), but we leave it in for a 
simpler interface.
"""
function effective_saturation(ν_eff::FT, ϑ_l::FT, θr::FT) where {FT}
    ϑ_l_safe = max(ϑ_l, θr + sqrt(eps(FT)))
    ν_eff_safe = max(ν_eff, θr + sqrt(eps(FT)))
    S_l = (ϑ_l_safe - θr) / (ν_eff_safe - θr)
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
    # effective saturation clips ν_eff and ϑ_l
    # as needed so that S_l ∈ (0,1].
    S_l_eff = effective_saturation(ν_eff, ϑ_l, θ_r)
    ϑ_l_safe = max(ϑ_l, θ_r + sqrt(eps(FT)))
    ν_eff_safe = max(ν_eff, θ_r + sqrt(eps(FT)))

    if S_l_eff <= FT(1.0)
        ψ = matric_potential(cm, S_l_eff)
    else
        ψ = (ϑ_l_safe - ν_eff_safe) / S_s
    end
    return ψ
end

"""
   dψdϑ(cm::vanGenuchten{FT}, ϑ, ν_eff, θ_r, S_s)

Computes and returns the derivative of the pressure head
with respect to ϑ for the van Genuchten formulation.
"""
function dψdϑ(cm::vanGenuchten{FT}, ϑ, ν_eff, θ_r, S_s) where {FT}
    # effective saturation clips ν_eff and ϑ
    # as needed so that S_l ∈ (0,1],
    # but we use ν_eff alone, so we clip that here.
    # The second clipping in effective saturation
    # will not change it further.
    ν_eff_safe = max(ν_eff, θ_r + sqrt(eps(FT)))
    S = effective_saturation(ν_eff_safe, ϑ, θ_r)
    (; α, m, n) = cm
    if S < 1.0
        return FT(
            1.0 / (α * m * n) / (ν_eff_safe - θ_r) *
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
   dψdϑ(cm::BrooksCorey{FT}, ϑ, ν_eff, θ_r, S_s)

Computes and returns the derivative of the pressure head
with respect to ϑ for the Brooks and Corey formulation.
"""
function dψdϑ(cm::BrooksCorey{FT}, ϑ, ν_eff, θ_r, S_s) where {FT}
    # effective saturation clips ν_eff and ϑ
    # as needed so that S_l ∈ (0,1],
    # but we use ν_eff alone, so we clip that here.
    # The second clipping in effective saturation
    # will not change it further.
    ν_eff_safe = max(ν_eff, θ_r + sqrt(eps(FT)))
    S = effective_saturation(ν_eff_safe, ϑ, θ_r)
    (; ψb, c) = cm
    if S < 1.0
        return -ψb / (c * (ν_eff_safe - θ_r)) * S^(-(1 + 1 / c))
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
    # effective saturation clips ν_eff and ϑ_l
    # as needed so that S_l ∈ (0,1].
    S_l_eff = effective_saturation(ν_eff, ϑ_l, θ_r)
    ϑ_l_safe = max(ϑ_l, θ_r + sqrt(eps(FT)))
    ν_eff_safe = max(ν_eff, θ_r + sqrt(eps(FT)))
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
