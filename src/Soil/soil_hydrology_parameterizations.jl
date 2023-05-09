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
   dψdϑ(ϑ, ν, θ_r, vg_α, vg_n, vg_m, S_s)

Computes and returns the derivative of the pressure head 
with respect to ϑ for the van Genuchten formulation.
"""
function dψdϑ(ϑ, ν, θ_r, vg_α, vg_n, vg_m, S_s)
    S = effective_saturation(ν, ϑ, θ_r)
    if S < 1.0
        return 1.0 / (vg_α * vg_m * vg_n) / (ν - θ_r) *
               (S^(-1 / vg_m) - 1)^(1 / vg_n - 1) *
               S^(-1 / vg_m - 1)
    else
        return 1.0 / S_s
    end
end

"""
    volumetric_liquid_fraction(ϑ_l::FT, ν_eff::FT) where {FT}

A pointwise function returning the volumetric liquid fraction
given the augmented liquid fraction and the effective porosity.
"""
function volumetric_liquid_fraction(ϑ_l::FT, ν_eff::FT) where {FT}
    if ϑ_l < ν_eff
        θ_l = ϑ_l
    else
        θ_l = ν_eff
    end
    return θ_l
end


"""
     matric_potential(α::FT, n::FT, m::FT, S::FT) where {FT}

A point-wise function returning the matric potential, using the
van Genuchten formulation.
"""
function matric_potential(α::FT, n::FT, m::FT, S::FT) where {FT}
    ψ = -((S^(-FT(1) / m) - FT(1)) * α^(-n))^(FT(1) / n)
    return ψ
end

"""
     inverse_matric_potential(α::FT, n::FT, m::FT, ψ::FT) where {FT}

A point-wise function returning the effective saturation, given
the matric potential, using the
van Genuchten formulation.
"""
function inverse_matric_potential(α::FT, n::FT, m::FT, ψ::FT) where {FT}
    ψ > 0 && error("Matric potential is positive")
    S = (FT(1) + (α * abs(ψ))^n)^(-m)
    return S
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
    pressure_head(
        α::FT,
        n::FT,
        m::FT,
        θ_r::FT,
        ϑ_l::FT,
        ν_eff::FT,
        S_s::FT,
    ) where {FT}

A point-wise function returning the pressure head in
variably saturated soil, using the van Genuchten formulation
for matric potential if the soil is not saturated, and
an approximation of the positive pressure in the soil if the
soil is saturated.
"""
function pressure_head(
    α::FT,
    n::FT,
    m::FT,
    θ_r::FT,
    ϑ_l::FT,
    ν_eff::FT,
    S_s::FT,
) where {FT}
    S_l_eff = effective_saturation(ν_eff, ϑ_l, θ_r)
    if S_l_eff <= FT(1.0)
        ψ = matric_potential(α, n, m, S_l_eff)
    else
        ψ = (ϑ_l - ν_eff) / S_s
    end
    return ψ
end

"""
     hydraulic_conductivity(K_sat::FT, m::FT, S::FT) where {FT}

A point-wise function returning the hydraulic conductivity, using the
van Genuchten formulation.
"""
function hydraulic_conductivity(K_sat::FT, m::FT, S::FT) where {FT}
    if S < FT(1)
        K = sqrt(S) * (FT(1) - (FT(1) - S^(FT(1) / m))^m)^FT(2)
    else
        K = FT(1)
    end
    return K * K_sat
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
