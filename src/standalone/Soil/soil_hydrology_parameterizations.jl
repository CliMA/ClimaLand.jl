export volumetric_liquid_fraction,
    effective_saturation,
    matric_potential,
    inverse_matric_potential,
    pressure_head,
    hydraulic_conductivity,
    impedance_factor,
    viscosity_factor,
    dψdϑ,
    update_albedo!,
    ice_fraction,
    soil_tortuosity,
    soil_resistance,
    dry_soil_layer_thickness

"""
    albedo_from_moisture(surface_eff_sat::FT, albedo_dry::FT, albedo_wet::FT)

Calculates pointwise albedo for any band as a function of soil effective saturation given
the dry and wet albedo values for that band.
"""
function albedo_from_moisture(
    surface_eff_sat::FT,
    albedo_dry::FT,
    albedo_wet::FT,
) where {FT}
    return albedo_dry * (1 - surface_eff_sat) + albedo_wet * surface_eff_sat
end


"""
    update_albedo!(bc::AtmosDrivenFluxBC, p, soil_domain, model_parameters)

Calculates and updates PAR and NIR albedo as a function of volumetric soil water content at
the top of the soil. If the soil layers are larger than the specified `albedo_calc_top_thickness`,
the water content of the top layer is used in the calclulation. For the PAR and NIR bands,

α_band = α_{band,dry} * (1 - S_e) +  α_{band,wet} * (S_e)

where S_e is the relative soil wetness above some depth, `albedo_calc_top_thickness`. This
is a modified version of Equation (1) of:

Braghiere, R. K., Wang, Y., Gagné-Landmann, A., Brodrick, P. G., Bloom, A. A., Norton,
A. J., et al. (2023). The importance of hyperspectral soil albedo information for improving
Earth system model projections. AGU Advances, 4, e2023AV000910. https://doi.org/10.1029/2023AV000910

where effective saturation is used in place of volumetric soil water content.The dry and wet
albedo values come from a global soil color map and soil color to albedo map from CLM.

CLM reference: Lawrence, P.J., and Chase, T.N. 2007. Representing a MODIS consistent land surface in the Community Land Model
(CLM 3.0). J. Geophys. Res. 112:G01023. DOI:10.1029/2006JG000168.
"""
function update_albedo!(bc::AtmosDrivenFluxBC, p, soil_domain, model_parameters)
    (;
        ν,
        θ_r,
        PAR_albedo_dry,
        NIR_albedo_dry,
        PAR_albedo_wet,
        NIR_albedo_wet,
        albedo_calc_top_thickness,
    ) = model_parameters
    S_sfc = p.soil.sfc_S_e
    FT = eltype(soil_domain.fields.Δz_top)
    # checks if there is at least 1 layer centered within the top soil depth
    if minimum(soil_domain.fields.Δz_top) < albedo_calc_top_thickness
        # ∫H_S_e_dz is the integral of effective saturation from (surface-albedo_calc_top_thickness) to surface
        ∫H_S_e_dz = p.soil.sfc_S_e
        # ∫H_dz is integral of 1 from (surface-albedo_calc_top_thickness) to surface
        ∫H_dz = p.soil.sfc_scratch
        # zero all centers lower than boundary, set everything above to one
        @. p.soil.sub_sfc_scratch = ClimaLand.heaviside(
            albedo_calc_top_thickness + sqrt(eps(FT)),
            soil_domain.fields.z_sfc - soil_domain.fields.z,
        )
        ClimaCore.Operators.column_integral_definite!(
            ∫H_dz,
            p.soil.sub_sfc_scratch,
        )
        # zeros all effective saturation at levels centered lower than boundary
        @. p.soil.sub_sfc_scratch =
            ClimaLand.heaviside(
                albedo_calc_top_thickness + sqrt(eps(FT)),
                soil_domain.fields.z_sfc - soil_domain.fields.z,
            ) * effective_saturation(ν, p.soil.θ_l, θ_r)
        ClimaCore.Operators.column_integral_definite!(
            ∫H_S_e_dz,
            p.soil.sub_sfc_scratch,
        )
        @. S_sfc = ∫H_S_e_dz / ∫H_dz
    else
        # in the case where no layer is centered above boundary, use the values of the top layer
        S_sfc .= ClimaLand.Domains.top_center_to_surface(
            effective_saturation.(ν, p.soil.θ_l, θ_r),
        )
    end
    @. p.soil.PAR_albedo =
        albedo_from_moisture(S_sfc, PAR_albedo_dry, PAR_albedo_wet)
    @. p.soil.NIR_albedo =
        albedo_from_moisture(S_sfc, NIR_albedo_dry, NIR_albedo_wet)
end

"""
update_albedo!(bc::AbstractWaterBC, p, soil_domain, model_parameters)

    Does nothing for boundary conditions where albedo is not used
"""
function update_albedo!(
    bc::BC,
    p,
    soil_domain,
    model_parameters,
) where {BC <: AbstractEnergyHydrologyBC} end
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
    ice_fraction(θ_l::FT, θ_i::FT, ν::FT, θ_r::FT, ρ_l::FT, ρ_i::FT)::FT 
where {FT}

Computes and returns the
fraction of the humidity in the pore space 
that is due to sublimation.

f = θ_iρ_i/(θ_iρ_i+(θ_l-θ_r)ρ_l)
"""
function ice_fraction(
    θ_l::FT,
    θ_i::FT,
    ν::FT,
    θ_r::FT,
    ρ_l::FT,
    ρ_i::FT,
)::FT where {FT}
    m_l = (θ_l - θ_r) * ρ_l
    m_i = θ_i * ρ_i
    f = m_i / (m_i + m_l)
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
                    ice_frac::FT,
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
    ice_frac::FT,
    hydrology_cm::C,
    ν::FT,
    θ_r::FT,
    d_ds::FT,
    earth_param_set::EP,
) where {FT, EP, C}
    (; S_c, α) = hydrology_cm
    _D_vapor = LP.D_vapor(earth_param_set)
    _ρ_liq = LP.ρ_cloud_liq(earth_param_set)
    _grav = LP.grav(earth_param_set)
    _σ = FT(7.2e-2) # need to add to CP, surface tension of water N/m

    S_w = effective_saturation(ν, θ_l + θ_i, θ_r)
    τ_a = soil_tortuosity(θ_l, θ_i, ν)
    dsl::FT = dry_soil_layer_thickness(S_w, S_c, d_ds)
    r_pore::FT = 2 * _σ * α / _ρ_liq / _grav
    θ_safe = max(eps(FT), (θ_i + θ_l - θ_r))
    r_shell::FT = r_pore / _D_vapor / (4 * θ_safe) * (π - 2 * (θ_safe)^(1 / 2))
    # This factor just damps r_shell to zero more quickly as x -> 0. numerical.
    x = θ_safe / FT(0.001)
    factor = 1 / (1 - exp(-x))^2
    r_soil = dsl / (_D_vapor * τ_a) + factor * r_shell# [s\m]
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
