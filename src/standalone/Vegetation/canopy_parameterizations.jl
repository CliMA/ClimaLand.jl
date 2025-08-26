using ..ClimaLand.Canopy
export canopy_sw_rt_beer_lambert, # Radiative transfer
    canopy_sw_rt_two_stream,
    extinction_coeff,
    compute_G,
    moisture_stress,
    # Conductance
    medlyn_term,
    medlyn_conductance,
    penman_monteith,
    quadratic_soil_moisture_stress

# 1. Radiative transfer

"""
    compute_G(
        G::ConstantGFunction,
        _,
    )

Returns the constant leaf angle distribution value for the given G function.
Takes in an arbitrary value for the cosine of the solar zenith angle, which is not used.
"""
function compute_G(G::ConstantGFunction, _)
    return G.ld
end

"""
    compute_G(
        G::CLMGFunction,
        cosθs,
    )

Returns the leaf angle distribution value for CLM G function as a function of the
cosine of the solar zenith angle and the leaf orientation index.
See section 3.1 of https://www2.cesm.ucar.edu/models/cesm2/land/CLM50_Tech_Note.pdf.

Note that the zenith angle is defined ∈ [0,2π), so to prevent a negative value
of G when the sun is below the horizon, we clip cosθs >= 0.
"""
function compute_G(G::CLMGFunction, cosθs::FT) where {FT}
    χl = G.χl
    ϕ1 = 0.5 - 0.633 * χl - 0.33 * χl^2
    ϕ2 = 0.877 * (1 - 2 * ϕ1)
    return FT(ϕ1 + ϕ2 * max(cosθs, 0))
end

"""
    compute_fractional_absorbances!(
        p,
        RT::BeerLambertModel{FT},
        LAI,
        α_soil_PAR,
        α_soil_NIR,
    )

Computes the PAR and NIR fractional absorbances, reflectances, and tranmittances
for a canopy in the case of the
Beer-Lambert model. The absorbances are a function of the radiative transfer
model, as well as the leaf area index, the clumping index,
the cosine of the zenith angle, the leaf angle distribution,
the extinction coefficient, and the
soil albedo in the PAR and NIR bands. Returns a
NamedTuple of NamedTuple, of the form:
(; par = (; refl = , trans = , abs = ),  nir = (; refl = , trans = , abs = ))
"""
function compute_fractional_absorbances!(
    p,
    RT::BeerLambertModel{FT},
    LAI,
    α_soil_PAR,
    α_soil_NIR,
) where {FT}
    RTP = RT.parameters
    cosθs = p.drivers.cosθs
    @. p.canopy.radiative_transfer.par = canopy_sw_rt_beer_lambert(
        RTP.G_Function,
        cosθs,
        RTP.Ω,
        RTP.α_PAR_leaf,
        LAI,
        α_soil_PAR,
    )
    @. p.canopy.radiative_transfer.nir = canopy_sw_rt_beer_lambert(
        RTP.G_Function,
        cosθs,
        RTP.Ω,
        RTP.α_NIR_leaf,
        LAI,
        α_soil_NIR,
    )
end

"""
    compute_fractional_absorbances!(p,
        RT::TwoStreamModel{FT},
        LAI,
        α_soil_PAR,
        α_soil_NIR,
    )

Computes the PAR and NIR fractional absorbances, reflectances, and tranmittances
for a canopy in the case of the
Two-stream model. The absorbances are a function of the radiative transfer
model, as well as the leaf area index, the clumping index,
the cosine of the zenith angle, the leaf angle distribution,
the extinction coefficient, and the
soil albedo in the PAR and NIR bands.

This model also depends on the diffuse fraction.
Returns a
NamedTuple of NamedTuple, of the form:
(; par = (; refl = , trans = , abs = ),  nir = (; refl = , trans = , abs = ))
"""
function compute_fractional_absorbances!(
    p,
    RT::TwoStreamModel{FT},
    LAI,
    α_soil_PAR,
    α_soil_NIR,
) where {FT}
    RTP = RT.parameters
    cosθs = p.drivers.cosθs
    frac_diff = p.drivers.frac_diff
    @. p.canopy.radiative_transfer.par = canopy_sw_rt_two_stream(
        RTP.G_Function,
        RTP.Ω,
        RTP.n_layers,
        RTP.α_PAR_leaf,
        RTP.τ_PAR_leaf,
        LAI,
        cosθs,
        α_soil_PAR,
        frac_diff,
    )
    @. p.canopy.radiative_transfer.nir = canopy_sw_rt_two_stream(
        RTP.G_Function,
        RTP.Ω,
        RTP.n_layers,
        RTP.α_NIR_leaf,
        RTP.τ_NIR_leaf,
        LAI,
        cosθs,
        α_soil_NIR,
        frac_diff,
    )
end

"""
    canopy_sw_rt_beer_lambert(
        G_Function,
        cosθs::FT,
        Ω::FT,
        α_leaf::FT,
        LAI::FT,
        α_soil::FT,
    )

Computes the absorbed, reflected, and transmitted flux fractions by radiation band.

This applies the Beer-Lambert law, which is a function of leaf reflectance
(`α_leaf`), the leaf angle distribution and zenith angle (defined via `G_Function`, and `cosθs`), leaf area index (`LAI`),
and the albedo of the soil (`α_soil`).

Returns a tuple of reflected, absorbed, and transmitted radiation fractions.
"""
function canopy_sw_rt_beer_lambert(
    G_Function,
    cosθs::FT,
    Ω::FT,
    α_leaf::FT,
    LAI::FT,
    α_soil::FT,
) where {FT}
    K = extinction_coeff(G_Function, cosθs)
    AR = (1 - α_leaf) * (1 - exp(-K * LAI * Ω)) * (1 - α_soil)
    TR = exp(-K * LAI * Ω)
    RR = FT(1) - AR - TR * (1 - α_soil)
    return (; abs = AR, refl = RR, trans = TR)
end

"""
    canopy_sw_rt_two_stream(
        G_Function,
        Ω::FT,
        n_layers::UInt64,
        SW_d::FT,
        α_leaf::FT,
        τ_leaf::FT,
        LAI::FT,
        cosθs::FT,
        α_soil::FT,
        frac_diff::FT,
    )

Computes the absorbed, reflected, and transmitted flux fractions by radiation band.

This applies the two-stream radiative transfer solution which takes into account
the impacts of scattering within the canopy. The function takes in all
parameters from the parameter struct of a TwoStreamModel, along with the
incident radiation, LAI, extinction coefficient K, soil albedo from the
canopy soil_driver, the cosine of the solar zenith angle, and τ.

Returns a tuple of reflected, absorbed, and transmitted radiation fractions.
"""
function canopy_sw_rt_two_stream(
    G_Function,
    Ω::FT,
    n_layers::UInt64,
    α_leaf::FT,
    τ_leaf::FT,
    LAI::FT,
    cosθs::FT,
    α_soil::FT,
    frac_diff::FT,
) where {FT}
    α_soil = max(eps(FT), α_soil) # this prevents division by zero, below.
    cosθs = max(eps(FT), cosθs) # The insolations package returns θs > π/2 (nighttime), but this assumes cosθs >0
    G = compute_G(G_Function, cosθs)
    K = extinction_coeff(G_Function, cosθs)
    # Compute μ̄, the average inverse diffuse optical length per LAI
    μ̄ = 1 / max(2G, eps(FT))

    # Clip this to prevent dividing by zero; the sum must also be < 1. Note that using eps(FT)
    # as the threshold leads to numerical errors, so we use 1e-4
    ω = min(max(α_leaf + τ_leaf, FT(1e-4)), 1 - FT(1e-4))

    # Compute aₛ, the single scattering albedo
    aₛ =
        0.5 * ω * (1 - cosθs * log((abs(cosθs) + 1) / max(abs(cosθs), eps(FT))))

    # Compute β₀, the direct upscattering parameter
    β₀ = (1 / ω) * aₛ * (1 + μ̄ * K) / (μ̄ * K)

    # Compute β, the diffuse upscattering parameter
    diff = α_leaf - τ_leaf
    # With uniform distribution, Dickinson integral becomes following:
    c²θ̄ = pi * G / 4
    β = 0.5 * (ω + diff * c²θ̄) / ω

    # Compute coefficients for two-stream solution
    b = 1 - ω + ω * β
    c = ω * β
    d = ω * β₀ * μ̄ * K
    f = ω * μ̄ * K * (1 - β₀)
    h = √(b^2 - c^2) / μ̄
    σ = (μ̄ * K)^2 + c^2 - b^2

    u₁ = b - c / α_soil
    u₂ = b - c * α_soil
    u₃ = f + c * α_soil

    s₁ = exp(-h * LAI * Ω)
    s₂ = exp(-K * LAI * Ω)

    p₁ = b + μ̄ * h
    p₂ = b - μ̄ * h
    p₃ = b + μ̄ * K
    p₄ = b - μ̄ * K

    d₁ = p₁ * (u₁ - μ̄ * h) / s₁ - p₂ * (u₁ + μ̄ * h) * s₁
    d₂ = (u₂ + μ̄ * h) / s₁ - (u₂ - μ̄ * h) * s₁

    # h coefficients for direct upward flux
    h₁ = -d * p₄ - c * f
    h₂ =
        1 / d₁ * (
            (d - h₁ / σ * p₃) * (u₁ - μ̄ * h) / s₁ -
            p₂ * s₂ * (d - c - h₁ / σ * (u₁ + μ̄ * K))
        )
    h₃ =
        -1 / d₁ * (
            (d - h₁ / σ * p₃) * (u₁ + μ̄ * h) * s₁ -
            p₁ * s₂ * (d - c - h₁ / σ * (u₁ + μ̄ * K))
        )

    # h coefficients for direct downward flux
    h₄ = -f * p₃ - c * d
    h₅ =
        -1 / d₂ *
        (h₄ * (u₂ + μ̄ * h) / (σ * s₁) + (u₃ - h₄ / σ * (u₂ - μ̄ * K)) * s₂)
    h₆ =
        1 / d₂ *
        (h₄ / σ * (u₂ - μ̄ * h) * s₁ + (u₃ - h₄ / σ * (u₂ - μ̄ * K)) * s₂)

    # h coefficients for diffuse upward flux
    h₇ = c * (u₁ - μ̄ * h) / (d₁ * s₁)
    h₈ = -c * s₁ * (u₁ + μ̄ * h) / d₁

    # h coefficients for diffuse downward flux
    h₉ = (u₂ + μ̄ * h) / (d₂ * s₁)
    h₁₀ = -s₁ * (u₂ - μ̄ * h) / d₂

    # Compute the LAI per layer for this canopy
    Lₗ = LAI / n_layers

    # Initialize the fraction absorbed value and layer counter
    F_abs = 0
    i = 0

    # Total light reflected from top of canopy
    F_refl = 0

    # Total light transmitted through the canopy on the downward pass
    F_trans = 0

    # Intialize vars to save computed fluxes from each layer for the next layer
    I_dir_up_prev = 0
    I_dir_dn_prev = 0
    I_dif_up_prev = 0
    I_dif_dn_prev = 0


    # Compute F_abs in each canopy layer
    while i <= n_layers

        # Compute cumulative LAI at this layer
        L = i * Lₗ

        # Compute the direct fluxes into/out of the layer
        I_dir_up =
            h₁ * exp(-K * L * Ω) / σ +
            h₂ * exp(-h * L * Ω) +
            h₃ * exp(h * L * Ω)
        I_dir_dn =
            h₄ * exp(-K * L * Ω) / σ +
            h₅ * exp(-h * L * Ω) +
            h₆ * exp(h * L * Ω)

        # Add collimated radiation to downward flux
        I_dir_dn += exp(-K * L * Ω)

        # Compute the diffuse fluxes into/out of the layer
        I_dif_up = h₇ * exp(-h * L * Ω) + h₈ * exp(h * L * Ω)
        I_dif_dn = h₉ * exp(-h * L * Ω) + h₁₀ * exp(h * L * Ω)

        # Energy balance giving radiation absorbed in the layer
        if i == 0
            I_dir_abs = 0
            I_dif_abs = 0
        else
            I_dir_abs = I_dir_up - I_dir_up_prev - I_dir_dn + I_dir_dn_prev
            I_dif_abs = I_dif_up - I_dif_up_prev - I_dif_dn + I_dif_dn_prev
        end

        if i == 1 # prev = top of layer
            F_refl =
                (1 - frac_diff) * I_dir_up_prev + (frac_diff) * I_dif_up_prev
        end
        if i == n_layers # not prev = bottom of layer
            F_trans = (1 - frac_diff) * I_dir_dn + (frac_diff) * I_dif_dn
        end


        # Add radiation absorbed in the layer to total absorbed radiation
        F_abs += (1 - frac_diff) * I_dir_abs + (frac_diff) * I_dif_abs

        # Save input/output values to compute energy balance of next layer
        I_dir_up_prev = I_dir_up
        I_dir_dn_prev = I_dir_dn
        I_dif_up_prev = I_dif_up
        I_dif_dn_prev = I_dif_dn

        # Move on to the next layer
        i += 1
    end

    # Reflected is reflected from the land surface. F_abs
    # refers to radiation absorbed by the canopy on either pass.
    # F_trans refers to the downward pass only. Note that
    # (1-α_soil)*F_trans + F_abs = total absorbed fraction by land, so the following
    # must hold
    # @assert (1 - α_soil) * FT(F_trans) + FT(F_abs) + FT(F_refl) ≈ 1
    # This is tested in test/standalone/Vegetation/test_two_stream.jl
    return (; abs = FT(F_abs), refl = FT(F_refl), trans = FT(F_trans))
end

"""
    extinction_coeff(G_Function,
                     cosθs::FT) where {FT}

Computes the vegetation extinction coefficient (`K`), as a function
of the cosine of the sun zenith angle (`cosθs`),
and the leaf angle distribution function(`G_Function`).

In the two-stream scheme, values of K ~ 1/epsilon can lead to numerical issues.
Here we clip it to 1e6.
"""
function extinction_coeff(G_Function, cosθs::FT) where {FT}
    G = compute_G(G_Function, cosθs)
    K = min(G / max(cosθs, eps(FT)), FT(1e6))
    return K
end

"""
    moisture_stress(pl::FT,
                    sc::FT,
                    pc::FT) where {FT}

Computes the moisture stress factor (`β`), which is unitless,
 as a function of
a constant (`sc`, 1/Pa), a reference pressure (`pc`, Pa), and
the leaf water pressure (`pl`, Pa) .

See Eqn 12.57 of G. Bonan's textbook,
Climate Change and Terrestrial Ecosystem Modeling (2019).
"""
function moisture_stress(pl::FT, sc::FT, pc::FT) where {FT}
    β = min(FT(1), (1 + exp(sc * pc)) / (1 + exp(sc * (pc - pl))))
    return β
end

"""
    conductance_molar_flux_to_m_per_s(gs::FT, T::FT, R::FT, P::FT) where {FT}

This currently takes a conductance (moles per leaf area per second)
and converts it to m/s.
"""
function conductance_molar_flux_to_m_per_s(
    gs::FT,
    T::FT,
    R::FT,
    P::FT,
) where {FT}
    return gs * (R * T) / P # convert to m s-1
end


# 3. Stomatal conductance model
"""
    medlyn_term(g1::FT, T_air::FT, P_air::FT, q_air::FT, thermo_params) where {FT}

Computes the Medlyn term, equal to `1+g1/sqrt(VPD)`,
by first computing the `VPD`,
where `VPD` is the vapor pressure deficit in the atmosphere
(Pa), and `g_1` is a constant with units of `sqrt(Pa)`.

`thermo_params` is the Thermodynamics.jl parameter set.

Note that in supersaturated conditions the vapor pressure deficit will be negative,
which leads to an imaginary Medlyn term `m`. Clipping to zero solves this, but this leads
to division by zero, so we regularize the division by adding a small quantity.

An alternative to consider in the future is to compute the inverse of this quantity
and stomatal resistance instead of conductance.
"""
function medlyn_term(
    g1::FT,
    T_air::FT,
    P_air::FT,
    q_air::FT,
    thermo_params,
) where {FT}
    VPD = max(
        ClimaLand.vapor_pressure_deficit(T_air, P_air, q_air, thermo_params),
        FT(0),
    ) # clip negative values of VPD to zero
    return 1 + g1 / sqrt(VPD + sqrt(eps(FT))) # regularize values of VPD which are smaller than sqrt(eps(FT))
end


"""
    medlyn_conductance(g0::FT,
                       Drel::FT,
                       medlyn_term::FT,
                       An::FT,
                       ca::FT) where {FT}

Computes the stomatal conductance according to Medlyn, as a function of
the minimum stomatal conductance (`g0`),
the relative diffusivity of water vapor with respect to CO2 (`Drel`),
the Medlyn term (unitless), the biochemical demand for CO2 (`An`), and the
atmospheric concentration of CO2 (`ca`).

This returns the conductance in units of mol/m^2/s. It must be converted to
m/s using the molar density of water prior to use in SurfaceFluxes.jl.
"""
function medlyn_conductance(
    g0::FT,
    Drel::FT,
    medlyn_term::FT,
    An::FT,
    ca::FT,
) where {FT}
    gs = g0 + Drel * medlyn_term * (An / ca)
    return gs
end

"""
    penman_monteith(
        Δ::FT, # Rate of change of saturation vapor pressure with air temperature. (Pa K−1)
        Rn::FT, # Net irradiance (W m−2)
        G::FT, # Ground heat flux (W m−2)
        ρa::FT, # Dry air density (kg m−3)
        cp::FT, # Specific heat capacity of air (J kg−1 K−1)
        VPD::FT, # vapor pressure deficit (Pa)
        ga::FT, # atmospheric conductance (m s−1)
        γ::FT, # Psychrometric constant (γ ≈ 66 Pa K−1)
        gs::FT, # surface or stomatal conductance (m s−1)
        Lv::FT, # Volumetric latent heat of vaporization (J m-3)
        ) where {FT}

Computes the evapotranspiration in m/s using the Penman-Monteith equation.
"""
function penman_monteith(
    Δ::FT,
    Rn::FT,
    G::FT,
    ρa::FT,
    cp::FT,
    VPD::FT,
    ga::FT,
    γ::FT,
    gs::FT,
    Lv::FT,
) where {FT}
    ET = (Δ * (Rn - G) + ρa * cp * VPD * ga) / ((Δ + γ * (1 + ga / gs)) * Lv)
    return ET
end

"""
    quadratic_soil_moisture_stress(
        θ::FT,
        meanalpha::FT = FT(1.0),
        a_hat::FT = FT(0.0),
        b_hat::FT = FT(0.685),
        θ0::FT = FT(0.0),
        θ1::FT = FT(0.6)
    ) where {FT}

Computes an empirical soil moisture stress factor (`β`) according to the quadratic
functional form of Stocker et al. (2020), Eq 21, using the soil moisture (`θ`),
AET/PET (`meanalpha`), and parameters `a_hat`, `b_hat`, `θ0`, and `θ1`. Default values
are from the original paper.
"""
function quadratic_soil_moisture_stress(
    θ::FT,
    meanalpha::FT = FT(1.0),
    a_hat::FT = FT(0.0),
    b_hat::FT = FT(0.733),
    θ0::FT = FT(0.0),
    θ1::FT = FT(0.6),
) where {FT}
    β0 = a_hat + b_hat * meanalpha
    q = (FT(1.0) - β0) / (θ0 - θ1)^2
    β = FT(1.0) - q * (θ - θ1)^2

    # Bound to [0,1]
    β = clamp(β, FT(0.0), FT(1.0))

    # Force 1.0 above the soil-moisture threshold θ1
    β = ifelse(θ > θ1, FT(1.0), β)

    return β
end
