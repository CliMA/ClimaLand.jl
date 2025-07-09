using ..ClimaLand.Canopy
export canopy_sw_rt_beer_lambert,
    canopy_sw_rt_two_stream,
    extinction_coeff,
    compute_G,
    arrhenius_function,
    intercellular_co2,
    co2_compensation,
    rubisco_assimilation,
    light_assimilation,
    max_electron_transport,
    electron_transport,
    optimality_max_photosynthetic_rates,
    net_photosynthesis,
    moisture_stress,
    dark_respiration,
    compute_GPP,
    MM_Kc,
    MM_Ko,
    compute_Vcmax,
    medlyn_term,
    medlyn_conductance,
    upscale_leaf_conductance,
    penman_monteith,
    nitrogen_content,
    plant_respiration_maintenance,
    plant_respiration_growth,
    enforce_albedo_constraint,
    intrinsic_quantum_yield,
    compute_viscosity_ratio,
    compute_Kmm,
    optimal_co2_ratio_c3,
    pmodel_vcmax,
    compute_LUE,
    compute_mj_with_jmax_limitation,
    co2_compensation_p,
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
See section 3.1 of https://www2.cesm.ucar.edu/models/cesm2/land/CLM50_Tech_Note.pdf
"""
function compute_G(G::CLMGFunction, cosθs::FT) where {FT}
    χl = G.χl
    ϕ1 = 0.5 - 0.633 * χl - 0.33 * χl^2
    ϕ2 = 0.877 * (1 - 2 * ϕ1)
    return FT(ϕ1 + ϕ2 * cosθs)
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

# 2. Photosynthesis, Farquhar model


"""
    intercellular_co2(ca::FT, Γstar::FT, medlyn_factor::FT) where{FT}

Computes the intercellular CO2 concentration (mol/mol) given the
atmospheric concentration (`ca`, mol/mol), the CO2 compensation (`Γstar`,
 mol/mol), and the Medlyn factor (unitless).
"""
function intercellular_co2(ca::FT, Γstar::FT, medlyn_term::FT) where {FT}
    c_i = max(ca * (1 - 1 / medlyn_term), Γstar)
    return c_i
end

"""
    co2_compensation(Γstar25::FT,
                     ΔHΓstar::FT,
                     T::FT,
                     To::FT,
                     R::FT) where {FT}

Computes the CO2 compensation point (`Γstar`),
in units of mol/mol,
as a function of its value at 25 °C (`Γstar25`),
a constant energy of activation (`ΔHΓstar`), a standard temperature (`To`),
the unversal gas constant (`R`), and the temperature (`T`).

See Table 11.5 of G. Bonan's textbook, Climate Change and Terrestrial Ecosystem Modeling (2019).
"""
function co2_compensation(
    Γstar25::FT,
    ΔHΓstar::FT,
    T::FT,
    To::FT,
    R::FT,
) where {FT}
    Γstar = Γstar25 * arrhenius_function(T, To, R, ΔHΓstar)
    return Γstar
end

"""
    rubisco_assimilation(is_c3::AbstractFloat, args...)

Calls the correct rubisco assimilation function based on the `is_c3`.

A `is_c3` value of 1.0 corresponds to C3 photosynthesis and calls
`c3_rubisco_assimilation`, while 0.0 corresponds to C4 photsynthesis and calls
`c4_rubisco_assimilation`.
"""
function rubisco_assimilation(is_c3::AbstractFloat, args...)
    is_c3 > 0.5 ? c3_rubisco_assimilation(args...) :
    c4_rubisco_assimilation(args...)
end

"""
    c3_rubisco_assimilation(Vcmax::FT,
                         ci::FT,
                         Γstar::FT,
                         Kc::FT,
                         Ko::FT,
                         oi::FT) where {FT}

Computes the Rubisco limiting rate of photosynthesis for C3 plants (`Ac`),
in units of moles CO2/m^2/s,
as a function of the maximum rate of carboxylation of Rubisco (`Vcmax`),
the leaf internal carbon dioxide partial pressure (`ci`),
the CO2 compensation point (`Γstar`), and Michaelis-Menten parameters
for CO2 and O2, respectively, (`Kc`) and (`Ko`).

The empirical parameter oi is equal to 0.209 (mol/mol).
See Table 11.5 of G. Bonan's textbook, Climate Change and Terrestrial Ecosystem Modeling (2019).
"""
function c3_rubisco_assimilation(
    Vcmax::FT,
    ci::FT,
    Γstar::FT,
    Kc::FT,
    Ko::FT,
    oi::FT,
) where {FT}
    Ac = Vcmax * (ci - Γstar) / (ci + Kc * (1 + oi / Ko))
    return Ac
end


function c3_rubisco_assimilation(
    Vcmax::FT,
    ci::FT,
    Γstar::FT,
    Kmm::FT
) where {FT}
    Ac = Vcmax * (ci - Γstar) / (ci + Kmm) 
    return Ac
end 

"""
    c4_rubisco_assimilation(Vcmax::FT,_...) where {FT}

Computes the Rubisco limiting rate of photosynthesis for C4 plants (`Ac`)
in units of moles CO2/m^2/s,
as equal to the maximum rate of carboxylation of Rubisco (`Vcmax`).
"""
function c4_rubisco_assimilation(Vcmax::FT, _...) where {FT}
    Ac = Vcmax
    return Ac
end

"""
    light_assimilation(is_c3::AbstractFloat, args...)

Calls the correct light assimilation function based on the `is_c3`.

A `is_c3` value of 1.0 corresponds to C3 photosynthesis and calls
`c3_light_assimilation`, while 0.0 corresponds to C4 photsynthesis and calls
`c4_light_assimilation`.
"""
function light_assimilation(is_c3::AbstractFloat, args...)
    is_c3 > 0.5 ? c3_light_assimilation(args...) :
    c4_light_assimilation(args...)
end
"""
    c3_light_assimilation(
                       J::FT,
                       ci::FT,
                       Γstar::FT,
                       ::FT,
                       ::FT) where {FT}

Computes the electron transport limiting rate (`Aj`),
in units of moles CO2/m^2/s.

For C3 plants, this is a function of
the rate of electron transport (`J`), the leaf internal carbon dioxide partial pressure (`ci`),
and the CO2 compensation point (`Γstar`).
See Table 11.5 of G. Bonan's textbook, Climate Change and Terrestrial Ecosystem Modeling (2019).
"""
function c3_light_assimilation(J::FT, ci::FT, Γstar::FT, _...) where {FT}
    Aj = J * (ci - Γstar) / (4 * (ci + 2 * Γstar))
    return Aj
end

"""
    light_assimilation(::FT, ::FT, ::FT, APAR::FT, E::FT) where {FT}

Computes the electron transport limiting rate (`Aj`),
in units of moles CO2/m^2/s.

For C4 plants, this is a function of APAR and a efficiency parameter E, see Equation 11.70 of
 G. Bonan's textbook, Climate Change and Terrestrial Ecosystem Modeling (2019).
"""
function c4_light_assimilation(::FT, ::FT, ::FT, APAR::FT, E::FT) where {FT}
    Aj = APAR * E
    return Aj
end

"""
    max_electron_transport(Vcmax::FT) where {FT}

Computes the maximum potential rate of electron transport (`Jmax`),
in units of mol/m^2/s,
as a function of Vcmax at 25 °C (`Vcmax25`),
a constant (`ΔHJmax`), a standard temperature (`To`),
the unversal gas constant (`R`), and the temperature (`T`).

See Table 11.5 of G. Bonan's textbook,
Climate Change and Terrestrial Ecosystem Modeling (2019).
"""
function max_electron_transport(
    Vcmax25::FT,
    ΔHJmax::FT,
    T::FT,
    To::FT,
    R::FT,
) where {FT}
    Jmax25 = Vcmax25 * FT(exp(1))
    Jmax = Jmax25 * arrhenius_function(T, To, R, ΔHJmax)
    return Jmax
end

"""
    electron_transport(APAR::FT,
                       Jmax::FT,
                       θj::FT,
                       ϕ::FT) where {FT}

Computes the rate of electron transport (`J`),
in units of mol/m^2/s, as a function of
the maximum potential rate of electron transport (`Jmax`),
absorbed photosynthetically active radiation (`APAR`),
an empirical "curvature parameter" (`θj`; Bonan Eqn 11.21)
and the quantum yield of photosystem II (`ϕ`).

See Ch 11, G. Bonan's textbook, Climate Change and Terrestrial Ecosystem Modeling (2019).
"""
function electron_transport(APAR::FT, Jmax::FT, θj::FT, ϕ::FT) where {FT}
    # Light utilization of APAR
    IPSII = ϕ * APAR / 2
    # This is a solution to a quadratic equation
    # θj *J^2 - (IPSII+Jmax)*J+IPSII*Jmax = 0, Equation 11.21
    J =
        (IPSII + Jmax - sqrt((IPSII + Jmax)^2 - 4 * θj * IPSII * Jmax)) /
        (2 * θj)
    return J
end

"""
optimality_max_photosynthetic_rates(APAR::FT,  θj::FT, ϕ::FT, oi::FT, ci::FT, Γstar::FT, Kc::FT, Ko::FT)

Computes the photosynthesis rates Vcmax and Jmax in
mol/m^2/s given absorbed photosynthetically active radiation (`APAR`),
an empirical "curvature parameter" (`θj`; Bonan Eqn 11.21)
 the quantum yield of photosystem II (`ϕ`), the intercellular
o2 content (`oi`), the intercellular CO2 concentration (ci),
Γstar, and Kc and Ko.

See Smith et al. 2019.
"""
function optimality_max_photosynthetic_rates(
    APAR::FT,
    θj::FT,
    ϕ::FT,
    oi::FT,
    ci::FT,
    Γstar::FT,
    Kc::FT,
    Ko::FT,
    c::FT,
) where {FT}
    # Light utilization of APAR
    IPSII = ϕ * APAR / 2

    mc = (ci - Γstar) / (ci + Kc * (1 + oi / Ko))
    m = (ci - Γstar) / (ci + 2 * Γstar)

    # Corrected form of ω, see https://github.com/SmithEcophysLab/optimal_vcmax_R/issues/3
    if ((θj < 1) & (8 * c > m) & (4 * c < m) & (m / c < 8 * θj)) |
       ((θj > 1) & (4 * c > m))
        ω =
            (
                -2 + 4 * θj - sqrt(
                    ((-1 + θj) * (m - 8 * c * θj)^2) / (c * (-m + 4 * c * θj)),
                )
            ) / 2
    elseif ((θj < 1) & (8 * c < m)) |
           ((m / c > 8 * θj) & (8 * c > m) & (4 * c < m))
        ω =
            (
                -2 +
                4 * θj +
                sqrt(((-1 + θj) * (m - 8 * c * θj)^2) / (c * (-m + 4 * c * θj)))
            ) / 2
    else
        ω = FT(0)
    end
    ωstar = 1 + ω - sqrt((1 + ω)^2 - 4 * θj * ω)
    Jmax = IPSII * ω
    Vcmax = mc > eps(FT) ? IPSII * (m / mc) * ωstar / (8 * θj) : FT(0)
    return Jmax, Vcmax
end

"""
    net_photosynthesis(Ac::FT,
                       Aj::FT,
                       Rd::FT,
                       β::FT) where {FT}

Computes the total net carbon assimilation (`An`),
in units of mol CO2/m^2/s, as a function of
the Rubisco limiting factor (`Ac`), the electron transport limiting rate (`Aj`),
dark respiration (`Rd`), and the moisture stress factor (`β`).

See Table 11.5 of G. Bonan's textbook, Climate Change and Terrestrial Ecosystem Modeling (2019).
"""
function net_photosynthesis(Ac::FT, Aj::FT, Rd::FT, β::FT) where {FT}
    An = max(0, min(Ac, Aj) * β - Rd)
    return An
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
    dark_respiration(is_c3::AbstractFloat, args...)

Calls the correct dark respiration function based on `is_c3`.

A `is_c3` value of 1.0 corresponds to C3 photosynthesis and calls
`c3_dark_respiration`, while 0.0 corresponds to C4 photsynthesis and calls
`c4_dark_respiration`.
"""
function dark_respiration(is_c3::AbstractFloat, args...)
    is_c3 > 0.5 ? c3_dark_respiration(args...) : c4_dark_respiration(args...)
end
"""
    c4_dark_respiration(VCmax25::FT,
                        β::FT,
                        T::FT,
                        R::FT
                        To::FT,
                        ::FT,
                        ::FT,
                        Q10::FT,
                        s5::FT,
                        s6::FT,
                        fC4::FT) where {FT}

Computes dark respiration (`Rd`),
in units of mol CO2/m^2/s, as a function of
 the moisture stress factor (`β`),
the unversal gas constant (`R`), the temperature (`T`),
Vcmax25, and
other parameters.

See Equation 11.73 of G. Bonan's textbook,
Climate Change and Terrestrial Ecosystem Modeling (2019).
"""
function c4_dark_respiration(
    Vcmax25::FT,
    β::FT,
    T::FT,
    R::FT,
    To::FT,
    ::FT,
    ::FT,
    Q10::FT,
    s5::FT,
    s6::FT,
    fC4::FT,
) where {FT}
    Rd = fC4 * Vcmax25 * β * Q10^((T - To) / 10) / (1 + exp(s5 * (T - s6)))
    return Rd
end

"""
    c3_dark_respiration(Vcmax25::FT, β::FT,
                        T::FT,
                        R::FT,
                        To::FT,
                        fC3::FT,
                        ΔHRd::FT,) where {FT}

Computes dark respiration (`Rd`),
in units of mol CO2/m^2/s, as a function of
 the moisture stress factor (`β`),
the unversal gas constant (`R`), and the temperature (`T`),
Vcmax25,  and
other parameters.

See Table 11.5 of G. Bonan's textbook,
Climate Change and Terrestrial Ecosystem Modeling (2019).
"""
function c3_dark_respiration(
    Vcmax25::FT,
    β::FT,
    T::FT,
    R::FT,
    To::FT,
    fC3::FT,
    ΔHRd::FT,
    _...,
) where {FT}
    Rd = fC3 * Vcmax25 * β * arrhenius_function(T, To, R, ΔHRd)
    return Rd
end

"""
    compute_GPP(An::FT,
             K::FT,
             LAI::FT,
             Ω::FT) where {FT}

Computes the total canopy photosynthesis (`GPP`) as a function of
the total net carbon assimilation (`An`), the extinction coefficient (`K`),
leaf area index (`LAI`) and the clumping index (`Ω`).
"""
function compute_GPP(An::FT, K::FT, LAI::FT, Ω::FT) where {FT}
    GPP = An * (1 - exp(-K * LAI * Ω)) / (K * Ω)
    return GPP
end


"""
    upscale_leaf_conductance(gs::FT, LAI::FT, T::FT, R::FT, P::FT) where {FT}

This currently takes a leaf conductance (moles per leaf area per second)
and (1) converts it to m/s, (2) upscales to the entire canopy, by assuming
the leaves in the canopy are in parallel and hence multiplying
by LAI.
"""
function upscale_leaf_conductance(
    gs::FT,
    LAI::FT,
    T::FT,
    R::FT,
    P::FT,
) where {FT}
    # TODO: Check what CLM does, and check if we can use the same function
    #  for GPP from An, and make more general.
    canopy_conductance = gs * LAI * (R * T) / P # convert to m s-1
    return canopy_conductance
end

"""
    arrhenius_function(T::FT, To::FT, R::FT, ΔH::FT)

Computes the Arrhenius function at temperature `T` given
the reference temperature `To=298.15K`, the universal
gas constant `R`, and the energy activation `ΔH`.

See Table 11.5 of G. Bonan's textbook,
Climate Change and Terrestrial Ecosystem Modeling (2019).
"""
function arrhenius_function(T::FT, To::FT, R::FT, ΔH::FT) where {FT}
    return exp(ΔH * (T - To) / (To * R * T))
end

"""
    MM_Kc(Kc25::FT,
          ΔHkc::FT,
          T::FT,
          To::FT,
          R::FT) where {FT}

Computes the Michaelis-Menten coefficient for CO2 (`Kc`),
in units of mol/mol,
as a function of its value at 25 °C (`Kc25`),
a constant (`ΔHkc`), a standard temperature (`To`),
the unversal gas constant (`R`), and the temperature (`T`).

See Table 11.5 of G. Bonan's textbook,
Climate Change and Terrestrial Ecosystem Modeling (2019).
"""
function MM_Kc(Kc25::FT, ΔHkc::FT, T::FT, To::FT, R::FT) where {FT}
    Kc = Kc25 * arrhenius_function(T, To, R, ΔHkc)
    return Kc
end

"""
    MM_Ko(Ko25::FT,
          ΔHko::FT,
          T::FT,
          To::FT,
          R::FT) where {FT}

Computes the Michaelis-Menten coefficient for O2 (`Ko`),
in units of mol/mol,
as a function of its value at 25 °C (`Ko25`),
a constant (`ΔHko`), a standard temperature (`To`),
the universal gas constant (`R`), and the temperature (`T`).

See Table 11.5 of G. Bonan's textbook,
Climate Change and Terrestrial Ecosystem Modeling (2019).
"""
function MM_Ko(Ko25::FT, ΔHko::FT, T::FT, To::FT, R::FT) where {FT}
    Ko = Ko25 * arrhenius_function(T, To, R, ΔHko)
    return Ko
end

"""
    compute_Vcmax(is_c3::AbstractFloat, args...)

Calls the correct Vcmax function based on `is_c3`.

A `is_c3` value of 1.0 corresponds to C3 photosynthesis and calls
`c3_compute_Vcmax`, while 0.0 corresponds to C4 photsynthesis and calls
`c4_compute_Vcmax`.
"""
function compute_Vcmax(is_c3::AbstractFloat, args...)
    is_c3 > 0.5 ? c3_compute_Vcmax(args...) : c4_compute_Vcmax(args...)
end

"""
    c4_compute_Vcmax(Vcmax25::FT, T::FT, R::FT, To::FT, ::FT, Q10::FT, s1::FT, s2::FT,
                     s3::FT, s4::FT) where {FT}

Computes the maximum rate of carboxylation of Rubisco (`Vcmax`),
in units of mol/m^2/s,
as a function of temperature (`T`), the universal
gas constant `R`, Vcmax25,  and other parameters.

For C4 photosynthesis, this uses Equation 11.73 from G. Bonan
Climate Change and Terrestrial Ecosystem Modeling (2019).
"""
function c4_compute_Vcmax(
    Vcmax25::FT,
    T::FT,
    R::FT,
    To::FT,
    ::FT,
    Q10::FT,
    s1::FT,
    s2::FT,
    s3::FT,
    s4::FT,
) where {FT}
    Vcmax =
        Vcmax25 * Q10^((T - To) / 10) / (1 + exp(s1 * (T - s2))) /
        (1 + exp(s3 * (s4 - T)))
    return Vcmax
end

"""
    c3_compute_Vcmax(Vcmax25::FT, T::FT, R::FT, To::FT, ΔHVcmax::FT) where {FT}

Computes the maximum rate of carboxylation of Rubisco (`Vcmax`),
in units of mol/m^2/s,
as a function of temperature (`T`), the universal
gas constant `R`, Vcmax25, and other parameters.

For C3 photosynthesis, this uses Table 11.5 from G. Bonan:
Climate Change and Terrestrial Ecosystem Modeling (2019).
"""
function c3_compute_Vcmax(
    Vcmax25::FT,
    T::FT,
    R::FT,
    To::FT,
    ΔHVcmax::FT,
    args...,
) where {FT}
    Vcmax = Vcmax25 * arrhenius_function(T, To, R, ΔHVcmax)
    return Vcmax
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
    nitrogen_content(
                     ne::FT, # Mean leaf nitrogen concentration (kg N (kg C)-1)
                     Vcmax25::FT, #
                     LAI::FT, # Leaf area index
                     SAI::FT,
                     RAI::FT,
                     ηsl::FT, # live stem  wood coefficient (kg C m-3)
                     h::FT, # canopy height (m)
                     σl::FT # Specific leaf density (kg C m-2 [leaf])
                     μr::FT, # Ratio root nitrogen to top leaf nitrogen (-), typical value 1.0
                     μs::FT, # Ratio stem nitrogen to top leaf nitrogen (-), typical value 0.1
                    ) where {FT}

Computes the nitrogen content of leafs (Nl), roots (Nr) and stems (Ns).
"""
function nitrogen_content(
    ne::FT, # Mean leaf nitrogen concentration (kg N (kg C)-1)
    Vcmax25::FT, #
    LAI::FT, # Leaf area index
    SAI::FT,
    RAI::FT,
    ηsl::FT, # live stem  wood coefficient (kg C m-3)
    h::FT, # canopy height (m)
    σl::FT, # Specific leaf density (kg C m-2 [leaf])
    μr::FT, # Ratio root nitrogen to top leaf nitrogen (-), typical value 1.0
    μs::FT, # Ratio stem nitrogen to top leaf nitrogen (-), typical value 0.1
) where {FT}
    Sc = ηsl * h * LAI * ClimaLand.heaviside(SAI)
    Rc = σl * RAI
    nm = Vcmax25 / ne
    Nl = nm * σl * LAI
    Nr = μr * nm * Rc
    Ns = μs * nm * Sc
    return Nl, Nr, Ns
end

"""
    plant_respiration_maintenance(
        Rd::FT, # Dark respiration
        β::FT, # Soil moisture factor
        Nl::FT, # Nitrogen content of leafs
        Nr::FT, # Nitrogen content of roots
        Ns::FT, # Nitrogen content of stems
        ) where {FT}

Computes plant maintenance respiration as a function of dark respiration (Rd),
the nitrogen content of leafs (Nl), roots (Nr) and stems (Ns),
and the soil moisture factor (β).
"""
function plant_respiration_maintenance(
    Rd::FT, # Dark respiration
    β::FT, # Soil moisture factor
    Nl::FT, # Nitrogen content of leafs
    Nr::FT, # Nitrogen content of roots
    Ns::FT, # Nitrogen content of stems
) where {FT}
    # When LAI is zero, Nl = 0
    Rpm = Rd * (β + (Nr + Ns) / max(Nl, eps(FT)))
    return Rpm
end

"""
    plant_respiration_growth(
        Rel::FT, # Factor of relative contribution
        An::FT, # Net photosynthesis
        Rpm::FT # Plant maintenance respiration
        ) where {FT}

Computes plant growth respiration as a function of net photosynthesis (An),
plant maintenance respiration (Rpm), and a relative contribution factor, Rel.
"""
function plant_respiration_growth(Rel::FT, An::FT, Rpm::FT) where {FT}
    Rg = Rel * (An - Rpm)
    return Rg
end

"""
    enforce_albedo_constraint(α, τ)
A function which enforces α+τ <= 1.
"""
enforce_albedo_constraint(α, τ) = 1 - α - τ > 0 ? α : 1 - τ


# P-model 
"""
    intrinsic_quantum_yield(
        T::FT, 
        c::FT,
        ϕa0::FT,
        ϕa1::FT,
        ϕa2::FT
    ) where {FT}

    Computes the intrinsic quantum yield of photosynthesis ϕ (mol/mol) 
    as a function of temperature T (K) and a calibratable parameter c (unitless). 
    The functional form given in Bernacchi et al (2003) and used in Stocker 
    et al. (2020) is a second order polynomial in T (deg C) with coefficients ϕa0, 
    ϕa1, and ϕa2.
"""
function intrinsic_quantum_yield(T::FT, c::FT, ϕa0::FT, ϕa1::FT, ϕa2::FT) where {FT}
    # convert to C
    T = T - FT(273.15)
    ϕ = c * (ϕa0 + ϕa1 * T + ϕa2 * T^2)
    return max(ϕ, FT(0)) # Ensure non-negative quantum yield
end


"""
    density_h20(
        T::FT, 
        p::FT
    ) where {FT}

    Computes the density of water in kg/m^3 given temperature T (K) and pressure p (Pa)
    according to F.H. Fisher and O.E Dial, Jr. (1975) Equation of state of pure water and 
    sea water, Tech. Rept., Marine Physical Laboratory, San Diego, CA.

    (can consider removing the higher order terms?) 
"""
function density_h2o(T::FT, p::FT) where {FT}
    # Convert temperature to Celsius
    T = T - FT(273.15) 
    # Convert pressure to bar (1 bar = 1e5 Pa)
    pbar = FT(1e-5) * p

    # λ(T): bar cm3/g
    λ = FT(1788.316) +
        FT(21.55053) * T -
        FT(0.4695911) * T^2 +
        FT(3.096363e-3) * T^3 -
        FT(7.341182e-6) * T^4

    # p0(T): bar
    p0 = FT(5918.499) +
         FT(58.05267) * T -
         FT(1.1253317) * T^2 +
         FT(6.6123869e-3) * T^3 -
         FT(1.4661625e-5) * T^4

    # v_inf(T): cm3/g
    T_powers = (T .^ (0:9))  # T^0 to T^9
    coeffs_vinf = FT.([
        0.6980547,
       -7.435626e-4,
        3.704258e-5,
       -6.315724e-7,
        9.829576e-9,
       -1.197269e-10,
        1.005461e-12,
       -5.437898e-15,
        1.69946e-17,
       -2.295063e-20
    ])
    v_inf = dot(coeffs_vinf, T_powers)

    # Specific volume [cm3/g]
    v = v_inf + λ / (p0 + pbar)

    # Convert to density [kg/m3]
    return FT(1e3) / v
end


"""
    viscosity_h20(
        T::FT, 
        p::FT,
        constant_density::Bool
    ) where {FT}

    Computes the viscosity of water in Pa s given temperature T (K) and pressure p (Pa)
    according to Huber et al. (2009) [https://doi.org/10.1063/1.3088050]. 

    Can consider simplifying if this level of precision is not needed
"""
function viscosity_h2o(T::FT, p::FT, constant_density::Bool) where {FT}
    # Reference constants
    tk_ast  = FT(647.096)    # K
    ρ_ast = FT(322.0)      # kg/m^3
    μ_ast  = FT(1e-6)       # Pa s

    # Get density of water [kg/m^3]
    earth_param_set = LP.LandParameters(FT)
    ρ0 = LP.ρ_cloud_liq(earth_param_set)
    constant_density ? ρ = ρ0 : ρ = density_h2o(T, p)

    # Dimensionless variables
    tbar  = T / tk_ast
    tbarx = sqrt(tbar)
    tbar2 = tbar^2
    tbar3 = tbar^3
    ρbar  = ρ / ρ_ast
    
    # Calculate μ0 (Eq. 11 & Table 2)
    μ0 = FT(1.67752) + FT(2.20462) / tbar + FT(0.6366564) / tbar2 - FT(0.241605) / tbar3
    μ0 = FT(1e2) * tbarx / μ0

    # Coefficients h_array from Table 3
    h_array = FT.([
        0.520094    0.0850895   -1.08374   -0.289555   0.0         0.0;
        0.222531    0.999115     1.88797    1.26613    0.0         0.120573;
       -0.281378   -0.906851    -0.772479  -0.489837  -0.257040    0.0;
        0.161913    0.257399     0.0        0.0        0.0         0.0;
       -0.0325372   0.0          0.0        0.0698452  0.0         0.0;
        0.0         0.0          0.0        0.0        0.00872102  0.0;
        0.0         0.0          0.0       -0.00435673 0.0        -0.000593264
    ])

    # Compute μ1 (Eq. 12 & Table 3)
    μ1 = FT(0.0)
    ctbar = (FT(1.0) / tbar) - FT(1.0)
    for i in 1:6
        coef1 = ctbar^(i - 1)
        coef2 = FT(0.0)
        for j in 1:7
            coef2 += h_array[j, i] * (ρbar - FT(1.0))^(j - 1)
        end
        μ1 += coef1 * coef2
    end
    μ1 = exp(ρbar * μ1)

    μ_bar = μ0 * μ1       # Eq. 2
    μ = μ_bar * μ_ast     # Eq. 1
    return FT(μ)          # Pa s
end


"""
    compute_viscosity_ratio(
        T::FT, 
        p::FT,
        constant_density::Bool
    ) where {FT}

    Computes η*, the ratio of the viscosity of water at temperature T and pressure p
    to the viscosity of water at STP. If `constant_density` is true, the density of water
    is taken to be constant (1000.0 kg/m^3). Otherwise we use an EOS to compute the density
    at the given temperature and pressure. 
"""
function compute_viscosity_ratio(T::FT, p::FT, constant_density::Bool) where {FT} 
    η25 = viscosity_h2o(FT(298.15), FT(101325.0), constant_density)
    ηstar = viscosity_h2o(T, p, constant_density) / η25
    return FT(ηstar)
end


"""
    po2(
        P_air::FT, 
        oi::FT
    ) where {FT}

    Computes the partial pressure of O2 in the air (Pa) given atmospheric pressure (`P_air`)
    and a constant mixing ratio of O2 (`oi`), typically 0.209. 
"""
function po2(P_air::FT, oi::FT) where {FT} 
    return oi * P_air 
end

"""
    co2_compensation_p(
        T::FT,
        To::FT,
        p::FT,
        R::FT, 
        ΔHΓstar::FT
        Γstar25::FT
    ) where {FT}

    Computes the CO2 compensation point (`Γstar`), in units Pa, as a function of temperature T (K)
    and pressure p (Pa). See Equation B5 of Stocker et al. (2020). 
"""
function co2_compensation_p(
    T::FT,
    To::FT,
    p::FT,
    R::FT, 
    ΔHΓstar::FT,
    Γstar25::FT 
) where {FT}
    Γstar = Γstar25 * p / FT(101325.0) * arrhenius_function(T, To, R, ΔHΓstar)
    return Γstar
end

"""
    compute_Kmm(
        T::FT, 
        p::FT, 
        Kc25::FT, 
        Ko25::FT, 
        ΔHkc::FT, 
        ΔHko::FT, 
        To::FT, 
        R::FT,
        oi::FT
    ) where {FT}

    Computes the effective Michaelis-Menten coefficient for Rubisco-limited photosynthesis (`Kmm`),
    in units Pa, as a function of temperature T (K), atmospheric pressure p (Pa), and constants:
    Kc25 (Michaelis-Menten coefficient for CO2 at 25 °C), Ko25 (Michaelis-Menten coefficient for O2 at 25 °C),
    ΔHkc (effective enthalpy of activation for Kc), ΔHko (effective enthalpy of activation for Ko),
    To (reference temperature, typically 298.15 K), R (universal gas constant), and oi (O2 mixing ratio, 
    typically 0.209).
"""
function compute_Kmm(
    T::FT, 
    p::FT, 
    Kc25::FT, 
    Ko25::FT, 
    ΔHkc::FT, 
    ΔHko::FT, 
    To::FT, 
    R::FT,
    oi::FT
) where {FT}
    Kc = MM_Kc(Kc25, ΔHkc, T, To, R)
    Ko = MM_Ko(Ko25, ΔHko, T, To, R)

    return Kc * (1 + po2(p, oi) / Ko) 
end

"""
    optimal_co2_ratio_c3(
        Kmm::FT,
        Γstar::FT, 
        ηstar::FT, 
        ca::FT, 
        VPD::FT, 
        β::FT,
        Drel::FT 
    ) where {FT}

    The p-model assumptions, that 1) plants optimize the relative costs of transpiration per unit
    carbon assimlated and costs of maintaining carboxylation capacity per unit carbon assimilated;
    2) coordination hypothesis (assimilation is limited simultaneously by both light and Rubisco) 
    are applied to compute the optimal ratio of intercellular to ambient CO2 concentration (`χ`)
    and auxiliary variables ξ, mj, and mc. mj and mc represent capacities for light and Rubisco-
    limited photosynthesis, respectively. 

    Parameters: Kmm (effective Michaelis-Menten coefficient for Rubisco-limited photosynthesis, Pa),
    Γstar (CO2 compensation point, Pa), ηstar (viscosity ratio), ca (ambient CO2 partial pressure, Pa),
    VPD (vapor pressure deficit, Pa), β (moisture stress factor, unitless), Drel = 1.6 (relative 
    diffusivity of water vapor with respect to CO2, unitless).
"""
function optimal_co2_ratio_c3(
    Kmm::FT,
    Γstar::FT, 
    ηstar::FT, 
    ca::FT, 
    VPD::FT, 
    β::FT,
    Drel::FT 
) where {FT}
    ξ = sqrt(β * (Kmm + Γstar) / (Drel * ηstar))
    χ = Γstar / ca + (1 - Γstar / ca) * ξ / (ξ + sqrt(VPD)) 

    # define some auxiliary variables 
    γ = Γstar / ca 
    κ = Kmm / ca 

    mj = (χ - γ) / (χ + 2 * γ) # eqn 11 in Stocker et al. (2020)
    mc = (χ - γ) / (χ + κ) # eqn 7 in Stocker et al. (2020)

    return χ, ξ, mj, mc
end 

function optimal_ξ_c3(
    Kmm::FT,
    Γstar::FT, 
    ηstar::FT, 
    β::FT,
    Drel::FT 
) where {FT}
    return sqrt(β * (Kmm + Γstar) / (Drel * ηstar))
end 

function compute_mj(
    ξ::FT,
    Γstar::FT,
    ca::FT,
    VPD::FT
) where {FT}
    γ = Γstar / ca 
    χ = Γstar / ca + (1 - Γstar / ca) * ξ / (ξ + sqrt(VPD)) 
    return (χ - γ) / (χ + 2 * γ)
end

function compute_mc(
    ξ::FT,
    Kmm::FT,
    Γstar::FT,
    ca::FT
) where {FT}
    γ = Γstar / ca 
    κ = Kmm / ca 
    χ = Γstar / ca + (1 - Γstar / ca) * ξ / (ξ + sqrt(VPD)) 
    return (χ - γ) / (χ + κ)
end

function compute_ci(
    ξ::FT,
    ca::FT, 
    Γstar::FT,
    VPD::FT
) where {FT}
    return (ca * ξ + Γstar * sqrt(VPD)) / (ξ + sqrt(VPD))
end

"""
    pmodel_gs(
        χ::FT, 
        ca::FT,
        A::FT
    ) where {FT}

    Computes the stomatal conductance of CO2 (`gs`), in units of mol CO2/m^2/s
    via Fick's law. Parameters are the ratio of intercellular to ambient CO2 
    concentration (`χ`), the ambient CO2 partial pressure (`ca`, in Pa), and the 
    assimilation rate (`A`). This is related to the conductance of H2O by a 
    factor Drel = 1.6. 
"""
function pmodel_gs(
    χ::FT, 
    ca::FT,
    A::FT
) where {FT}
    return A / (ca * (1 - χ)) 
end

"""
    compute_mj_with_jmax_limitation(
        mj::FT,
        cstar::FT
    ) where {FT}

    Computes m' such that Aj = ϕ0 I_abs * m' (a LUE model) by assuming that dA/dJmax = c
    is constant. cstar is defined as 4c, a free parameter. Wang etal (2017) derive cstar = 0.412
    at STP and using Vcmax/Jmax = 1.88. 
"""
# TODO: test if cstar should be made a free parameter? 
function compute_mj_with_jmax_limitation(
    mj::FT, 
    cstar::FT
) where {FT}
    arg = 1 - (cstar / mj)^(FT(2/3))
    sqrt_arg = ifelse(arg < 0, FT(0.0), sqrt(arg)) # avoid complex numbers
    return FT(mj * sqrt_arg)
end 


"""
    compute_LUE(
        ϕ0::FT, 
        β::FT,
        mprime::FT,
        Mc::FT
    ) where {FT} 

    Computes light use efficiency (LUE) in kg C/mol from intrinsic quantum yield (`ϕ0`),
    moisture stress factor (`β`), and a Jmax modified capacity (`mprime`); see Eqn 17 and 19
    in Stocker et al. (2020). Mc is the molar mass of carbon (kg/mol) = 0.0120107 kg/mol.
"""
function compute_LUE(
    ϕ0::FT, 
    β::FT,
    mprime::FT,
    Mc::FT 
) where {FT} 
    return ϕ0 * β * mprime * Mc 
end 


"""
    pmodel_vcmax(
        ϕ0::FT,
        I_abs::FT,
        mprime::FT,
        mc::FT
        βm::FT
    ) where {FT}

    Computes the maximum rate of carboxylation assuming optimality and Aj = Ac using 
    the intrinsic quantum yield (`ϕ0`), absorbed radiation (`I_abs`), Jmax-adjusted capacity
    (`mprime`), a Rubisco-limited capacity (`mc`), and empirical soil moisture stress factor
    (`βm`). See Eqns 16 and 6 in Stocker et al. (2020). 
"""
function pmodel_vcmax(
    ϕ0::FT, 
    I_abs::FT,
    mprime::FT,
    mc::FT,
    βm::FT
) where {FT}
    Vcmax = βm * ϕ0 * I_abs * mprime / mc 
    return Vcmax
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
    θ1::FT = FT(0.6) 
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

"""
    electron_transport_pmodel(
        ϕ0::FT,
        I_abs::FT,
        Jmax::FT 
    ) where {FT}
    
    Computes the rate of electron transport (`J`) in mol electrons/m^2/s for the pmodel.
"""
function electron_transport_pmodel(
    ϕ0::FT,
    I_abs::FT,
    Jmax::FT 
) where {FT}
    J = FT(4) * ϕ0 * I_abs / sqrt(FT(1) + (FT(4) * ϕ0 * I_abs / max(Jmax, eps(FT)))^2) 
    return J
end



"""
inst_temp_scaling(
    T_canopy::FT;
    T_acclim::FT = T_canopy,
    To::FT,
    Ha::FT,
    Hd::FT,
    aS::FT,
    bS::FT,
    R::FT
) where {FT}

Given Vcmax or Jmax that have acclimated according to T_acclim, this function computes
the instantaneous temperature scaling factor f ∈ [0, ∞) for these maximum rates at the 
instantaneous current temperature T_canopy. To is a reference temperature for the constants
and should be set to 298.15 K (25 °C). By default we assume that T_acclim = T_canopy. 

The parameters (`Ha`, `Hd`, `aS`, `bS`) come from Kattge & Knorr (2007)

| Quantity | Ha (J/mol) | Hd (J/mol) | aS (J/mol/K) | bS (J/mol/K^2) |
|----------|------------|------------|--------------|----------------|
| Vcmax    |   71 513   | 200 000    |   668.39     | 1.07           |
| Jmax     |   49 884   | 200 000    |   659.70     | 0.75           |

"""
function inst_temp_scaling(
    T_canopy::FT,
    T_acclim::FT,
    To::FT,
    Ha::FT,
    Hd::FT,
    aS::FT,
    bS::FT,
    R::FT
 ) where {FT}
    T_acclim = T_acclim - FT(273.15)    # °C for ΔS(T)
    ΔS        = aS - bS * T_acclim      # entropy term (J mol⁻¹ K⁻¹)

    # Arrhenius-type activation scaling factor
    f_act = arrhenius_function(T_canopy, To, R, Ha)

    # high temperature deactivation scaling factor
    num = 1 + exp( (To * ΔS - Hd) / (R * To) )
    den = 1 + exp( (T_canopy * ΔS - Hd) / (R * T_canopy) )
    f_deact = num / den

    return f_act * f_deact 
end


"""
inst_temp_scaling_rd(
    T_canopy::FT,
    To::FT, 
    aRd::FT, 
    bRd::FT
) where {FT}

Computes the instantaneous temperature scaling factor for dark respiration (Rd) 
at canopy temperature `T_canopy` given reference temperature `To`, the first order
coefficient `aRd`, and the second order coefficient `bRd`. 

Usees the log-quadratic functional form of Heskel et al. (2016) 
https://www.pnas.org/doi/full/10.1073/pnas.1520282113
"""
function inst_temp_scaling_rd(
    T_canopy::FT,
    To::FT, 
    aRd::FT, 
    bRd::FT
) where {FT}
    return exp(
        aRd * (T_canopy - To) + bRd * ((T_canopy - FT(273.15))^2 - (To - FT(273.15))^2)
    )
end
