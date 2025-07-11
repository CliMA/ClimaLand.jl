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
    enforce_albedo_constraint

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
    @. p.canopy.radiative_transfer.par = canopy_sw_rt_two_stream_analytical(
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
    @. p.canopy.radiative_transfer.nir = canopy_sw_rt_two_stream_analytical(
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

# todo: docstring
function canopy_sw_rt_two_stream_analytical(
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
    aₛ =
        FT(0.5) *
        ω *
        (1 - cosθs * log((abs(cosθs) + 1) / max(abs(cosθs), eps(FT))))
    β₀ = (1 / ω) * aₛ * (1 + μ̄ * K) / (μ̄ * K)
    diff = α_leaf - τ_leaf
    # With uniform distribution, Dickinson integral becomes following:
    c²θ̄ = FT(pi) * G / FT(4)
    β = FT(0.5) * (ω + diff * c²θ̄) / ω
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
        FT(1) / d₁ * (
            (d - h₁ / σ * p₃) * (u₁ - μ̄ * h) / s₁ -
            p₂ * s₂ * (d - c - h₁ / σ * (u₁ + μ̄ * K))
        )
    h₃ =
        FT(-1) / d₁ * (
            (d - h₁ / σ * p₃) * (u₁ + μ̄ * h) * s₁ -
            p₁ * s₂ * (d - c - h₁ / σ * (u₁ + μ̄ * K))
        )

    # h coefficients for direct downward flux
    h₄ = -f * p₃ - c * d
    h₅ =
        FT(-1) / d₂ *
        (h₄ * (u₂ + μ̄ * h) / (σ * s₁) + (u₃ - h₄ / σ * (u₂ - μ̄ * K)) * s₂)
    h₆ =
        FT(1) / d₂ *
        (h₄ / σ * (u₂ - μ̄ * h) * s₁ + (u₃ - h₄ / σ * (u₂ - μ̄ * K)) * s₂)
    # h coefficients for diffuse upward flux
    h₇ = c * (u₁ - μ̄ * h) / (d₁ * s₁)
    h₈ = -c * s₁ * (u₁ + μ̄ * h) / d₁
    h₉ = (u₂ + μ̄ * h) / (d₂ * s₁)
    h₁₀ = -s₁ * (u₂ - μ̄ * h) / d₂
    L = LAI
    I_dir_up =
        h₁ * exp(-K * L * Ω) / σ + h₂ * exp(-h * L * Ω) + h₃ * exp(h * L * Ω)
    I_dir_dn =
        h₄ * exp(-K * L * Ω) / σ + h₅ * exp(-h * L * Ω) + h₆ * exp(h * L * Ω)
    I_dir_dn += exp(-K * L * Ω)

    # Compute the diffuse fluxes into/out of the layer
    I_dif_up = h₇ * exp(-h * L * Ω) + h₈ * exp(h * L * Ω)
    I_dif_dn = h₉ * exp(-h * L * Ω) + h₁₀ * exp(h * L * Ω)
    I_dir_abs = I_dir_up + I_dir_dn
    I_dif_abs = I_dif_up + I_dif_dn
    @assert I_dir_abs >= 0
    @assert I_dif_abs >= 0
    F_abs = (1 - frac_diff) * I_dir_abs + (frac_diff) * I_dif_abs
    F_refl = (1 - frac_diff) * I_dir_up + (frac_diff) * I_dif_up
    F_trans = (1 - frac_diff) * I_dir_dn + (frac_diff) * I_dif_dn
    # Main.@infiltrate
    return (; abs = F_abs, refl = F_refl, trans = F_trans)
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

TODO: Check what CLM does, and check if we can use the same function
for GPP from An, and make more general.
"""
function upscale_leaf_conductance(
    gs::FT,
    LAI::FT,
    T::FT,
    R::FT,
    P::FT,
) where {FT}
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
"""
function medlyn_term(
    g1::FT,
    T_air::FT,
    P_air::FT,
    q_air::FT,
    thermo_params,
) where {FT}
    VPD = ClimaLand.vapor_pressure_deficit(T_air, P_air, q_air, thermo_params)
    return 1 + g1 / sqrt(VPD)
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
