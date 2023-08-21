using ..ClimaLSM.Canopy
export plant_absorbed_pfd,
    extinction_coeff,
    arrhenius_function,
    intercellular_co2,
    co2_compensation,
    rubisco_assimilation,
    light_assimilation,
    max_electron_transport,
    electron_transport,
    net_photosynthesis,
    moisture_stress,
    dark_respiration,
    compute_GPP,
    MM_Kc,
    MM_Ko,
    compute_Vcmax,
    medlyn_term,
    medlyn_conductance,
    canopy_surface_fluxes,
    upscale_leaf_conductance,
    penman_monteith,
    nitrogen_content,
    plant_respiration_maintenance,
    plant_respiration_growth

# 1. Radiative transfer

"""
    compute_absorbances(
        RT::BeerLambertModel{FT},
        PAR,
        NIR,
        LAI,
        K,
        α_soil_PAR,
        α_soil_NIR,
        _,
        _,
    )

Computes the PAR and NIR absorbances, reflectances, and tranmittances
for a canopy in the case of the 
Beer-Lambert model. The absorbances are a function of the radiative transfer 
model, as well as the magnitude of incident PAR and NIR radiation in moles of 
photons, the leaf area index, the extinction coefficient, and the
soil albedo in the PAR and NIR bands. Returns a
NamedTuple of NamedTuple, of the form:
(; par = (; refl = , trans = , abs = ),  nir = (; refl = , trans = , abs = ))
"""
function compute_absorbances(
    RT::BeerLambertModel{FT},
    PAR,
    NIR,
    LAI,
    K,
    α_soil_PAR,
    α_soil_NIR,
    _...,
) where {FT}
    RTP = RT.parameters
    canopy_par =
        @. plant_absorbed_pfd(RT, PAR, RTP.α_PAR_leaf, LAI, K, α_soil_PAR)
    canopy_nir =
        @. plant_absorbed_pfd(RT, NIR, RTP.α_NIR_leaf, LAI, K, α_soil_NIR)
    return (; par = canopy_par, nir = canopy_nir)
end

"""
    compute_absorbances(
        RT::TwoStreamModel{FT},
        PAR,
        NIR,
        LAI,
        K,
        α_soil_PAR,
        α_soil_NIR,
        θs,
        frac_diff,
    )

Computes the PAR and NIR absorbances, reflectances, and tranmittances
for a canopy in the case of the 
Beer-Lambert model. The absorbances are a function of the radiative transfer 
model, as well as the magnitude of incident PAR and NIR radiation in moles of 
photons, the leaf area index, the extinction coefficient, and the
soil albedo in the PAR and NIR bands. 

This model also depends on the diffuse fraction and the zenith angle.
Returns a
NamedTuple of NamedTuple, of the form:
(; par = (; refl = , trans = , abs = ),  nir = (; refl = , trans = , abs = ))
"""
function compute_absorbances(
    RT::TwoStreamModel{FT},
    PAR,
    NIR,
    LAI,
    K,
    α_soil_PAR,
    α_soil_NIR,
    θs,
    frac_diff,
) where {FT}
    RTP = RT.parameters
    canopy_par = @. plant_absorbed_pfd(
        RT,
        PAR,
        RTP.α_PAR_leaf,
        RTP.τ_PAR_leaf,
        LAI,
        K,
        θs,
        α_soil_PAR,
        frac_diff,
    )
    canopy_nir = @. plant_absorbed_pfd(
        RT,
        NIR,
        RTP.α_NIR_leaf,
        RTP.τ_NIR_leaf,
        LAI,
        K,
        θs,
        α_soil_NIR,
        frac_diff,
    )
    return (; par = canopy_par, nir = canopy_nir)
end

"""
    plant_absorbed_pfd(
        RT::BeerLambertModel{FT},
        SW_IN:FT,
        α_leaf::FT,
        LAI::FT,
        K::FT,
        α_soil::FT
    )

Computes the absorbed, reflected, and transmitted photon flux density 
in terms of mol photons per m^2 per 
second for a radiation band.

This applies the Beer-Lambert law, which is a function of incident 
radiation (`SW_IN`; moles of photons/m^2/), leaf reflectance
(`α_leaf`), the extinction coefficient (`K`), leaf area index (`LAI`),
and the albedo of the soil (`α_soil`). 

Returns a tuple of reflected, absorbed, and transmitted radiation in
mol photons/m^2/s.
"""
function plant_absorbed_pfd(
    RT::BeerLambertModel{FT},
    SW_IN::FT,
    α_leaf::FT,
    LAI::FT,
    K::FT,
    α_soil::FT,
) where {FT}
    RTP = RT.parameters
    AR = SW_IN * (1 - α_leaf) * (1 - exp(-K * LAI * RTP.Ω)) * (1 - α_soil)
    TR = SW_IN * exp(-K * LAI * RTP.Ω)
    RR = SW_IN - AR - TR * (1 - α_soil)
    return (; abs = AR, refl = RR, trans = TR)
end

"""
    plant_absorbed_pfd(
        RT::TwoStreamModel{FT},
        α_leaf,
        SW_IN::FT,
        LAI::FT,
        K::FT,
        τ_leaf,
        θs::FT,
        α_soil::FT,
    )

Computes the absorbed, transmitted, and reflected  photon flux density 
in terms of mol photons  per m^2 per second for a radiation band. 

This applies the two-stream radiative transfer solution which takes into account
the impacts of scattering within the canopy. The function takes in all 
parameters from the parameter struct of a TwoStreamModel, along with the 
incident radiation, LAI, extinction coefficient K, soil albedo from the 
canopy soil_driver, solar zenith angle, and τ.

Returns a tuple of reflected, absorbed, and transmitted radiation in
mol photons/m^2/s.
"""
function plant_absorbed_pfd(
    RT::TwoStreamModel{FT},
    SW_IN::FT,
    α_leaf::FT,
    τ_leaf::FT,
    LAI::FT,
    K::FT,
    θs::FT,
    α_soil::FT,
    frac_diff::FT,
) where {FT}

    (; ld, Ω, n_layers) = RT.parameters

    # Compute μ̄, the average inverse diffuse optical length per LAI
    μ̄ = 1 / (2 * ld)

    ω = α_leaf + τ_leaf

    # Compute aₛ, the single scattering albedo
    aₛ = 0.5 * ω * (1 - cos(θs) * log((abs(cos(θs)) + 1) / abs(cos(θs))))

    # Compute β₀, the direct upscattering parameter
    β₀ = (1 / ω) * aₛ * (1 + μ̄ * K) / (μ̄ * K)

    # Compute β, the diffuse upscattering parameter
    diff = α_leaf - τ_leaf
    # With uniform distribution, Dickinson integral becomes following:
    c²θ̄ = pi * ld / 4
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

    # Total light reflected form top of canopy
    F_refl = 0

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

        # Add collimated radiation to downard flux
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

        if i == 1
            F_refl = (1 - frac_diff) * I_dir_up + (frac_diff) * I_dif_up
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

    # Convert fractional absorption into absorption and return
    # Ensure floating point precision is correct (it may be different for PAR)
    F_trans = (1 - F_abs - F_refl) / (1 - α_soil)
    return (;
        abs = FT(SW_IN * F_abs),
        refl = FT(SW_IN * F_refl),
        trans = FT(SW_IN * F_trans),
    )
end

"""
    diffuse_fraction(td::FT, T::FT, SW_IN::FT, RH::FT, θs::FT) where {FT}

Computes the fraction of diffuse radiation (`diff_frac`) as a function
of the solar zenith angle (`θs`), the total surface incident shortwave radiation (`SW_IN`),
the air temperature (`T`), the relative humidity (`RH`), and the day of the year
(`td`). 

See Appendix A of Braghiere, "Evaluation of turbulent fluxes of CO2, sensible heat,
and latent heat as a function of aerosol optical depth over the course of deforestation
in the Brazilian Amazon" 2013.

Note that cos(θs) is equal to zero when θs = π/2, and this is a coefficient
of k₀, which we divide by in this expression. This can amplify small errors
when θs is near π/2. 

This formula is empirical and can yied negative numbers depending on the 
input, which, when dividing by something very near zero, 
can become large negative numbers.

Because of that, we cap the returned value to lie within [0,1].
"""
function diffuse_fraction(td, T::FT, SW_IN::FT, RH::FT, θs::FT) where {FT}
    k₀ = FT(1370 * (1 + 0.033 * cos(2π * td / 365))) * cos(θs)
    kₜ = SW_IN / k₀
    if kₜ ≤ 0.3
        diff_frac = FT(
            kₜ *
            (1 - 0.232 * kₜ + 0.0239 * cos(θs) - 6.82e-4 * T + 0.0195 * RH),
        )
    elseif kₜ ≤ 0.78
        diff_frac = FT(
            kₜ * (
                1.329 - 1.716 * kₜ + 0.267 * cos(θs) - 3.57e-3 * T + 0.106 * RH
            ),
        )
    else
        diff_frac =
            FT(kₜ * (0.426 * kₜ + 0.256 * cos(θs) - 3.49e-3 * T + 0.0734 * RH))
    end
    return max(min(diff_frac, FT(1)), FT(0))
end


"""
    extinction_coeff(ld::FT,
                     θs::FT) where {FT}

Computes the vegetation extinction coefficient (`K`), as a function
of the sun zenith angle (`θs`), and the leaf angle distribution (`ld`).
"""
function extinction_coeff(ld::FT, θs::FT) where {FT}
    K = ld / max(cos(θs), eps(FT))
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
    rubisco_assimilation(::C3,
                         Vcmax::FT,
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
function rubisco_assimilation(
    ::C3,
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
    rubisco_assimilation(::C4, Vcmax::FT,_...) where {FT}

Computes the Rubisco limiting rate of photosynthesis for C4 plants (`Ac`)
in units of moles CO2/m^2/s,
as equal to the maximum rate of carboxylation of Rubisco (`Vcmax`).
"""
function rubisco_assimilation(::C4, Vcmax::FT, _...) where {FT}
    Ac = Vcmax
    return Ac
end

"""
    light_assimilation(::C3,
                       J::FT,
                       ci::FT,
                       Γstar::FT) where {FT}

Computes the electron transport limiting rate (`Aj`),
in units of moles CO2/m^2/s, for C3 plants as a function of
the rate of electron transport (`J`), the leaf internal carbon dioxide partial pressure (`ci`),
and the CO2 compensation point (`Γstar`).

See Table 11.5 of G. Bonan's textbook, Climate Change and Terrestrial Ecosystem Modeling (2019).
"""
function light_assimilation(::C3, J::FT, ci::FT, Γstar::FT) where {FT}
    Aj = J * (ci - Γstar) / (4 * (ci + 2 * Γstar))
    return Aj
end

"""
    light_assimilation(::C4, J::FT, _...) where {FT}

Computes the electron transport limiting rate (`Aj`),
in units of moles CO2/m^2/s, for C4 plants, as equal to
the rate of electron transport (`J`).
"""
function light_assimilation(::C4, J::FT, _...) where {FT}
    Aj = J
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
    dark_respiration(Vcmax25::FT,
                     β::FT,
                     f::FT,
                     ΔHkc::FT,
                     T::FT,
                     To::FT,
                     R::FT) where {FT}

Computes dark respiration (`Rd`),
in units of mol CO2/m^2/s, as a function of the maximum rate of carboxylation of Rubisco (`Vcmax25`),
and the moisture stress factor (`β`), an empirical factor `f` is equal to 0.015,
a constant (`ΔHRd`), a standard temperature (`To`),
the unversal gas constant (`R`), and the temperature (`T`).

See Table 11.5 of G. Bonan's textbook, 
Climate Change and Terrestrial Ecosystem Modeling (2019).
"""
function dark_respiration(
    Vcmax25::FT,
    β::FT,
    f::FT,
    ΔHRd::FT,
    T::FT,
    To::FT,
    R::FT,
) where {FT}
    Rd = f * Vcmax25 * β * arrhenius_function(T, To, R, ΔHRd)
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
    compute_Vcmax(Vcmax25::FT,
           T::FT,
           To::FT,
           R::FT,
           ep5::FT) where {FT}

Computes the maximum rate of carboxylation of Rubisco (`Vcmax`),
in units of mol/m^2/s, 
as a function of temperature (`T`), Vcmax at the reference temperature 25 °C (`Vcmax25`),
the universal gas constant (`R`), and the reference temperature (`To`).

See Table 11.5 of G. Bonan's textbook, 
Climate Change and Terrestrial Ecosystem Modeling (2019).
"""
function compute_Vcmax(
    Vcmax25::FT,
    T::FT,
    To::FT,
    R::FT,
    ΔHVcmax::FT,
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
    VPD = ClimaLSM.vapor_pressure_deficit(T_air, P_air, q_air, thermo_params)
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
                     RAI::FT,
                     ηsl::FT, # live stem  wood coefficient (kg C m-3) 
                     h::FT, # canopy height (m)
                     σl::FT # Specific leaf density (kg C m-2 [leaf])
                     μr::FT, # Ratio root nitrogen to top leaf nitrogen (-), typical value 1.0
                     μs::FT, # Ratio stem nitrogen to top leaf nitrogen (-), typical value 0.1 
                    ) where {FT}

Computes the nitrogen content of leafs (Nl), roots (Nr) and stems (Ns)
as a function of leaf area index (LAI), specific leaf density (σl),
the carbon content of roots (Rc), the carbon content of stems (Rs), 
and mean leaf nitrogen concentration (nm).
"""
function nitrogen_content(
    ne::FT, # Mean leaf nitrogen concentration (kg N (kg C)-1)
    Vcmax25::FT, #
    LAI::FT, # Leaf area index
    RAI::FT,
    ηsl::FT, # live stem  wood coefficient (kg C m-3) 
    h::FT, # canopy height (m)
    σl::FT, # Specific leaf density (kg C m-2 [leaf])
    μr::FT, # Ratio root nitrogen to top leaf nitrogen (-), typical value 1.0
    μs::FT, # Ratio stem nitrogen to top leaf nitrogen (-), typical value 0.1 
) where {FT}
    Sc = ηsl * h * LAI
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
        f::FT # Factor to convert from mol CO2 to kg C
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
    f1::FT, # Factor to convert from mol CO2 to kg C
) where {FT}
    Rpm = f1 * Rd * (β + (Nr + Ns) / Nl)
    return Rpm
end

"""
    plant_respiration_growth(
        f::FT, # Factor of relative contribution
        GPP::FT, # Gross primary productivity
        Rpm::FT # Plant maintenance respiration
        ) where {FT}

Computes plant growth respiration as a function of gross primary productivity (GPP),
plant maintenance respiration (Rpm), and a relative contribution factor, f.
"""
function plant_respiration_growth(
    f2::FT, # Factor of relative contribution, usually 0.25
    GPP::FT, # Gross primary productivity
    Rpm::FT, # Plant maintenance respiration
) where {FT}
    Rg = f2 * (GPP - Rpm)
    return Rg
end
