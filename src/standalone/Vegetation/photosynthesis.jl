abstract type AbstractPhotosynthesisModel{FT} <: AbstractCanopyComponent{FT} end
export get_Vcmax25_leaf,
    get_Rd_leaf,
    get_An_leaf,
    get_J_over_Jmax,
    get_GPP_canopy,
    rubisco_assimilation,
    light_assimilation,
    net_photosynthesis,
    gross_photosynthesis,
    arrhenius_function
ClimaLand.name(model::AbstractPhotosynthesisModel) = :photosynthesis
# Photosynthesis models must add auxiliary variables and a method for  
# update_photosynthesis!(p,Y,photosynthesis_model, canopy_model)
# which updatesthose auxiliary variables in place.
# Downstream SIF, autotrophic respiration, and diagnostics require Vcmax, J/Jmax, and Rd.
ClimaLand.auxiliary_vars(model::AbstractPhotosynthesisModel) = ()
ClimaLand.auxiliary_types(model::AbstractPhotosynthesisModel{FT}) where {FT} =
    ()
ClimaLand.auxiliary_domain_names(::AbstractPhotosynthesisModel) = ()
function update_photosynthesis!(
    p,
    Y,
    photosynthesis_model::AbstractPhotosynthesisModel,
    canopy_model,
)
    nothing
end
get_GPP_canopy(p, photosynthesis_model::AbstractPhotosynthesisModel) = nothing # used by diagnostics
get_Vcmax_leaf(p, photosynthesis_model::AbstractPhotosynthesisModel) = nothing # used by solar induced fluorescence, diagnostics, autotrophic respiration
get_J_over_Jmax(
    Y,
    p,
    canopy_model,
    photosynthesis_model::AbstractPhotosynthesisModel,
) = nothing # used by solar induced fluorescence
get_Rd_leaf(p, photosynthesis_model::AbstractPhotosynthesisModel) = nothing # used with autotrophic respiration and diagnostics
get_An_leaf(p, photosynthesis_model::AbstractPhotosynthesisModel) = nothing # used with Medlyn conductance, autotrophic respiration and diagnostics
# Farquhar photosynthesis computations which are common to multiple Photosynthesis models

"""
    net_photosynthesis(A::FT,
                       Rd::FT) where {FT}

Computes the total net carbon assimilation (`An`),
in units of mol CO2/m^2/s, as a function of
the gross assimilation (`A`) and 
dark respiration (`Rd`).

Note that if  moisture stress is factored into the model, it must be applied correctly to Rd
and A prior to this function.

Both A and Rd must either be leaf level (standard Farquhar), or 
canopy level (p-model), but not mixed. The returned quantity is 
similarily either leaf or canopy level.

See Table 11.5 of G. Bonan's textbook, Climate Change and Terrestrial Ecosystem Modeling (2019).
"""
function net_photosynthesis(A::FT, Rd::FT) where {FT}
    An = max(0, A - Rd)
    return An
end


"""
    gross_photosynthesis(Ac::FT,
                         Aj::FT) where {FT}

Computes the gross carbon assimilation (`A`),
in units of mol CO2/m^2/s, as a function of
the Rubisco limiting factor (`Ac`), the electron transport limiting rate (`Aj`).

Note that if  moisture stress is factored into the model, it must be applied correctly to Ac
and Aj prior to this function.

Both Ac and Aj must either be leaf level (standard Farquhar), or 
canopy level (p-model), but not mixed. The returned quantity is 
similarily either leaf or canopy level.

See Table 11.5 of G. Bonan's textbook, Climate Change and Terrestrial Ecosystem Modeling (2019).
"""
function gross_photosynthesis(Ac::FT, Aj::FT) where {FT}
    return min(Ac, Aj)
end

"""
    rubisco_assimilation(is_c3::AbstractFloat, args...)

Calls the correct rubisco assimilation function based on the `is_c3`;
the returned quantity may be canopy-level (p-model) or leaf-level (standard
Farquhar).

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
                         Kmm::FT) where {FT}

Computes the Rubisco limiting rate of photosynthesis for C3 plants (`Ac`),
in units of moles CO2/m^2/s,
as a function of the maximum rate of carboxylation of Rubisco (`Vcmax`),
the leaf internal carbon dioxide partial pressure (`ci`),
the CO2 compensation point (`Γstar`), and the effective Michaelis-Menten parameter
`K_mm`.

The units of ci, Γstar, and Kmm may be Pa or mol/mol, but they must be consistent. The
units of Vcmax must be mol CO2/m^2.s. Similarily, the assimilation rate
may be leaf-level (standard Farquhar) or canopy level (p-model), 
depending on the units of Vcmax.

See Table 11.5 of G. Bonan's textbook, Climate Change and Terrestrial Ecosystem Modeling (2019).
"""
function c3_rubisco_assimilation(
    Vcmax::FT,
    ci::FT,
    Γstar::FT,
    Kmm::FT,
) where {FT}
    Ac = Vcmax * (ci - Γstar) / (ci + Kmm)
    return Ac
end

"""
    c4_rubisco_assimilation(Vcmax::FT,_...) where {FT}

Computes the Rubisco limiting rate of photosynthesis for C4 plants (`Ac`)
in units of moles CO2/m^2/s,
as equal to the maximum rate of carboxylation of Rubisco (`Vcmax`).

The assimilation rate
may be leaf-level (standard Farquhar) or canopy level (p-model), 
depending on the units of Vcmax.
"""
function c4_rubisco_assimilation(Vcmax::FT, _...) where {FT}
    Ac = Vcmax
    return Ac
end

"""
    light_assimilation(is_c3::AbstractFloat, args...)

Calls the correct light assimilation function based on the `is_c3`;
the returned quantity may be canopy-level (p-model) or leaf-level (standard
Farquhar).

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

The returned assimilation rate
may be leaf-level (standard Farquhar) or canopy level (p-model), 
depending on the units of J.

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

The returned assimilation rate
may be leaf-level (standard Farquhar) or canopy level (p-model), 
depending on the units of APAR (per unit leaf area, or per ground area).

For C4 plants, this is a function of APAR and a efficiency parameter E, see Equation 11.70 of
 G. Bonan's textbook, Climate Change and Terrestrial Ecosystem Modeling (2019).
"""
function c4_light_assimilation(::FT, ::FT, ::FT, APAR::FT, E::FT) where {FT}
    Aj = APAR * E
    return Aj
end


"""
    MM_Kc(Kc25::FT,
          ΔHkc::FT,
          T::FT,
          To::FT,
          R::FT) where {FT}

Computes the Michaelis-Menten coefficient for CO2 (`Kc`),
in units of mol/mol or Pa (depending on the units of `Kc25`),
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
in units of mol/mol or Pa (depending on the units of `Ko25`),
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
    compute_APAR_canopy_moles(
        f_abs::FT,
        par_d::FT,
        λ_γ_PAR::FT,
        lightspeed::FT,
        planck_h::FT,
        N_a::FT
    ) where {FT}

Computes the absorbed photosynthetically active radiation for the canopy in mol photons m^-2 s^-1,
given the fraction of absorbed PAR (`f_abs`), the PAR downwelling flux (`par_d`, in W m^-2),
and the wavelength of PAR (`λ_γ_PAR`, in m), and the physical constants necessary to compute
the energy per mol PAR photons.
"""
function compute_APAR_canopy_moles(
    f_abs::FT,
    par_d::FT,
    λ_γ_PAR::FT,
    lightspeed::FT,
    planck_h::FT,
    N_a::FT,
) where {FT}
    energy_per_mole_photon_par = planck_h * lightspeed * N_a / λ_γ_PAR
    return f_abs * par_d / energy_per_mole_photon_par
end

"""
    compute_APAR_leaf_moles(
        f_abs::FT,
        par_d::FT,
        λ_γ_PAR::FT,
        lightspeed::FT,
        planck_h::FT,
        N_a::FT,
        LAI::FT
    ) where {FT}

Computes the absorbed photosynthetically active radiation for the leaf in mol photons m^-2 s^-1,
given the fraction of absorbed PAR (`f_abs`), the PAR downwelling flux (`par_d`, in W m^-2),
and the wavelength of PAR (`λ_γ_PAR`, in m), the leaf area index `LAI`, and the physical constants necessary to compute
the energy per mol PAR photons.
"""
function compute_APAR_leaf_moles(
    f_abs::FT,
    par_d::FT,
    λ_γ_PAR::FT,
    lightspeed::FT,
    planck_h::FT,
    N_a::FT,
    LAI::FT,
) where {FT}
    APAR_canopy = compute_APAR_canopy_moles(
        f_abs,
        par_d,
        λ_γ_PAR,
        lightspeed,
        planck_h,
        N_a,
    )
    return APAR_canopy / max(LAI, sqrt(eps(FT)))
end
