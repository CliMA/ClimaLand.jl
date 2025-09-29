# NOTE: A future to-do is to unite the pmodel and Farquhar functions where
# possible. For example, the rubisco and light assimilation functions
# do not need to be different, but we currently model the c4 light
# assimilation slightly different between the two. Additionally, we should
# use the same leaf-level representation throughout both, rather than have
# one be canopy level (Pmodel), and one leaf-level. Finally, we should
# use the same units throughout both.
abstract type AbstractPhotosynthesisModel{FT} <: AbstractCanopyComponent{FT} end
export get_Vcmax25_leaf,
    get_Rd_leaf,
    get_An_leaf,
    get_J_over_Jmax,
    get_GPP_canopy,
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
