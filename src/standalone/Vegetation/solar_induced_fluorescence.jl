export SIFParameters, Lee2015SIFModel

abstract type AbstractSIFModel{FT} <: AbstractCanopyComponent{FT} end

"""
    SIFParameters{FT<:AbstractFloat}

The required parameters for the SIF parameterisation 
Lee et al, 2015. Global Change Biology 21, 3469-3477, doi:10.1111/gcb.12948.
$(DocStringExtensions.FIELDS)
"""
@kwdef struct SIFParameters{FT <: AbstractFloat}
    "The rate coefficient for florescence, unitless"
    kf::FT
    "Parameter used to compute the rate coefficient for heat loss in dark-adapted conditions, Tol et al. 2014, unitless"
    kd_p1::FT
    "Parameter used to compute the rate coefficient for heat loss in dark-adapted conditions, Tol et al. 2014, unitless"
    kd_p2::FT
    "Parameter used to compute the rate coefficient for heat loss in dark-adapted conditions, Tol et al. 2014, unitless"
    min_kd::FT
    "Parameter used to compute the rate coefficient for heat loss in light-adapted conditions, Lee et al 2013 (unitless)"
    kn_p1::FT
    "Parameter used to compute the rate coefficient for heat loss in light-adapted conditions, Lee et al 2013 (unitless)"
    kn_p2::FT
    "Rate coefficient for photochemical quenching"
    kp::FT
    "Slope of line relating leaf-level fluorescence to spectrometer-observed fluorescence as a function of Vcmax 25 (leaf level). Lee et al 2015."
    kappa_p1::FT
    "Intercept of line relating leaf-level fluorescence to spectrometer-observed fluorescence as a function of Vcmax 25 (leaf level).  Lee et al 2015."
    kappa_p2::FT
end

Base.eltype(::SIFParameters{FT}) where {FT} = FT

"""
    Lee2015SIFModel{FT, SP <: SIFParameters{FT}} <: AbstractSIFModel{FT}

`Lee2015SIFModel` struct containing the parameters needed to compute
the SIF parameterization described in:
Lee et al, 2015. Global Change Biology 21, 3469-3477, doi:10.1111/gcb.12948.
"""
struct Lee2015SIFModel{FT, SP <: SIFParameters{FT}} <: AbstractSIFModel{FT}
    parameters::SP
end

ClimaLand.name(model::AbstractSIFModel) = :sif
ClimaLand.auxiliary_vars(model::Lee2015SIFModel) = (:SIF,)
ClimaLand.auxiliary_types(model::Lee2015SIFModel{FT}) where {FT} = (FT,)
ClimaLand.auxiliary_domain_names(::Lee2015SIFModel) = (:surface,)


# 4 Solar Induced Fluorescence (SIF)

# call function below inside photosynthesis.jl p

"""
    update_SIF!(p, Y, sif_model::Lee2015SIFModel, canopy)

Updates observed SIF at 755 nm in W/m^2. Note that Tc is in Kelvin, and photo
synthetic rates are in mol/m^2/s, and APAR is in PPFD.
Lee et al, 2015. Global Change Biology 21, 3469-3477, doi:10.1111/gcb.12948
"""
function update_SIF!(p, Y, sif_model::Lee2015SIFModel, canopy)
    SIF = p.canopy.sif.SIF
    earth_param_set = canopy.parameters.earth_param_set

    # Compute APAR
    f_abs_par = p.canopy.radiative_transfer.par.abs
    par_d = p.canopy.radiative_transfer.par_d
    (; λ_γ_PAR,) = canopy.radiative_transfer.parameters
    c = LP.light_speed(earth_param_set)
    planck_h = LP.planck_constant(earth_param_set)
    N_a = LP.avogadro_constant(earth_param_set)
    APAR_canopy_moles = @. lazy(
        compute_APAR_canopy_moles(f_abs_par, par_d, λ_γ_PAR, c, planck_h, N_a),
    )

    # Get max photosynthesis rates
    T_canopy = canopy_temperature(canopy.energy, canopy, Y, p)
    Vcmax25_leaf = get_Vcmax25_leaf(p, canopy.photosynthesis)
    J_over_Jmax = get_J_over_Jmax(Y, p, canopy, canopy.photosynthesis)

    T_freeze = LP.T_freeze(earth_param_set)
    sif_parameters = sif_model.parameters

    @. SIF = compute_SIF_at_a_point(
        APAR_canopy_moles,
        T_canopy,
        Vcmax25_leaf,
        J_over_Jmax,
        T_freeze,
        sif_parameters,
    )
end
function lazy_SIF(p, Y, sif_model::Lee2015SIFModel, canopy)
    SIF = p.canopy.sif.SIF
    earth_param_set = canopy.parameters.earth_param_set

    # Compute APAR
    f_abs_par = p.canopy.radiative_transfer.par.abs
    par_d = p.canopy.radiative_transfer.par_d
    (; λ_γ_PAR,) = canopy.radiative_transfer.parameters
    c = LP.light_speed(earth_param_set)
    planck_h = LP.planck_constant(earth_param_set)
    N_a = LP.avogadro_constant(earth_param_set)
    APAR_canopy_moles = @. lazy(
        compute_APAR_canopy_moles(f_abs_par, par_d, λ_γ_PAR, c, planck_h, N_a),
    )

    # Get max photosynthesis rates
    T_canopy = canopy_temperature(canopy.energy, canopy, Y, p)
    Vcmax25_leaf = get_Vcmax25_leaf(p, canopy.photosynthesis)
    J_over_Jmax = get_J_over_Jmax(Y, p, canopy, canopy.photosynthesis)

    T_freeze = LP.T_freeze(earth_param_set)
    sif_parameters = sif_model.parameters

    return @. lazy(compute_SIF_at_a_point(
        APAR_canopy_moles,
        T_canopy,
        Vcmax25_leaf,
        J_over_Jmax,
        T_freeze,
        sif_parameters,
    ))
end

Base.broadcastable(m::SIFParameters) = tuple(m)


"""
    compute_SIF_at_a_point(
        APAR_canopy_moles::FT,
        Tc::FT,
        Vcmax25_leaf::FT,
        J_over_Jmax::FT,
        T_freeze::FT,
        sif_parameters::SIFParameters{FT},
    ) where {FT}
    
Computes the Solar Induced Fluorescence (SIF) at 755 nm in W/m^2 using the Lee et al. 2015 model.
This takes as parameters `APAR` (absorbed photosynthetically active radiation, mol/m^2/s), `Tc`
(canopy temperature, K), `Vcmax25` (maximum carboxylation rate at 25 °C, mol/m^2/s), `Jmax`
(electron transport rate, mol/m^2/s), `J` (electron transport rate, mol/m^2/s), `T_freeze` 
(freezing temperature, K), `sif_parameters` (SIF parameters). 
"""
function compute_SIF_at_a_point(
    APAR_canopy_moles::FT,
    Tc::FT,
    Vcmax25_leaf::FT,
    J_over_Jmax::FT,
    T_freeze::FT,
    sif_parameters::SIFParameters{FT},
) where {FT}
    (; kf, kd_p1, kd_p2, min_kd, kn_p1, kn_p2, kp, kappa_p1, kappa_p2) =
        sif_parameters
    kd = max(kd_p1 * (Tc - T_freeze) + kd_p2, min_kd)
    x = 1 - J_over_Jmax
    kn = (kn_p1 * x - kn_p2) * x
    ϕp0 = kp / max(kf + kp + kn, eps(FT))
    ϕp = J_over_Jmax * ϕp0
    ϕf = kf / max(kf + kd + kn, eps(FT)) * (1 - ϕp)
    κ = kappa_p1 * Vcmax25_leaf * FT(1e6) + kappa_p2 # formula expects Vcmax25 in μmol/m^2/s
    F = APAR_canopy_moles * ϕf
    SIF_755 = F / max(κ, eps(FT))

    return SIF_755
end
