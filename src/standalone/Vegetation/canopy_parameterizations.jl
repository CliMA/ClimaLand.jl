using ..ClimaLand.Canopy
export canopy_sw_rt_beer_lambert, # Radiative transfer
    canopy_sw_rt_two_stream,
    extinction_coeff,
    compute_G,
    # Conductance
    medlyn_term,
    medlyn_conductance,
    penman_monteith,
    quadratic_soil_moisture_stress

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
        Thermodynamics.vapor_pressure_deficit(
            thermo_params,
            T_air,
            P_air,
            q_air,
        ),
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
