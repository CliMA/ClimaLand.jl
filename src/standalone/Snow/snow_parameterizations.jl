export specific_heat_capacity,
    snow_thermal_conductivity,
    snow_bulk_temperature,
    liquid_mass_fraction,
    maximum_liquid_mass_fraction,
    runoff_timescale,
    compute_water_runoff,
    energy_from_q_l_and_swe,
    energy_from_T_and_swe,
    update_snow_cover_fraction!,
    phase_change_flux,
    update_snow_albedo!,
    energy_flux_falling_snow,
    energy_flux_falling_rain

"""
    update_snow_albedo!(α, m::ConstantAlbedoModel, Y, p, t, earth_param_set)

Updates the snow albedo `α` in place with the current albedo,
according to the ConstantAlbedoModel.
"""
function update_snow_albedo!(
    α,
    m::ConstantAlbedoModel,
    Y,
    p,
    t,
    earth_param_set,
)
    @. α = m.α
end

"""
    update_snow_albedo!(α, m::ZenithAngleAlbedoModel, Y, p, t, earth_param_set)

Updates the snow albedo `α` in place with the current albedo,
according to the ZenithAngleAlbedoModel.
"""
function update_snow_albedo!(
    α,
    m::ZenithAngleAlbedoModel,
    Y,
    p,
    t,
    earth_param_set,
)
    FT = eltype(earth_param_set)
    _ρ_liq = LP.ρ_cloud_liq(earth_param_set)
    @. α =
        min(1 - m.β * (p.snow.ρ_snow / _ρ_liq - m.x0), 1) *
        (m.α_0 + m.Δα * exp(-m.k * max(p.drivers.cosθs, eps(FT))))
    @. α = max(min(α, 1), 0)
end

"""
    update_snow_cover_fraction!(x::FT; z0 = FT(1e-1), β_scf = FT(2))::FT where {FT}

Returns the snow cover fraction from snow depth `z`, from
Wu, Tongwen, and Guoxiong Wu. "An empirical formula to compute
snow cover fraction in GCMs." Advances in Atmospheric Sciences
21 (2004): 529-535.
"""
function update_snow_cover_fraction!(
    scf,
    m::WuWuSnowCoverFractionModel,
    Y,
    p,
    t,
    earth_param_set,
)
    z = p.snow.z_snow
    @. scf = min(m.β_scf * (z / m.z0) / (z / m.z0 + 1), 1)
end

"""
    ClimaLand.surface_height(
        model::SnowModel{FT},
        Y,
        p,
    ) where {FT}

Returns the surface height of the `Snow` model; surface (ground) elevation
and ignores snow depth (CLM).

Once topography or land ice is incorporated, this will need to change to
z_sfc + land_ice_depth. Note that land ice can
be ~1-3 km thick on Greenland/

In order to compute surface fluxes, this cannot be larger than the
height of the atmosphere measurement location (z_atmos > z_land always).

This assumes that the surface elevation is zero.
"""
function ClimaLand.surface_height(model::SnowModel{FT}, Y, p) where {FT}
    return FT(0)
end

"""
    surface_albedo(model::SnowModel, Y, p)

A helper function which computes and returns the snow albedo.
"""
function ClimaLand.surface_albedo(model::SnowModel, Y, p)
    return p.snow.α_snow
end

"""
    surface_emissivity(model::SnowModel, Y, p)

A helper function which computes and returns the snow emissivity.
"""
function ClimaLand.surface_emissivity(model::SnowModel, _...)
    return model.parameters.ϵ_snow
end

"""
    ClimaLand.component_temperature(model::SnowModel, Y, p)

a helper function which returns the surface temperature for the snow
model, which is stored in the aux state.
"""
function ClimaLand.component_temperature(model::SnowModel, Y, p)
    return p.snow.T_sfc
end

"""
    ClimaLand.component_specific_humidity(model::SnowModel, Y, p)

Returns the precomputed specific humidity over snow as a weighted
fraction of the saturated specific humidity over liquid and frozen
water.

This uses the atmospheric T, P, q from p.drivers.
"""
function ClimaLand.component_specific_humidity(model::SnowModel, Y, p)
    h_sfc = ClimaLand.surface_height(model, Y, p)
    h_air = model.boundary_conditions.atmos.h
    @. p.snow.q_sfc = snow_surface_specific_humidity(
        p.snow.T_sfc,
        p.snow.q_l,
        p.drivers.T,
        p.drivers.P,
        p.drivers.q,
        h_air - h_sfc,
        model.parameters,
    )

    return p.snow.q_sfc
end

"""
    ClimaLand.surface_roughness_model(model::SnowModel, Y, p)

a helper function which returns the surface roughness model for the snow
model.
"""
function ClimaLand.surface_roughness_model(
    model::SnowModel{FT},
    Y,
    p,
) where {FT}
    return SurfaceFluxes.ConstantRoughnessParams{FT}(
        model.parameters.z_0m,
        model.parameters.z_0b,
    )
end

"""
    snow_surface_specific_humidity(T_sfc::FT, q_l::FT, T_air::FT, P_air::FT, q_air::FT, Δz::FT, parameters) where {FT}

Computes the snow surface specific humidity at a point, assuming a weighted averaged (by mass fraction)
of the saturated specific humidity over ice and over liquid, at temperature T_sfc.
"""
function snow_surface_specific_humidity(
    T_sfc::FT,
    q_l::FT,
    T_air::FT,
    P_air::FT,
    q_air::FT,
    Δz::FT,
    parameters,
) where {FT}
    surface_flux_params =
        LP.surface_fluxes_parameters(parameters.earth_param_set)
    thermo_params = LP.thermodynamic_parameters(parameters.earth_param_set)
    ρ_sfc = ClimaLand.compute_ρ_sfc(
        surface_flux_params,
        T_air,
        P_air,
        q_air,
        Δz,
        T_sfc,
    )
    qsat_over_ice = Thermodynamics.q_vap_saturation(
        thermo_params,
        T_sfc,
        ρ_sfc,
        Thermodynamics.Ice(),
    )
    qsat_over_liq = Thermodynamics.q_vap_saturation(
        thermo_params,
        T_sfc,
        ρ_sfc,
        Thermodynamics.Liquid(),
    )
    return qsat_over_ice * (1 - q_l) + q_l * (qsat_over_liq)
end

"""
    specific_heat_capacity(q_l::FT,
                           parameters::SnowParameters{FT}
                           ) where {FT}

Computes the specific heat capacity of the snow, neglecting
any contribution from air in the pore spaces, given
the liquid water mass fraction q_l and other parameters.
"""
function specific_heat_capacity(
    q_l::FT,
    parameters::SnowParameters{FT},
) where {FT}
    q_i = 1 - q_l
    _cp_i = FT(LP.cp_i(parameters.earth_param_set))

    _cp_l = FT(LP.cp_l(parameters.earth_param_set))

    cp_s = q_l * _cp_l + q_i * _cp_i
    return cp_s
end

"""
    snow_thermal_conductivity(ρ_snow::FT,
                         parameters::SnowParameters{FT},
                         ) where {FT}

Computes the thermal conductivity, given the density
of the snow, according to Equation
5.33 from Bonan's textbook, which in turn is taken from
Jordan (1991).

We have adjusted the original equation to make the coefficients
non-dimensional by multiplying by the first by x = ρ\\_ice/ρ\\_ice
and the second by x², with ρ\\_ice in kg/m³.

When ρ\\_snow = ρ\\_ice, we recover κ\\_snow = κ\\_ice.
"""
function snow_thermal_conductivity(
    ρ_snow::FT,
    parameters::SnowParameters{FT},
) where {FT}
    _κ_air = FT(LP.K_therm(parameters.earth_param_set))
    _ρ_ice = FT(LP.ρ_cloud_ice(parameters.earth_param_set))
    κ_ice = parameters.κ_ice
    return _κ_air +
           (FT(0.07) * (ρ_snow / _ρ_ice) + FT(0.93) * (ρ_snow / _ρ_ice)^2) *
           (κ_ice - _κ_air)
end

"""
    diurnal_damping_depth(κ::FT, ρ::FT, parameters::SnowParameters{FT})::FT where {FT}

Provides the characteristic depth for diurnal variations in temperature in snow, used in the
calculation of the surface tempearture. Formula is replicated from the Utah Energy Balance (UEB) model.
A paper describing the formula can be found at https://hess.copernicus.org/articles/18/5061/2014/.
"""
function diurnal_damping_depth(
    κ::FT,
    ρ::FT,
    parameters::SnowParameters{FT},
)::FT where {FT}
    _cp_i = FT(LP.cp_i(parameters.earth_param_set))
    _DT_ = FT(parameters.earth_param_set.insol_params.day)
    return sqrt(κ * _DT_ / (pi * _cp_i * ρ))
end

"""
    surface_temp_scaling_depth(κ::FT, ρ::FT, parameters::SnowParameters{FT})::FT where {FT}

Provides the scaling length scale for calculation of the surface temperature, which is stable
in the small-snowpack limit.
"""
function surface_temp_scaling_length(
    κ::FT,
    ρ::FT,
    z::FT,
    parameters::SnowParameters{FT},
)::FT where {FT}
    d0 = diurnal_damping_depth(κ, ρ, parameters)
    return d0 * (1 - exp(-z / d0))
end

"""
    liquid_mass_fraction(S::FT, S_l::FT)::FT where {FT}

Diagnoses the liquid mass fraction from the prognostic variables S_l and S,
while preventing from division by zero and enforcing that q_l = S_l/S must
lie between 0 and 1.
"""
function liquid_mass_fraction(S::FT, S_l::FT)::FT where {FT}
    S_safe = max(S, eps(FT))
    S_l_safe = min(max(S_l, FT(0)), S_safe)
    return S_l_safe / S_safe
end

"""
    snow_bulk_temperature(U::FT,
                          S::FT,
                          q_l::FT,
                          parameters::SnowParameters{FT}) where {FT}

Computes the bulk snow temperature from the snow water equivalent S,
energy per unit area U, liquid water fraction q_l, and specific heat
capacity c_s, along with other needed parameters.

If there is no snow (U = S = q_l = 0), the bulk temperature is the reference temperature,
which is 273.16K.
"""
function snow_bulk_temperature(
    U::FT,
    S::FT,
    q_l::FT,
    parameters::SnowParameters{FT},
) where {FT}
    S_safe = max(S, FT(0))
    _ρ_l = FT(LP.ρ_cloud_liq(parameters.earth_param_set))
    _T_ref = FT(LP.T_0(parameters.earth_param_set))
    _LH_f0 = FT(LP.LH_f0(parameters.earth_param_set))
    cp_s = specific_heat_capacity(q_l, parameters)
    _ΔS = parameters.ΔS
    _T_freeze = FT(LP.T_freeze(parameters.earth_param_set))
    return _T_ref +
           (U + _ρ_l * _LH_f0 * S_safe * (1 - q_l)) /
           (_ρ_l * cp_s * (S_safe + _ΔS))
end

"""
    snow_bulk_density(SWE::FT, z::FT, parameters::SnowParameters{FT}) where {FT}

Returns the snow density given the current model state when depth and SWE are available.
"""
function snow_bulk_density(
    SWE::FT,
    z::FT,
    parameters::SnowParameters{FT},
)::FT where {FT}
    _ρ_l = LP.ρ_cloud_liq(parameters.earth_param_set)
    ε = eps(FT) #for preventing dividing by zero
    #return SWE/z * ρ_l but tend to ρ_l as SWE → 0
    #also handle instabilities when z, SWE both near machine precision
    return max(SWE, ε) / max(z, SWE, ε) * _ρ_l
end

"""
    maximum_liquid_mass_fraction(ρ_snow::FT, T::FT, parameters::SnowParameters{FT}) where {FT}

Computes the maximum liquid water mass fraction, given
the density of the snow ρ_snow and other parameters.
"""
function maximum_liquid_mass_fraction(
    ρ_snow::FT,
    T::FT,
    parameters::SnowParameters{FT},
) where {FT}
    _ρ_l = LP.ρ_cloud_liq(parameters.earth_param_set)
    _T_freeze = LP.T_freeze(parameters.earth_param_set)
    if T > _T_freeze
        return FT(0)
    else
        parameters.θ_r * _ρ_l / ρ_snow
    end
end

"""
    runoff_timescale(z::FT, Ksat::FT, Δt::FT) where {FT}

Computes the timescale for liquid water to percolate and leave the snowpack,
given the depth of the snowpack z and the hydraulic conductivity Ksat.
"""
function runoff_timescale(z::FT, Ksat::FT, Δt::FT) where {FT}
    τ = max(Δt, z / Ksat)
    return τ
end

"""
    volumetric_internal_energy_liq(T, parameters)

Computes the volumetric internal energy of the liquid water
in the snowpack at a point at temperature T.
"""
function volumetric_internal_energy_liq(T, parameters)
    _ρ_l = LP.ρ_cloud_liq(parameters.earth_param_set)
    _T_freeze = LP.T_freeze(parameters.earth_param_set)
    _cp_l = LP.cp_l(parameters.earth_param_set)
    _T_ref = LP.T_0(parameters.earth_param_set)

    I_liq = _ρ_l * _cp_l * (T .- _T_ref)
    return I_liq
end

"""
    compute_energy_runoff(S::FT, S_l::FT, T::FT, parameters) where {FT}

Computes the rate of change in the snow water equivalent S due to loss of
liquid water (runoff) from the snowpack.

Runoff occurs as the snow melts and exceeds the water holding capacity.
"""
function compute_water_runoff(
    S::FT,
    S_l::FT,
    T::FT,
    ρ_snow::FT,
    z::FT,
    parameters,
) where {FT}
    τ = runoff_timescale(z, parameters.Ksat, parameters.Δt)
    q_l_max::FT = maximum_liquid_mass_fraction(ρ_snow, T, parameters)
    S_safe = max(S, FT(0))
    return -(S_l - q_l_max * S_safe) / τ * heaviside(S_l - q_l_max * S_safe)
end

"""
     phase_change_flux(U::FT, S::FT, q_l::FT, energy_flux::FT, parameters) where {FT}

Computes the volume flux of liquid water undergoing phase change, given the
applied energy flux and current state of U,S,q_l.
"""
function phase_change_flux(
    U::FT,
    S::FT,
    q_l::FT,
    energy_flux::FT,
    parameters,
) where {FT}
    S_safe = max(S, FT(0))

    energy_at_T_freeze = energy_from_q_l_and_swe(S_safe, q_l, parameters)
    Upred = U - energy_flux * parameters.Δt
    energy_excess = Upred - energy_at_T_freeze

    _LH_f0 = LP.LH_f0(parameters.earth_param_set)
    _ρ_liq = LP.ρ_cloud_liq(parameters.earth_param_set)
    _cp_i = LP.cp_i(parameters.earth_param_set)
    _cp_l = LP.cp_l(parameters.earth_param_set)
    _T_ref = LP.T_0(parameters.earth_param_set)
    _T_freeze = LP.T_freeze(parameters.earth_param_set)
    if energy_excess > 0 || (energy_excess < 0 && q_l > 0)
        return -energy_excess / parameters.Δt / _ρ_liq /
               ((_cp_l - _cp_i) * (_T_freeze - _T_ref) + _LH_f0)
    else
        return FT(0)
    end
end

"""
    energy_from_q_l_and_swe(S::FT, q_l::FT, parameters) where {FT}

A helper function for compute the snow energy per unit area, given snow
water equivalent S, liquid fraction q_l, and snow model parameters.

This assumes that the snow is at the freezing point.
"""
function energy_from_q_l_and_swe(S::FT, q_l::FT, parameters) where {FT}
    _T_freeze = FT(LP.T_freeze(parameters.earth_param_set))
    _ρ_l = FT(LP.ρ_cloud_liq(parameters.earth_param_set))
    _T_ref = FT(LP.T_0(parameters.earth_param_set))
    _LH_f0 = FT(LP.LH_f0(parameters.earth_param_set))
    _ΔS = parameters.ΔS

    c_snow = specific_heat_capacity(q_l, parameters)
    return _ρ_l * (S + _ΔS) * c_snow * (_T_freeze - _T_ref) -
           _ρ_l * S * (1 - q_l) * _LH_f0
end

"""
    energy_from_T_and_swe(S::FT, T::FT, parameters) where {FT}

A helper function for compute the snow energy per unit area, given snow
water equivalent S, bulk temperature T, and snow model parameters.

The liquid mass fraction is assumed to be zero if T<=T_freeze, and 1 otherwise.
"""
function energy_from_T_and_swe(S::FT, T::FT, parameters) where {FT}
    _ρ_l = FT(LP.ρ_cloud_liq(parameters.earth_param_set))
    _T_ref = FT(LP.T_0(parameters.earth_param_set))
    _T_freeze = FT(LP.T_freeze(parameters.earth_param_set))
    _LH_f0 = FT(LP.LH_f0(parameters.earth_param_set))
    _cp_i = FT(LP.cp_i(parameters.earth_param_set))
    _cp_l = FT(LP.cp_l(parameters.earth_param_set))
    _ΔS = parameters.ΔS

    if T <= _T_freeze
        return _ρ_l * _cp_i * (S + _ΔS) * (T - _T_ref) - _ρ_l * S * _LH_f0
    else
        return _ρ_l * (S + _ΔS) * _cp_l * (T - _T_ref)
    end

end

"""
    energy_flux_falling_snow(atmos, p, parameters)

Returns the energy flux of falling snow for a PrescribedAtmosphere,
approximated as ρe_snow * P_snow, where ρe_snow = -LH_f0 * _ρ_liq.
This is a negative internal energy, due the to negative contribution of
the latent heat of melting to the energy of the snow,
 and it neglects the sensible heat portion of the snow. The energy
 is per unit volume of liquid water, and P_snow is expressed as
the volume flux of liquid water resulting from the snow.

This method can be extended to coupled simulations, where atmos is of type
CoupledAtmosphere, and the energy flux of the falling snow is passed in the
cache `p`. In that case, this should specify `atmos::PrescribedAtmosphere`.
"""
function energy_flux_falling_snow(atmos, p, parameters)
    _LH_f0 = LP.LH_f0(parameters.earth_param_set)
    _ρ_liq = LP.ρ_cloud_liq(parameters.earth_param_set)
    ρe_snow = -_LH_f0 * _ρ_liq
    return @. lazy(ρe_snow * p.drivers.P_snow)
end

"""
    energy_flux_falling_rain(atmos, p, parameters)

Returns the energy flux of falling rain for a PrescribedAtmosphere,
approximated as ρ_l e_l(T_atmos) * P_liq. The energy is per unit volume of liquid water,
and P_liq is expressed as the volume flux of liquid water resulting from the rain.

This method can be extended to coupled simulations, where atmos is of type
CoupledAtmosphere, and the energy flux of the falling rain is passed in the
cache `p`.  In that case, this should specify `atmos::PrescribedAtmosphere`.
"""
function energy_flux_falling_rain(atmos, p, parameters)
    return @. lazy(
        volumetric_internal_energy_liq(p.drivers.T, parameters) *
        p.drivers.P_liq,
    )
end

"""
    update_density_and_depth!(ρ_snow, z_snow, density::MinimumDensityModel, Y, p, params::SnowParameters)

Extends the update_density_and_depth! function for the MinimumDensityModel type; updates the snow density and depth in place.
"""
function update_density_and_depth!(
    ρ_snow,
    z_snow,
    density::MinimumDensityModel{FT},
    Y,
    p,
    params::SnowParameters{FT},
) where {FT}
    _ρ_l = LP.ρ_cloud_liq(params.earth_param_set)
    @. ρ_snow = density.ρ_min * (1 - p.snow.q_l) + _ρ_l * p.snow.q_l
    @. z_snow = _ρ_l * Y.snow.S / ρ_snow
end

"""
    update_density_prog!(density::MinimumDensityModel, model::SnowModel, Y, p)

Updates all prognostic variables associated with density/depth given the current model state.
This is the default method for the constant minimum density model,
 which has no prognostic variables.
"""
function update_density_prog!(
    density::MinimumDensityModel,
    model::SnowModel,
    dY,
    Y,
    p,
)
    return nothing
end

"""
    flux_balance(
        T_sfc_guess::FT,
        T_bulk::FT,
        z::FT,
        κ::FT,
        ρ_snow::FT,
        ϵ_snow::FT,
        SW_net::FT,
        LW_d::FT,
        q_l::FT,
        h_sfc::FT,
        displ::FT,
        P_atmos::FT,
        T_atmos::FT,
        q_atmos::FT,
        u_atmos,
        roughness_model,
        atmos_h::FT,
        gustiness::FT,
        parameters::SnowParameters{FT},
    )

Returns the balance (difference) in surface energy fluxes between a conductive flux for
a proposed surface temperature, as well as input forcing surface fluxes. Evaluated point-wise
to enable the broadcasting of the numerical root-solving algorithms over the resulting ClimaFields.
"""
function flux_balance(
    T_sfc_guess::FT,
    T_bulk::FT,
    z::FT,
    κ::FT,
    ρ_snow::FT,
    ϵ_snow::FT,
    SW_net::FT,
    LW_d::FT,
    q_l::FT,
    h_sfc::FT,
    displ::FT,
    P_atmos::FT,
    T_atmos::FT,
    q_atmos::FT,
    u_atmos,
    roughness_model,
    atmos_h::FT,
    gustiness,
    parameters::SnowParameters{FT},
)::FT where {FT}

    # For some snowpack states, the SecantMethod can cause a nonpositive T_sfc to
    # be attempted as part of the iteration trajectory - however, the algorithm
    # can sometimes recover and still converge to the right value with the branch
    # below. Negative T_sfc cause domain errors in SurfaceFluxes, which would prevent such recovery.
    if T_sfc_guess <= eps(FT)
        return κ * (T_sfc_guess - T_bulk)
    end

    _σ = LP.Stefan(parameters.earth_param_set)
    LW_net = -ϵ_snow * (LW_d - _σ * T_sfc_guess^4) #match sign convention in ./shared_utilities/drivers.jl
    d = surface_temp_scaling_length(κ, ρ_snow, z, parameters)

    q_sfc = snow_surface_specific_humidity(
        T_sfc_guess,
        q_l,
        T_atmos,
        P_atmos,
        q_atmos,
        atmos_h - h_sfc,
        parameters,
    )

    #we only need the lhf, shf, so the derivative functions
    #can be specified as a simple computation:
    blank_deriv(args...) = FT(1)

    latent_flux, sensible_flux, _, _, _, _, _ =
        ClimaLand.compute_turbulent_fluxes_at_a_point(
            P_atmos,
            T_atmos,
            q_atmos,
            u_atmos,
            atmos_h,
            T_sfc_guess,
            q_sfc,
            roughness_model,
            nothing,
            nothing,
            h_sfc,
            displ,
            blank_deriv,
            blank_deriv,
            gustiness,
            parameters.earth_param_set,
        )

    #For numerical stability, the scaling length of conductive flux
    #is multiplied on the LHS, instead of dividing the temperature difference:
    #(This function output is not needed in flux units, at present)
    return FT(
        d * (SW_net + LW_net + latent_flux + sensible_flux) +
        κ * (T_sfc_guess - T_bulk),
    )
end

"""
    solve_for_surface_temp_at_a_point(
        T_bulk::FT,
        z_snow::FT,
        κ_snow::FT,
        ρ_snow::FT,
        ϵ_snow::FT,
        SW_net::FT,
        LW_d::FT,
        q_l::FT,
        h_sfc::FT,
        displ::FT,
        P_atmos::FT,
        T_atmos::FT,
        q_atmos::FT,
        u_atmos,
        roughness_model,
        atmos_h::FT,
        gustiness,
        parameters::SnowParameters{FT},
    )

Determines the surface temperature from the surface energy balance, at a point location.
Makes use of a root-solving algorithm (SecantMethod), which is broadcasted over the input ClimaFields.
"""
function solve_for_surface_temp_at_a_point(
    T_bulk::FT,
    z_snow::FT,
    κ_snow::FT,
    ρ_snow::FT,
    ϵ_snow::FT,
    SW_net::FT,
    LW_d::FT,
    q_l::FT,
    h_sfc::FT,
    displ::FT,
    P_atmos::FT,
    T_atmos::FT,
    q_atmos::FT,
    u_atmos,
    roughness_model,
    atmos_h::FT,
    gustiness,
    parameters::SnowParameters{FT},
)::FT where {FT}

    flux_balance_closure(T_sfc_guess::FT)::FT = flux_balance(
        T_sfc_guess,
        T_bulk,
        z_snow,
        κ_snow,
        ρ_snow,
        ϵ_snow,
        SW_net,
        LW_d,
        q_l,
        h_sfc,
        displ,
        P_atmos,
        T_atmos,
        q_atmos,
        u_atmos,
        roughness_model,
        atmos_h,
        gustiness,
        parameters,
    )

    tol = parameters.surf_temp.tol
    N_iters = parameters.surf_temp.N_iters

    #Starting points picked to minimize number of failed convergences:
    method = SecantMethod(T_atmos, T_atmos - FT(1))

    sol = RootSolvers.find_zero(
        T -> flux_balance_closure(T),
        method,
        CompactSolution(),
        ResidualTolerance(tol),
        N_iters,
    )

    return ifelse(sol.converged, sol.root, T_bulk)
end

"""
    surface_residual_flux(
        T_sfc_root::FT,
        κ::FT,
        ρ::FT,
        z::FT,
        parameters::SnowParameters{FT}
    )

Determines the residual surface energy flux to dump into the bulk energy whenever the correct surface
temperature to satisfy Fourier's Law of Conduction at the surface is greater than T_freeze.
"""
function surface_residual_flux(
    T_sfc_root::FT,
    κ::FT,
    ρ::FT,
    z::FT,
    parameters::SnowParameters{FT},
)::FT where {FT}
    _T_freeze = FT(LP.T_freeze(parameters.earth_param_set))
    d_safe = max(surface_temp_scaling_length(κ, ρ, z, parameters), eps(FT))
    return T_sfc_root > _T_freeze ? FT(-κ * (T_sfc_root - _T_freeze) / d_safe) :
           FT(0)
end

"""
    update_surf_temp!(model::SnowModel, surf_temp::EquilibriumGradientTemperatureModel, SW_net, LW_down, Y, p, t)

Updates the surface temperature variable, storing residual energy flux in p.snow.surf_residual_flux to be added to the bulk.
"""
function update_surf_temp!(
    model::SnowModel,
    surf_temp::EquilibriumGradientTemperatureModel,
    SW_net,
    LW_d,
    Y,
    p,
    t,
)
    bc = model.boundary_conditions
    _T_freeze = LP.T_freeze(model.parameters.earth_param_set)

    #If roughness length is updated during the Monin-Obukhov iterations, we might need to move this inside of the root solve.
    h_sfc = ClimaLand.surface_height(model, Y, p)
    roughness_model = ClimaLand.surface_roughness_model(model, Y, p)
    displ = ClimaLand.surface_displacement_height(model, Y, p)
    #might need to update this call as gustiness models change:
    gustiness = SurfaceFluxes.ConstantGustinessSpec(bc.atmos.gustiness)

    #get surf_temp values, even if they are > T_freeze:
    p.snow.T_sfc .=
        solve_for_surface_temp_at_a_point.(
            p.snow.T,
            p.snow.z_snow,
            p.snow.κ,
            p.snow.ρ_snow,
            model.parameters.ϵ_snow,
            SW_net, #passing radiation as args allows different radiation situations, i.e. open ground vs canopy-shielded
            LW_d, #passing radiation as args enables different radiation situations, i.e. open ground vs canopy-shielded
            p.snow.q_l,
            h_sfc,
            displ,
            p.drivers.P,
            p.drivers.T,
            p.drivers.q,
            p.drivers.u,
            roughness_model,
            bc.atmos.h,
            gustiness,
            model.parameters,
        )

    #set residual flux using values if T_sfc > T_freeze:
    p.snow.surf_residual_flux .=
        surface_residual_flux.(
            p.snow.T_sfc,
            p.snow.κ,
            p.snow.ρ_snow,
            p.snow.z_snow,
            model.parameters,
        )

    #reset T_sfc accordingly:
    p.snow.T_sfc .= min.(_T_freeze, p.snow.T_sfc)
    return nothing
end


"""
    update_surf_temp!(model::SnowModel, surf_temp::BulkSurfaceTemperatureModel, Y, p, t)

Updates the surface temperature variable so that it matches the bulk temperature.
Note, this update function is not called for the integrated land model - the update to the surface temperature happens via a different
function that handles the radiation sent from the canopy (see `lsm_radiant_energy_fluxes!()`)
"""
function update_surf_temp!(
    model::SnowModel,
    surf_temp::BulkSurfaceTemperatureModel,
    SW_net,
    LW_d,
    Y,
    p,
    t,
)
    p.snow.T_sfc .= p.snow.T
    return nothing
end

"""
    get_residual_surface_flux(surf_temp::BulkSurfaceTemperatureModel, Y, p)

Returns any residual surface flux as defined by the surface temperature parameterization choice.
"""
function get_residual_surface_flux(surf_temp::BulkSurfaceTemperatureModel, Y, p)
    return FTfromY(Y)(0)
end

"""
    get_residual_surface_flux(surf_temp::EquilibriumGradientTemperatureModel, Y, p)

Returns any residual surface flux as defined by the surface temperature parameterization choice.
"""
function get_residual_surface_flux(
    surf_temp::EquilibriumGradientTemperatureModel,
    Y,
    p,
)
    return p.snow.surf_residual_flux
end
