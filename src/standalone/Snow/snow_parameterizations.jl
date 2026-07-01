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

Returns the snow cover fraction from ground-area snow depth `z`, from
Wu, Tongwen, and Guoxiong Wu. "An empirical formula to compute
snow cover fraction in GCMs." Advances in Atmospheric Sciences
21 (2004): 529-535.
"""
function update_snow_cover_fraction!(
    p,
    m::WuWuSnowCoverFractionModel,
    Y,
    t,
    earth_param_set,
    prognostic_land_components,
)
    scf = p.snow.snow_cover_fraction
    z = p.snow.z_snow #ground-area snow depth in this case
    @. scf = min(m.β_scf * (z / m.z0) / (z / m.z0 + 1), 1)
    maximum_snow_cover_fraction!(p, Val(prognostic_land_components))
end

"""
    maximum_snow_cover_fraction!(p, prognostic_land_components)

Default method that clips `p.snow.snow_cover_fraction` against any
upper bound implied by the set of prognostic land components; the
fallback does nothing.

Currently, this only acts when an inland water (lake) component is
present in an integrated `LandModel`, in which case the snow cover
fraction is restricted to the non-lake area.
"""
maximum_snow_cover_fraction!(p, plc_val) = nothing

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
    surface_flux_params = LP.surface_fluxes_parameters(model.parameters.earth_param_set)
    thermo_params = LP.thermodynamic_parameters(model.parameters.earth_param_set)

    @. p.snow.q_sfc = snow_surface_specific_humidity(
        p.snow.T_sfc,
        p.snow.q_l,
        p.drivers.T,
        p.drivers.P,
        p.drivers.q,
        h_air - h_sfc,
        surface_flux_params,
        thermo_params
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
    snow_surface_specific_humidity(T_sfc::FT, q_l::FT, T_air::FT, P_air::FT, q_air::FT, Δz::FT, surface_flux_params, thermo_params) where {FT}

Computes the snow surface specific humidity at a point, assuming a weighted averaged (by mass fraction)
of the saturated specific humidity over ice and over liquid, at temperature T_sfc.

Be aware that if this function changes you must also change the internals of `update_T_sfc_scheme`.
"""
function snow_surface_specific_humidity(
    T_sfc::FT,
    q_l::FT,
    T_air::FT,
    P_air::FT,
    q_air::FT,
    Δz::FT,
    surface_flux_params,
    thermo_params
) where {FT}
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
                           earth_param_set
                           ) where {FT}

Computes the specific heat capacity of the snow, neglecting
any contribution from air in the pore spaces, given
the liquid water mass fraction q_l and other parameters.
"""
function specific_heat_capacity(q_l::FT, earth_param_set) where {FT}
    q_i = 1 - q_l
    _cp_i = FT(LP.cp_i(earth_param_set))
    _cp_l = FT(LP.cp_l(earth_param_set))
    cp_s = q_l * _cp_l + q_i * _cp_i
    return cp_s
end

"""
    snow_thermal_conductivity(
        scheme::JordanSnowConductivityModel,
        ρ_snow::FT,
        earth_param_set,
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
    scheme::JordanSnowConductivityModel,
    ρ_snow::FT,
    earth_param_set,
) where {FT}
    κ_ice = scheme.κ_ice
    _κ_air = LP.K_therm(earth_param_set)
    _ρ_ice = LP.ρ_cloud_ice(earth_param_set)
    return _κ_air +
           (
        scheme.linear_coeff * (ρ_snow / _ρ_ice) +
        scheme.quadratic_coeff * (ρ_snow / _ρ_ice)^2
    ) * (κ_ice - _κ_air)
end

"""
    snow_thermal_conductivity(
        scheme::SturmSnowConductivityModel,
        ρ_snow::FT,
        earth_param_set,
    ) where {FT}

Computes the thermal conductivity, given the density
of the snow, according to Sturm, M., Holmgren, J., König, M., 
and Morris, K.: The thermal conductivity of seasonal snow, J. Glaciol., 43, 26–41,
https://doi.org/10.3189/S0022143000002781, 1997.

We have adjusted the original equation to make the coefficients
non-dimensional by multiplying and dividing by ρ\\_ice as appropriate.
"""
function snow_thermal_conductivity(
    scheme::SturmSnowConductivityModel,
    ρ_snow::FT,
    earth_param_set,
) where {FT}
    _ρ_ice = LP.ρ_cloud_ice(earth_param_set)

    x = min(ρ_snow / _ρ_ice, scheme.max / _ρ_ice)
    if x < scheme.threshold
        return scheme.b1 + scheme.m1 * x
    else
        return scheme.b2 + scheme.m2 * x + scheme.q2 * x^2
    end
end

"""
    diurnal_damping_depth(κ::FT, ρ::FT, earth_param_set)::FT where {FT}

Provides the characteristic depth for diurnal variations in temperature in snow, used in the
calculation of the surface tempearture. Formula is replicated from the Utah Energy Balance (UEB) model.
A paper describing the formula can be found at https://hess.copernicus.org/articles/18/5061/2014/.
"""
function diurnal_damping_depth(κ::FT, ρ::FT, earth_param_set)::FT where {FT}
    _cp_i = FT(LP.cp_i(earth_param_set))
    _DT_ = FT(earth_param_set.insol_params.day)
    return sqrt(κ * _DT_ / (pi * _cp_i * ρ))
end

"""
    surface_temp_scaling_depth(κ::FT, ρ::FT, earth_param_set)::FT where {FT}

Provides the scaling length scale for calculation of the surface temperature, which is stable
in the small-snowpack limit (it becomes (output) -> (input z) as (input z) -> 0).
We recomment inputting the z-per-ground-area (i.e. averaged over the grid by the snow cover fraction) to this function,
since 0 < z-per-ground-area <= z-per-snow-area, which means the reference depth will always remain inside the snowpack as (true z) -> 0.
"""
function surface_temp_scaling_length(
    κ::FT,
    ρ::FT,
    z::FT,
    earth_param_set,
)::FT where {FT}
    d0 = diurnal_damping_depth(κ, ρ, earth_param_set)
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
                          _ΔS::FT,
                          earth_param_set) where {FT}

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
    _ΔS::FT,
    earth_param_set,
) where {FT}
    S_safe = max(S, FT(0))
    _ρ_l = FT(LP.ρ_cloud_liq(earth_param_set))
    _T_ref = FT(LP.T_0(earth_param_set))
    _LH_f0 = FT(LP.LH_f0(earth_param_set))
    cp_s = specific_heat_capacity(q_l, earth_param_set)
    return _T_ref +
           (U + _ρ_l * _LH_f0 * S_safe * (1 - q_l)) /
           (_ρ_l * cp_s * (S_safe + _ΔS))
end

"""
    snow_bulk_density(SWE::FT, z::FT, earth_param_set) where {FT}

Returns the snow density given the current model state when depth and SWE are available.
Ensure the passed values are consistent in whether both have been multiplied by the snow-cover-fraction
or not (both must be per-snow-area or per-ground-area.)
"""
function snow_bulk_density(
    SWE_area::FT,
    z_area::FT,
    earth_param_set,
)::FT where {FT}
    _ρ_l = LP.ρ_cloud_liq(earth_param_set)
    ε = eps(FT) #for preventing dividing by zero
    #return SWE/z * ρ_l but tend to ρ_l as SWE → 0
    #also handle instabilities when z, SWE both near machine precision
    return max(SWE_area, ε) / max(z_area, SWE_area, ε) * _ρ_l
end

"""
    maximum_liquid_mass_fraction(ρ_snow::FT, T::FT, θ_r::FT, earth_param_set) where {FT}

Computes the maximum liquid water mass fraction, given
the density of the snow ρ_snow and other parameters.
"""
function maximum_liquid_mass_fraction(
    ρ_snow::FT,
    T::FT,
    θ_r::FT,
    earth_param_set,
) where {FT}
    _ρ_l = LP.ρ_cloud_liq(earth_param_set)
    _T_freeze = LP.T_freeze(earth_param_set)
    if T > _T_freeze
        return FT(0)
    else
        θ_r * _ρ_l / ρ_snow
    end
end

"""
    runoff_timescale(z::FT, Ksat::FT, Δt::FT) where {FT}

Computes the timescale for liquid water to percolate and leave the snowpack,
given the depth of the snowpack z and the hydraulic conductivity Ksat.
the passed "z" should be p.snow.z; that of the snow depth averaged over
the ground area (not the snow area), as it is in a ratio with Y.snow.S.
"""
function runoff_timescale(z::FT, Ksat::FT, Δt::FT) where {FT}
    τ = max(Δt, z / Ksat)
    return τ
end

"""
    volumetric_internal_energy_liq(T, earth_param_set)

Computes the volumetric internal energy of the liquid water
in the snowpack at a point at temperature T.
"""
function volumetric_internal_energy_liq(T::FT, earth_param_set) where {FT}
    _ρ_l = LP.ρ_cloud_liq(earth_param_set)
    _cp_l = LP.cp_l(earth_param_set)
    _T_ref = LP.T_0(earth_param_set)

    I_liq = _ρ_l * _cp_l * (T .- _T_ref)
    return I_liq
end

"""
    compute_water_runoff(
        S::FT,
        S_l::FT,
        T::FT,
        ρ_snow::FT,
        z::FT,
        Ksat::FT,
        Δt::FT,
        θ_r::FT,
        earth_param_set,
    )

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
    Ksat::FT,
    Δt::FT,
    θ_r::FT,
    earth_param_set,
) where {FT}
    τ = runoff_timescale(z, Ksat, Δt)
    q_l_max::FT = maximum_liquid_mass_fraction(ρ_snow, T, θ_r, earth_param_set)
    S_safe = max(S, FT(0))
    return -(S_l - q_l_max * S_safe) / τ * heaviside(S_l - q_l_max * S_safe)
end

"""
     phase_change_flux(U::FT, S::FT, q_l::FT, energy_flux::FT, Δt::FT, ΔS::FT, earth_param_set) where {FT}

Computes the volume flux of liquid water undergoing phase change, given the
applied energy flux and current state of U,S,q_l.
"""
function phase_change_flux(
    U::FT,
    S::FT,
    q_l::FT,
    energy_flux::FT,
    Δt::FT,
    ΔS::FT,
    earth_param_set,
) where {FT}
    S_safe = max(S, FT(0))

    energy_at_T_freeze =
        energy_from_q_l_and_swe(S_safe, q_l, ΔS, earth_param_set)
    Upred = U - energy_flux * Δt
    energy_excess = Upred - energy_at_T_freeze

    _LH_f0 = LP.LH_f0(earth_param_set)
    _ρ_liq = LP.ρ_cloud_liq(earth_param_set)
    _cp_i = LP.cp_i(earth_param_set)
    _cp_l = LP.cp_l(earth_param_set)
    _T_ref = LP.T_0(earth_param_set)
    _T_freeze = LP.T_freeze(earth_param_set)
    if energy_excess > 0 || (energy_excess < 0 && q_l > 0)
        return -energy_excess / Δt / _ρ_liq /
               ((_cp_l - _cp_i) * (_T_freeze - _T_ref) + _LH_f0)
    else
        return FT(0)
    end
end

"""
    energy_from_q_l_and_swe(S::FT, q_l::FT, ΔS::FT, earth_param_set) where {FT}

A helper function for compute the snow energy per unit area, given snow
water equivalent S, liquid fraction q_l, and snow model parameters.

This assumes that the snow is at the freezing point.
"""
function energy_from_q_l_and_swe(
    S::FT,
    q_l::FT,
    _ΔS::FT,
    earth_param_set,
) where {FT}
    _T_freeze = FT(LP.T_freeze(earth_param_set))
    _ρ_l = FT(LP.ρ_cloud_liq(earth_param_set))
    _T_ref = FT(LP.T_0(earth_param_set))
    _LH_f0 = FT(LP.LH_f0(earth_param_set))
    c_snow = specific_heat_capacity(q_l, earth_param_set)
    return _ρ_l * (S + _ΔS) * c_snow * (_T_freeze - _T_ref) -
           _ρ_l * S * (1 - q_l) * _LH_f0
end

"""
    energy_from_T_and_swe(S::FT, T::FT, _ΔS::FT, earth_param_set) where {FT}

A helper function for compute the snow energy per unit area, given snow
water equivalent S, bulk temperature T, and snow model parameters.

The liquid mass fraction is assumed to be zero if T<=T_freeze, and 1 otherwise.
"""
function energy_from_T_and_swe(
    S::FT,
    T::FT,
    _ΔS::FT,
    earth_param_set,
) where {FT}
    _ρ_l = FT(LP.ρ_cloud_liq(earth_param_set))
    _T_ref = FT(LP.T_0(earth_param_set))
    _T_freeze = FT(LP.T_freeze(earth_param_set))
    _LH_f0 = FT(LP.LH_f0(earth_param_set))
    _cp_i = FT(LP.cp_i(earth_param_set))
    _cp_l = FT(LP.cp_l(earth_param_set))

    if T <= _T_freeze
        return _ρ_l * _cp_i * (S + _ΔS) * (T - _T_ref) - _ρ_l * S * _LH_f0
    else
        return _ρ_l * (S + _ΔS) * _cp_l * (T - _T_ref)
    end

end

"""
    energy_flux_falling_snow(atmos, p, earth_param_set)

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
function energy_flux_falling_snow(atmos, p, earth_param_set)
    _LH_f0 = LP.LH_f0(earth_param_set)
    _ρ_liq = LP.ρ_cloud_liq(earth_param_set)
    ρe_snow = -_LH_f0 * _ρ_liq
    return @. lazy(ρe_snow * p.drivers.P_snow)
end

"""
    energy_flux_falling_rain(atmos, p, earth_param_set)

Returns the energy flux of falling rain for a PrescribedAtmosphere,
approximated as ρ_l e_l(T_atmos) * P_liq. The energy is per unit volume of liquid water,
and P_liq is expressed as the volume flux of liquid water resulting from the rain.

This method can be extended to coupled simulations, where atmos is of type
CoupledAtmosphere, and the energy flux of the falling rain is passed in the
cache `p`.  In that case, this should specify `atmos::PrescribedAtmosphere`.
"""
function energy_flux_falling_rain(atmos, p, earth_param_set)
    return @. lazy(
        volumetric_internal_energy_liq(p.drivers.T, earth_param_set) *
        p.drivers.P_liq,
    )
end

"""
    update_density_and_depth!(ρ_snow, z_snow, density::MinimumDensityModel, Y, p, earth_param_set)

Extends the update_density_and_depth! function for the MinimumDensityModel type; updates the snow density and depth in place.
"""
function update_density_and_depth!(
    ρ_snow,
    z_snow,
    density::MinimumDensityModel{FT},
    Y,
    p,
    earth_param_set,
) where {FT}
    _ρ_l = LP.ρ_cloud_liq(earth_param_set)
    @. ρ_snow = density.ρ_min * (1 - p.snow.q_l) + _ρ_l * p.snow.q_l
    @. z_snow = _ρ_l * Y.snow.S / ρ_snow #z_snow is per ground area, matching Y.snow.S
end

"""
    compute_extra_prog_tendency!(
        parameterization::Union{MinimumDensityModel, ConstantAlbedoModel, ZenithAngleAlbedoModel},
        model::SnowModel,
        dY,
        Y,
        p,
        t
        )

Updates all prognostic variables associated with the albedo and/or density models given the current model state.
This is the default method for the constant minimum density model and ConstantAlbedoModel or ZenithAngleAlbedoModel,
 which have no additional prognostic variables.
"""
function compute_extra_prog_tendency!(
    parameterization::Union{
        MinimumDensityModel,
        ConstantAlbedoModel,
        ZenithAngleAlbedoModel,
    },
    model::SnowModel,
    dY,
    Y,
    p,
    t,
)
    return nothing
end

"""
    update_q_vap_sfc_scheme(
        ζ,
        param_set,
        thermo_params,
        inputs,
        scheme,
        T_sfc,
        u_star,
        z_0m,
        z_0b,
        q_l
)

A helper function for the Surface Fluxes computation which, given a 
required and fixed set of positional arguments (ζ, `param_set`, `thermo_params`
inputs, scheme, `T_sfc`, `u_star`, `z_0m`, `z_0b`) and the snow
liquid fraction `q_l`, computes the snow surface specific humidity.

This function is required because with the EquilibriumGradient parameterization
for surface temperature, we must recompute q.
"""
function update_q_vap_sfc_scheme(
    ζ,
    param_set,
    thermo_params,
    inputs,
    scheme,
    T_sfc,
    u_star,
    z_0m,
    z_0b,
    q_l,
)
    q_atmos = inputs.q_tot_int
    ρ_atmos = inputs.ρ_int
    Δz = inputs.Δz
    T_atmos = inputs.T_int
    P_atmos = Thermodynamics.air_pressure(thermo_params, T_atmos, ρ_atmos, q_atmos)

     q_sfc = snow_surface_specific_humidity(
         T_sfc,
         q_l,
         T_atmos,
         P_atmos,
         q_atmos,
         Δz,
         param_set,
         thermo_params
    )
    return q_sfc
end

"""
    update_T_sfc_scheme(
        ζ,
        param_set,
        thermo_params,
        inputs,
        scheme,
        u_star,
        z_0m,
        z_0b,
        q_l,
        T_bulk, 
        κ,
        d,
        σ,
        ϵ,
        SW_n,
        LW_d,
)

A helper function for the Surface Fluxes computation which, given a 
required and fixed set of positional arguments (ζ, `param_set`, `thermo_params`
inputs, scheme, `u_star`, `z_0m`, `z_0b`) and the snow
liquid fraction `q_l`, bulk temperature `T_bulk`, snow conductivity κ, effective
depth of the top layer d, Stefan Boltzman constant σ, emissivity ϵ, net shortwave
radiation `SW_n`, and long wave radiation downwards `LW_d`, computes an updated
snow surface temperature.

It makes this estimate by incrementing the initial guess for snow surface temperature
`T_0` (stored in `inputs`) by the Newton update ΔT, where `ΔT = -f(T_0)/f'(T_0)` and
f(T) = SW_n + LW_n(T) + H(T) + L(T) +κ(T-T̄)/d = 0.

Be aware that if the snow surface specific humidity parameterization changes, 
we must also change the internals of this function.

"""
function update_T_sfc_scheme(
    ζ,
    param_set,
    thermo_params,
    inputs,
    scheme,
    u_star,
    z_0m,
    z_0b,
    q_l,
    T_bulk,
    κ,
    d,
    σ,
    ϵ,
    SW_n,
    LW_d,
)
    T_sfc = inputs.T_sfc_guess
    T_atmos = inputs.T_int
    ρ_atmos = inputs.ρ_int
    q_atmos = inputs.q_tot_int
    Δz = inputs.Δz
    P_atmos =
        Thermodynamics.air_pressure(thermo_params, T_atmos, ρ_atmos, q_atmos)
    ρ_sfc =
        ClimaLand.compute_ρ_sfc(param_set, T_atmos, P_atmos, q_atmos, Δz, T_sfc)
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
    q_sfc = qsat_over_ice * (1 - q_l) + q_l * (qsat_over_liq)
    ∂q∂T =
        Thermodynamics.∂q_vap_sat_∂T_from_L(
            thermo_params,
            qsat_over_ice,
            Thermodynamics.latent_heat_sublim(thermo_params, T_sfc),
            T_sfc,
        ) * (1 - q_l) +
        q_l * Thermodynamics.∂q_vap_sat_∂T_from_L(
            thermo_params,
            qsat_over_liq,
            Thermodynamics.latent_heat_vapor(thermo_params, T_sfc),
            T_sfc,
        )

    g_h = SurfaceFluxes.heat_conductance(
        param_set,
        ζ,
        u_star,
        inputs,
        z_0m,
        z_0b,
        scheme,
    )
    b_flux = SurfaceFluxes.buoyancy_flux(param_set, ζ, u_star, inputs)
    E = SurfaceFluxes.evaporation(
        param_set,
        inputs,
        g_h,
        q_atmos,
        q_sfc,
        ρ_sfc,
        inputs.moisture_model,
    )
    L = SurfaceFluxes.latent_heat_flux(
        param_set,
        inputs,
        E,
        inputs.moisture_model,
    )
    H = SurfaceFluxes.sensible_heat_flux(
        param_set,
        inputs,
        g_h,
        T_atmos,
        T_sfc,
        ρ_sfc,
        E,
    )
    _LH_v0 = Thermodynamics.Parameters.LH_v0(thermo_params)
    cp_d = Thermodynamics.Parameters.cp_d(thermo_params)
    ∂L∂T = ρ_sfc * g_h * _LH_v0 * ∂q∂T
    ∂H∂T = ρ_sfc * g_h * cp_d
    LW_n = -ϵ * (LW_d - σ * T_sfc^4)
    ∂LW_n∂T = 4 * ϵ * σ * T_sfc^3
    ΔT =
        -(d * (SW_n + LW_n + L + H) + κ * (T_sfc - T_bulk)) /
        (d * (∂LW_n∂T + ∂L∂T + ∂H∂T) + κ)
    return T_sfc + ΔT
end

"""
    solve_for_surface_temp_at_a_point(
        T_initial_guess::FT,
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
        earth_param_set,
        surf_temp::EquilibriumGradientTemperatureModel,
    )::FT where {FT}

Solves for T satisfying:
(A) f(T) = SW_n + LW_n(T) + H(T) + L(T) +κ(T-T̄)/d = 0

by
(1) first defining the update_T(T; ...) and update_q(T; ...) functions, which create
        T_{i+1} = ΔT(T_i) + T_i; ΔT = -f(T_i)/f'(T_i)
        q_{i+1} = q_sat(T_{i+1}; liq)* q_l + (1-q_l) q_sat(T_{i+1}; ice)
(2) Solving for the root of f(T) by evaluating these update functions each iteration of
    the surface fluxes solve.

Please note that we cannot use `ClimaLand.turbulent_fluxes!` and that functionality directly,
because the actual surface temperature may be different from the value found in the root solve.
If the value found in the root solve is above the freezing temp, we convert the excess fluxes into
melting, and recompute the surface fluxes using the freezing temperature of water.
"""
function solve_for_surface_temp_at_a_point(
    T_initial_guess::FT,
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
    earth_param_set,
    surf_temp::EquilibriumGradientTemperatureModel,
)::FT where {FT}
    config = SurfaceFluxes.SurfaceFluxConfig(roughness_model, gustiness)
    positional_default_args = (
        scheme = SurfaceFluxes.PointValueScheme(),
        solver_opts = nothing,
        flux_specs = nothing,
    )
    # u is already a vector when we get it from a coupled atmosphere, otherwise we need to make it one
    if u_atmos isa FT
        u = (u_atmos, FT(0))
    else
        u = u_atmos
    end

    thermo_params = LP.thermodynamic_parameters(earth_param_set)
    surface_flux_params = LP.surface_fluxes_parameters(earth_param_set)
    _grav = LP.grav(earth_param_set)
    _σ = LP.Stefan(earth_param_set)
    d = surface_temp_scaling_length(κ_snow, ρ_snow, z_snow, earth_param_set)
    ρ_atmos =
        Thermodynamics.air_density(thermo_params, T_atmos, P_atmos, q_atmos)
    update_q(args...) = update_q_vap_sfc_scheme(args..., q_l)
    update_T(args...) =
        update_T_sfc_scheme(args..., q_l, T_bulk, κ_snow, d, _σ, ϵ_snow, SW_net, LW_d)
    q_sfc = snow_surface_specific_humidity(
        T_initial_guess,
        q_l,
        T_atmos,
        P_atmos,
        q_atmos,
        atmos_h - h_sfc,
        surface_flux_params,
        thermo_params
    )

    output = SurfaceFluxes.surface_fluxes(
        surface_flux_params,
        T_atmos,
        q_atmos,
        FT(0),#phase_partition_atmos.liq,
        FT(0),#,phase_partition_atmos.ice,
        ρ_atmos,
        T_initial_guess,
        q_sfc,
        _grav * h_sfc,
        atmos_h - h_sfc,
        displ,
        u,
        (FT(0), FT(0)), # u_sfc
        nothing, # roughness inputs
        config,
        positional_default_args...,
        update_T,
        update_q,
    )
    return output.T_sfc
end

"""
    surface_residual_flux(
        T_sfc_root::FT,
        κ::FT,
        ρ::FT,
        z::FT,
        earth_param_set
    )

Determines the residual surface energy flux to dump into the bulk energy whenever the correct surface
temperature to satisfy Fourier's Law of Conduction at the surface is greater than T_freeze.
"""
function surface_residual_flux(
    T_sfc_root::FT,
    κ::FT,
    ρ::FT,
    z::FT,
    earth_param_set,
)::FT where {FT}
    _T_freeze = FT(LP.T_freeze(earth_param_set))
    d_safe = max(surface_temp_scaling_length(κ, ρ, z, earth_param_set), eps(FT))
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
            p.snow.T_sfc,
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
            model.parameters.earth_param_set,
            surf_temp,
        )

    #set residual flux using values if T_sfc > T_freeze:
    p.snow.surf_residual_flux .=
        surface_residual_flux.(
            p.snow.T_sfc,
            p.snow.κ,
            p.snow.ρ_snow,
            p.snow.z_snow,
            model.parameters.earth_param_set,
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
