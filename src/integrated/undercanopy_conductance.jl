#=
Under-canopy resistance for the turbulent fluxes of the soil surface.

In the integrated land models, the turbulent fluxes of the soil are computed by
SurfaceFluxes.jl against the atmospheric reference state, with the roughness of
bare soil. Under a canopy, however, turbulent exchange between the ground and
the air above is strongly reduced (e.g., Sellers et al., 1996; Zeng et al.,
2005; Sakaguchi and Zeng, 2009). Following CLM5 (`CanopyFluxesMod`), we
represent this with a parallel-path ("mosaic") conductance: a fraction
`W = exp(-(LAI + SAI))` of the ground (the canopy gap fraction) exchanges
directly with the atmosphere with the bare-soil conductance `g_h`, and the
remaining fraction exchanges through the under-canopy resistance

    r' = (1 + γ min(max(Ri, 0), 10)) / (C_s u*_c)

between the ground and the canopy air, in series with the aerodynamic
resistance between the canopy air and the atmosphere:

    g_eff = W g_h + (1 - W) / (1/g_h + r').

Here `C_s` is the dense-canopy ground transfer coefficient (0.004 in CLM5),
`u*_c` is the friction velocity above the canopy, and the stability correction
with `γ = 0.5` and the bulk Richardson number of the canopy air space
`Ri = g h_c (T_af - T_g) / (T_af u*_c²)` applies when the canopy air
(temperature `T_af`) is warmer than the ground (temperature `T_g`) (Sakaguchi
and Zeng, 2009). The aerodynamic resistance between the canopy air and the
atmosphere is approximated by the bare-soil resistance `1/g_h` from
SurfaceFluxes.jl, so that it includes the effects of atmospheric stability
consistently with the bare-soil path and `g_eff ≤ g_h` (see
`Soil.soil_conductance_ratio`); this overestimates the resistance above a
rough canopy, which is typically smaller than `r'`.

The friction velocity above the canopy is estimated from neutral similarity
with the canopy roughness length and displacement height, and the canopy air
temperature from the conductance-weighted mean of the air and leaf
temperatures. The ground temperature in the Richardson number is the soil
surface temperature at the time of evaluation (the top layer temperature, as
the skin temperature is solved for afterwards).

The resulting `W` and `r'` are stored in `p.soil.W_gap` and
`p.soil.r_undercanopy`, and are used in the soil skin temperature solve and in
the soil turbulent fluxes (see
`src/standalone/Soil/soil_surface_temperature.jl`). As `W → 1` (no
vegetation), `g_eff → g_h`, recovering bare soil.
=#

"""
    undercanopy_resistance_at_a_point(
        LAI, h_c, z_0m_c, z_0b_c, displ_c, leaf_Cd,
        T_canopy, T_ground, T_air, u_air, Δz_ref, gustiness, C_s,
        earth_param_set,
    )

Returns the resistance (s/m) between the ground and the canopy air in a dense
canopy, given the atmospheric state at the reference height `Δz_ref` above the
ground. See the discussion at the top of this file.
"""
function undercanopy_resistance_at_a_point(
    LAI::FT,
    h_c::FT,
    z_0m_c::FT,
    z_0b_c::FT,
    displ_c::FT,
    leaf_Cd::FT,
    T_canopy::FT,
    T_ground::FT,
    T_air::FT,
    u_air,
    Δz_ref::FT,
    gustiness::FT,
    C_s::FT,
    earth_param_set,
) where {FT}
    surface_flux_params = LP.surface_fluxes_parameters(earth_param_set)
    thermo_params = LP.thermodynamic_parameters(earth_param_set)
    κ = SurfaceFluxes.Parameters.von_karman_const(surface_flux_params)
    grav = LP.grav(earth_param_set)
    cp_d = Thermodynamics.Parameters.cp_d(thermo_params)
    U = max(u_air isa FT ? abs(u_air) : hypot(u_air[1], u_air[2]), gustiness)
    # Neutral similarity above the canopy. When the reference height is close
    # to the canopy top (or below it), the logarithmic profile does not apply;
    # we then use the log-law value at the height z_0m e above the
    # displacement height, which bounds u*_c by κU.
    z_eff = max(Δz_ref - displ_c, z_0m_c * FT(ℯ))
    u_star_c = κ * U / log(z_eff / z_0m_c)
    g_ac = κ * u_star_c / log(z_eff / z_0b_c)
    # Canopy air temperature: conductance-weighted mean of the air (brought
    # dry-adiabatically to the surface) and leaf temperatures
    T_a = T_air + grav * Δz_ref / cp_d
    g_leaf = leaf_Cd * u_star_c * max(LAI, FT(0))
    T_af = (g_ac * T_a + g_leaf * T_canopy) / (g_ac + g_leaf)
    # Stability correction (Sakaguchi and Zeng, 2009)
    Ri = grav * h_c * (T_af - T_ground) / (T_af * u_star_c^2)
    stability_factor = 1 + FT(0.5) * min(max(Ri, FT(0)), FT(10))
    return stability_factor / (C_s * u_star_c)
end

"""
    update_undercanopy_conductance!(p, soil::Soil.EnergyHydrology, canopy, Y)

Computes the canopy gap fraction `p.soil.W_gap` and the under-canopy
resistance `p.soil.r_undercanopy` used in the soil turbulent fluxes of
integrated land models with a canopy. This must be called after the auxiliary
state is updated (which resets these to their bare soil values) and before the
soil skin temperature and turbulent fluxes are computed. See the discussion at
the top of this file.
"""
function update_undercanopy_conductance!(
    p,
    soil::Soil.EnergyHydrology,
    canopy,
    Y,
)
    bc = soil.boundary_conditions.top
    bc isa Soil.AtmosDrivenFluxBC || return nothing
    sfp = canopy.boundary_conditions.turbulent_flux_parameterization
    sfp isa Canopy.MoninObukhovCanopyFluxes || return nothing
    atmos = bc.atmos
    earth_param_set = soil.parameters.earth_param_set
    area_index = p.canopy.biomass.area_index
    T_canopy = Canopy.canopy_temperature(canopy.energy, canopy, Y, p)
    T_ground = ClimaLand.Domains.top_center_to_surface(p.soil.T)
    h_sfc = ClimaLand.surface_height(soil, Y, p)
    C_s = soil.parameters.C_s_undercanopy
    @. p.soil.W_gap = exp(-(max(area_index.leaf, 0) + max(area_index.stem, 0)))
    @. p.soil.r_undercanopy = undercanopy_resistance_at_a_point(
        area_index.leaf,
        canopy.biomass.height,
        sfp.z_0m,
        sfp.z_0b,
        sfp.displ,
        sfp.Cd,
        T_canopy,
        T_ground,
        p.drivers.T,
        p.drivers.u,
        atmos.h - h_sfc,
        atmos.gustiness,
        C_s,
        earth_param_set,
    )
    return nothing
end
