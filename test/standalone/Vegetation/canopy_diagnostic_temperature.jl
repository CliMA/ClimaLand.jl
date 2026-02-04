using SurfaceFluxes
using Thermodynamics
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaParams
FT = Float64
toml_dict = LP.create_toml_dict(FT)
earth_param_set = LP.LandParameters(toml_dict)
thermo_params = LP.thermodynamic_parameters(earth_param_set)
z_0m = FT(0.13)
z_0b = FT(0.013)
roughness_model = SurfaceFluxes.ConstantRoughnessParams{FT}(z_0m,z_0b)
gustiness = FT(1)
displ = FT(0.0)
LAI = FT(1.0)
r_stomata_canopy = FT(80)
leaf_Cd = FT(0.08)
SW_n = FT(300.0)
LW_d = FT(200.0)
T_air = FT(298)
P_air = FT(101325)
q_tot_air = FT(0.006)
h_air = FT(10);
u_air = FT(3);
h_sfc = FT(0)
q_vap_sfc_guess = Thermodynamics.q_vap_saturation(
    thermo_params,
    T_sfc,
    Thermodynamics.air_density(thermo_params, T_air, P_air, q_air),
    Thermodynamics.Liquid(),
)
T_sfc_guess = T_air
ϵ_canopy = FT(0.98)
ϵ_ground = FT(0.97)
T_ground = T_air
function LW(T_canopy; T_ground, ϵ_canopy, LW_d, _σ, ϵ_ground)
    LW_d_canopy = (1 - ϵ_canopy) * LW_d + ϵ_canopy * _σ * T_canopy^4
    LW_u_ground = ϵ_ground * _σ * T_ground^4 + (1 - ϵ_ground) * LW_d_canopy
    return ϵ_canopy * LW_d - 2 * ϵ_canopy * _σ * T_canopy^4 + ϵ_canopy * LW_u_ground
end

function update_T_sfc_at_a_point(
    ζ,
    param_set,
    thermo_params,
    inputs,
    scheme,
    u_star,
    z_0m,
    z_0b,
    leaf_Cd,
    AI,
    SW_d,
    LW_d,
    ϵ_ground,
    T_ground,
    ϵ_canopy
    _σ
    )
    function findzero(T_canopy; ζ,
                      param_set,
                      thermo_params,
                      inputs,
                      scheme,
                      u_star,
                      z_0m,
                      z_0b,
                      leaf_Cd,
                      AI,
                      SW_n,
                      LW_d,
                      ϵ_ground,
                      T_ground,
                      ϵ_canopy
                      _σ
                      )
        FT = eltype(param_set)
        Φ_sfc = SurfaceFluxes.surface_geopotential(inputs)
        Φ_int = SurfaceFluxes.interior_geopotential(param_set, inputs)
        T_int = inputs.T_int
        
        
        b_flux = buoyancy_flux(param_set, ζ, u_star, inputs)
        Δz_eff = effective_height(inputs)
        ΔU = windspeed(inputs, param_set, b_flux)
        ΔU_safe = max(ΔU, eps(FT))
        Cd = drag_coefficient(param_set, ζ, z_0m, Δz_eff, scheme)
        Ch = heat_exchange_coefficient(param_set, ζ, z_0m, z_0h, Δz_eff, scheme)
        g_h = SurfaceFluxes.heat_conductance(
            param_set,
            ζ,
            u_star,
            inputs,
            z_0m,
            z_0b,
            scheme,
    )
        ws = SurfaceFluxes.windspeed(param_set, ζ, u_star, inputs)
        g_land = leaf_Cd * u_star * AI
        
        ΔΦ = Φ_int - Φ_sfc
        cp_d = Thermodynamics.Parameters.cp_d(thermo_params)
        T_sfc =
            (T_int + T_canopy * g_land / g_h + ΔΦ / cp_d) / (1 + g_land / g_h)
        
        (shf, lhf, E, ρτxz, ρτyz) = compute_flux_components(
        param_set, inputs, Ch, Cd, T_sfc, inputs.q_vap_sfc, inputs.ρ_int, b_flux,
        )
        LW_d_canopy = (1 - ϵ_canopy) * LW_d + ϵ_canopy * _σ * T_canopy^4
        LW_u_ground = ϵ_ground * _σ * T_ground^4 + (1 - ϵ_ground) * LW_d_canopy
        LW_n = ϵ_canopy * LW_d - 2 * ϵ_canopy * _σ * T_canopy^4 + ϵ_canopy * LW_u_ground

        return SW_n+LW_n - shf - lhf
    end
    
    
end

function update_q_vap_sfc_at_a_point(
    ζ,
    param_set,
    thermo_params,
    inputs,
    scheme,
    T_sfc,
    u_star,
    z_0m,
    z_0b,
    leaf_Cd,
    LAI,
    r_stomata_canopy,
)
    FT = eltype(param_set)
    g_leaf = leaf_Cd * u_star * LAI
    g_stomata = 1 / r_stomata_canopy
    g_land = g_stomata * g_leaf / (g_leaf + g_stomata)
    g_h = SurfaceFluxes.heat_conductance(
        param_set,
        ζ,
        u_star,
        inputs,
        z_0m,
        z_0b,
        scheme,
    )
    
    q_vap_int = inputs.q_tot_int - inputs.q_liq_int - inputs.q_ice_int
    T_canopy = inputs.T_sfc_guess # UPDATE
    q_canopy = Thermodynamics.q_vap_saturation(
        thermo_params,
        T_canopy,
        inputs.ρ_int,
        Thermodynamics.Liquid(),
    )
    q_new = (g_land / g_h * q_canopy + q_vap_int) / (1 + g_land / g_h)
    return q_new
end
function dfuncT(u_star, g_h, earth_param_set)
    FT = eltype(earth_param_set)
    return FT(1)
end
function dfuncq(u_star,
                g_h,
                q_sat,
                T_sfc,
                earth_param_set,
                )
    FT = eltype(earth_param_set)
    return FT(1)
end
compute_turbulent_fluxes_at_a_point(P_air, T_air, q_tot_air, u_air, h_air, T_sfc_guess, q_vap_sfc_guess, roughness_model, update_T_sfc, update_q_vap_sfc, h_sfc, displ, dfuncT, dfuncq, gustiness, earth_param_set)
