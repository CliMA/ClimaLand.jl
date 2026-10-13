export TwoComponentResNetLAIParameters,
    compute_resnet_carbon_allocation,
    compute_aridity_index,
    compute_seasonality_index,
    compute_woody_fraction,
    compute_two_component_partitioning,
    compute_infiltration_pulse_hysteresis,
    compute_conifer_cold_buffer,
    compute_two_component_lai,
    compute_dryland_phenology_index,
    compute_thermal_freeze_gate,
    compute_drought_dormancy_gate,
    compute_dual_trigger_dormancy,
    compute_hydro_optical_greenness,
    compute_drought_accelerated_woody_turnover,
    compute_structural_canopy_floor

"""
    TwoComponentResNetLAIParameters{FT<:AbstractFloat}

Parameter container for the Two-Component Canopy (slow woody layer + fast herbaceous pulse layer)
prognostic LAI model with ResNet-PINN dynamic carbon allocation and Cycle 86 Grand Hydro-Thermal
Dynamic Equilibrium architecture.

# Theoretical Background
The model advances beyond standard single-pool models by representing the canopy as two coupled
dynamical layers:
1. **Slow Woody Structural Layer** (`L_wood`): Governed by slow turnover (\\alpha_w \\approx 0.010 \\text{ day}^{-1},
   turnover time \\tau_w \\approx 100 \\text{ days}), capturing dense structural wood and forest canopies,
   with continuous drought-deciduous senescence acceleration under root-zone moisture deficits.
2. **Fast Herbaceous Pulse Layer** (`L_herb`): Governed by rapid turnover (\\alpha_h \\approx 0.592 \\text{ day}^{-1},
   turnover time \\tau_h \\approx 1.7 \\text{ days}), capturing opportunistic grasses, savannas, and shrub green-up.

Dynamic carbon allocation is modulated via a 2-layer Physics-Informed Residual Network (ResNet-PINN)
with Swish activations mapping environmental state vectors [A_0, \\text{VPD}, S_{\\text{soil}}] to allocation efficiency:
```math
h_1 = \\operatorname{swish}(\\mathbf{W}_1 \\mathbf{x} + \\mathbf{b}_1)
h_2 = \\operatorname{swish}(\\mathbf{W}_2 h_1 + \\mathbf{b}_2) + h_1
p_{\\text{alloc}} = \\operatorname{clamp}\\left(0.90 + 0.50 \\tanh(\\mathbf{W}_3 h_2 + \\mathbf{b}_3), 0.60, 1.50\\right)
```

Dual-trigger dormancy couples sub-zero thermal freezing collapse (f_freeze) and
dryland soil moisture drought collapse (f_drought):
```math
f_{\\text{dormant}} = f_{\\text{freeze}} \\cdot f_{\\text{drought}}
```

# References
- Zhou et al. (2025), "A General Model for the Seasonal to Decadal Dynamics of Leaf Area", Global Change Biology. https://onlinelibrary.wiley.com/doi/pdf/10.1111/gcb.70125

\$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct TwoComponentResNetLAIParameters{FT <: AbstractFloat}
    "Light extinction coefficient (dimensionless), default 0.5"
    k::FT = FT(0.5)
    "Leaf construction/maintenance cost (mol m^-2 yr^-1), default 12.227"
    z::FT = FT(12.227)
    "Base departure parameter from square-wave dynamics, default 0.771"
    sigma::FT = FT(0.771)
    "Base smoothing turnover rate (day^-1), default 0.0305"
    alpha_base::FT = FT(0.0305)
    "Water-limitation turnover acceleration gain, default 0.1154"
    alpha_water::FT = FT(0.1154)
    "Slow woody pool turnover rate multiplier (day^-1), default 0.010"
    alpha_wood::FT = FT(0.010)
    "Fast herbaceous pool turnover rate multiplier (day^-1), default 0.59225"
    alpha_herb::FT = FT(0.59225)
    "Herbaceous pulse flush gain, default 1.13265"
    herb_pulse::FT = FT(1.13265)
    "ResNet layer 1 weight for A0 (GPP), default 0.51342"
    w1_a0::FT = FT(0.51342)
    "ResNet layer 1 weight for VPD, default -0.78294"
    w1_vpd::FT = FT(-0.78294)
    "ResNet layer 1 weight for soil moisture, default 0.65428"
    w1_sw::FT = FT(0.65428)
    "ResNet layer 1 bias, default -0.05279"
    b1::FT = FT(-0.05279)
    "ResNet layer 2 hidden weight, default 0.50581"
    w2_h::FT = FT(0.50581)
    "ResNet layer 2 bias, default 0.18767"
    b2::FT = FT(0.18767)
    "ResNet output layer weight, default 0.60059"
    w3_out::FT = FT(0.60059)
    "ResNet output layer bias, default 0.46032"
    b3_out::FT = FT(0.46032)
    "Autumn senescence acceleration factor, default 0.73843"
    autumn_damp::FT = FT(0.73843)
    "High-latitude conifer winter cold buffer gain, default 0.39384"
    enf_cold_gain::FT = FT(0.39384)
    "Topsoil pulse infiltration activation threshold (mm), default 3.12612"
    pulse_thresh::FT = FT(3.12612)
    "Optimum co-limitation smoothing exponent, default 7.5440"
    gamma_opt::FT = FT(7.5440)
    "Directional-to-diffuse radiation extinction ratio, default 1.2595"
    k_dir_scale::FT = FT(1.2595)
    "Diffuse radiation enhancement coefficient under cloud cover, default 0.1288"
    beta_diffuse::FT = FT(0.1288)
    "Canopy self-shading saturation exponent, default 1.2278"
    shade_sat_exp::FT = FT(1.2278)
    
    # --- Cycle 86 Grand Hydro-Thermal Dynamic Equilibrium Parameters ---
    "Drought-accelerated woody pool turnover senescence gain, default 0.5000"
    gamma_drought::FT = FT(0.5000)
    "Hydro-optical greenness leaf curing weight, default 0.1531"
    hydro_weight::FT = FT(0.1531)
    "Dual-trigger drought dormancy steepness exponent, default 3.5000"
    beta_dry::FT = FT(3.5000)
    "Minimum optical greenness fraction during full dormancy/curing, default 0.6500"
    green_min::FT = FT(0.6500)
    "Sub-zero thermal freezing collapse threshold (Kelvin), default 273.15 K"
    theta_freeze::FT = FT(273.15)
    "Thermal chilling collapse steepness (1/K), default 0.80"
    kappa_cold::FT = FT(0.80)
    "Aridity structural floor scaling threshold, default 1.8000"
    theta_ai_scale::FT = FT(1.8000)
    "Soil thermal thaw memory timescale (days), default 1.0000"
    tau_thaw::FT = FT(1.0000)
    "Cryogenic photoinhibition critical threshold, default 0.79989"
    thaw_crit::FT = FT(0.7998920800344707)
    "Cryogenic photoinhibition steepness, default 14.65825"
    kappa_thaw::FT = FT(14.658247142245926)
    "High-latitude conifer retention baseline fraction, default 0.1500"
    enf_retention_scale::FT = FT(0.1500)
end

Base.eltype(::TwoComponentResNetLAIParameters{FT}) where {FT} = FT
Base.broadcastable(x::TwoComponentResNetLAIParameters) = tuple(x)

"""
    compute_resnet_carbon_allocation(a0_norm, vpd_norm, sw_norm, params)

Computes the dynamic carbon allocation factor using a 2-layer ResNet with Swish activation and residual skip connection.

# Arguments
- `a0_norm::FT`: Normalized potential GPP (A0 / A0_peak).
- `vpd_norm::FT`: Normalized vapor pressure deficit (VPD / 20 hPa).
- `sw_norm::FT`: Normalized root-zone soil moisture (S_soil / 120 mm).
- `params::TwoComponentResNetLAIParameters{FT}`: Parameter container.

# Returns
- `p_alloc::FT`: Carbon allocation multiplier bounded in [0.60, 1.50].
"""
@inline function compute_resnet_carbon_allocation(
    a0_norm::FT,
    vpd_norm::FT,
    sw_norm::FT,
    params::TwoComponentResNetLAIParameters{FT},
) where {FT <: AbstractFloat}
    lin1 =
        params.w1_a0 * a0_norm +
        params.w1_vpd * vpd_norm +
        params.w1_sw * sw_norm +
        params.b1
    h1 = lin1 / (FT(1) + exp(-lin1)) # Swish activation
    lin2 = params.w2_h * h1 + params.b2
    h2 = (lin2 / (FT(1) + exp(-lin2))) + h1 # Residual Skip Connection
    return clamp(
        FT(0.90) + FT(0.50) * tanh(params.w3_out * h2 + params.b3_out),
        FT(0.60),
        FT(1.50),
    )
end

@inline function compute_resnet_carbon_allocation(
    params::TwoComponentResNetLAIParameters{FT},
    a0_norm::FT,
    vpd_norm::FT,
    sw_norm::FT,
) where {FT <: AbstractFloat}
    return compute_resnet_carbon_allocation(a0_norm, vpd_norm, sw_norm, params)
end

"""
    compute_aridity_index(precip_annual, A0_annual, vpd_gs)

Computes the dimensionless aridity index (supply-to-demand ratio proxy) from mean annual precipitation
(mol H2O m^-2 yr^-1), annual potential productivity (mol CO2 m^-2 yr^-1), and mean growing season VPD (Pa).
"""
@inline function compute_aridity_index(
    precip_annual::FT,
    A0_annual::FT,
    vpd_gs::FT = FT(1000.0),
) where {FT <: AbstractFloat}
    demand_proxy = max(A0_annual * (vpd_gs / FT(100.0)) * FT(25.0), eps(FT))
    return clamp(precip_annual / demand_proxy, FT(0.05), FT(5.0))
end

"""
    compute_seasonality_index(GSL)

Computes the dimensionless seasonality index from growing season length (GSL in days).
Tropical evergreen systems (GSL ~ 365) yield SI ~ 0, whereas strongly seasonal deciduous
or boreal systems (GSL ~ 120-180) yield SI ~ 0.5-0.7.
"""
@inline function compute_seasonality_index(GSL::FT) where {FT <: AbstractFloat}
    return clamp(FT(1.0) - GSL / FT(365.25), FT(0.0), FT(1.0))
end

"""
    compute_woody_fraction(aridity_index, seasonality_index)

Computes the continuous structural woody fraction `f_wood` from eco-climatic metrics.
"""
@inline function compute_woody_fraction(
    aridity_index::FT,
    seasonality_index::FT,
) where {FT <: AbstractFloat}
    # Humid and low-seasonality regions favor structural woody biomass; arid/seasonal favor herbaceous flush
    return clamp(
        FT(0.50) + FT(0.40) * tanh(aridity_index - FT(1.0)) -
        FT(0.30) * tanh(seasonality_index - FT(0.5)),
        FT(0.05),
        FT(0.95),
    )
end

"""
    compute_two_component_partitioning(aridity_index, seasonality_index)

Computes the continuous structural woody fraction `f_wood` and herbaceous fraction `f_herb = 1 - f_wood`
from eco-climatic metrics without static biome classifications.
"""
@inline function compute_two_component_partitioning(
    aridity_index::FT,
    seasonality_index::FT,
) where {FT <: AbstractFloat}
    f_wood = compute_woody_fraction(aridity_index, seasonality_index)
    return f_wood, FT(1) - f_wood
end

"""
    compute_infiltration_pulse_hysteresis(s_top, precip_flux, pulse_thresh, dt)

Updates the shallow topsoil moisture reservoir capturing post-rain opportunistic green-up.
"""
@inline function compute_infiltration_pulse_hysteresis(
    s_top::FT,
    precip_flux::FT,
    pulse_thresh::FT,
    dt::FT,
) where {FT <: AbstractFloat}
    s_in = max(FT(0), precip_flux * dt)
    s_new = (s_top + s_in) * exp(-dt / FT(86400.0 * 3.0)) # 3-day topsoil memory
    pulse_activity = clamp(
        (s_new - pulse_thresh) / max(pulse_thresh, eps(FT)),
        FT(0),
        FT(2.0),
    )
    return s_new, pulse_activity
end

"""
    compute_conifer_cold_buffer(t_air_k, enf_gain)

Applies optical cold attenuation for evergreen conifer needles during sub-freezing winter conditions.
"""
@inline function compute_conifer_cold_buffer(
    t_air_k::FT,
    enf_gain::FT,
) where {FT <: AbstractFloat}
    t_celsius = t_air_k - FT(273.15)
    cold_factor = FT(1) / (FT(1) + exp(-FT(0.25) * (t_celsius + FT(5.0))))
    return FT(1) - enf_gain * (FT(1) - cold_factor)
end

"""
    compute_dryland_phenology_index(rainforest_idx, evergreen_idx, p_mm_ann)

Computes continuous dryland/savanna sensitivity index isolating ecosystems subject to seasonal moisture drought.
"""
@inline function compute_dryland_phenology_index(
    rainforest_idx::FT,
    evergreen_idx::FT,
    p_mm_ann::FT,
) where {FT <: AbstractFloat}
    return clamp(
        (FT(1.0) - rainforest_idx) *
        (FT(1.0) - evergreen_idx) *
        max(FT(0.0), FT(1.0) - p_mm_ann / FT(850.0)),
        FT(0.0),
        FT(1.0),
    )
end

"""
    compute_thermal_freeze_gate(t_air_k, theta_freeze, kappa_cold)

Computes sub-zero thermal freezing collapse factor for structural overwintering floors.
Operates on physical temperature in Kelvin (default theta_freeze = 273.15 K, kappa_cold = 0.80 1/K).
"""
@inline function compute_thermal_freeze_gate(
    t_air_k::FT,
    theta_freeze::FT = FT(273.15),
    kappa_cold::FT = FT(0.80),
) where {FT <: AbstractFloat}
    return FT(1) / (FT(1) + exp(-kappa_cold * (t_air_k - theta_freeze)))
end

"""
    compute_drought_dormancy_gate(soil_stress, dryland_pheno_idx, beta_dry)

Computes drought dormancy gate collapsing structural biomass floors during prolonged severe dry seasons.
"""
@inline function compute_drought_dormancy_gate(
    soil_stress::FT,
    dryland_pheno_idx::FT,
    beta_dry::FT = FT(3.5),
) where {FT <: AbstractFloat}
    return FT(1.0) - dryland_pheno_idx * (FT(1.0) - soil_stress)^beta_dry
end

"""
    compute_dual_trigger_dormancy(f_freeze, f_drought)

Combines cryogenic freeze dormancy and hydro-drought dormancy multiplicatively.
"""
@inline function compute_dual_trigger_dormancy(
    f_freeze::FT,
    f_drought::FT,
) where {FT <: AbstractFloat}
    return f_freeze * f_drought
end

"""
    compute_hydro_optical_greenness(deciduous_idx, evergreen_idx, dryland_pheno_idx, pheno_gate, cryo_gate, f_freeze, eff_soil_stress, green_min, hydro_weight)

Computes the active optical green leaf area fraction, accounting for photoperiod, cryo-thermal photoinhibition,
and drought-induced leaf curing.
"""
@inline function compute_hydro_optical_greenness(
    deciduous_idx::FT,
    evergreen_idx::FT,
    dryland_pheno_idx::FT,
    pheno_gate::FT,
    cryo_gate::FT,
    f_freeze::FT,
    eff_soil_stress::FT,
    green_min::FT = FT(0.65),
    hydro_weight::FT = FT(0.1531),
) where {FT <: AbstractFloat}
    hydro_curing = hydro_weight * dryland_pheno_idx * (FT(1.0) - eff_soil_stress)
    f_green =
        FT(1.0) -
        (FT(1.0) - green_min) * (
            deciduous_idx * (FT(1.0) - pheno_gate) +
            FT(0.50) * evergreen_idx * (FT(1.0) - cryo_gate) * (FT(1.0) - f_freeze) +
            hydro_curing
        )
    return clamp(f_green, FT(0.05), FT(1.0))
end

"""
    compute_drought_accelerated_woody_turnover(alpha_w_base, gamma_drought, dryland_pheno_idx, soil_stress)

Computes woody pool senescence rate with continuous drought acceleration in savannas and drylands.
"""
@inline function compute_drought_accelerated_woody_turnover(
    alpha_w_base::FT,
    gamma_drought::FT,
    dryland_pheno_idx::FT,
    soil_stress::FT,
) where {FT <: AbstractFloat}
    drought_senescence =
        FT(1.0) + gamma_drought * dryland_pheno_idx * (FT(1.0) - soil_stress)
    return alpha_w_base * drought_senescence
end

"""
    compute_two_component_lai(l_wood, l_herb, l_steady_wood, l_steady_herb, alpha_w, alpha_h, dt)

Updates slow woody and fast herbaceous canopy pools towards their respective steady-state targets
via analytical exponential relaxation over time step `dt`.
"""
@inline function compute_two_component_lai(
    l_wood::FT,
    l_herb::FT,
    l_steady_wood::FT,
    l_steady_herb::FT,
    alpha_w::FT,
    alpha_h::FT,
    dt::FT,
) where {FT <: AbstractFloat}
    new_wood = l_wood + (FT(1) - exp(-alpha_w * dt)) * (l_steady_wood - l_wood)
    new_herb = l_herb + (FT(1) - exp(-alpha_h * dt)) * (l_steady_herb - l_herb)
    return new_wood, new_herb
end

"""
    compute_structural_canopy_floor(rainforest_idx, evergreen_idx, deciduous_idx, ai, f_dormant, theta_ai_scale, enf_retention_scale, enf_cold_gain, t_air_k)

Computes the continuous biophysical structural canopy floor preventing unrealistic vegetation collapse while respecting drought and thermal dormancy.
"""
@inline function compute_structural_canopy_floor(
    rainforest_idx::FT,
    evergreen_idx::FT,
    deciduous_idx::FT,
    ai::FT,
    f_dormant::FT,
    theta_ai_scale::FT = FT(1.8),
    enf_retention_scale::FT = FT(0.15),
    enf_cold_gain::FT = FT(0.39384),
    t_air_k::FT = FT(293.15),
) where {FT <: AbstractFloat}
    arid_floor_scale = clamp(tanh(ai / max(FT(0.1), theta_ai_scale))^FT(1.5), FT(0.02), FT(1.0))
    rf_floor = FT(2.05) * rainforest_idx
    decid_floor = (FT(0.28) + FT(0.20) * tanh(ai - FT(1.0))) * deciduous_idx * arid_floor_scale * f_dormant
    t_celsius = t_air_k - FT(273.15)
    cold_factor = FT(1) / (FT(1) + exp(-FT(0.25) * (t_celsius + FT(5.0))))
    winter_cold_buffer = enf_cold_gain * evergreen_idx * (FT(1.0) - deciduous_idx) * (FT(1.0) - cold_factor)
    conifer_retention = (enf_retention_scale - winter_cold_buffer) * evergreen_idx * (FT(1.0) - deciduous_idx) * tanh(ai / FT(1.5)) * arid_floor_scale * f_dormant
    struct_floor = rf_floor + max(FT(0.15) * evergreen_idx * arid_floor_scale * f_dormant, conifer_retention) + decid_floor
    return max(FT(0.05), struct_floor)
end
