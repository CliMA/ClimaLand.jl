"""
    get_domain_info(FT, ::Val{:US_NR1}; dz_bottom = FT(1.25),
        dz_top = FT(0.05), nelements = 20, zmin = FT(-10), zmax = FT(0))

Gets and returns primary domain information for the US-NR1 (Colorado Niwot Ridge)
Fluxnet site.

The data source comes from: Unknown.
"""
function get_domain_info(
    FT,
    ::Val{:US_NR1};
    dz_bottom = FT(1.25),
    dz_top = FT(0.05),
    nelements = 20,
    zmin = FT(-10),
    zmax = FT(0),
)

    dz_tuple = (dz_bottom, dz_top)

    return (
        dz_tuple = dz_tuple,
        nelements = nelements,
        zmin = zmin,
        zmax = zmax,
    )
end

"""
    get_parameters(::Val{:US_NR1}; kwargs...)

Gets parameters for the Fluxnet site US-NR1 (Colorado Niwot Ridge)
and returns them as a Named Tuple.

Data sources:

Atmosphere height:
    - Metzger, Stefan & Burba, George & Burns, Sean & Blanken, Peter & Li, 
    Jiahong & Luo, Hongyan & Zulueta, Rommel. (2016). https://doi.org/10.5194/amt-9-1341-2016
Conductance parameters:
    - Wang et al. 2021 https://doi.org/10.5194/gmd-14-6741-2021
Photosynthesis parameters:
    - Holtzman, Nataniel, et al. 2023, https://doi.org/10.1029/2023WR035481
Hydraulics parameters:
    - Holtzman, Nataniel, et al. 2023, https://doi.org/10.1029/2023WR035481
"""
function get_parameters(
    ::Val{:US_NR1};
    time_offset = 7,
    lat = FT(40.0329),
    long = FT(-105.5464),
    atmos_h = FT(21.5),
    soil_ν = FT(0.45),
    soil_K_sat = FT(4e-7),
    soil_S_s = FT(1e-3),
    soil_vg_n = FT(2.05),
    soil_vg_α = FT(0.04),
    θ_r = FT(0.0),
    ν_ss_quartz = FT(0.1),
    ν_ss_om = FT(0.1),
    ν_ss_gravel = FT(0.0),
    z_0m_soil = FT(0.1),
    z_0b_soil = FT(0.1),
    soil_ϵ = FT(0.98),
    soil_α_PAR = FT(0.2),
    soil_α_NIR = FT(0.2),
    Ω = FT(0.71),
    ld = FT(0.5),
    G_Function = ConstantGFunction(ld),
    α_PAR_leaf = FT(0.1),
    λ_γ_PAR = FT(5e-7),
    τ_PAR_leaf = FT(0.05),
    α_NIR_leaf = FT(0.35),
    τ_NIR_leaf = FT(0.25),
    ϵ_canopy = FT(0.97),
    ac_canopy = FT(3e3),
    g1 = FT(141),
    Drel = FT(1.6),
    g0 = FT(1e-4),
    SAI = FT(1.0),
    f_root_to_shoot = FT(3.5),
    K_sat_plant = 5e-9,
    ψ63 = FT(-4 / 0.0098),
    Weibull_param = FT(4),
    a = FT(0.05 * 0.0098),
    conductivity_model = PlantHydraulics.Weibull{FT}(
        K_sat_plant,
        ψ63,
        Weibull_param,
    ),
    retention_model = PlantHydraulics.LinearRetentionCurve{FT}(a),
    plant_ν = FT(8.06e-4),
    plant_S_s = FT(1e-2 * 0.0098),
    rooting_depth = FT(1.0),
    n_stem = Int64(1),
    h_leaf = FT(6.5),
    h_stem = FT(7.5),
    h_canopy = h_leaf + h_stem,
    z0_m = FT(0.13) * h_canopy,
    z0_b = FT(0.1) * z0_m,
) where {FT}

    return (;
        time_offset,
        lat,
        long,
        atmos_h,
        soil_ν,
        soil_K_sat,
        soil_S_s,
        soil_vg_n,
        soil_vg_α,
        θ_r,
        ν_ss_quartz,
        ν_ss_om,
        ν_ss_gravel,
        z_0m_soil,
        z_0b_soil,
        soil_ϵ,
        soil_α_PAR,
        soil_α_NIR,
        Ω,
        ld,
        G_Function,
        α_PAR_leaf,
        λ_γ_PAR,
        τ_PAR_leaf,
        α_NIR_leaf,
        τ_NIR_leaf,
        ϵ_canopy,
        ac_canopy,
        g1,
        Drel,
        g0,
        Vcmax25,
        SAI,
        f_root_to_shoot,
        K_sat_plant,
        ψ63,
        Weibull_param,
        a,
        conductivity_model,
        retention_model,
        plant_ν,
        plant_S_s,
        rooting_depth,
        n_stem,
        n_leaf,
        h_leaf,
        h_stem,
        h_canopy,
        z0_m,
        z0_b,
    )

end
