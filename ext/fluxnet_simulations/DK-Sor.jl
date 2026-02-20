#=
This file contains domain specifications, parameters, and other site-level information for running
ClimaLand at the DK-Sor (Soroe, Denmark) fluxtower site.
DK-Sor is a temperate deciduous broadleaf (beech) forest, FLUXNET2015 Tier 1 site.
Citation: Pilegaard, K., et al. (2011), Increasing net CO2 uptake by a Danish beech forest during
the period from 1996 to 2009, Agricultural and Forest Meteorology, 151(7), 934-946.
=#

"""
    get_domain_info(FT, ::Val{:DK_Sor}; kwargs...)

Gets and returns primary domain information for the DK-Sor (Soroe, Denmark) site,
which is a temperate deciduous broadleaf (beech) forest.
"""
function FluxnetSimulations.get_domain_info(
    FT,
    ::Val{:DK_Sor};
    dz_bottom = FT(2.0),
    dz_top = FT(0.05),
    nelements = 20,
    zmin = FT(-10.0),
    zmax = FT(0),
)
    dz_tuple = (dz_bottom, dz_top)
    return (; dz_tuple, nelements, zmin, zmax)
end

"""
    get_location(FT, ::Val{:DK_Sor}; kwargs...)

Returns geographical information for DK-Sor (Soroe, Denmark) site.
The `time_offset` is the difference from UTC in hours
and excludes daylight savings time, following Fluxnet convention.
For this site, the local time is UTC+1 (CET).
"""
function FluxnetSimulations.get_location(
    FT,
    ::Val{:DK_Sor};
    time_offset = 1,
    lat = FT(55.48587),
    long = FT(11.64464),
)
    return (; time_offset, lat, long)
end

"""
    get_fluxtower_height(FT, ::Val{:DK_Sor}; kwargs...)

Returns atmosphere height for DK-Sor site.
Reference height = 57 m, canopy height = 25 m (from NetCDF metadata).
"""
function FluxnetSimulations.get_fluxtower_height(
    FT,
    ::Val{:DK_Sor};
    atmos_h = FT(57),
)
    return (; atmos_h,)
end

"""
    get_parameters(FT, ::Val{:DK_Sor}; kwargs...)

Gets parameters for the DK-Sor (Soroe, Denmark) beech forest site.

Soil parameters from SoilGrids/literature for temperate forest on loamy soil.
Plant hydraulics parameters for European beech (Fagus sylvatica).
"""
function FluxnetSimulations.get_parameters(
    FT,
    ::Val{:DK_Sor};
    soil_ν = FT(0.45),
    soil_K_sat = FT(4.42e-6),  # ~1.6 cm/hr, loamy soil
    soil_S_s = FT(1e-3),
    soil_vg_n = FT(1.56),
    soil_vg_α = FT(3.6),
    θ_r = FT(0.067),
    ν_ss_quartz = FT(0.25),
    ν_ss_om = FT(0.04),        # higher OM for forest soil
    ν_ss_gravel = FT(0.0),
    z_0m_soil = FT(0.01),
    z_0b_soil = FT(0.001),
    soil_ϵ = FT(0.98),
    soil_α_PAR = FT(0.2),
    soil_α_NIR = FT(0.2),
    Ω = FT(0.69),
    χl = FT(0.1),              # slightly planophile beech leaves
    G_Function = ConstantGFunction(FT(0.5)),
    α_PAR_leaf = FT(0.1),
    λ_γ_PAR = FT(5e-7),
    τ_PAR_leaf = FT(0.05),
    α_NIR_leaf = FT(0.45),
    τ_NIR_leaf = FT(0.25),
    ϵ_canopy = FT(0.97),
    g1 = FT(141),              # Medlyn g1 for temperate deciduous broadleaf
    Drel = FT(1.6),
    g0 = FT(1e-4),
    Vcmax25 = FT(6e-5),        # ~60 μmol/m²/s for beech
    ac_canopy = FT(2500),
    pc = FT(-2e6),
    sc = FT(5e-6),
    SAI = FT(1.0),             # stem area index for forest
    f_root_to_shoot = FT(0.25),
    K_sat_plant = FT(1.8e-8),
    ψ63 = FT(-4.0 / 0.0098),  # ~-4 MPa P50 for beech
    Weibull_param = FT(3),
    a = FT(0.05 * 0.0098),
    conductivity_model = PlantHydraulics.Weibull{FT}(
        K_sat_plant,
        ψ63,
        Weibull_param,
    ),
    retention_model = PlantHydraulics.LinearRetentionCurve{FT}(a),
    plant_ν = FT(0.005),
    plant_S_s = FT(1e-2 * 0.0098),
    rooting_depth = FT(1.5),   # deep roots for beech forest
    n_stem = Int64(1),
    n_leaf = Int64(1),
    h_leaf = FT(1.0),
    h_stem = FT(24.0),         # ~25 m canopy height
    h_canopy = h_leaf + h_stem,
)
    return (;
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
        χl,
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
    )
end
