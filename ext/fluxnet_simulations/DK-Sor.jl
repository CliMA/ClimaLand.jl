#=
This file contains domain specifications, parameters, and other site-level information for running
ClimaLand at the DK-Sor (Denmark Sorø) CalMIP site.

Site information:
- Location: Sorø, Denmark (55.486°N, 11.6446°E)
- Vegetation type: Deciduous Broadleaf Forest (Beech - Fagus sylvatica)
- CalMIP Phase 1a test calibration site
- Reference: https://github.com/callmip-org/Phase1
=#

"""
    FluxnetSimulations.get_domain_info(FT, ::Val{:DK_Sor}; dz_bottom = FT(1.5), dz_top = FT(0.1),
        nelements = 20, zmin = FT(-10), zmax = FT(0))

Gets and returns primary domain information for the DK-Sor (Denmark Sorø) CalMIP site.
Default parameters are provided and can be overriden using keyword arguments.
"""
function FluxnetSimulations.get_domain_info(
    FT,
    ::Val{:DK_Sor};
    dz_bottom = FT(1.5),
    dz_top = FT(0.1),
    nelements = 20,
    zmin = FT(-10),
    zmax = FT(0),
)
    dz_tuple = (dz_bottom, dz_top)
    return (; dz_tuple, nelements, zmin, zmax)
end

"""
    get_location(FT, ::Val{:DK_Sor}; kwargs)

Returns geographical information for DK-Sor (Denmark Sorø) CalMIP site.
The values are provided as defaults, and can be overwritten by passing the
corresponding keyword arguments to this function.

The `time_offset` is the difference from UTC in hours
and excludes daylight savings time, following Fluxnet convention.
For this site, Denmark is in the Central European Time zone (CET), which is UTC+1.
"""
function FluxnetSimulations.get_location(
    FT,
    ::Val{:DK_Sor};
    time_offset = 1,
    lat = FT(55.486),
    long = FT(11.6446),
)
    return (; time_offset, lat, long)
end

"""
    get_fluxtower_height(FT, ::Val{:DK_Sor}; kwargs...)

Returns atmosphere height for DK-Sor (Denmark Sorø) CalMIP site.
The values are provided as defaults, and can be overwritten by passing the
corresponding keyword arguments to this function.
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

Gets parameters for the CalMIP site DK-Sor (Denmark Sorø), which is a Deciduous Broadleaf Forest
dominated by beech trees (Fagus sylvatica), and returns them as a Named Tuple. 
The values are provided as defaults, and can be overwritten by passing the 
corresponding keyword arguments to this function.

Note: These parameters are initial estimates and should be calibrated using CalMIP protocol.
For calibration targets and observations, see the CalMIP Phase 1 protocol at:
https://github.com/callmip-org/Phase1

Data sources:
- Site characteristics from CalMIP Phase 1 protocol
- Initial parameter estimates based on deciduous broadleaf forest typical values
- Conductance parameters adapted from Wang et al. 2021 https://doi.org/10.5194/gmd-14-6741-2021
"""
function FluxnetSimulations.get_parameters(
    FT,
    ::Val{:DK_Sor};
    # Soil hydraulic and thermal parameters
    soil_ν = FT(0.50),
    soil_K_sat = FT(5e-7),
    soil_S_s = FT(1e-2),
    soil_vg_n = FT(1.8),
    soil_vg_α = FT(0.04),
    θ_r = FT(0.05),
    # Soil composition
    ν_ss_quartz = FT(0.15),
    ν_ss_om = FT(0.12),
    ν_ss_gravel = FT(0.0),
    # Soil surface properties
    z_0m_soil = FT(0.01),
    z_0b_soil = FT(0.01),
    soil_ϵ = FT(0.98),
    soil_α_PAR = FT(0.15),
    soil_α_NIR = FT(0.15),
    # Canopy structure parameters
    Ω = FT(0.70),
    χl = FT(0.1),
    G_Function = ConstantGFunction(FT(0.5)),
    # Leaf optical properties
    α_PAR_leaf = FT(0.11),
    λ_γ_PAR = FT(5e-7),
    τ_PAR_leaf = FT(0.05),
    α_NIR_leaf = FT(0.40),
    τ_NIR_leaf = FT(0.25),
    # Canopy radiative properties
    ϵ_canopy = FT(0.97),
    ac_canopy = FT(5e2),
    # Photosynthesis and stomatal conductance parameters
    g1 = FT(150),
    Drel = FT(1.6),
    g0 = FT(1e-4),
    Vcmax25 = FT(6.5e-5),
    # Water stress parameters
    pc = FT(-2.0e6),
    sc = FT(5e-6),
    # Structural parameters
    SAI = FT(1.0),
    f_root_to_shoot = FT(3.0),
    # Plant hydraulics parameters
    K_sat_plant = 8e-8,
    ψ63 = FT(-4 / 0.0098),
    Weibull_param = FT(4),
    a = FT(0.1 * 0.0098),
    conductivity_model = PlantHydraulics.Weibull{FT}(
        K_sat_plant,
        ψ63,
        Weibull_param,
    ),
    retention_model = PlantHydraulics.LinearRetentionCurve{FT}(a),
    plant_ν = FT(1.44e-4),
    plant_S_s = FT(1e-2 * 0.0098),
    # Root and canopy height parameters
    rooting_depth = FT(1.0),
    n_stem = Int64(1),
    n_leaf = Int64(1),
    h_stem = FT(25),
    h_leaf = FT(5),
    h_canopy = h_stem + h_leaf,
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
