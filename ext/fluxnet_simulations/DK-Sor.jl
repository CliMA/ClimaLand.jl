#=
This file contains domain specifications, parameters, and other site-level information for running
ClimaLand at the DK-Sor fluxtower site.
Citation: Sorø, Denmark beech forest FLUXNET site metadata used in the CalLMIP Phase 1a workflow.
=#

"""
    get_domain_info(FT, ::Val{:DK_Sor}; dz_bottom = FT(1.5), dz_top = FT(0.025),
        nelements = 20, zmin = FT(-9), zmax = FT(0))

Gets and returns primary domain information for the DK-Sor (Sorø, Denmark)
Fluxnet site. The values are provided as defaults, and can be overwritten by
passing the corresponding keyword arguments to this function.

The CalLMIP Phase 1a DK-Sor workflow uses a 20-layer stretched soil column
extending to approximately -9 m.
"""
function FluxnetSimulations.get_domain_info(
    FT,
    ::Val{:DK_Sor};
    dz_bottom = FT(1.5),
    dz_top = FT(0.025),
    nelements = 20,
    zmin = FT(-9),
    zmax = FT(0),
)
    dz_tuple = (dz_bottom, dz_top)
    return (; dz_tuple, nelements, zmin, zmax)
end

"""
    get_location(FT, ::Val{:DK_Sor}; kwargs...)

Returns geographical information for DK-Sor (Sorø, Denmark) Fluxnet site.
The values are provided as defaults, and can be overwritten by passing the
corresponding keyword arguments to this function.

The `time_offset` is the difference from UTC in hours and excludes daylight
savings time, following Fluxnet convention. For this site, the local standard
time is UTC+1 (Central European Time).
"""
function FluxnetSimulations.get_location(
    FT,
    ::Val{:DK_Sor};
    time_offset = 1,
    lat = FT(55.49),
    long = FT(11.64),
)
    return (; time_offset, lat, long)
end

"""
    get_fluxtower_height(FT, ::Val{:DK_Sor}; kwargs...)

Returns atmosphere height for DK-Sor (Sorø, Denmark) Fluxnet site.
The values are provided as defaults, and can be overwritten by passing the
corresponding keyword arguments to this function.
"""
function FluxnetSimulations.get_fluxtower_height(
    FT,
    ::Val{:DK_Sor};
    atmos_h = FT(32),
)
    return (; atmos_h,)
end

"""
    get_parameters(FT, ::Val{:DK_Sor}; kwargs...)

Gets parameters for the Fluxnet site DK-Sor (Sorø, Denmark), which is a
temperate deciduous broadleaf forest, and returns them as a NamedTuple.
Default parameters are provided and can be overwritten using keyword arguments.

These defaults are aligned with the DK-Sor single-site calibration and CalLMIP
forward-model interfaces in this repository.
"""
function FluxnetSimulations.get_parameters(
    FT,
    ::Val{:DK_Sor};
    soil_ν = FT(0.45),
    soil_K_sat = FT(1e-5),
    soil_S_s = FT(1e-3),
    soil_vg_n = FT(1.6),
    soil_vg_α = FT(1.6),
    θ_r = FT(0.07),
    ν_ss_quartz = FT(0.47),
    ν_ss_om = FT(0.03),
    ν_ss_gravel = FT(0.12),
    z_0m_soil = FT(0.01),
    z_0b_soil = FT(0.001),
    soil_ϵ = FT(0.98),
    soil_α_PAR = FT(0.2),
    soil_α_NIR = FT(0.2),
    Ω = FT(1),
    χl = FT(0.25),
    G_Function = ConstantGFunction(FT(0.5)),
    α_PAR_leaf = FT(0.1),
    λ_γ_PAR = FT(5e-7),
    τ_PAR_leaf = FT(0.05),
    α_NIR_leaf = FT(0.45),
    τ_NIR_leaf = FT(0.25),
    ϵ_canopy = FT(0.97),
    ac_canopy = FT(2500),
    g1 = FT(141),
    Drel = FT(1.6),
    g0 = FT(1e-4),
    Vcmax25 = FT(9e-5),
    SAI = FT(1.5),
    f_root_to_shoot = FT(3.5),
    K_sat_plant = FT(5e-9),
    ψ63 = FT(-4 / 0.0098),
    Weibull_param = FT(4),
    a = FT(0.05 * 0.0098),
    conductivity_model = Canopy.Weibull{FT}(K_sat_plant, ψ63, Weibull_param),
    retention_model = Canopy.LinearRetentionCurve{FT}(a),
    plant_ν = FT(2.46e-4),
    plant_S_s = FT(1e-2 * 0.0098),
    rooting_depth = FT(0.3),
    h_canopy = FT(25),
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
        h_canopy,
    )
end
