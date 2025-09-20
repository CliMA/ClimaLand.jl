#=
This file contains domain specifications, parameters, and other site-level information for running
ClimaLand at the US-NR1 fluxtower site.
Citation: Peter D. Blanken, Russel K. Monson, Sean P. Burns, David R. Bowling, Andrew A. Turnipseed (2022), AmeriFlux FLUXNET-1F US-NR1 Niwot Ridge Forest (LTER NWT1), Ver. 3-5, AmeriFlux AMP, (Dataset). https://doi.org/10.17190/AMF/1871141
=#

"""
    get_domain_info(FT, ::Val{:US_NR1}; dz_bottom = FT(1.25), dz_top = FT(0.05),
        nelements = 20, zmin = FT(-10), zmax = FT(0))

Gets and returns primary domain information for the US-NR1 (Colorado Niwot Ridge)
Fluxnet site. The values are provided as defaults, and can be overwritten by passing the
corresponding keyword arguments to this function.
"""
function FluxnetSimulations.get_domain_info(
    FT,
    ::Val{:US_NR1};
    dz_bottom = FT(1.25),
    dz_top = FT(0.05),
    nelements = 20,
    zmin = FT(-10),
    zmax = FT(0),
)

    dz_tuple = (dz_bottom, dz_top)
    return (; dz_tuple, nelements, zmin, zmax)
end

"""
    get_location(FT, ::Val{:US_NR1}; kwargs...)

Returns geographical information for US-NR1 (Colorado Niwot Ridge) Fluxnet site.
The values are provided as defaults, and can be overwritten by passing the
corresponding keyword arguments to this function.

The `time_offset` is the difference from UTC in hours
and excludes daylight savings time, following Fluxnet convention.
For this site, the local time is UTC-7 for Mountain Standard Time (MST).
"""
function FluxnetSimulations.get_location(
    FT,
    ::Val{:US_NR1};
    time_offset = -7,
    lat = FT(40.0329),
    long = FT(-105.5464),
    atmos_h = FT(21.5),
)
    return (; time_offset, lat, long, atmos_h)
end

"""
    get_parameters(FT, ::Val{:US_NR1}; kwargs...)

Gets parameters for the Fluxnet site US-NR1 (Colorado Niwot Ridge),
which is an Evergreen Needleleaf Forest, and returns them as a Named Tuple.
The values are provided as defaults, and can be overwritten by passing the
corresponding keyword arguments to this function.

Data sources:

Conductance parameters:
    - Wang et al. 2021 https://doi.org/10.5194/gmd-14-6741-2021
Photosynthesis parameters:
    - Holtzman, Nataniel, et al. 2023, https://doi.org/10.1029/2023WR035481
Hydraulics parameters:
    - Holtzman, Nataniel, et al. 2023, https://doi.org/10.1029/2023WR035481
"""
function FluxnetSimulations.get_parameters(
    FT,
    ::Val{:US_NR1};
    soil_ν = FT(0.45),
    soil_K_sat = FT(4e-7),
    soil_S_s = FT(1e-3),
    soil_hydrology_cm = vanGenuchten{FT}(; α = FT(0.04), n = FT(2.05)),
    θ_r = FT(0.0),
    ν_ss_quartz = FT(0.1),
    ν_ss_om = FT(0.1),
    ν_ss_gravel = FT(0.0),
    z_0m_soil = FT(0.1),
    z_0b_soil = FT(0.1),
    soil_ϵ = FT(0.98),
    soil_albedo = Soil.ConstantTwoBandSoilAlbedo{FT}(;
        PAR_albedo = FT(0.2),
        NIR_albedo = FT(0.2),
    ),
    Ω = FT(0.71),
    χl = FT(0.5),
    G_Function = ConstantGFunction(FT(0.5)),
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
    Vcmax25 = FT(9e-5),
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
    n_leaf = Int64(1),
    n_stem = Int64(1),
    h_leaf = FT(6.5),
    h_stem = FT(7.5),
    h_canopy = h_leaf + h_stem,
    z_0m = FT(0.13) * h_canopy,
    z_0b = FT(0.1) * z_0m,
)
    return (;
        soil_ν,
        soil_K_sat,
        soil_S_s,
        soil_hydrology_cm,
        θ_r,
        ν_ss_quartz,
        ν_ss_om,
        ν_ss_gravel,
        z_0m_soil,
        z_0b_soil,
        soil_ϵ,
        soil_albedo,
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
        z_0m,
        z_0b,
    )

end
