#=
This file contains domain specifications, parameters, and other site-level information for running
ClimaLand at the US-Var fluxtower site.
Citation: Siyan Ma, Liukang Xu, Joseph Verfaillie, Dennis Baldocchi (2023), AmeriFlux FLUXNET-1F US-Var Vaira Ranch- Ione, Ver. 3-5, AmeriFlux AMP, (Dataset). https://doi.org/10.17190/AMF/1993904
=#

"""
    get_domain_info(FT, ::Val{:US_Var}; dz_tuple = nothing,
        nelements = 14, zmin = FT(-0.5), zmax = FT(0))

Gets and returns primary domain information for the US-Var (California Vaira Ranch Ione), which is a grassland,
Fluxnet site. The values are provided as defaults, and can be overwritten by passing the corresponding
keyword arguments to this function.

Domain depth source: Xu and Baldocchi, 2003.
"""
function FluxnetSimulations.get_domain_info(
    FT,
    ::Val{:US_Var};
    dz_tuple = FT.((0.05, 0.02)),
    nelements = 24,
    zmin = FT(-0.5),
    zmax = FT(0),
)
    return (; dz_tuple, nelements, zmin, zmax)
end

"""
    get_location(::Val{:US_Var}; kwargs...)

Returns geographical information for US-Var (California Vaira Ranch Ione) Fluxnet site.
The values are provided as defaults, and can be overwritten by passing the
corresponding keyword arguments to this function.

The `time_offset` is the difference from UTC in hours
and excludes daylight savings time, following Fluxnet convention.
For this site, the local time is UTC-8 for Pacific Standard Time (PST).
"""
function FluxnetSimulations.get_location(
    FT,
    ::Val{:US_Var};
    time_offset = -8,
    lat = FT(38.4133),
    long = FT(-120.9508),
    atmos_h = FT(2),
)
    return (; time_offset, lat, long, atmos_h)
end

"""
    get_parameters(FT, ::Val{:US_Var}; kwargs...)

Gets parameters for the Fluxnet site US-Var (California Vaira Ranch Ione)
and returns them as a Named Tuple. The values are provided as defaults, and can
be overwritten by passing the keyword arguments to this function.

Data sources:

Soil parameters:
    - Bonan, G. Climate change and terrestrial ecosystem modeling. Cambridge University Press, 2019.
    - Xu and Baldocchi 2003, doi:10.1016/j.agrformet.2003.10.004
Conductance parameters:
    - CLM 5.0 Tech Note: https://www2.cesm.ucar.edu/models/cesm2/land/CLM50_Tech_Note.pdf
Photosynthesis parameters:
    - CLM 5.0 Tech Note: https://www2.cesm.ucar.edu/models/cesm2/land/CLM50_Tech_Note.pdf
    - Slevin et al. 2015, https://doi.org/10.5194/gmd-8-295-2015
Hydraulics parameters:
    - Holtzman, Nataniel, et al. 2023, https://doi.org/10.1029/2023WR035481
    - Xu and Baldocchi 2003, doi:10.1016/j.agrformet.2003.10.004
"""
function FluxnetSimulations.get_parameters(
    FT,
    ::Val{:US_Var};
    soil_ν = FT(0.5),
    soil_K_sat = FT(0.45 / 3600 / 100),
    soil_S_s = FT(1e-3),
    soil_hydrology_cm = vanGenuchten{FT}(; α = FT(2.0), n = FT(2.0)),
    θ_r = FT(0.067),
    ν_ss_quartz = FT(0.3),
    ν_ss_om = FT(0.02),
    ν_ss_gravel = FT(0.0),
    z_0m_soil = FT(0.01),
    z_0b_soil = FT(0.01),
    soil_ϵ = FT(0.98),
    soil_albedo = Soil.ConstantTwoBandSoilAlbedo{FT}(;
        PAR_albedo = FT(0.2),
        NIR_albedo = FT(0.2),
    ),
    Ω = FT(1.0),
    χl = FT(0.5),
    G_Function = ConstantGFunction(χl),
    α_PAR_leaf = FT(0.11),
    λ_γ_PAR = FT(5e-7),
    τ_PAR_leaf = FT(0.05),
    α_NIR_leaf = FT(0.45),
    τ_NIR_leaf = FT(0.34),
    ϵ_canopy = FT(0.97),
    g1 = FT(166),
    Drel = FT(1.6),
    g0 = FT(1e-4),
    Vcmax25 = FT(2.5e-5),
    ac_canopy = FT(745),
    pc = FT(-3e5),
    sc = FT(4e-6),
    SAI = FT(0),
    f_root_to_shoot = FT(3.5),
    K_sat_plant = 2e-8,
    ψ63 = FT(-2.7 / 0.0098),
    Weibull_param = FT(4),
    a = FT(0.05 * 0.0098),
    conductivity_model = PlantHydraulics.Weibull{FT}(
        K_sat_plant,
        ψ63,
        Weibull_param,
    ),
    retention_model = PlantHydraulics.LinearRetentionCurve{FT}(a),
    plant_ν = FT(8.93e-3),
    plant_S_s = FT(1e-2 * 0.0098),
    rooting_depth = FT(0.3),
    n_stem = Int64(0),
    n_leaf = Int64(1),
    h_leaf = FT(0.5),
    h_stem = FT(0),
    h_canopy = h_leaf + h_stem,
    z0_m = FT(0.13) * h_canopy,
    z0_b = FT(0.1) * z0_m,
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
        z0_m,
        z0_b,
    )
end
