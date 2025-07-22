"""
    get_domain_info(FT, site_ID::Symbol; dz_bottom = FT(1.5), dz_top = FT(0.1),
        nelements = 20, zmin = FT(-10),  zmax = FT(0))

Gets and returns primary domain information for a generic Fluxnet site,
using autofilled values from US-MOz (Missouri Ozark) site.
"""
function get_domain_info(
    FT,
    site_ID::Symbol;
    dz_bottom = FT(1.5),
    dz_top = FT(0.1),
    nelements = 20,
    zmin = FT(-10),
    zmax = FT(0),
) where {FT}

    dz_tuple = (dz_bottom, dz_top)

    return (
        dz_tuple = dz_tuple,
        nelements = nelements,
        zmin = zmin,
        zmax = zmax,
    )
end


"""
    get_parameters(site_ID::Symbol; kwargs...)

Gets parameters for a generic Fluxnet site,
using autofilled values from US-MOz (Missouri Ozark) site.

Data sources:

Conductance parameters:
    - Wang et al. 2021 https://doi.org/10.5194/gmd-14-6741-2021
Hydraulics parameters:
    - Holtzman, Nataniel, et al. 2023, https://doi.org/10.1029/2023WR035481
"""
function get_parameters(
    site_ID::Symbol;
    time_offset = 7,
    atmos_h = FT(32),
    lat = FT(38.7441),
    long = FT(-92.2000),
    soil_ν = FT(0.55),
    soil_K_sat = FT(4e-7),
    soil_S_s = FT(1e-2),
    soil_vg_n = FT(2.0),
    soil_vg_α = FT(0.05),
    θ_r = FT(0.04),
    ν_ss_quartz = FT(0.1),
    ν_ss_om = FT(0.1),
    ν_ss_gravel = FT(0.0),
    z_0m_soil = FT(0.01),
    z_0b_soil = FT(0.01),
    soil_ϵ = FT(0.98),
    soil_α_PAR = FT(0.2),
    soil_α_NIR = FT(0.2),
    Ω = FT(0.69),
    χl = FT(0.1),
    G_Function = CLMGFunction(χl),
    α_PAR_leaf = FT(0.1),
    λ_γ_PAR = FT(5e-7),
    τ_PAR_leaf = FT(0.05),
    α_NIR_leaf = FT(0.45),
    τ_NIR_leaf = FT(0.25),
    ϵ_canopy = FT(0.97),
    ac_canopy = FT(5e2),
    g1 = FT(141),
    Drel = FT(1.6),
    g0 = FT(1e-4),
    Vcmax25 = FT(6e-5),
    pc = FT(-2.0e6),
    sc = FT(5e-6),
    SAI = FT(1.0),
    f_root_to_shoot = FT(3.5),
    K_sat_plant = 7e-8,
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
    rooting_depth = FT(0.5),
    n_stem = Int64(1),
    n_leaf = Int64(1),
    h_stem = FT(9),
    h_leaf = FT(9.5),
    h_canopy = h_stem + h_leaf,
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

###################################
#            UTILITIES            #
###################################

"""
    replace_hyphen(old_site_ID::String)

Replaces all instances of hyphens in a given site ID string with underscores
and returns a Symbol of the reformatted site ID to be used as a Val{} type.
"""
function replace_hyphen(old_site_ID::String)
    new_site_ID = replace(old_site_ID, "-" => "_")

    return Symbol(new_site_ID)
end
