###################################
#        MODULE FUNCTIONS         #
###################################

"""
    get_domain_info(FT; dz_bottom = FT(1.5), dz_top = FT(0.1),
        nelements = 20, zmin = FT(-10),  zmax = FT(0))

Gets and returns primary domain information for a generic Fluxnet site, using
default values corresponding to a 10m deep soil column with 20 layers, with
a resolution of 10 cm at the top of the domain and 1.5m at the bottom of the domain.
"""
function FluxnetSimulations.get_domain_info(
    FT;
    dz_bottom = FT(1.5),
    dz_top = FT(0.1),
    nelements = 20,
    zmin = FT(-10),
    zmax = FT(0),
)

    dz_tuple = (dz_bottom, dz_top)

    return (; dz_tuple, nelements, zmin, zmax)
end

function FluxnetSimulations.get_location(
    site_ID::String;
    site_info = FluxnetSimulationsExt.get_site_info(site_ID),
    lat = site_info[:lat],
    long = site_info[:long],
    time_offset = site_info[:time_offset],
)
    atmos_h = site_info[:atmospheric_sensor_height][1] # take the first height as default
    return (; time_offset, lat, long, atmos_h)
end

"""
    get_parameters(FT, site_ID, domain, pft::Pft; kwargs...)

Gets parameters for a generic Fluxnet site based on PFT (plant functional type),
supplemented with additional autofilled values from US-MOz (Missouri Ozark) site.

Data sources:

Hydraulics parameters:
    - Holtzman, Nataniel, et al. 2023, https://doi.org/10.1029/2023WR035481
"""
function FluxnetSimulations.get_parameters(
    FT,
    site_ID,
    domain,
    pft::Pft = ClimaLand.Canopy.clm_dominant_pft(domain.space.surface),
    ;
    soil_params_gupta = ClimaLand.Soil.soil_vangenuchten_parameters(
        domain.space.subsurface,
        FT,
    ),
    soil_params_grids = ClimaLand.Soil.soil_composition_parameters(
        domain.space.subsurface,
        FT,
    ),
<<<<<<< HEAD
<<<<<<< HEAD
    soil_ν = ClimaCore.Fields.field2array(soil_params_gupta[:ν])[1],
    soil_K_sat = ClimaCore.Fields.field2array(soil_params_gupta[:K_sat])[1],
    soil_S_s = FT(1e-2),
    soil_hydrology_cm = ClimaCore.MatrixFields.column_field2array(
        soil_params_gupta[:hydrology_cm],
    )[1],
<<<<<<< HEAD
    θ_r = ClimaCore.Fields.field2array(soil_params_gupta[:θ_r])[1],
    ν_ss_quartz = ClimaCore.Fields.field2array(
        soil_params_grids[:ν_ss_quartz],
    )[1],
    ν_ss_om = ClimaCore.Fields.field2array(soil_params_grids[:ν_ss_om])[1],
    ν_ss_gravel = ClimaCore.Fields.field2array(
        soil_params_grids[:ν_ss_gravel],
    )[1],
    z_0m_soil = FT(0.01),
    z_0b_soil = FT(0.01),
    soil_ϵ = FT(0.98),
    soil_albedo_params = ClimaLand.Soil.clm_soil_albedo_parameters(
        domain.space.surface,
    ),
    soil_albedo = ClimaLand.Soil.CLMTwoBandSoilAlbedo{FT}(;
        NamedTuple{keys(soil_albedo_params)}((
            ClimaCore.Fields.field2array(v)[1] for
            v in values(soil_albedo_params)
        ))...,
=======
    soil_ν = soil_params_gupta[:ν],
    soil_K_sat = soil_params_gupta[:K_sat],
=======
    soil_ν = ClimaCore.Fields.field2array(soil_params_gupta[:ν])[1],
    soil_K_sat = ClimaCore.Fields.field2array(soil_params_gupta[:K_sat])[1],
>>>>>>> a523e6aba (pain)
    soil_S_s = FT(1e-2),
    soil_hydrology_cm = ClimaCore.MatrixFields.column_field2array(soil_params_gupta[:hydrology_cm])[1],
=======
>>>>>>> a199d603d (cleaning up for a branch clone)
    θ_r = ClimaCore.Fields.field2array(soil_params_gupta[:θ_r])[1],
    ν_ss_quartz = ClimaCore.Fields.field2array(
        soil_params_grids[:ν_ss_quartz],
    )[1],
    ν_ss_om = ClimaCore.Fields.field2array(soil_params_grids[:ν_ss_om])[1],
    ν_ss_gravel = ClimaCore.Fields.field2array(
        soil_params_grids[:ν_ss_gravel],
    )[1],
    z_0m_soil = FT(0.01),
    z_0b_soil = FT(0.01),
    soil_ϵ = FT(0.98),
    soil_albedo_params = ClimaLand.Soil.clm_soil_albedo_parameters(
        domain.space.surface,
    ),
    soil_albedo = ClimaLand.Soil.CLMTwoBandSoilAlbedo{FT}(;
<<<<<<< HEAD
<<<<<<< HEAD
        ClimaLand.Soil.clm_soil_albedo_parameters(domain.space.surface)...,
>>>>>>> 4ccce0bef (for FluxnetSimulationsExt module, wrote access functions to get simulation info + parameters for 4 specific sites with hardcoded information and any other site with mapped info)
=======
         NamedTuple{keys(soil_albedo_params)}(
                  (ClimaCore.Fields.field2array(v)[1] for v in values(soil_albedo_params)))...,
>>>>>>> a523e6aba (pain)
=======
        NamedTuple{keys(soil_albedo_params)}((
            ClimaCore.Fields.field2array(v)[1] for
            v in values(soil_albedo_params)
        ))...,
>>>>>>> a199d603d (cleaning up for a branch clone)
    ),
    Ω = pft.parameters.Ω,
    χl = pft.parameters.χl,
    G_Function = CLMGFunction(χl),
    α_PAR_leaf = pft.parameters.α_PAR_leaf,
    λ_γ_PAR = FT(5e-7),
    α_NIR_leaf = pft.parameters.α_NIR_leaf,
    τ_PAR_leaf = pft.parameters.τ_PAR_leaf,
    τ_NIR_leaf = pft.parameters.τ_NIR_leaf,
    ϵ_canopy = pft.parameters.ϵ_canopy,
    ac_canopy = pft.parameters.ac_canopy,
    g1 = pft.parameters.g1,
    Drel = FT(1.6),
    g0 = ClimaCore.Fields.field2array(
        ClimaLand.Canopy.clm_medlyn_g0(domain.space.surface),
    )[1],
    Vcmax25 = pft.parameters.Vcmax25,
    pc = FT(-2.0e6),
    sc = FT(5e-6),
    SAI = FT(0),
    f_root_to_shoot = pft.parameters.f_root_to_shoot,
    K_sat_plant = pft.parameters.K_sat_plant,
    ψ63 = pft.parameters.ψ63,
    Weibull_param = FT(4),
    a = FT(0.1 * 0.0098),
    conductivity_model = PlantHydraulics.Weibull{FT}(
        K_sat_plant,
        ψ63,
        Weibull_param,
    ),
    retention_model = PlantHydraulics.LinearRetentionCurve{FT}(a),
    plant_ν = pft.parameters.plant_ν,
    plant_S_s = FT(1e-2 * 0.0098),
    rooting_depth = pft.parameters.rooting_depth,
    n_stem = Int64(0),
    n_leaf = Int64(1),
    h_leaf = ClimaCore.Fields.field2array(
        ClimaLand.Canopy.clm_z_top(domain.space.surface),
    )[1],
    h_canopy = FluxnetSimulationsExt.get_canopy_height(site_ID),
<<<<<<< HEAD
<<<<<<< HEAD
    h_stem = ((h_canopy - h_leaf)) > 0 ? h_canopy - h_leaf : FT(0.0),
=======
    h_stem = h_canopy - h_leaf,
>>>>>>> 4ccce0bef (for FluxnetSimulationsExt module, wrote access functions to get simulation info + parameters for 4 specific sites with hardcoded information and any other site with mapped info)
=======
    h_stem = ((h_canopy - h_leaf)) > 0 ? h_canopy - h_leaf : FT(0.0),
>>>>>>> fe513ef05 (run from fluxnet2015)
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

###################################
#            UTILITIES            #
###################################

"""
    replace_hyphen(old_site_ID::String)

Replaces all instances of hyphens in a given site ID string with underscores
and returns a Symbol of the reformatted site ID to be used as a Val{} type.

For example, an input string "US-MOz" would be output as "US_MOz".
"""
function FluxnetSimulations.replace_hyphen(old_site_ID::String)
    new_site_ID = replace(old_site_ID, "-" => "_")

    return Symbol(new_site_ID)
end
