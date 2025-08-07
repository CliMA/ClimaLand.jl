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
function get_domain_info(
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

function get_location(lat, long; time_offset = get_time_offset(lat, long))
    return (; time_offset, lat, long)
end

"""
    get_parameters(lat, long, pft::Pft, domain; kwargs...)

Gets parameters for a generic Fluxnet site based on PFT (plant functional type),
supplemented with additional autofilled values from US-MOz (Missouri Ozark) site.

Data sources:

Hydraulics parameters:
    - Holtzman, Nataniel, et al. 2023, https://doi.org/10.1029/2023WR035481
"""
function get_parameters(FT,
    lat,
    long,
    pft::Pft = get_pft(lat, long),
    domain = create_site_column(FT, lat, long);
    soil_params_gupta = 
        ClimaLand.Soil.soil_vangenuchten_parameters(domain.space.subsurface, FT),
    soil_params_grids = 
        ClimaLand.Soil.soil_composition_parameters(domain.space.subsurface, FT),
    time_offset = get_time_offset(long),
    atmos_h = FT(50), # take max of each site, height >= h_canopy + 10
    soil_ν = soil_params_gupta[:ν], # gupta
    soil_K_sat = soil_params_gupta[:K_sat], # gupta
    soil_S_s = FT(1e-2), # TODO, constant for now
    soil_hydrology_cm = soil_params_gupta[:hydrology_cm],
    θ_r = soil_params_gupta[:θ_r], # gupta, residual
    ν_ss_quartz = soil_params_grids[:ν_ss_quartz], # soil grids
    ν_ss_om = soil_params_grids[:ν_ss_om], # soil grids
    ν_ss_gravel = soil_params_grids[:ν_ss_gravel], # soil grids
    z_0m_soil = FT(0.01), # constant?
    z_0b_soil = FT(0.01), # constant?
    soil_ϵ = FT(0.98), # constant
    soil_albedo = 
        ClimaLand.Soil.CLMTwoBandSoilAlbedo{FT}(; ClimaLand.Soil.clm_soil_albedo_parameters(domain.space.surface)...),
    Ω = pft.parameters.Ω,
    χl = pft.parameters.χl,
    G_Function = CLMGFunction(χl),
    α_PAR_leaf = pft.parameters.α_PAR_leaf,
    λ_γ_PAR = FT(5e-7), # constant
    α_NIR_leaf = pft.parameters.α_NIR_leaf,
    τ_PAR_leaf = pft.parameters.τ_PAR_leaf,
    τ_NIR_leaf = pft.parameters.τ_NIR_leaf,
    ϵ_canopy = pft.parameters.ϵ_canopy,
    ac_canopy = pft.parameters.ac_canopy,
    g1 = pft.parameters.g1,
    Drel = FT(1.6), # constant, ratio
    g0 = ClimaCore.Fields.field2array(ClimaLand.Canopy.clm_medlyn_g0(domain.space.surface))[1], # CLM, vegetation_properties map
    Vcmax25 = pft.parameters.Vcmax25,
    pc = FT(-2.0e6), # climaparameters, water pressure
    sc = FT(5e-6), # climaparameters, low water sensitivity
    SAI = FT(0), # constant
    f_root_to_shoot = pft.parameters.f_root_to_shoot,
    K_sat_plant = pft.parameters.K_sat_plant,
    ψ63 = pft.parameters.ψ63,
    Weibull_param = FT(4), # constant
    a = FT(0.1 * 0.0098), # constant
    conductivity_model = PlantHydraulics.Weibull{FT}(
        K_sat_plant,
        ψ63,
        Weibull_param,
    ),
    retention_model = PlantHydraulics.LinearRetentionCurve{FT}(a),
    plant_ν = pft.parameters.plant_ν,
    plant_S_s = FT(1e-2 * 0.0098), # constant, note: incorrectly set as static across globe
    rooting_depth = pft.parameters.rooting_depth,
    n_stem = Int64(0),
    n_leaf = Int64(1),
    h_stem = FT(0),
    h_leaf = FT(9.5), # map, z_top
    h_canopy = h_stem + h_leaf,
    z0_m = FT(0.13) * h_canopy,
    z0_b = FT(0.1) * z0_m,
) 
    return (;
        time_offset,
        lat,
        long,
        atmos_h,
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

"""
    get_time_offset(longitude::Float64)

Given a longitude, returns the time offset from UTC in hours.
"""
function FluxnetSimulations.get_time_offset(longitude::Float64)
    hours = round(abs(longitude) / 15)

    if (hours > 11)
        hours = 11 - hours
    end

    return hours
end


"""
    create_site_column(FT, lat, long)

Creates and returns a Column domain for user-given latitude and longitude coordinates.
"""
function FluxnetSimulations.create_site_column( # might not be necessary?
    FT,
    lat,
    long;
    dz_bottom = FT(1.5),
    dz_top = FT(0.1),
    nelements = 20,
    zmin = FT(-10),
    zmax = FT(0),
)
    (; dz_tuple, nelements, zmin, zmax) = 
        get_domain_info(FT; dz_bottom, dz_top, nelements, zmin, zmax)

    zlim = (zmin, zmax)
    longlat = (long, lat)

    return Domains.Column(; zlim, nelements, longlat, dz_tuple)
end

"""
    get_pft(lat, long)

Given latitude and longitude, returns the dominant PFT at the coordinate.
"""
function FluxnetSimulations.get_pft(lat, long) # use ClimaUtilities, space-varying input,
    # get path of map from artifact
    dominant_PFT_map_path = joinpath(
        ClimaLand.Artifacts.clm_data_folder_path(),
        "dominant_PFT_map.nc",
    )

    pft_map = NCDataset(dominant_PFT_map_path)

    dominant_PFTs = pft_map["dominant_PFT"][:, :]
    lati_XY = pft_map["LATIXY"][:, :]
    long_XY = pft_map["LONGXY"][:, :]

    # find closest recorded point
    dist = (lati_XY .- lat) .^ 2 .+ (long_XY .- long) .^ 2
    closest_idx = argmin(dist)
    pft_idx = CartesianIndices(dist)[closest_idx]

    # convert to Pft()
    dominant_pft_idx = dominant_PFTs[pft_idx] - 1

    close(pft_map)

    if (dominant_pft_idx > 14)
        if (dominant_pft_idx == 7)
            dominant_pft_idx == 14
        else
            dominant_pft_idx == 13
        end
    end

    dominant_pft = ClimaLand.default_pfts[dominant_pft_idx]

    return dominant_pft
end
