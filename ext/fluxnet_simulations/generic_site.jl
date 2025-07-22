###################################
#        MODULE FUNCTIONS         #
###################################

"""
    get_domain_info(FT, site_ID::Symbol; dz_bottom = FT(1.5), dz_top = FT(0.1),
        nelements = 20, zmin = FT(-10),  zmax = FT(0))

Gets and returns primary domain information for a generic Fluxnet site, using
default values corresponding to a 10m deep soil column with 20 layers, with
a resolution of 10 cm at the top of the domain and 1.5m at the bottom of the domain.
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

function get_location(lat, long; time_offset = get_time_offset(lat, long))
    return (; time_offset, lat, long)
end

"""
    get_parameters(lat, long, pft::Pft; kwargs...)

Gets parameters for a generic Fluxnet site based on PFT (plant functional type),
supplemented with additional autofilled values from US-MOz (Missouri Ozark) site.

Data sources:

Hydraulics parameters:
    - Holtzman, Nataniel, et al. 2023, https://doi.org/10.1029/2023WR035481
"""
function get_parameters(
    lat,
    long;
    pft::Pft = get_pft(lat, long), soil_params_gupta = get_soil_params_gupta(lat, long),
    time_offset = get_time_offset(long), # calculate offset based on lat-long, DONE
    atmos_h = FT(50), # take max of each site, height >= h_canopy + 10
    soil_ν = FT(0.55), # gupta
    soil_K_sat = FT(4e-7), # gupta
    soil_S_s = FT(1e-2), # TODO, constant for now
    soil_vg_n = FT(2.0), # gupta
    soil_vg_α = FT(0.05), # gupta
    θ_r = FT(0.04), # gupta, residual
    ν_ss_quartz = FT(0.1), # soil grids
    ν_ss_om = FT(0.1), # soil grids
    ν_ss_gravel = FT(0.0), # soil grids
    z_0m_soil = FT(0.01), # constant?
    z_0b_soil = FT(0.01), # constant?
    soil_ϵ = FT(0.98), # constant
    Ω = pft.Ω,
    χl = pft.χl,
    G_Function = CLMGFunction(χl),
    α_PAR_leaf = pft.α_PAR_leaf,
    λ_γ_PAR = FT(5e-7), # constant
    α_NIR_leaf = pft.α_NIR_leaf,
    τ_PAR_leaf = pft.τ_PAR_leaf,
    τ_NIR_leaf = pft.τ_NIR_leaf,
    ϵ_canopy = pft.ϵ_canopy,
    ac_canopy = pft.ac_canopy,
    g1 = pft.g1,
    Drel = FT(1.6), # constant, ratio
    g0 = FT(1e-4), # CLM, vegetation_properties map
    Vcmax25 = pft.Vcmax25,
    pc = FT(-2.0e6), # climaparameters, water pressure
    sc = FT(5e-6), # climaparameters, low water sensitivity
    SAI = FT(0), # constant
    f_root_to_shoot = pft.f_root_to_shoot,
    K_sat_plant = pft.K_sat_plant,
    ψ63 = pft.ψ63,
    Weibull_param = FT(4), # CLM map
    a = FT(0.1 * 0.0098), # constant
    conductivity_model = PlantHydraulics.Weibull{FT}(
        K_sat_plant,
        ψ63,
        Weibull_param,
    ),
    retention_model = PlantHydraulics.LinearRetentionCurve{FT}(a),
    plant_ν = pft.plant_ν,
    plant_S_s = FT(1e-2 * 0.0098), # constant, note: incorrectly set as static across globe
    rooting_depth = pft.rooting_depth,
    n_stem = Int64(0),
    n_leaf = Int64(1),
    h_stem = FT(0),
    h_leaf = FT(9.5), # map, z_top
    h_canopy = h_stem + h_leaf,
    z0_m = FT(0.13) * h_canopy,
    z0_b = FT(0.1) * z0_m,
) where {FT}
    soil_α_PAR = FT(0.2) # function of clm_soil_albedo_parameters, soil moisture, can be calculated within
    soil_α_NIR = FT(0.2) # function of clm_soil_albedo_parameters, soil moisture
    # ://github.com/CliMA/ClimaLand.jl/blob/b3fce30e5e9b4e3024f0bbfd13e2e9654b541296/src/standalone/Soil/soil_albedo.jl#L193
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
function get_time_offset(longitude::Float64)
    hours = round(longitude / 15)

    if (hours > 11)
        hours = 11 - hours
    end

    return hours
end

"""
    get_pft(lat, long)

Given latitude and longitude, returns the dominant PFT at the coordinate.
"""
function get_pft(lat, long) # use ClimaUtilities, space-varying input,
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

function get_soil_params_gupta(lat, long; depth = 6) # soil_composition_parameters
    # note: latitude and longitudes values are consistent across maps

    # get soil porosity
    # long x lat x depth
    soil_porosity_map_path = joinpath(
        ClimaLand.Artifacts.soil_params_artifact_folder_path(),
        "porosity_map_gupta_etal2020_1.0x1.0x4.nc",
    )
    soil_porosity_map = NCDataset(soil_porosity_map_path)

    porosities = soil_porosity_map["ν"]
    lats = soil_porosity_map["lat"][:]
    lons = soil_porosity_map["lon"][:]

    lat_idx = findmin(abs.(lats .- lat))[2]
    lon_idx = findmin(abs.(lons .- long))[2]
    soil_ν = porosities[lat_idx, lon_idx, depth]

    # get soil ksat
    # long x lat x depth
    soil_ksat_map_path = joinpath(
        ClimaLand.Artifacts.soil_params_artifact_folder_path(),
        "ksat_map_gupta_etal2020_1.0x1.0x4.nc",
    )
    soil_ksat_map = NCDataset(soil_ksat_map_path)

    ksats = soil_ksat_map["Ksat"]
    soil_K_sat = k_sats[lat_idx, lon_idx, depth]

    # get soil residual
    # long x lat x depth
    soil_residual_map_path = joinpath(
        ClimaLand.Artifacts.soil_params_artifact_folder_path(),
        "residual_map_gupta_etal2020_1.0x1.0x4.nc",
    )
    soil_residual_map = NCDataset(soil_residual_map_path)

    residuals = soil_residual_map["θ_r"]
    θ_r = residuals[lat_idx, lon_idx, depth]

    # get soil vg α
    # long x lat x depth
    soil_vg_α_map_path = joinpath(
        ClimaLand.Artifacts.soil_params_artifact_folder_path(),
        "vGalpha_map_gupta_etal2020_1.0x1.0x4.nc",
    )
    soil_vg_α_map = NCDataset(soil_vg_α_map_path)

    soil_vg_α_values = soil_vg_α_map["α"]
    soil_vg_α = soil_vg_α_values[lat_idx, lon_idx, depth]

    # get soil vg n
    # long x lat x depth
    soil_vg_n_map_path = joinpath(
        ClimaLand.Artifacts.soil_params_artifact_folder_path(),
        "vGn_map_gupta_etal2020_1.0x1.0x4.nc",
    )
    soil_vg_n_map = NCDataset(soil_vg_n_map_path)

    soil_vg_n_values = soil_vg_n_map["n"]
    soil_vg_n = soil_vg_n_values[lat_idx, lon_idx, depth]

    return (; soil_ν, soil_K_sat, θ_r, soil_vg_α, soil_vg_n)
end

function get_soil_grid_params(lat, long)
    soil_grid_map_path = ClimaLand.Artifacts.soil_grids_params_artifact_path()
end
