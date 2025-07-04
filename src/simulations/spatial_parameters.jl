export clm_canopy_parameters
function mean(x::AbstractArray{T}) where {T}
    sum(x) / length(x)
end

"""
    use_lowres_clm(space)

Returns true if the space is closer in resolution to 0.9x1.25 degree lat/long than 0.125x0.125 degree lat/long.
This is used to determine which resolution of CLM data to use.
"""
function use_lowres_clm(
    surface_space::ClimaCore.Spaces.AbstractSpectralElementSpace,
)
    node_scale = ClimaCore.Spaces.node_horizontal_length_scale(surface_space)
    surface_mesh = ClimaCore.Spaces.topology(surface_space).mesh
    if surface_mesh isa ClimaCore.Meshes.AbstractCubedSphere
        # in this case, node_scale is in meters
        sphere_radius = surface_mesh.domain.radius
        horizontal_length_scale(lat_res, long_res) = sqrt(
            4 * pi * sphere_radius^2 / ((360 / long_res) * (180 / lat_res)),
        )
        highres_scale = horizontal_length_scale(0.125, 0.125)
        lowres_scale = horizontal_length_scale(0.9, 1.25)
    elseif surface_mesh isa ClimaCore.Meshes.RectilinearMesh
        # in this case, node_scale is in degrees
        highres_scale = 0.125
        lowres_scale = sqrt(1.25 * 0.9)
    else
        return false
    end
    return abs(lowres_scale - node_scale) < abs(highres_scale - node_scale)
end

# This method will likely never be called, but is included for completeness
use_lowres_clm(surface_space::ClimaCore.Spaces.AbstractSpace) = false

"""
    clm_canopy_parameters(
        surface_space;
        regridder_type = :InterpolationsRegridder,
        extrapolation_bc = (
            Interpolations.Periodic(),
            Interpolations.Flat(),
            Interpolations.Flat(),
        ),
        lowres=use_lowres_clm(surface_space),
    )

Reads spatially varying parameters for the canopy, from NetCDF files
based on CLM and MODIS data, and regrids them to the grid defined by the
`surface_space` of the Clima simulation. Returns a NamedTuple of ClimaCore
Fields.

In particular, this file returns a field for
- clumping index Ω
- albedo and transmissitivy in PAR and NIR bands
- leaf angle distribution G function parameter χl
- Medlyn g1
- C3 flag
- VCmax25

The values correspond to the value of the dominant PFT at each point.

The NetCDF files are stored in ClimaArtifacts and more detail on their origin
is provided there. The keyword arguments `regridder_type` and `extrapolation_bc`
affect the regridding by (1) changing how we interpolate to ClimaCore points which
are not in the data, and (2) changing how extrapolate to points beyond the range of the
data. The keyword argument lowres is a flag that determines if the 0.9x1.25 or 0.125x0.125
resolution CLM data artifact is used. If the lowres flag is not provided, the clm artifact
with the closest resolution to the surface_space is used.
"""
function clm_canopy_parameters(
    surface_space;
    regridder_type = :InterpolationsRegridder,
    extrapolation_bc = (
        Interpolations.Periodic(),
        Interpolations.Flat(),
        Interpolations.Flat(),
    ),
    lowres = use_lowres_clm(surface_space),
)
    context = ClimaComms.context(surface_space)
    clm_artifact_path =
        ClimaLand.Artifacts.clm_data_folder_path(; context, lowres)
    # Foliage clumping index data derived from MODIS
    modis_ci_artifact_path =
        ClimaLand.Artifacts.modis_ci_data_folder_path(; context)

    # TwoStreamModel parameters
    nans_to_one(x) = isnan(x) ? eltype(x)(1) : x
    Ω = SpaceVaryingInput(
        joinpath(modis_ci_artifact_path, "He_et_al_2012_1x1.nc"),
        "ci",
        surface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc,),
        file_reader_kwargs = (; preprocess_func = nans_to_one,),
    )
    χl = SpaceVaryingInput(
        joinpath(clm_artifact_path, "vegetation_properties_map.nc"),
        "xl",
        surface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc,),
    )
    G_Function = ClimaLand.Canopy.CLMGFunction.(χl)
    α_PAR_leaf = SpaceVaryingInput(
        joinpath(clm_artifact_path, "vegetation_properties_map.nc"),
        "rholvis",
        surface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc,),
    )
    τ_PAR_leaf = SpaceVaryingInput(
        joinpath(clm_artifact_path, "vegetation_properties_map.nc"),
        "taulvis",
        surface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc,),
    )
    α_NIR_leaf = SpaceVaryingInput(
        joinpath(clm_artifact_path, "vegetation_properties_map.nc"),
        "rholnir",
        surface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc,),
    )
    τ_NIR_leaf = SpaceVaryingInput(
        joinpath(clm_artifact_path, "vegetation_properties_map.nc"),
        "taulnir",
        surface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc,),
    )
    # Conductance Model
    # g1 is read in units of sqrt(kPa) and then converted to sqrt(Pa)
    g1 = SpaceVaryingInput(
        joinpath(clm_artifact_path, "vegetation_properties_map.nc"),
        "medlynslope",
        surface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc,),
        file_reader_kwargs = (; preprocess_func = (data) -> data * 10^(3 / 2),),
    )

    #Photosynthesis model
    # vcmax is read in units of umol CO2/m^2/s and then converted to mol CO2/m^2/s
    Vcmax25 = SpaceVaryingInput(
        joinpath(clm_artifact_path, "vegetation_properties_map.nc"),
        "vcmx25",
        surface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc,),
        file_reader_kwargs = (; preprocess_func = (data) -> data / 1_000_000,),
    )
    # photosynthesis mechanism is read as a float, where 1.0 indicates c3 and 0.0 c4
    is_c3 = SpaceVaryingInput(
        joinpath(clm_artifact_path, "vegetation_properties_map.nc"),
        "c3_dominant",
        surface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc,),
    )
    rooting_depth = SpaceVaryingInput(
        joinpath(clm_artifact_path, "vegetation_properties_map.nc"),
        "rooting_depth",
        surface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc,),
    )
    return (;
        Ω = Ω,
        rooting_depth = rooting_depth,
        is_c3 = is_c3,
        Vcmax25 = Vcmax25,
        g1 = g1,
        G_Function = G_Function,
        α_PAR_leaf = α_PAR_leaf,
        τ_PAR_leaf = τ_PAR_leaf,
        α_NIR_leaf = α_NIR_leaf,
        τ_NIR_leaf = τ_NIR_leaf,
    )
end

"""
    default_spatially_varying_soil_parameters(
        surface_space,
        subsurface_space,
        regridder_type = :InterpolationsRegridder,
        extrapolation_bc = (
            Interpolations.Periodic(),
            Interpolations.Flat(),
            Interpolations.Flat(),
        ),
        lowres=use_lowres_clm(surface_space),
    )

Reads spatially varying parameters for the soil model, from NetCDF files
based on SoilGrids and the van Genuchten data product from Gupta et al 2020,
 and regrids them to the grid defined by the
`subsurface_space` and `surface_space` of the Clima simulation, as appropriate.
Returns a NamedTuple of ClimaCore Fields. The albedo parameters are read from
CLM data.

In particular, this file returns a field for
- (α, n, m) (van Genuchten parameters)
- Ksat
- porosity
- residual water content
- various texture variables: volumetric fractions of organic matter, coarse fragments, and quartz.
- TOPMODEL parameters
- soil albedo parameters
- specific storativity.

The NetCDF files are stored in ClimaArtifacts and more detail on their origin
is provided there. The keyword arguments `regridder_type` and `extrapolation_bc`
affect the regridding by (1) changing how we interpolate to ClimaCore points which
are not in the data, and (2) changing how extrapolate to points beyond the range of the
data. The keyword argument lowres is a flag that determines if the 0.9x1.25 or 0.125x0.125
resolution CLM data artifact is used. If the lowres flag is not provided, the clm artifact
with the closest resolution to the surface_space is used.


Since these parameters are read from discretized data sets, 
they carry an inherent land/sea mask. This land/sea mask may not match the
underlying land sea mask of the simulation. While values over the ocean do
not matter, we need to ensure that values in the simulation are set to
something physical, even if they are not set in the data.
In the future, this should be handled by ClimaUtilities via extrpolation.
Here we set them manually.
"""
function default_spatially_varying_soil_parameters(
    subsurface_space,
    surface_space,
    FT;
    regridder_type = :InterpolationsRegridder,
    extrapolation_bc = (
        Interpolations.Periodic(),
        Interpolations.Flat(),
        Interpolations.Flat(),
    ),
    lowres = use_lowres_clm(surface_space),
)
    context = ClimaComms.context(surface_space)
    soil_params_artifact_path =
        ClimaLand.Artifacts.soil_params_artifact_folder_path(; context)
    vg_α = SpaceVaryingInput(
        joinpath(
            soil_params_artifact_path,
            "vGalpha_map_gupta_etal2020_1.0x1.0x4.nc",
        ),
        "α",
        subsurface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc,),
    )
    vg_n = SpaceVaryingInput(
        joinpath(
            soil_params_artifact_path,
            "vGn_map_gupta_etal2020_1.0x1.0x4.nc",
        ),
        "n",
        subsurface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc,),
    )
    # The soil parameters may not have been defined on the same land
    # sea mask. Here we make sure that soil parameters over land,
    # but not set in the soil data, are set to something reasonable
    soil_params_mask = SpaceVaryingInput(
        joinpath(
            soil_params_artifact_path,
            "vGalpha_map_gupta_etal2020_1.0x1.0x4.nc",
        ),
        "α",
        subsurface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc,),
        file_reader_kwargs = (; preprocess_func = (data) -> data > 0,),
    )
    # If the parameter mask =  0, set to physical value
    # (equal to the mean where we have data; the mean of α is in log space)
    # This function in applied in **simulation mask** aware
    # manner.
    # That is, we replace values in the simulation, but without data, with the mean
    # over the data.
    masked_to_value(field, mask, value) =
        mask == 1.0 ? field : eltype(field)(value)

    μ = FT(0.33)
    vg_α .= masked_to_value.(vg_α, soil_params_mask, 10.0^μ)
    μ = FT(1.74)
    vg_n .= masked_to_value.(vg_n, soil_params_mask, μ)

    vg_fields_to_hcm_field(α::FT, n::FT) where {FT} =
        ClimaLand.Soil.vanGenuchten{FT}(; @NamedTuple{α::FT, n::FT}((α, n))...)
    hydrology_cm = vg_fields_to_hcm_field.(vg_α, vg_n)

    θ_r = SpaceVaryingInput(
        joinpath(
            soil_params_artifact_path,
            "residual_map_gupta_etal2020_1.0x1.0x4.nc",
        ),
        "θ_r",
        subsurface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc,),
    )

    ν = SpaceVaryingInput(
        joinpath(
            soil_params_artifact_path,
            "porosity_map_gupta_etal2020_1.0x1.0x4.nc",
        ),
        "ν",
        subsurface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc,),
    )
    K_sat = SpaceVaryingInput(
        joinpath(
            soil_params_artifact_path,
            "ksat_map_gupta_etal2020_1.0x1.0x4.nc",
        ),
        "Ksat",
        subsurface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc,),
    )

    # Set missing values to the mean. For Ksat, we use the mean in log space.
    μ = FT(-5.08)
    K_sat .= masked_to_value.(K_sat, soil_params_mask, 10.0^μ)

    ν .= masked_to_value.(ν, soil_params_mask, 0.47)

    θ_r .= masked_to_value.(θ_r, soil_params_mask, 0.09)

    S_s = ClimaCore.Fields.zeros(subsurface_space) .+ FT(1e-3)

    soilgrids_artifact_path =
        ClimaLand.Artifacts.soil_grids_params_artifact_path(;
            lowres = true,
            context,
        )
    ν_ss_om = SpaceVaryingInput(
        soilgrids_artifact_path,
        "nu_ss_om",
        subsurface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc,),
    )

    ν_ss_quartz = SpaceVaryingInput(
        soilgrids_artifact_path,
        "nu_ss_sand",
        subsurface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc,),
    )

    ν_ss_gravel = SpaceVaryingInput(
        soilgrids_artifact_path,
        "nu_ss_cf",
        subsurface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc,),
    )
    # we require that the sum of these be less than 1 and equal to or bigger than zero.
    # The input should satisfy this almost exactly, but the regridded values may not.
    # Values of zero are OK here, so we dont need to apply any masking
    # if sum > 1, normalize by sum, else "normalize" by 1 (i.e., do not normalize)
    texture_norm = @. max(ν_ss_gravel + ν_ss_quartz + ν_ss_om, 1)
    @. ν_ss_gravel = ν_ss_gravel / texture_norm
    @. ν_ss_om = ν_ss_om / texture_norm
    @. ν_ss_quartz = ν_ss_quartz / texture_norm

    # Read in f_max data and topmodel params land sea mask
    infile_path = ClimaLand.Artifacts.topmodel_data_path()
    f_max =
        SpaceVaryingInput(infile_path, "fmax", surface_space; regridder_type)
    topmodel_params_mask = SpaceVaryingInput(
        infile_path,
        "landsea_mask",
        surface_space;
        regridder_type,
    )
    soil_params_mask_sfc =
        ClimaLand.Domains.top_center_to_surface(soil_params_mask)
    f_max = masked_to_value.(f_max, topmodel_params_mask, FT(0))
    f_max = masked_to_value.(f_max, soil_params_mask_sfc, FT(0))
    PAR_albedo_dry, NIR_albedo_dry, PAR_albedo_wet, NIR_albedo_wet = map(
        s -> SpaceVaryingInput(
            joinpath(
                ClimaLand.Artifacts.clm_data_folder_path(; context, lowres),
                "soil_properties_map.nc",
            ),
            s,
            surface_space;
            regridder_type,
            regridder_kwargs = (; extrapolation_bc,),
        ),
        (
            "PAR_albedo_dry",
            "NIR_albedo_dry",
            "PAR_albedo_wet",
            "NIR_albedo_wet",
        ),
    )

    return (;
        ν = ν,
        ν_ss_om = ν_ss_om,
        ν_ss_quartz = ν_ss_quartz,
        ν_ss_gravel = ν_ss_gravel,
        hydrology_cm = hydrology_cm,
        K_sat = K_sat,
        S_s = S_s,
        θ_r = θ_r,
        PAR_albedo_wet = PAR_albedo_wet,
        NIR_albedo_wet = NIR_albedo_wet,
        PAR_albedo_dry = PAR_albedo_dry,
        NIR_albedo_dry = NIR_albedo_dry,
        f_max = f_max,
    )
end
