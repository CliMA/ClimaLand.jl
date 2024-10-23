using ClimaComms
using ClimaUtilities.ClimaArtifacts
import Interpolations
import ClimaUtilities.SpaceVaryingInputs: SpaceVaryingInput
import ClimaUtilities.Regridders: InterpolationsRegridder
export clm_canopy_parameters
function mean(x::AbstractArray{T}) where {T}
    sum(x) / length(x)
end

"""
    clm_canopy_parameters(
        surface_space;
        regridder_type = :InterpolationsRegridder,
        extrapolation_bc = (
            Interpolations.Periodic(),
            Interpolations.Flat(),
            Interpolations.Flat(),
        ),
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
data.
"""
function clm_canopy_parameters(
    surface_space;
    regridder_type = :InterpolationsRegridder,
    extrapolation_bc = (
        Interpolations.Periodic(),
        Interpolations.Flat(),
        Interpolations.Flat(),
    ),
)
    context = ClimaComms.context(surface_space)
    clm_artifact_path = ClimaLand.Artifacts.clm_data_folder_path(; context)
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
    G_Function = CLMGFunction(χl)
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
        α_NIR_leaf = α_PAR_leaf,
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
    )

Reads spatially varying parameters for the soil model, from NetCDF files
based on SoilGrids and the van Genuchten data product from Gupta et al 2020,
 and regrids them to the grid defined by the
`subsurface_space` and `surface_space` of the Clima simulation, as appropriate. 
Returns a NamedTuple of ClimaCore Fields. 

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
data.
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
)
    context = ClimaComms.context(surface_space)
    soil_params_artifact_path =
        ClimaLand.Artifacts.soil_params_artifact_folder_path(; context)
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
    oceans_to_value(field, mask, value) =
        mask == 1.0 ? field : eltype(field)(value)

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
    x = parent(vg_α)
    μ = mean(log10.(x[x .> 0]))
    vg_α .= oceans_to_value.(vg_α, soil_params_mask, 10.0^μ)

    x = parent(vg_n)
    μ = mean(x[x .> 0])
    vg_n .= oceans_to_value.(vg_n, soil_params_mask, μ)

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

    x = parent(K_sat)
    μ = mean(log10.(x[x .> 0]))
    K_sat .= oceans_to_value.(K_sat, soil_params_mask, 10.0^μ)

    ν .= oceans_to_value.(ν, soil_params_mask, 1)

    θ_r .= oceans_to_value.(θ_r, soil_params_mask, 0)


    S_s =
        oceans_to_value.(
            ClimaCore.Fields.zeros(subsurface_space) .+ FT(1e-3),
            soil_params_mask,
            1,
        )
    ν_ss_om =
        oceans_to_value.(
            ClimaCore.Fields.zeros(subsurface_space) .+ FT(0.1),
            soil_params_mask,
            0,
        )
    ν_ss_quartz =
        oceans_to_value.(
            ClimaCore.Fields.zeros(subsurface_space) .+ FT(0.1),
            soil_params_mask,
            0,
        )
    ν_ss_gravel =
        oceans_to_value.(
            ClimaCore.Fields.zeros(subsurface_space) .+ FT(0.1),
            soil_params_mask,
            0,
        )
    PAR_albedo = ClimaCore.Fields.zeros(surface_space) .+ FT(0.2)
    NIR_albedo = ClimaCore.Fields.zeros(surface_space) .+ FT(0.2)

    soil_params_mask_sfc =
        ClimaLand.Domains.top_center_to_surface(soil_params_mask)

    # Read in f_max data and land sea mask
    infile_path = ClimaLand.Artifacts.topmodel_data_path()
    f_max =
        SpaceVaryingInput(infile_path, "fmax", surface_space; regridder_type)
    mask = SpaceVaryingInput(
        infile_path,
        "landsea_mask",
        surface_space;
        regridder_type,
    )
    # Unsure how to handle two masks
    f_max = oceans_to_value.(f_max, mask, FT(0.0))
    f_max = oceans_to_value.(f_max, soil_params_mask_sfc, FT(0.0))
    return (;
        ν = ν,
        ν_ss_om = ν_ss_om,
        ν_ss_quartz = ν_ss_quartz,
        ν_ss_gravel = ν_ss_gravel,
        hydrology_cm = hydrology_cm,
        K_sat = K_sat,
        S_s = S_s,
        θ_r = θ_r,
        PAR_albedo = PAR_albedo,
        NIR_albedo = NIR_albedo,
        f_max = f_max,
        mask = mask,
    )
end
