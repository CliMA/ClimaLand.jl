using ClimaComms
using ClimaCore
using Interpolations
import ClimaUtilities.Regridders: InterpolationsRegridder
import ClimaUtilities.SpaceVaryingInputs: SpaceVaryingInput
import ClimaUtilities.ClimaArtifacts: @clima_artifact
import ClimaLand: use_lowres_clm, Artifacts

regridder_type = :InterpolationsRegridder
extrapolation_bc =
    (Interpolations.Periodic(), Interpolations.Flat(), Interpolations.Flat())

"""
    clm_canopy_radiation_parameters(
        surface_space;
        regridder_type = :InterpolationsRegridder,
        extrapolation_bc = (
            Interpolations.Periodic(),
            Interpolations.Flat(),
            Interpolations.Flat(),
        ),
        lowres=use_lowres_clm(surface_space),
    )

Reads spatially varying parameters for the canopy radiative transfer schemes,
 from NetCDF files based on CLM and MODIS data, and regrids them to the grid defined by the
`surface_space` of the Clima simulation. Returns a NamedTuple of ClimaCore
Fields.

In particular, this file returns a field for
- clumping index Ω
- albedo and transmissitivy in PAR and NIR bands
- leaf angle distribution G function parameter χl

The values correspond to the value of the dominant PFT at each point.

The NetCDF files are stored in ClimaArtifacts and more detail on their origin
is provided there. The keyword arguments `regridder_type` and `extrapolation_bc`
affect the regridding by (1) changing how we interpolate to ClimaCore points which
are not in the data, and (2) changing how extrapolate to points beyond the range of the
data. The keyword argument lowres is a flag that determines if the 0.9x1.25 or 0.125x0.125
resolution CLM data artifact is used. If the lowres flag is not provided, the clm artifact
with the closest resolution to the surface_space is used.
"""
function clm_canopy_radiation_parameters(
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
    clm_artifact_path = Artifacts.clm_data_folder_path(; context, lowres)
    # Foliage clumping index data derived from MODIS
    modis_ci_artifact_path = Artifacts.modis_ci_data_folder_path(; context)

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
    return (;
        Ω = Ω,
        G_Function = G_Function,
        α_PAR_leaf = α_PAR_leaf,
        τ_PAR_leaf = τ_PAR_leaf,
        α_NIR_leaf = α_NIR_leaf,
        τ_NIR_leaf = τ_NIR_leaf,
    )
end

"""
    clm_photosynthesis_parameters(
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
based on CLM data, and regrids them to the grid defined by the
`surface_space` of the Clima simulation. Returns a NamedTuple of ClimaCore
Fields.

In particular, this file returns a field for
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
function clm_photosynthesis_parameters(
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
    clm_artifact_path = Artifacts.clm_data_folder_path(; context, lowres)
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
    return (; is_c3 = is_c3, Vcmax25 = Vcmax25)
end


"""
    clm_rooting_depth(
        surface_space;
        regridder_type = :InterpolationsRegridder,
        extrapolation_bc = (
            Interpolations.Periodic(),
            Interpolations.Flat(),
            Interpolations.Flat(),
        ),
        lowres=use_lowres_clm(surface_space),
    )

Reads spatially varying rooting depth for the canopy, from a NetCDF file
based on CLM data, and regrids it to the grid defined by the
`surface_space` of the Clima simulation. Returns a NamedTuple of ClimaCore
Fields.

The values correspond to the value of the dominant PFT at each point.

The NetCDF files are stored in ClimaArtifacts and more detail on their origin
is provided there. The keyword arguments `regridder_type` and `extrapolation_bc`
affect the regridding by (1) changing how we interpolate to ClimaCore points which
are not in the data, and (2) changing how extrapolate to points beyond the range of the
data. The keyword argument lowres is a flag that determines if the 0.9x1.25 or 0.125x0.125
resolution CLM data artifact is used. If the lowres flag is not provided, the clm artifact
with the closest resolution to the surface_space is used.
"""
function clm_rooting_depth(
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
    clm_artifact_path = Artifacts.clm_data_folder_path(; context, lowres)
    rooting_depth = SpaceVaryingInput(
        joinpath(clm_artifact_path, "vegetation_properties_map.nc"),
        "rooting_depth",
        surface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc,),
    )
    return rooting_depth
end


"""
    clm_medlyn_g1(
        surface_space;
        regridder_type = :InterpolationsRegridder,
        extrapolation_bc = (
            Interpolations.Periodic(),
            Interpolations.Flat(),
            Interpolations.Flat(),
        ),
        lowres=use_lowres_clm(surface_space),
    )

Reads spatially varying g1 for the canopy, from a NetCDF file
based on CLM data, and regrids it to the grid defined by the
`surface_space` of the Clima simulation. Returns a NamedTuple of ClimaCore
Fields.

The values correspond to the value of the dominant PFT at each point.

The NetCDF files are stored in ClimaArtifacts and more detail on their origin
is provided there. The keyword arguments `regridder_type` and `extrapolation_bc`
affect the regridding by (1) changing how we interpolate to ClimaCore points which
are not in the data, and (2) changing how extrapolate to points beyond the range of the
data. The keyword argument lowres is a flag that determines if the 0.9x1.25 or 0.125x0.125
resolution CLM data artifact is used. If the lowres flag is not provided, the clm artifact
with the closest resolution to the surface_space is used.
"""
function clm_medlyn_g1(
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
    clm_artifact_path = Artifacts.clm_data_folder_path(; context, lowres)

    # g1 is read in units of sqrt(kPa) and then converted to sqrt(Pa)
    g1 = SpaceVaryingInput(
        joinpath(clm_artifact_path, "vegetation_properties_map.nc"),
        "medlynslope",
        surface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc,),
        file_reader_kwargs = (; preprocess_func = (data) -> data * 10^(3 / 2),),
    )
    return g1
end
