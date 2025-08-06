using ClimaComms
using ClimaCore
using Interpolations
import ClimaUtilities.Regridders: InterpolationsRegridder
import ClimaUtilities.SpaceVaryingInputs: SpaceVaryingInput
import ClimaUtilities.ClimaArtifacts: @clima_artifact
using ClimaLand
import ClimaLand: use_lowres_clm, Artifacts

masked_to_value(field, mask, value) = mask == 1.0 ? field : eltype(field)(value)

"""
    clm_soil_albedo_parameters(
        surface_space,
        regridder_type = :InterpolationsRegridder,
        extrapolation_bc = (
            Interpolations.Periodic(),
            Interpolations.Flat(),
            Interpolations.Flat(),
        ),
        interpolation_method = Interpolations.Constant(),
        lowres=use_lowres_clm(surface_space),
    )

Reads spatially varying albedo parameters for the soil model, from NetCDF files
of CLM data for the PAR and NIR albedo of wet and dry soil.

The NetCDF files are stored in ClimaArtifacts and more detail on their origin
is provided there. The keyword arguments `regridder_type`, `extrapolation_bc`, and
`interpolation_method` 
affect the regridding by (1) changing how we interpolate to ClimaCore points which
are not in the data, and (2) changing how extrapolate to points beyond the range of the
data, and (3) changed the spatial interpolation method. 
The keyword argument lowres is a flag that determines if the 0.9x1.25 or 0.125x0.125
resolution CLM data artifact is used. If the lowres flag is not provided, the clm artifact
with the closest resolution to the surface_space is used.

Since these parameters are read from discretized data sets, 
they carry an inherent land/sea mask. This land/sea mask may not match the
underlying land sea mask of the simulation. 
"""
function clm_soil_albedo_parameters(
    surface_space;
    regridder_type = :InterpolationsRegridder,
    extrapolation_bc = (
        Interpolations.Periodic(),
        Interpolations.Flat(),
        Interpolations.Flat(),
    ),
    interpolation_method = Interpolations.Constant(),
    lowres = use_lowres_clm(surface_space),
)
    context = ClimaComms.context(surface_space)
    PAR_albedo_dry, NIR_albedo_dry, PAR_albedo_wet, NIR_albedo_wet = map(
        s -> SpaceVaryingInput(
            joinpath(
                Artifacts.clm_data_folder_path(; context, lowres),
                "soil_properties_map.nc",
            ),
            s,
            surface_space;
            regridder_type,
            regridder_kwargs = (; extrapolation_bc, interpolation_method),
        ),
        (
            "PAR_albedo_dry",
            "NIR_albedo_dry",
            "PAR_albedo_wet",
            "NIR_albedo_wet",
        ),
    )
    return (;
        PAR_albedo_wet = PAR_albedo_wet,
        NIR_albedo_wet = NIR_albedo_wet,
        PAR_albedo_dry = PAR_albedo_dry,
        NIR_albedo_dry = NIR_albedo_dry,
    )
end

"""
    soil_vangenuchten_parameters(
        surface_space,
        subsurface_space,
        regridder_type = :InterpolationsRegridder,
        extrapolation_bc = (
            Interpolations.Periodic(),
            Interpolations.Flat(),
            Interpolations.Flat(),
        ),
       interpolation_method = Interpolations.Linear(),
    )

Reads spatially varying van Genuchten parameters for the soil model, from NetCDF files
based the van Genuchten data product from Gupta et al 2020,
 and regrids them to the grid defined by the
`subsurface_space` and `surface_space` of the Clima simulation, as appropriate.
Returns a NamedTuple of ClimaCore Fields. 

In particular, this file returns a field for
- (α, n, m) (van Genuchten parameters)
- Ksat
- porosity
- residual water content

The NetCDF files are stored in ClimaArtifacts and more detail on their origin
is provided there. The keyword arguments `regridder_type`, `extrapolation_bc`, and
`interpolation_method` 
affect the regridding by (1) changing how we interpolate to ClimaCore points which
are not in the data, and (2) changing how extrapolate to points beyond the range of the
data, and (3) changed the spatial interpolation method.

Since these parameters are read from discretized data sets, 
they carry an inherent land/sea mask. This land/sea mask may not match the
underlying land sea mask of the simulation. While values over the ocean do
not matter, we need to ensure that values in the simulation are set to
something physical, even if they are not set in the data.
In the future, this should be handled by ClimaUtilities via extrpolation.
Here we set them manually.
"""
function soil_vangenuchten_parameters(
    subsurface_space,
    FT;
    regridder_type = :InterpolationsRegridder,
    extrapolation_bc = (
        Interpolations.Periodic(),
        Interpolations.Flat(),
        Interpolations.Flat(),
    ),
    interpolation_method = Interpolations.Linear(),
)
    context = ClimaComms.context(subsurface_space)
    soil_params_artifact_path =
        Artifacts.soil_params_artifact_folder_path(; context)
    vg_α = SpaceVaryingInput(
        joinpath(
            soil_params_artifact_path,
            "vGalpha_map_gupta_etal2020_1.0x1.0x4.nc",
        ),
        "α",
        subsurface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc, interpolation_method),
    )
    vg_n = SpaceVaryingInput(
        joinpath(
            soil_params_artifact_path,
            "vGn_map_gupta_etal2020_1.0x1.0x4.nc",
        ),
        "n",
        subsurface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc, interpolation_method),
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
        regridder_kwargs = (; extrapolation_bc, interpolation_method),
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
        regridder_kwargs = (; extrapolation_bc, interpolation_method),
    )

    ν = SpaceVaryingInput(
        joinpath(
            soil_params_artifact_path,
            "porosity_map_gupta_etal2020_1.0x1.0x4.nc",
        ),
        "ν",
        subsurface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc, interpolation_method),
    )
    K_sat = SpaceVaryingInput(
        joinpath(
            soil_params_artifact_path,
            "ksat_map_gupta_etal2020_1.0x1.0x4.nc",
        ),
        "Ksat",
        subsurface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc, interpolation_method),
    )

    # Set missing values to the mean. For Ksat, we use the mean in log space.
    μ = FT(-5.08)
    K_sat .= masked_to_value.(K_sat, soil_params_mask, 10.0^μ)
    K_sat .= max.(K_sat, sqrt(eps(FT)))

    ν .= masked_to_value.(ν, soil_params_mask, 0.47)

    θ_r .= masked_to_value.(θ_r, soil_params_mask, 0.09)

    return (; ν = ν, hydrology_cm = hydrology_cm, K_sat = K_sat, θ_r = θ_r)
end

"""
    soil_composition_parameters(
        surface_space,
        subsurface_space,
        regridder_type = :InterpolationsRegridder,
        extrapolation_bc = (
            Interpolations.Periodic(),
            Interpolations.Flat(),
            Interpolations.Flat(),
        ),
        interpolation_method = Interpolations.Linear(),
        path = Artifacts.soil_grids_params_artifact_path(;
                                                                   lowres = true,
                                                                   ClimaComms.context(subsurface_space),
                                                                   )
    )

Reads spatially varying parameters for the soil model, from NetCDF files
based on SoilGrids,
 and regrids them to the grid defined by the
`subsurface_space` and `surface_space` of the Clima simulation, as appropriate.
Returns a NamedTuple of ClimaCore Fields. 

In particular, this file returns a field for
- various texture variables: volumetric fractions of organic matter, coarse fragments, and quartz.

The NetCDF files are stored in ClimaArtifacts and more detail on their origin
is provided there. The keyword arguments `regridder_type`, `extrapolation_bc`, and
`interpolation_method` 
affect the regridding by (1) changing how we interpolate to ClimaCore points which
are not in the data, and (2) changing how extrapolate to points beyond the range of the
data, and (3) changed the spatial interpolation method.

Any `path` can be provided, but this assumes that the parameters are stored under the names
`nu_ss_om`, `nu_ss_sand`, and `nu_ss_cf`.
"""
function soil_composition_parameters(
    subsurface_space,
    FT;
    regridder_type = :InterpolationsRegridder,
    extrapolation_bc = (
        Interpolations.Periodic(),
        Interpolations.Flat(),
        Interpolations.Flat(),
    ),
    interpolation_method = Interpolations.Linear(),
    path = Artifacts.soil_grids_params_artifact_path(;
        lowres = true,
        context = ClimaComms.context(subsurface_space),
    ),
)
    ν_ss_om = SpaceVaryingInput(
        path,
        "nu_ss_om",
        subsurface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc, interpolation_method),
    )

    ν_ss_quartz = SpaceVaryingInput(
        path,
        "nu_ss_sand",
        subsurface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc, interpolation_method),
    )

    ν_ss_gravel = SpaceVaryingInput(
        path,
        "nu_ss_cf",
        subsurface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc, interpolation_method),
    )
    # we require that the sum of these be less than 1 and equal to or bigger than zero.
    # The input should satisfy this almost exactly, but the regridded values may not (if linear
    # interpolation is used).
    # Values of zero are OK here, so we dont need to apply any masking
    # if sum > 1, normalize by sum, else "normalize" by 1 (i.e., do not normalize)
    texture_norm = @. max(ν_ss_gravel + ν_ss_quartz + ν_ss_om, 1)
    @. ν_ss_gravel = ν_ss_gravel / texture_norm
    @. ν_ss_om = ν_ss_om / texture_norm
    @. ν_ss_quartz = ν_ss_quartz / texture_norm

    return (;
        ν_ss_om = ν_ss_om,
        ν_ss_quartz = ν_ss_quartz,
        ν_ss_gravel = ν_ss_gravel,
    )
end

"""
    topmodel_fmax(
        surface_space,
        FT;
        regridder_type = :InterpolationsRegridder,
        extrapolation_bc = (
            Interpolations.Periodic(),
            Interpolations.Flat(),
            Interpolations.Flat(),
        ),
        interpolation_method = Interpolations.Linear(),
    )

Returns the spatially varying parameter `fmax` of the TOPMODEL
parameterization; this is read from a nc file and then regridded
to the simulation grid.

The keyword arguments `regridder_type`, `extrapolation_bc`, and
`interpolation_method` 
affect the regridding by (1) changing how we interpolate to ClimaCore points which
are not in the data, and (2) changing how extrapolate to points beyond the range of the
data, and (3) changed the spatial interpolation method.
"""
function topmodel_fmax(
    surface_space,
    FT;
    regridder_type = :InterpolationsRegridder,
    extrapolation_bc = (
        Interpolations.Periodic(),
        Interpolations.Flat(),
        Interpolations.Flat(),
    ),
    interpolation_method = Interpolations.Linear(),
)
    # Read in f_max data and topmodel data land sea mask
    infile_path = Artifacts.topmodel_data_path()
    f_max =
        SpaceVaryingInput(infile_path, "fmax", surface_space; regridder_type)
    topmodel_params_mask = SpaceVaryingInput(
        infile_path,
        "landsea_mask",
        surface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc, interpolation_method),
    )
    f_max = masked_to_value.(f_max, topmodel_params_mask, FT(0))
    return f_max
end
