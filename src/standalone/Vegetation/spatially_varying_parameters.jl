using ClimaComms
using ClimaCore
import Interpolations
import ClimaUtilities.Regridders: InterpolationsRegridder
import ClimaUtilities.SpaceVaryingInputs: SpaceVaryingInput
import ClimaUtilities.ClimaArtifacts: @clima_artifact
import ClimaLand: Artifacts

"""
    clm_canopy_radiation_parameters(
        surface_space;
        regridder_type = :InterpolationsRegridder,
        extrapolation_bc = (
            Interpolations.Periodic(),
            Interpolations.Flat(),
            Interpolations.Flat(),
        ),
        interpolation_method = Interpolations.Linear(),
        lowres = ClimaLand.Domains.use_lowres_clm(surface_space),
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
is provided there. The keyword arguments `regridder_type`, `extrapolation_bc`, and
`regridder_kwargs`
affect the regridding by (1) changing how we interpolate to ClimaCore points which
are not in the data, and (2) changing how extrapolate to points beyond the range of the
data, and (3) changed the spatial interpolation method.

The keyword argument lowres is a flag that determines if the 0.9x1.25 or 0.125x0.125
resolution CLM data artifact is used. If the lowres flag is not provided, the clm artifact
with the closest resolution to the surface_space is used.

By default linear interpolation is used. This can be changed to nearest neighbor by passing
`interpolation_method = Interpolations.Constant()`, but linear interpolation is recommended.
"""
function clm_canopy_radiation_parameters(
    surface_space;
    regridder_type = :InterpolationsRegridder,
    extrapolation_bc = (
        Interpolations.Periodic(),
        Interpolations.Flat(),
        Interpolations.Flat(),
    ),
    interpolation_method = Interpolations.Linear(),
    lowres = ClimaLand.Domains.use_lowres_clm(surface_space),
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
        regridder_kwargs = (; extrapolation_bc, interpolation_method),
        file_reader_kwargs = (; preprocess_func = nans_to_one,),
    )
    χl = SpaceVaryingInput(
        joinpath(clm_artifact_path, "vegetation_properties_map.nc"),
        "xl",
        surface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc, interpolation_method),
    )
    G_Function = ClimaLand.Canopy.CLMGFunction.(χl)
    α_PAR_leaf = SpaceVaryingInput(
        joinpath(clm_artifact_path, "vegetation_properties_map.nc"),
        "rholvis",
        surface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc, interpolation_method),
    )
    τ_PAR_leaf = SpaceVaryingInput(
        joinpath(clm_artifact_path, "vegetation_properties_map.nc"),
        "taulvis",
        surface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc, interpolation_method),
    )
    α_NIR_leaf = SpaceVaryingInput(
        joinpath(clm_artifact_path, "vegetation_properties_map.nc"),
        "rholnir",
        surface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc, interpolation_method),
    )
    τ_NIR_leaf = SpaceVaryingInput(
        joinpath(clm_artifact_path, "vegetation_properties_map.nc"),
        "taulnir",
        surface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc, interpolation_method),
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
        interpolation_method = Interpolations.Linear(),
        lowres = ClimaLand.Domains.use_lowres_clm(surface_space),
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
is provided there. The keyword arguments `regridder_type`, `extrapolation_bc`, and
`regridder_kwargs`
affect the regridding by (1) changing how we interpolate to ClimaCore points which
are not in the data, and (2) changing how extrapolate to points beyond the range of the
data, and (3) changed the spatial interpolation method.

The keyword argument lowres is a flag that determines if the 0.9x1.25 or 0.125x0.125
resolution CLM data artifact is used. If the lowres flag is not provided, the clm artifact
with the closest resolution to the surface_space is used.

By default linear interpolation is used. This can be changed to nearest neighbor by passing
`interpolation_method = Interpolations.Constant()`, but linear interpolation is recommended.
"""
function clm_photosynthesis_parameters(
    surface_space;
    regridder_type = :InterpolationsRegridder,
    extrapolation_bc = (
        Interpolations.Periodic(),
        Interpolations.Flat(),
        Interpolations.Flat(),
    ),
    interpolation_method = Interpolations.Linear(),
    lowres = ClimaLand.Domains.use_lowres_clm(surface_space),
)
    context = ClimaComms.context(surface_space)
    clm_artifact_path = Artifacts.clm_data_folder_path(; context, lowres)
    # vcmax is read in units of umol CO2/m^2/s and then converted to mol CO2/m^2/s
    Vcmax25 = SpaceVaryingInput(
        joinpath(clm_artifact_path, "vegetation_properties_map.nc"),
        "vcmx25",
        surface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc, interpolation_method),
        file_reader_kwargs = (; preprocess_func = (data) -> data / 1_000_000,),
    )
    # photosynthesis mechanism is read as a float, where 1.0 indicates c3 and 0.0 c4
    is_c3 = SpaceVaryingInput(
        joinpath(clm_artifact_path, "vegetation_properties_map.nc"),
        "c3_dominant",
        surface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc, interpolation_method),
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
        interpolation_method = Interpolations.Linear(),
        lowres = ClimaLand.Domains.use_lowres_clm(surface_space),
    )

Reads spatially varying rooting depth for the canopy, from a NetCDF file
based on CLM data, and regrids it to the grid defined by the
`surface_space` of the Clima simulation. Returns a NamedTuple of ClimaCore
Fields.

The values correspond to the value of the dominant PFT at each point.

The NetCDF files are stored in ClimaArtifacts and more detail on their origin
is provided there. The keyword arguments `regridder_type`, `extrapolation_bc`, and
`regridder_kwargs`
affect the regridding by (1) changing how we interpolate to ClimaCore points which
are not in the data, and (2) changing how extrapolate to points beyond the range of the
data, and (3) changed the spatial interpolation method.

The keyword argument lowres is a flag that determines if the 0.9x1.25 or 0.125x0.125
resolution CLM data artifact is used. If the lowres flag is not provided, the clm artifact
with the closest resolution to the surface_space is used.

By default linear interpolation is used. This can be changed to nearest neighbor by passing
`interpolation_method = Interpolations.Constant()`, but linear interpolation is recommended.
"""
function clm_rooting_depth(
    surface_space;
    regridder_type = :InterpolationsRegridder,
    extrapolation_bc = (
        Interpolations.Periodic(),
        Interpolations.Flat(),
        Interpolations.Flat(),
    ),
    interpolation_method = Interpolations.Linear(),
    lowres = ClimaLand.Domains.use_lowres_clm(surface_space),
)
    context = ClimaComms.context(surface_space)
    clm_artifact_path = Artifacts.clm_data_folder_path(; context, lowres)
    rooting_depth = SpaceVaryingInput(
        joinpath(clm_artifact_path, "vegetation_properties_map.nc"),
        "rooting_depth",
        surface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc, interpolation_method),
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
        interpolation_method = Interpolations.Linear(),
        lowres = ClimaLand.Domains.use_lowres_clm(surface_space),
    )

Reads spatially varying g1 for the canopy, from a NetCDF file
based on CLM data, and regrids it to the grid defined by the
`surface_space` of the Clima simulation. Returns a NamedTuple of ClimaCore
Fields.

The values correspond to the value of the dominant PFT at each point.

The NetCDF files are stored in ClimaArtifacts and more detail on their origin
is provided there. The keyword arguments `regridder_type`, `extrapolation_bc`, and
`regridder_kwargs`
affect the regridding by (1) changing how we interpolate to ClimaCore points which
are not in the data, and (2) changing how extrapolate to points beyond the range of the
data, and (3) changed the spatial interpolation method.

The keyword argument lowres is a flag that determines if the 0.9x1.25 or 0.125x0.125
resolution CLM data artifact is used. If the lowres flag is not provided, the clm artifact
with the closest resolution to the surface_space is used.

By default linear interpolation is used. This can be changed to nearest neighbor by passing
`interpolation_method = Interpolations.Constant()`, but linear interpolation is recommended.
"""
function clm_medlyn_g1(
    surface_space;
    regridder_type = :InterpolationsRegridder,
    extrapolation_bc = (
        Interpolations.Periodic(),
        Interpolations.Flat(),
        Interpolations.Flat(),
    ),
    interpolation_method = Interpolations.Linear(),
    lowres = ClimaLand.Domains.use_lowres_clm(surface_space),
)
    context = ClimaComms.context(surface_space)
    clm_artifact_path = Artifacts.clm_data_folder_path(; context, lowres)

    # g1 is read in units of sqrt(kPa) and then converted to sqrt(Pa)
    g1 = SpaceVaryingInput(
        joinpath(clm_artifact_path, "vegetation_properties_map.nc"),
        "medlynslope",
        surface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc, interpolation_method),
        file_reader_kwargs = (; preprocess_func = (data) -> data * 10^(3 / 2),),
    )
    return g1
end

export clm_canopy_height, effective_canopy_height

"""
    clm_canopy_height(
        surface_space;
        regridder_type = :InterpolationsRegridder,
        extrapolation_bc = (
            Interpolations.Periodic(),
            Interpolations.Flat(),
            Interpolations.Flat(),
        ),
        interpolation_method = Interpolations.Linear(),
        lowres = ClimaLand.Domains.use_lowres_clm(surface_space),
    )

Read spatially-varying canopy height data from CLM vegetation properties.

Returns a Field of raw canopy heights from the CLM dataset without any modifications.
For heights that respect atmospheric coupling constraints, use `effective_canopy_height`
to cap the heights appropriately.

# Arguments
- `surface_space`: The surface space for regridding

# Optional Arguments
- `regridder_type`: Type of regridder (default: `:InterpolationsRegridder`)
- `extrapolation_bc`: Boundary conditions for extrapolation
- `interpolation_method`: Spatial interpolation method (default: Linear)
- `lowres`: Use low-resolution CLM data

# Returns
- Field of canopy heights in meters

# Example
```julia
raw_height = clm_canopy_height(surface_space)
effective_height = effective_canopy_height(raw_height, 10.0; buffer=2.0)
```
"""
function clm_canopy_height(
    surface_space;
    regridder_type = :InterpolationsRegridder,
    extrapolation_bc = (
        Interpolations.Periodic(),
        Interpolations.Flat(),
        Interpolations.Flat(),
    ),
    interpolation_method = Interpolations.Linear(),
    lowres = ClimaLand.Domains.use_lowres_clm(surface_space),
)
    context = ClimaComms.context(surface_space)
    clm_artifact_path = Artifacts.clm_data_folder_path(; context, lowres)

    canopy_height = SpaceVaryingInput(
        joinpath(clm_artifact_path, "vegetation_properties_map.nc"),
        "z_top",
        surface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc, interpolation_method),
    )
    return canopy_height
end

"""
    effective_canopy_height(
        canopy_height::ClimaCore.Fields.Field,
        z_atm::FT;
        buffer::FT = FT(2.0)
    ) where {FT}

Caps canopy heights to ensure they remain below the atmospheric reference height.
This is necessary when atmospheric forcing data (e.g., ERA5) is provided at a fixed
reference height (typically 10m). The canopy height must be less than z_atm to ensure
proper flux calculations and atmospheric coupling.

# Arguments
- `canopy_height`: Field of canopy heights read from CLM data (meters)
- `z_atm`: Atmospheric reference height for forcing data (meters)
- `buffer`: Minimum distance between canopy top and atmospheric reference (default: 2m)

# Returns
- Field of effective canopy heights, capped at `z_atm - buffer`

# Example
```julia
raw_height = clm_canopy_height(surface_space)
z_atm = FT(10.0)  # ERA5 forcing at 10m
height = effective_canopy_height(raw_height, z_atm)
```

# Notes
The function will log a warning message indicating how many grid cells had their
canopy height capped and provide statistics about the original height distribution.
This typically affects ~7% of global cells, primarily in regions with tall forests
(e.g., tropical rainforests, temperate forests in Patagonia and New Zealand).
"""
function effective_canopy_height(
    canopy_height::ClimaCore.Fields.Field,
    z_atm::FT;
    buffer::FT = FT(2.0),
) where {FT}
    max_height = z_atm - buffer

    # Compute statistics on CPU to avoid GPU Boolean Field issues
    max_original = maximum(canopy_height)

    # Only warn if we actually have heights that need capping
    if max_original >= max_height
        # Count how many cells exceed the threshold by summing a float-converted mask
        # This avoids creating Boolean Fields on GPU
        n_capped_field = @. ifelse(canopy_height >= max_height, FT(1), FT(0))
        n_capped_sum = sum(n_capped_field)
        n_total = length(canopy_height)
        pct_capped = 100.0 * n_capped_sum / n_total
        @warn "Capping canopy heights: $(round(n_capped_sum))/$n_total cells ($(round(pct_capped, digits=2))%) exceed max_height=$max_height m. Original max height: $(round(max_original, digits=2)) m. This is expected for tall forests when using atmospheric forcing at z_atm=$(z_atm) m."
    end

    # Apply cap: min(height, max_height) at each point
    return min.(canopy_height, max_height)
end
