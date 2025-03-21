using ClimaAnalysis
using ClimaLand
using ClimaLand.Artifacts

"""
    make_training_locations(path_to_grid)

Creates a list of geographic training locations (longitude, latitude) for model calibration.

# Arguments
- `path_to_grid`: Path to the grid data used for resampling the ERA5 data

# Notes
- Uses ERA5 monthly averaged surface data from 1979-2024
- The function applies a land mask to identify suitable training locations
"""
function make_training_locations(path_to_grid)
    era5_path = joinpath(
        ClimaLand.Artifacts.era5_surface_data_1979_2024_path(),
        "era5_monthly_averages_surface_single_level_197901-202410.nc",
    )

    lhf_era5 = OutputVar(era5_path, "mslhf")

    simdir = SimDir(path_to_grid)

    lhf_cl = get(simdir, "lhf")

    land_mask_var = slice(
        ClimaAnalysis.apply_landmask(
            ClimaAnalysis.resampled_as(lhf_era5, lhf_cl),
        ),
        time = 0,
    )

    training_locations = [
        (lon, lat) for
        (i, lon) in enumerate(ClimaAnalysis.longitudes(land_mask_var)) for
        (j, lat) in enumerate(ClimaAnalysis.latitudes(land_mask_var)) if
        isnan(land_mask_var.data[i, j])
    ]
    return training_locations
end
