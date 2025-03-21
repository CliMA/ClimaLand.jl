using ClimaAnalysis
using ClimaLand
using ClimaLand.Artifacts
using ClimaComms
ClimaComms.@import_required_backends

"""
    make_training_locations(path_to_grid)

Creates a list of geographic training locations (longitude, latitude) for model calibration.

It assumes that the output is on the default grid, as determined by
`ClimaLand.default_num_points`. It also assumes that the domain is the full
globe.

# Notes
- Uses ERA5 monthly averaged surface data from 1979-2024
- The function applies a land mask to identify suitable training locations
"""
function make_training_locations(nelements)
    radius = 0.1 # These don't matter
    depth = 0.1 # These don't matter

    domain = ClimaLand.Domains.SphericalShell(; radius, depth, nelements)
    num_long, num_lat, _ = ClimaLand.Diagnostics.default_diagnostic_num_points(domain)
    longs = collect(range(-180., 180., length=num_long))
    lats = collect(range(-90., 90., length=num_lat))

    var = ClimaAnalysis.OutputVar(
        Dict("long" => longs,
             "lat" => lats),
        zeros(num_lat, num_long)
    )

    # Apply_oceanmask applies the mask over the ocean
    vars_in_land = ClimaAnalysis.apply_oceanmask(var)

    training_locations = [
        (lon, lat) for
        (j, lon) in enumerate(ClimaAnalysis.longitudes(vars_in_land)) for
        (i, lat) in enumerate(ClimaAnalysis.latitudes(vars_in_land)) if
        !isnan(vars_in_land.data[i, j])
    ]
    return training_locations
end
