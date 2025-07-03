using ClimaAnalysis
using ClimaLand
using ClimaLand.Artifacts
using ClimaComms
using ClimaCore
ClimaComms.@import_required_backends


"""
    diagnostics_lat_lon(nelements)

Return latitudes and longitude for the diagnostics of the land
model run at a given resolution (`nelements`).

It assumes that the output is on the default grid, as determined by
`ClimaLand.default_num_points`. It also assumes that the domain is the full
globe.
"""
function diagnostics_lat_lon(nelements)
    cpu_comms_ctx = ClimaComms.context(ClimaComms.CPUSingleThreaded())
    domain = ClimaLand.Domains.SphericalShell(;
        radius = 0.1,
        depth = 0.1,
        nelements,
        comms_ctx = cpu_comms_ctx,
    ) # values of radius and depth (set to 0.1) do not matter; we just need information from the domain in the latitude/longitude directions.
    # If the default number of diagnostic latitude and longitude points is not used when running the simulations, this will need to be changed.
    num_long, num_lat, _ =
        ClimaLand.Diagnostics.default_diagnostic_num_points(domain)
    longs = collect(range(-180.0, 180.0, length = num_long))
    lats = collect(range(-90.0, 90.0, length = num_lat))
    return lats, longs
end


"""
    make_training_locations(nelements)

Create a list of geographic training locations (longitude, latitude) for model calibration.

# Notes
- ClimaLand simulates only locations where its mask == 1 (based on a threshold set in the simulations). Since parameters are not defined over the ocean, this can lead to unphysical parameter combinations right at the coastline.
- The interpolated diagnostics are not mask aware; this again means that grid points treated as land by the simulation are combined with grid points that are not acted on by the simulation (typically set to zero, but possibly undefined). Again the result is that grid points on the coastline are not guaranteed to be physically realistic.
- Consequently, here we define the mask used to select calibration points using a more stringent criterion (essentially, that the grid point in the simulation is only comprised of land, and the point in the diagnostic is only comprised of land)
"""
function make_training_locations(nelements)
    lats, longs = diagnostics_lat_lon(nelements)
    cpu_comms_ctx = ClimaComms.context(ClimaComms.CPUSingleThreaded())
    domain = ClimaLand.Domains.SphericalShell(;
        radius = 0.1,
        depth = 0.1,
        nelements,
        comms_ctx = cpu_comms_ctx,
    )

    # We need the mask in order to determine which points were treated as land by the simulation.
    # We set the threshold of fractional area of land to be 0.99 here - cells with > 1% ocean are ignored in calibration.
    mask = ClimaLand.landsea_mask(domain; threshold = 0.99)

    # Interpolate the ClimaLand land sea mask (Field) to the diagnostics grid (matrix)
    target_hcoords = [
        ClimaCore.Geometry.LatLongPoint(lat, lon) for lat in lats, lon in longs
    ]
    # Select only locations where the interpolated mask is equal to 1
    interpolated_mask =
        Array(ClimaCore.Remapping.interpolate(mask; target_hcoords))
	# skip the exact pole at -90, 90 in latitude.
    training_locations = [
        (lon, lat) for (j, lon) in enumerate(longs) for
        (i, lat) in enumerate(lats[2:end-1]) if !iszero(interpolated_mask[i+1, j])
    ]
    @show training_locations

    return training_locations
end
