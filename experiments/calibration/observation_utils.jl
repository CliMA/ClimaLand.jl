import ClimaLand
import ClimaCore
import ClimaAnalysis

"""
    get_lat_lon_from_resolution(nelements)

Return latitudes and longitude for the diagnostics of the land model run at a
given resolution (`nelements`).

It assumes that the output is on the default grid, as determined by
`ClimaLand.default_num_points`. It also assumes that the domain is the full
globe.
"""
function get_lat_lon_from_resolution(nelements)
    # Values of radius and depth (set to 0.1) do not matter; we just need
    # information from the domain in the latitude/longitude directions.
    # We do not want to load the CUDA backend or construct this on GPU which is
    # wasteful, so we construct this on CPU instead
    domain = ClimaLand.Domains.SphericalShell(;
        radius = 0.1,
        depth = 0.1,
        nelements,
        context = ClimaComms.SingletonCommsContext{ClimaComms.CPUSingleThreaded}(
            ClimaComms.CPUSingleThreaded(),
        ),
    )
    # If the default number of diagnostic latitude and longitude points is not
    # used when running the simulations, this will need to be changed.
    num_long, num_lat, _ =
        ClimaLand.Diagnostics.default_diagnostic_num_points(domain)
    Δ_long = 360.0 / num_long
    Δ_lat = 180.0 / num_lat
    longs = collect(range(-180.0; length = num_long, step = Δ_long))
    lats = collect(range(-90.0; length = num_lat, step = Δ_lat))
    return lats, longs
end

"""
    make_ocean_mask(nelements)

Make an ocean mask.

# Notes
- ClimaLand simulates only locations where its mask == 1 (based on a threshold
  set in the simulations). Since parameters are not defined over the ocean, this
  can lead to unphysical parameter combinations right at the coastline.
- The interpolated diagnostics are not mask aware; this again means that grid
  points treated as land by the simulation are combined with grid points that
  are not acted on by the simulation (typically set to zero, but possibly
  undefined). Again the result is that grid points on the coastline are not
  guaranteed to be physically realistic.
- Consequently, here we define the mask used to select calibration points using
  a more stringent criterion (essentially, that the grid point in the simulation
  is only comprised of land, and the point in the diagnostic is only comprised
  of land)
"""
function make_ocean_mask(nelements)
    lats, longs = get_lat_lon_from_resolution(nelements)

    # We do not want to load the CUDA backend or construct this on GPU which is
    # wasteful, so we construct this on CPU instead
    domain = ClimaLand.Domains.SphericalShell(;
        radius = 0.1,
        depth = 0.1,
        nelements,
        context = ClimaComms.SingletonCommsContext{ClimaComms.CPUSingleThreaded}(
            ClimaComms.CPUSingleThreaded(),
        ),
    )

    # We need the mask in order to determine which points were treated as land
    # by the simulation. We set the threshold of fractional area of land to be
    # 0.99 here - cells with > 1% ocean are ignored in calibration.
    mask = ClimaLand.Domains.landsea_mask(domain; threshold = 0.99)

    # Interpolate the ClimaLand land sea mask (Field) to the diagnostics grid
    # (matrix)
    target_hcoords = [
        ClimaCore.Geometry.LatLongPoint(lat, lon) for lat in lats, lon in longs
    ]
    # Select only locations where the interpolated mask is equal to 1
    interpolated_mask =
        Array(ClimaCore.Remapping.interpolate(mask; target_hcoords))

    # Make a OutputVar to make a mask from
    # TODO: This touches on the internals of ClimaAnalysis which should not be
    # done.
    attribs =
        Dict("short_name" => "landsea_mask", "long_name" => "Landsea mask")
    dim_attribs = Dict(
        "latitude" => Dict("units" => "degrees_north"),
        "longitude" => Dict("units" => "degrees_east"),
    )
    dims = Dict("latitude" => lats, "longitude" => longs)
    landsea_var =
        ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, interpolated_mask)

    # For ERA5 data, it shouldn't matter what the threshold is since we already
    # resampled onto the diagnostics grid
    return ClimaAnalysis.generate_lonlat_mask(
        landsea_var,
        NaN,
        1.0;
        threshold = 0.5,
    )
end
