# The purpose of this script is to produce NetCDF files compatible with ILAMB from
# the diagnostics produced by the land simulation.
#
# To run this script, you can run
# ```
# julia --project=.buildkite experiments/ilamb/ilamb_conversion.jl <output_dir> <save_dir>
# ```
# in the terminal.
# Or, directly call `make_compatible_with_ILAMB(output_dir, save_dir)`
# at the end of a simulation script after including this file.
#
# For more information about ILAMB, see the ILAMB and ILAMB-Data repo:
# https://github.com/rubisco-sfa/ILAMB
# https://github.com/rubisco-sfa/ILAMB-Data
#
# For more informaiton about the variable naming for ILAMB, see
# https://clipc-services.ceda.ac.uk/dreq/index/var.html

import ClimaAnalysis
import NCDatasets
import Dates

"""
    struct ILAMBMapping

Keep track of the short name, long name, conversion factor, and units for
making a NetCDF file produced by ClimaDiagnostics compatible with ILAMB.
"""
struct ILAMBMapping
    """Short name of the variable in the NetCDF file for ILAMB"""
    short_name::String

    """Long name of the variable in the NetCDF file for ILAMB"""
    long_name::String

    """Conversion factor from CliMA to ILAMB"""
    conversion_factor::Float64

    """Units of variable after applying conversion factor"""
    ilamb_units::String

    """Units of variable before conversion"""
    clima_units::String
end

# Map between short name in NetCDF files produced to ClimaDiagnostics to
# an ILAMBMapping
const ILAMB_VARIABLES = Dict(
    "gpp" => ILAMBMapping(
        "gpp",
        "Gross Primary Production",
        0.01201,
        "kg m-2 s-1", # kg C
        "mol CO2 m^-2 s^-1",
    ),
    "nee" => ILAMBMapping(
        "nee",
        "Net Ecosystem Exchange",
        0.01201,
        "kg m-2 s-1", # kg C
        "mol CO2 m^-2 s^-1",
    ),
    "er" => ILAMBMapping(
        "reco",
        "Ecosystem Respiration",
        0.01201,
        "kg m-2 s-1", # kg C
        "mol CO2 m^-2 s^-1",
    ),
    "ra" => ILAMBMapping(
        "ra",
        "Autotrophic Respiration",
        0.01201,
        "kg m-2 s-1", # kg C
        "mol CO2 m^-2 s^-1",
    ),
    "hr" => ILAMBMapping(
        "rh",
        "Heterotrophic Respiration",
        0.01201,
        "kg m-2 s-1", # kg C
        "mol CO2 m^-2 s^-1",
    ),
    "soc" => ILAMBMapping(
        "cSoil",
        "Soil Carbon",
        1.0,
        "kg m-2", # kg C
        "kg C m^-3",
    ),
    "cveg" => ILAMBMapping(
        "cVeg",
        "Vegetation Carbon",
        1.0,
        "kg m-2", # kg C
        "kg C m^-2",
    ),
    "et" => ILAMBMapping(
        "evspsbl",
        "Evapotranspiration",
        1.0,
        "kg m-2 s-1",
        "kg m^-2 s^-1",
    ),
    "lhf" =>
        ILAMBMapping("hfls", "Latent Heat Flux", 1.0, "W m-2", "W m^-2"),
    "shf" =>
        ILAMBMapping("hfss", "Sensible Heat Flux", 1.0, "W m-2", "W m^-2"),
    "lwu" => ILAMBMapping(
        "rlus",
        "Upward Longwave Radiation",
        1.0,
        "W m-2",
        "W m^-2",
    ),
    "lwd" => ILAMBMapping(
        "rlds",
        "Downward Longwave Radiation",
        1.0,
        "W m-2",
        "W m^-2",
    ),
    "swd" => ILAMBMapping(
        "rsds",
        "Downward Shortwave Radiation",
        1.0,
        "W m-2",
        "W m^-2",
    ),
    "swu" => ILAMBMapping(
        "rsus",
        "Upward Shortwave Radiation",
        1.0,
        "W m-2",
        "W m^-2",
    ),
    "lai" => ILAMBMapping("lai", "Leaf Area Index", 1.0, "1", "m^2 m^-2"),
    # Assuming the native range is 0–100 %
    "snowc" => ILAMBMapping("scf", "Snow Cover Fraction", 1.0, "1", ""),
    # Fixed: CF conventions state mrros is surface runoff, not mrro
    "sr" => ILAMBMapping(
        "mrros",
        "Surface Runoff",
        1000.0,
        "kg m-2 s-1",
        "m s^-1",
    ),
    # Fixed: Subsurface runoff uses ssro (no standard CF convention name)
    "ssr" => ILAMBMapping(
        "ssro",
        "Subsurface Runoff",
        1000.0,
        "kg m-2 s-1",
        "m s^-1",
    ),
    # Total runoff (surface + subsurface)
    "tr" => ILAMBMapping(
        "mrro",
        "Total Runoff",
        1000.0,
        "kg m-2 s-1",
        "m s^-1",
    ),
    "swe" =>
        ILAMBMapping("swe", "Snow Water Equivalent", 1000.0, "kg m-2", "m"),
    "iwc" => ILAMBMapping(
        "mrsos",
        "Soil Moisture (integrated over top 10cm)",
        1.0,
        "kg m-2",
        "kg/m^2",
    ),
    "msf" => ILAMBMapping("msf", "Moisture Stress Factor", 1.0, "1", ""),
    "trans" => ILAMBMapping(
        "tran",
        "Transpiration",
        1000.0,
        "kg m-2 s-1",
        "m s^-1",
    ),
    "precip" => ILAMBMapping(
        "pr",
        "Precipitation (includes snow)", # check if this meant to include both snow and rain
        -1.0,
        "kg m-2 s-1",
        "kg m^-2 s^-1",
    ),
    "tair" => ILAMBMapping("tas", "Air Temperature", 1.0, "K", "K"),
    "tsoil" => ILAMBMapping("tsl", "Temperature of Soil", 1.0, "K", "K"),
)

"""
    find_netcdf_files_for_ilamb(output_dir)

Given `output_dir`, a directory with NetCDF files produced by CliMA diagnostics,
find all the paths of the NetCDF files that can be converted to ILAMB and return
them as a vector of strings.

This function assumes the NetCDF files that have monthly averages and end with
"_1M_average.nc".
"""
function find_netcdf_files_for_ilamb(output_dir)
    all_filepaths = readdir(output_dir, join = true)
    filter!(filepath -> endswith(filepath, "_1M_average.nc"), all_filepaths)
    clima_short_names = keys(ILAMB_VARIABLES)
    filter!(
        filepath -> find_short_name(filepath) in clima_short_names,
        all_filepaths,
    )
    return all_filepaths
end

"""
    make_ilamb_netcdf_file(filepath, save_dir)

Given a `filepath` to a NetCDF file produced by ClimaDiagnostics, create a
NetCDF file that is compatible with ILAMB.
"""
function make_ilamb_netcdf_file(filepath, save_dir)
    clima_ds = NCDatasets.NCDataset(filepath)
    clima_short_name = find_short_name(filepath)

    # Add 14 days, since ClimaDiagnostics save dates at the start of the period
    # and ILAMB wants the days to be in the middle of the months
    dates = clima_ds["date"] .+ Dates.Day(14)
    start_date = first(dates)
    end_date = last(dates)

    ilamb_file_name =
        make_ilamb_file_name(clima_short_name, start_date, end_date)
    NCDatasets.NCDataset(joinpath(save_dir, ilamb_file_name), "c") do ilamb_ds
        add_temporal_dimension(ilamb_ds, clima_ds, dates)
        add_nontemporal_dimensions(ilamb_ds, clima_ds, clima_short_name)
        define_var(ilamb_ds, clima_ds, clima_short_name)
    end
    close(clima_ds)
end

"""
    add_nontemporal_dimensions(ilamb_ds, clima_ds, clima_short_name)

Add the nontemporal dimensions to `ilamb_ds` using `clima_ds`, where both are
open NetCDF files.

The nontemporal dimensions are longitude, latitude, and depth (for depth-dependent variables).
"""
function add_nontemporal_dimensions(ilamb_ds, clima_ds, clima_short_name)
    non_temporal_dimnames = ["lon", "lat"]
    # Add depth dimension for tsoil (map z -> depth for ILAMB/CMIP6 conventions)
    if clima_short_name == "tsoil" && haskey(clima_ds, "z")
        push!(non_temporal_dimnames, "z")
    end
    for dimname in non_temporal_dimnames
        # Map z coordinate to depth for ILAMB/CMIP6 conventions
        output_dimname = (dimname == "z") ? "depth" : dimname
        ilamb_ds.dim[output_dimname] = clima_ds.dim[dimname]

        # For depth, convert from negative z (depth below surface) to positive depth
        if dimname == "z"
            # ClimaLand z goes deepest-first (most negative) to shallowest (least negative)
            # Reverse so depth is ascending from surface (shallowest first) per CMIP6 convention
            depth_values = reverse(-Array(clima_ds[dimname]))
            attrib = Dict(
                "standard_name" => "depth",
                "long_name" => "depth",
                "units" => "m",
                "positive" => "down",
                "axis" => "Z",
                "bounds" => "depth_bnds",
            )
            NCDatasets.defVar(
                ilamb_ds,
                output_dimname,
                depth_values,
                (output_dimname,),
                attrib = attrib,
            )

            # Compute CF-compliant depth_bnds (layer edges, positive downward)
            # Per CF conventions: depth_bnds[i,0] = top of layer i, depth_bnds[i,1] = bottom
            # The top of the first soil layer must be exactly 0.0 m (the surface).
            # ILAMB's auto-generated bounds extrapolate and give a negative value for
            # the first layer's top bound, which is physically wrong.
            # Note: Julia NCDatasets uses column-major (Fortran) order, so dimensions are
            # reversed vs Python/NetCDF C-order.
            # Julia ("bnds", "depth") -> Python reads ("depth", "bnds") as required by CF.
            # Data must be (2, n) in Julia to be read as (n, 2) in Python.
            if !haskey(ilamb_ds.dim, "bnds")
                ilamb_ds.dim["bnds"] = 2
            end
            n = length(depth_values)
            depth_bnds = if haskey(clima_ds, "z_bnds")
                # Use source bounds if available: z_bnds is (n, 2) in Julia
                # Negate (z->depth), transpose to (2, n), reverse depth axis
                z_bnds = Array(clima_ds["z_bnds"])
                reverse(-z_bnds', dims = 2)
            else
                # Compute from layer centers. CF rule: top of first layer = 0.0 m (surface).
                # Store as (2, n): row 1 = top edges, row 2 = bottom edges
                bnds = zeros(eltype(depth_values), 2, n)
                bnds[1, 1] = 0.0  # Top of first layer = ground surface
                for i in 1:(n - 1)
                    midpoint = 0.5 * (depth_values[i] + depth_values[i + 1])
                    bnds[2, i] = midpoint      # Bottom of layer i
                    bnds[1, i + 1] = midpoint  # Top of layer i+1
                end
                # Extrapolate bottom of last layer
                bnds[2, n] =
                    depth_values[n] +
                    0.5 * (depth_values[n] - depth_values[n - 1])
                bnds
            end
            NCDatasets.defVar(
                ilamb_ds,
                "depth_bnds",
                depth_bnds,
                ("bnds", "depth"),
                attrib = Dict("units" => "m"),
            )
        else
            NCDatasets.defVar(
                ilamb_ds,
                output_dimname,
                Array(clima_ds[dimname]),
                NCDatasets.dimnames(clima_ds[dimname]),
                attrib = clima_ds[dimname].attrib,
            )
        end
    end
    return nothing
end

"""
    add_temporal_dimension(ilamb_ds, clima_ds, dates)

Add the temporal dimensions to `ilamb_ds` using `clima_ds`, where both are
open NetCDF files.

The temporal dimensions are time and time bounds.
"""
function add_temporal_dimension(ilamb_ds, clima_ds, dates)
    # Define length of dimensions
    ilamb_ds.dim["time"] = clima_ds.dim["time"]
    ilamb_ds.dim["bnds"] = 2
    NCDatasets.defVar(
        ilamb_ds,
        "time",
        dates,
        ("time",),
        attrib = Dict(
            "standard_name" => "time",
            "long_name" => "time",
            "units" => "days since 1850-01-01",
            "calendar" => "proleptic_gregorian",
            "axis" => "T",
            "bounds" => "time_bnds",
        ),
    )
    NCDatasets.defVar(
        ilamb_ds,
        "time_bnds",
        # Even though the dates need to be shifted by one month, the
        # date/time bounds are correct
        Array(clima_ds["date_bnds"]),
        ("bnds", "time"),
        attrib = Dict("units" => "days since 1850-01-01"),
    )
    return nothing
end

"""
    define_var(ilamb_ds, clima_ds, clima_short_name)

Define the variable with the name `clima_short_name` in `ilamb_ds` using
`clima_ds`, where both are open NetCDF files.
"""
function define_var(ilamb_ds, clima_ds, clima_short_name)
    # Use Float32 instead of Float64 to save space
    FT = Float32
    ilamb_mapping = ILAMB_VARIABLES[clima_short_name]
    (; short_name, long_name, conversion_factor, ilamb_units, clima_units) =
        ilamb_mapping

    # Use CF standard names where they differ from variable names
    cf_standard_name = if short_name == "tsl"
        "soil_temperature"
    else
        short_name
    end

    var_attribs = Dict(
        "_FillValue" => FT(1.0e20),
        "standard_name" => cf_standard_name,
        "long_name" => long_name,
        "units" => ilamb_units,
    )
    # For ILAMB, the order of dimensions should be lon, lat, and time
    # Note: Julia's NCDatasets uses column-major (Fortran) order, so dimension
    # ordering is reversed compared to Python/NetCDF C-order convention.
    # Julia ("lon", "lat", "time") -> Python reads ('time', 'lat', 'lon')
    # For tsoil: Julia ("lon", "lat", "depth", "time") -> Python reads ('time', 'depth', 'lat', 'lon')
    if clima_short_name == "tsoil"
        expected_dims = ["lon", "lat", "depth", "time"]
    else
        expected_dims = ["lon", "lat", "time"]
    end
    data = get_data(
        clima_ds,
        clima_short_name,
        clima_units,
        conversion_factor,
        expected_dims,
        FT,
    )
    NCDatasets.defVar(
        ilamb_ds,
        short_name,
        data,
        Tuple(expected_dims),
        attrib = var_attribs,
    )
    return nothing
end

"""
    get_data(
        clima_ds,
        clima_short_name,
        clima_units,
        conversion_factor,
        expected_dims
    )

Get the data from `clima_ds`, an open NetCDF file, with the name
`clima_short_name`.

Additional preprocessing includes converting to units that ILAMB expects with
`conversion_factor` and permuting the dimensions with `expected_dims`.
"""
function get_data(
    clima_ds,
    clima_short_name,
    clima_units,
    conversion_factor,
    expected_dims,
    FT,
)
    data = Array(clima_ds[clima_short_name])

    # For tsoil, reverse the depth (z) dimension so depth is ascending from surface
    if clima_short_name == "tsoil" && haskey(clima_ds, "z")
        all_dimnames = NCDatasets.dimnames(clima_ds[clima_short_name])
        z_dim_idx = findfirst(==("z"), all_dimnames)
        if !isnothing(z_dim_idx)
            data = reverse(data, dims = z_dim_idx)
        end
    end

    # Permute dimensions to match ILAMB
    # Map z -> depth for dimension matching
    dimnames_raw = NCDatasets.dimnames(clima_ds[clima_short_name])
    dimnames = [d == "z" ? "depth" : d for d in dimnames_raw]
    dimnames = filter(dimname -> dimname in expected_dims, dimnames)
    perm = indexin(expected_dims, collect(dimnames))
    clima_ds[clima_short_name].attrib["units"] != clima_units && error(
        "Units in the simulation-output NetCDF do not match the expected CliMA units for $clima_short_name. Update the conversion factor and units in the mapping and run this script again",
    )
    data .*= conversion_factor
    data = permutedims(data, perm)
    return FT.(data)
end

"""
    function make_ilamb_file_name(
        clima_short_name,
        start_date::Dates.DateTime,
        end_date::Dates.DateTime,
    )

Make the file name for the NetCDF file for use in the ILAMB leaderboard.

The file name will be
"[variable_name]_Lmon_CliMA_historical_r1i1p1f1_gn_[START_YEAR_MONTH]_[END_YEAR_MONTH].nc".

# Example
For example, for a simulation that started on `2000-03-01T00:00:00` and ended on
`2019-02-01T00:00:00` and variable name is `precip`, then the file name is
`pr_Lmon_CliMA_historical_r1i1p1f1_gn_200003-201902.nc`.
"""
function make_ilamb_file_name(
    clima_short_name,
    start_date::Dates.DateTime,
    end_date::Dates.DateTime,
)
    ilamb_short_name = ILAMB_VARIABLES[clima_short_name].short_name

    formatted_start_date = Dates.format(start_date, "YYYYmm")
    formatted_end_date = Dates.format(end_date, "YYYYmm")
    return "$(ilamb_short_name)_Lmon_CliMA_historical_r1i1p1f1_gn_$(formatted_start_date)-$(formatted_end_date).nc"
end

"""
    find_short_name(filepath)

Given a `filepath` to a NetCDF file produced by ClimaDiagnostics, return the
short name of the variable stored in the NetCDF file.
"""
function find_short_name(filepath)
    return first(ClimaAnalysis.Utils.match_nc_filename(basename(filepath)))
end

"""
    make_compatible_with_ILAMB(output_dir, save_dir)

Given the directory `output_dir` with NetCDF files produced by ClimaDiagnostics,
create NetCDF files compatible with ILAMB in the directory `save_dir`.
"""
function make_compatible_with_ILAMB(output_dir, save_dir)
    !ispath(save_dir) && mkpath(save_dir)
    nc_filepaths = find_netcdf_files_for_ilamb(output_dir)
    for nc_filepath in nc_filepaths
        make_ilamb_netcdf_file(nc_filepath, save_dir)
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) != 2
        error(
            "Usage: julia --project=.buildkite ilamb_conversion.jl <output_dir> <save_dir>",
        )
    end
    output_dir = ARGS[1]
    save_dir = ARGS[2]
    make_compatible_with_ILAMB(output_dir, save_dir)
end
