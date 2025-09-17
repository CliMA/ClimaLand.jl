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
    units::String
end

# Map between short name in NetCDF files produced to ClimaDiagnostics to
# an ILAMBMapping
const ILAMB_VARIABLES = Dict(
    "gpp" => ILAMBMapping(
        "gpp",
        "Gross Primary Production",
        0.01201,
        "kg m-2 s-1",
    ),
    "et" =>
        ILAMBMapping("evspsbl", "Evapotranspiration", 1.0, "kg m-2 s-1"),
    "lhf" => ILAMBMapping("hfls", "Latent Heat Flux", 1.0, "W m-2"),
    "shf" => ILAMBMapping("hfss", "Sensible Heat Flux", 1.0, "W m-2"),
    "lwu" =>
        ILAMBMapping("rlus", "Upward Longwave Radiation", 1.0, "W m-2"),
    "lwd" =>
        ILAMBMapping("rlds", "Downward Longwave Radiation", 1.0, "W m-2"),
    "swd" =>
        ILAMBMapping("rsds", "Downward Shortwave Radiation", 1.0, "W m-2"),
    "swu" =>
        ILAMBMapping("rsus", "Upward Shortwave Radiation", 1.0, "W m-2"),
    "lai" => ILAMBMapping("lai", "Leaf Area Index", 1.0, "1"),
    # Assuming the native range is 0–100 %
    "snowc" => ILAMBMapping("scf", "Snow Cover Fraction", 1.0, "1"),
    "sr" => ILAMBMapping(
        "mrro",
        "Surface Runoff",
        1000.0 / 2.63e6, # 1 month ≈ 2,629,440 seconds = 2.63 × 10e6 seconds
        "kg m-2 s-1",
    ),
    "ssr" =>
        ILAMBMapping("mrros", "Subsurface Runoff", 1000.0, "kg m-2 s-1"),
    "swe" => ILAMBMapping("swe", "Snow Water Equivalent", 1000.0, "kg m-2"),
    "iwc" =>
        ILAMBMapping("mrsos", "Soil Moisture (top layer)", 1.0, "kg m-2"),
    "msf" => ILAMBMapping("msf", "Moisture Stress Factor", 1.0, "1"),
    "trans" => ILAMBMapping("tran", "Transpiration", 1000.0, "kg m-2 s-1"),
    "swc" => ILAMBMapping(
        "swc",
        "Soil Water Content (top 10cm)",
        1000.0,
        "kg m-2",
    ),
    "tsoil" => ILAMBMapping("tsl", "Soil Temperature (top 10cm)", 1.0, "K"),
    "precip" => ILAMBMapping("pr", "Precipitation", -1.0, "kg m-2 s-1"),
    "tair" => ILAMBMapping("tas", "Air Temperature", 1.0, "K"),
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

    # Shift dates by one month and add 14 days, since ClimaDiagnostics save
    # dates at the end of the period instead of at the beginning and ILAMB wants
    # the days to be in the middle of the months
    dates = clima_ds["date"] .- Dates.Month(1) .+ Dates.Day(14)
    start_date = first(dates)
    end_date = last(dates)

    ilamb_file_name =
        make_ilamb_file_name(clima_short_name, start_date, end_date)
    NCDatasets.NCDataset(joinpath(save_dir, ilamb_file_name), "c") do ilamb_ds
        add_temporal_dimension(ilamb_ds, clima_ds, dates)
        add_nontemporal_dimensions(ilamb_ds, clima_ds)
        define_var(ilamb_ds, clima_ds, clima_short_name)
    end
    close(clima_ds)
end

"""
    add_nontemporal_dimensions(ilamb_ds, clima_ds)

Add the nontemporal dimensions to `ilamb_ds` using `clima_ds`, where both are
open NetCDF files.

The nontemporal dimensions are longitude and latitude.
"""
function add_nontemporal_dimensions(ilamb_ds, clima_ds)
    non_temporal_dimnames = ["lon", "lat"]
    for dimname in non_temporal_dimnames
        ilamb_ds.dim[dimname] = clima_ds.dim[dimname]
        NCDatasets.defVar(
            ilamb_ds,
            dimname,
            Array(clima_ds[dimname]),
            NCDatasets.dimnames(clima_ds[dimname]),
            attrib = clima_ds[dimname].attrib,
        )
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
    (; short_name, long_name, conversion_factor, units) = ilamb_mapping
    var_attribs = Dict(
        "_FillValue" => FT(1.0e20),
        "standard_name" => short_name,
        "long_name" => long_name,
        "units" => ilamb_units,
    )
    # For ILAMB, the order of dimensions should be lon, lat, and time
    expected_dims = ["lon", "lat", "time"]
    data = get_data(
        clima_ds,
        clima_short_name,
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
    get_data(clima_ds, clima_short_name, conversion_factor, expected_dims)

Get the data from `clima_ds`, an open NetCDF file, with the name
`clima_short_name`.

Additional preprocessing includes converting to units that ILAMB expects with
`conversion_factor` and permuting the dimensions with `expected_dims`.
"""
function get_data(
    clima_ds,
    clima_short_name,
    conversion_factor,
    expected_dims,
    FT,
)
    data = Array(clima_ds[clima_short_name])

    # Special cases for tsoil and swc to extract top 10 cm (z >= -0.1)
    if clima_short_name in ("tsoil", "swc")
        z = Array(clima_ds["z"])
        all_dimnames = NCDatasets.dimnames(clima_ds[clima_short_name])
        top_layer_indices = findall(val -> val >= -0.1, z)
        length(top_layer_indices) >= 1 ||
            error("No valid depths above -0.1m found for $(clima_short_name)")
        z_idx = last(top_layer_indices)

        index_tuple =
            (dimname == "z" ? z_idx : Colon() for dimname in all_dimnames)
        data = data[index_tuple...]
    end

    # Permute dimensions to match ILAMB
    dimnames = filter(
        dimname -> dimname in expected_dims,
        NCDatasets.dimnames(clima_ds[clima_short_name]),
    )
    perm = indexin(expected_dims, collect(dimnames))
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
