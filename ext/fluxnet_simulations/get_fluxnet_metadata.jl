###################################
#       FLUXNET2015 metadata      #
###################################
#
# These helpers read site metadata (lat/lon, time offset, sensor heights,
# canopy height) from the FLUXNET2015 metadata CSV bundled with the
# `fluxnet2015` artifact. The artifact is HPC-only by default — locally,
# call `generic_site_simulation` with explicit coordinate kwargs instead.

"""
    parse_array_field(field) -> Vector{Float64}

Parse a metadata cell that may hold a semicolon-separated list of floats
(e.g. multiple atmospheric sensor heights). Returns an empty vector for
empty/NaN cells.
"""
function parse_array_field(field)
    if field == "" || field == "NaN"
        return Float64[]
    elseif field isa Float64
        return isnan(field) ? Float64[] : [field]
    else
        return parse.(Float64, split(String(field), ";"))
    end
end

"""
    FluxnetSimulations.get_site_info(site_ID; fluxnet2015_metadata_path = nothing)

Look up site metadata for `site_ID` from the FLUXNET2015 `metadata_DD_clean.csv`.

Returns a NamedTuple:
- `lat::Float64` — latitude (deg)
- `long::Float64` — longitude (deg)
- `time_offset::Int` — local time offset from UTC (hours), using the ClimaLand
  convention `UTC = local_time + time_offset` (i.e. negated relative to the
  metadata's `utc_offset` column).
- `atmospheric_sensor_height::Vector{Float64}` — heights (m) of atmospheric
  sensors at the site; may be empty if the metadata cell is missing.

Missing fields produce a `@warn` and a `NaN`-valued entry rather than an error,
so a partial metadata row still mostly works.

The `fluxnet2015` artifact (containing the metadata CSV) is HPC-only — see
`ClimaLand.Artifacts.fluxnet2015_data_path`.
"""
function FluxnetSimulations.get_site_info(
    site_ID;
    fluxnet2015_metadata_path = nothing,
)
    if fluxnet2015_metadata_path === nothing
        fluxnet2015_metadata_path = joinpath(
            ClimaLand.Artifacts.fluxnet2015_data_path(),
            "metadata_DD_clean.csv",
        )
    end
    raw = readdlm(fluxnet2015_metadata_path, ',', Any)
    header = raw[1, :]
    data = raw[2:end, :]

    row_idx = findfirst(row -> row[1] == site_ID, eachrow(data))
    isnothing(row_idx) && error("Site ID $site_ID not found in metadata.")
    site_metadata = data[row_idx, :]

    varnames =
        ["latitude", "longitude", "utc_offset", "atmospheric_sensor_heights"]
    column_name_map = Dict(
        varname => findfirst(header[:] .== varname) for varname in varnames
    )

    nan_if_empty(x) =
        (x isa AbstractString && isempty(x)) ? NaN :
        (x isa AbstractFloat && isnan(x)) ? NaN : x

    for (varname, column) in column_name_map
        field = site_metadata[column]
        if (field isa AbstractString && isempty(field)) ||
           (field isa AbstractFloat && isnan(field))
            @warn "Field $(varname) is missing for site $site_ID"
        end
    end

    return (;
        lat = nan_if_empty(site_metadata[column_name_map["latitude"]]),
        long = nan_if_empty(site_metadata[column_name_map["longitude"]]),
        time_offset = Int(
            -1 * nan_if_empty(site_metadata[column_name_map["utc_offset"]]),
        ),
        atmospheric_sensor_height = parse_array_field(
            site_metadata[column_name_map["atmospheric_sensor_heights"]],
        ),
    )
end

"""
    FluxnetSimulations.get_canopy_height(site_ID; fluxnet2015_metadata_path = nothing)

Return canopy height (m) for `site_ID` from the FLUXNET2015 metadata CSV.
Same artifact requirements as `get_site_info`. Returns `NaN` (with a warning)
if the field is missing.
"""
function FluxnetSimulations.get_canopy_height(
    site_ID;
    fluxnet2015_metadata_path = nothing,
)
    if fluxnet2015_metadata_path === nothing
        fluxnet2015_metadata_path = joinpath(
            ClimaLand.Artifacts.fluxnet2015_data_path(),
            "metadata_DD_clean.csv",
        )
    end
    raw = readdlm(fluxnet2015_metadata_path, ',', Any)
    header = raw[1, :]
    data = raw[2:end, :]

    row_idx = findfirst(row -> row[1] == site_ID, eachrow(data))
    isnothing(row_idx) && error("Site ID $site_ID not found in metadata.")
    site_metadata = data[row_idx, :]

    col_idx = findfirst(header[:] .== "canopy_height")
    canopy_height = site_metadata[col_idx]

    if (canopy_height isa AbstractString && isempty(canopy_height)) ||
       (canopy_height isa AbstractFloat && isnan(canopy_height))
        @warn "Field canopy_height is missing for site $site_ID"
        return NaN
    end
    return Float64(canopy_height)
end

"""
    FluxnetSimulations.get_site_igbp(site_ID; fluxnet2015_metadata_path = nothing)

Return the IGBP land cover class string (e.g. "EBF", "DBF", "GRA", "SAV") for
`site_ID` from the FLUXNET2015 metadata CSV. Returns `"NA"` (with a warning)
if neither the column nor the cell is present, so callers can keep working
on partial metadata.

Tries a list of candidate column names, in order, since the column has
historically been named differently across snapshots of the metadata CSV.
"""
function FluxnetSimulations.get_site_igbp(
    site_ID;
    fluxnet2015_metadata_path = nothing,
)
    if fluxnet2015_metadata_path === nothing
        fluxnet2015_metadata_path = joinpath(
            ClimaLand.Artifacts.fluxnet2015_data_path(),
            "metadata_DD_clean.csv",
        )
    end
    raw = readdlm(fluxnet2015_metadata_path, ',', Any)
    header = raw[1, :]
    data = raw[2:end, :]

    row_idx = findfirst(row -> row[1] == site_ID, eachrow(data))
    isnothing(row_idx) && error("Site ID $site_ID not found in metadata.")
    site_metadata = data[row_idx, :]

    candidates = (
        "igbp_class",
        "igbp_landcover_class",
        "IGBP",
        "biome_code",
        "IGBP_class",
    )
    col_idx = nothing
    for cand in candidates
        col_idx = findfirst(header[:] .== cand)
        col_idx !== nothing && break
    end
    if col_idx === nothing
        @warn "No IGBP column in metadata CSV (tried $(candidates)) — returning NA for $site_ID"
        return "NA"
    end

    igbp = site_metadata[col_idx]
    if (igbp isa AbstractString && isempty(igbp)) ||
       (igbp isa AbstractFloat && isnan(igbp))
        @warn "IGBP class missing for $site_ID"
        return "NA"
    end
    return String(strip(string(igbp)))
end
