###################################
#       FLUXNET2015 metadata      #
###################################
#
# These helpers read site metadata (lat/lon, time offset, sensor heights,
# canopy height) from the metadata table of the `fluxnet2015` artifact. That
# artifact is too large to be downloadable; see the `fluxnet2015` entry in
# ClimaArtifacts for how to obtain it.

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

is_empty_field(x) =
    (x isa AbstractString && isempty(x)) || (x isa AbstractFloat && isnan(x))

"""
    read_site_metadata(site_ID, fluxnet2015_metadata_path)

Return `(header, row)` of the FLUXNET2015 metadata CSV for `site_ID`. If
`fluxnet2015_metadata_path` is `nothing`, the table of the `fluxnet2015`
artifact is used.
"""
function read_site_metadata(site_ID, fluxnet2015_metadata_path)
    path =
        isnothing(fluxnet2015_metadata_path) ?
        joinpath(
            ClimaLand.Artifacts.fluxnet2015_data_path(),
            "metadata_DD_clean.csv",
        ) : fluxnet2015_metadata_path
    raw = readdlm(path, ',', Any)
    row_idx = findfirst(==(site_ID), raw[2:end, 1])
    isnothing(row_idx) &&
        error("Site ID $site_ID not found in the FLUXNET2015 metadata.")
    return raw[1, :], raw[row_idx + 1, :]
end

"""
    FluxnetSimulations.get_site_info(site_ID; fluxnet2015_metadata_path = nothing)

Look up site metadata for `site_ID` from the FLUXNET2015 `metadata_DD_clean.csv`.

Returns a NamedTuple:
- `lat::Float64` — latitude (deg)
- `long::Float64` — longitude (deg)
- `time_offset` — local standard time minus UTC, in hours (e.g. -6 for US-MOz),
  the same convention as the hardcoded `get_location` methods. An `Int`, or a
  `Float64` for fractional offsets (e.g. 9.5 for sites in central Australia).
- `atmospheric_sensor_height::Vector{Float64}` — heights (m) of atmospheric
  sensors at the site, sorted; may be empty if the metadata cell is missing.

Missing fields produce a `@warn` and a `NaN`-valued entry rather than an error,
so a partial metadata row still mostly works.
"""
function FluxnetSimulations.get_site_info(
    site_ID;
    fluxnet2015_metadata_path = nothing,
)
    header, site_metadata =
        read_site_metadata(site_ID, fluxnet2015_metadata_path)
    field(varname) = site_metadata[findfirst(==(varname), header)]

    for varname in
        ("latitude", "longitude", "utc_offset", "atmospheric_sensor_heights")
        is_empty_field(field(varname)) &&
            @warn "Field $(varname) is missing for site $site_ID"
    end
    nan_if_empty(x) = is_empty_field(x) ? NaN : x

    utc_offset = nan_if_empty(field("utc_offset"))
    return (;
        lat = nan_if_empty(field("latitude")),
        long = nan_if_empty(field("longitude")),
        time_offset = isinteger(utc_offset) ? Int(utc_offset) : utc_offset,
        atmospheric_sensor_height = parse_array_field(
            field("atmospheric_sensor_heights"),
        ),
    )
end

"""
    site_ID_from_val(site_ID_val::Symbol)

Inverse of [`replace_hyphen`](@ref): `:US_MOz` becomes `"US-MOz"`.
"""
site_ID_from_val(site_ID_val::Symbol) = replace(String(site_ID_val), "_" => "-")

"""
    get_location(FT, ::Val{site_ID_val}; fluxnet2015_metadata_path = nothing)

Returns `(; time_offset, lat, long)` for a FLUXNET2015 site without a
hardcoded configuration, read from the FLUXNET2015 metadata with
[`get_site_info`](@ref). This requires the `fluxnet2015` artifact. Errors if
any of the three is missing from the metadata.
"""
function FluxnetSimulations.get_location(
    FT,
    ::Val{site_ID_val};
    fluxnet2015_metadata_path = nothing,
) where {site_ID_val}
    site_ID = site_ID_from_val(site_ID_val)
    (; time_offset, lat, long) =
        FluxnetSimulations.get_site_info(site_ID; fluxnet2015_metadata_path)
    any(isnan, (time_offset, lat, long)) && error(
        "The FLUXNET2015 metadata is missing the latitude, longitude or UTC \
         offset of $site_ID; set `lat`, `long` and `time_offset` directly.",
    )
    return (; time_offset, lat = FT(lat), long = FT(long))
end

"""
    get_fluxtower_height(FT, ::Val{site_ID_val}; fluxnet2015_metadata_path = nothing)

Returns `(; atmos_h)` for a FLUXNET2015 site without a hardcoded configuration:
the highest atmospheric sensor height in the FLUXNET2015 metadata. This
requires the `fluxnet2015` artifact.
"""
function FluxnetSimulations.get_fluxtower_height(
    FT,
    ::Val{site_ID_val};
    fluxnet2015_metadata_path = nothing,
) where {site_ID_val}
    site_ID = site_ID_from_val(site_ID_val)
    (; atmospheric_sensor_height) =
        FluxnetSimulations.get_site_info(site_ID; fluxnet2015_metadata_path)
    isempty(atmospheric_sensor_height) && error(
        "The FLUXNET2015 metadata has no atmospheric sensor height for \
         $site_ID; pass `atmos_h` to `prescribed_forcing_fluxnet` directly.",
    )
    return (; atmos_h = FT(maximum(atmospheric_sensor_height)))
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
    header, site_metadata =
        read_site_metadata(site_ID, fluxnet2015_metadata_path)
    canopy_height = site_metadata[findfirst(==("canopy_height"), header)]
    if is_empty_field(canopy_height)
        @warn "Field canopy_height is missing for site $site_ID"
        return NaN
    end
    return Float64(canopy_height)
end
