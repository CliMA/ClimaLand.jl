# Convert ClimaLand diagnostics into TRENDY submission files.
#
# The monthly NetCDF files written by ClimaDiagnostics (`<short_name>_1M_average.nc`)
# are renamed, converted to the TRENDY units and regridded to the regular TRENDY grid,
# and written as `<MODEL>/<SIM>/<MODEL>_<SIM>_<var>.nc`, following the TRENDY v15 /
# Global Carbon Budget 2026 protocol (variable list `trendy_listofvariables_GCP2026.xlsx`).
# The zonal NBP text file `<MODEL>/<MODEL>_zonalNBP_<SIM>.txt` and a model README
# template are written too.
#
# Usage, from the repository root:
# ```
# julia --project=.buildkite experiments/trendy/trendy_output.jl <diagnostics_dir> <output_root> <MODEL> <SIM> [resolution]
# ```
# or include this file and call `write_trendy_output(diagnostics_dir, output_root; ...)`.
#
# Writing the diagnostics directly on the TRENDY grid avoids a second regridding:
# pass `NetCDFWriter(space, outdir; start_date, horizontal_pts = trendy_grid(1.0))`
# as the output writer of `ClimaLand.default_diagnostics`. The native 1° grid of
# `global_box_domain` already has these points.

import NCDatasets
import Dates

const TRENDY_FILL_VALUE = -99999.0f0
"""Kilograms of carbon per mole of CO₂ (ClimaLand carbon fluxes are in mol CO₂ m⁻² s⁻¹)"""
const KG_C_PER_MOL_CO2 = 0.0120107
const ρ_LIQUID_WATER = 1000.0 # kg m⁻³
const ρ_ICE = 916.7 # kg m⁻³
const EARTH_RADIUS = 6.371e6 # m
const DAYS_PER_MONTH_NOLEAP = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
const TRENDY_TIME_UNITS = "days since 1700-01-01 00:00:00"

"""
    ClimaVariable

A ClimaLand diagnostic read from a monthly file: `data` has dimensions
(longitude, latitude, time) or (longitude, latitude, depth, time), with NaN
where there is no data, on the grid `lon`, `lat`, with depths `z` (m, negative
below the surface) or `nothing`, and `months` the (year, month) of each record.
"""
struct ClimaVariable
    data::Array{Float64}
    lon::Vector{Float64}
    lat::Vector{Float64}
    z::Union{Nothing, Vector{Float64}}
    months::Vector{Tuple{Int, Int}}
end

"""
    TRENDYVariable

How to compute a TRENDY variable from ClimaLand diagnostics: TRENDY `name`,
`long_name`, `units` and `priority` (from the TRENDY variable list), the ClimaLand
diagnostic short names it needs (`inputs`), and `compute(vars)`, which takes a
`Dict` of `ClimaVariable`s on the TRENDY grid and returns the TRENDY values with
the dimensions of the inputs. `layer_dim` names the soil layer dimension of
depth-resolved variables, and `description` says how the variable is obtained.
"""
struct TRENDYVariable
    name::String
    long_name::String
    units::String
    priority::Int
    inputs::Vector{String}
    compute::Function
    layer_dim::Union{Nothing, String}
    description::String
end

TRENDYVariable(name, long_name, units, priority, inputs, compute, description) =
    TRENDYVariable(
        name,
        long_name,
        units,
        priority,
        inputs,
        compute,
        nothing,
        description,
    )

"""
    layer_thickness(z; soil_depth = nothing)

Thickness (m) of the soil layers whose centres are at depths `z` (negative
below the surface): the layer boundaries are the surface, the midpoints between
centres, and the bottom of the soil at `-soil_depth`, or, if `soil_depth` is
`nothing`, the bottom centre extended by half of the last spacing.
"""
function layer_thickness(z; soil_depth = nothing)
    faces = layer_faces(z; soil_depth)
    return faces[1:(end - 1)] .- faces[2:end]
end

"""
    layer_faces(z; soil_depth = nothing)

Depths (m, negative below the surface) of the boundaries of the soil layers whose
centres are at `z`, from the surface (0) to the bottom of the soil; see
`layer_thickness`.
"""
function layer_faces(z; soil_depth = nothing)
    zs = sort(z; rev = true)
    zs == z || error("Expected depths ordered from the surface downward")
    n = length(z)
    n == 1 && return [0.0, isnothing(soil_depth) ? 2z[1] : -soil_depth]
    faces = zeros(n + 1)
    for i in 2:n
        faces[i] = (z[i - 1] + z[i]) / 2
    end
    faces[n + 1] =
        isnothing(soil_depth) ? z[n] - (z[n - 1] - z[n]) / 2 : -soil_depth
    return faces
end

"""
    top_layer_weights(z, depth; soil_depth = nothing)

Fraction of each soil layer (centres `z`) that lies within `depth` (m) of the
surface.
"""
function top_layer_weights(z, depth; soil_depth = nothing)
    faces = layer_faces(z; soil_depth)
    top = -faces[1:(end - 1)]
    bottom = -faces[2:end]
    return clamp.((min.(bottom, depth) .- top) ./ (bottom .- top), 0, 1)
end

"""Depth (m) of the top soil layer used for `mrsos`, as in CMIP6."""
const MRSOS_DEPTH = 0.1

"""
    column_integral(v::ClimaVariable; soil_depth = nothing)

Integrate `v.data` over the soil layers (third dimension), weighting by the
layer thickness.
"""
function column_integral(v::ClimaVariable; soil_depth = nothing)
    dz = layer_thickness(v.z; soil_depth)
    return dropdims(sum(v.data .* reshape(dz, 1, 1, :, 1); dims = 3); dims = 3)
end

carbon_flux(name) = vars -> vars[name].data .* KG_C_PER_MOL_CO2
copy_of(name; factor = 1.0) = vars -> vars[name].data .* factor

"""
    trendy_variables(; soil_depth = nothing)

Return the TRENDY variables that ClimaLand can provide, in the order of the
TRENDY variable list. `soil_depth` (m) is the depth of the soil column, used for
the thickness of the bottom layer in depth integrals.
"""
function trendy_variables(; soil_depth = nothing)
    layer_mass(vars) =
        (vars["swc"].data .* ρ_LIQUID_WATER .+ vars["si"].data .* ρ_ICE) .*
        reshape(layer_thickness(vars["swc"].z; soil_depth), 1, 1, :, 1)
    return [
        TRENDYVariable(
            "tas",
            "Near-Surface Air Temperature",
            "K",
            1,
            ["tair"],
            copy_of("tair"),
            "`tair`, the air temperature of the forcing",
        ),
        TRENDYVariable(
            "pr",
            "Precipitation",
            "kg m-2 s-1",
            1,
            ["precip"],
            copy_of("precip"; factor = -1.0),
            "`-precip` (ClimaLand precipitation is negative downward)",
        ),
        TRENDYVariable(
            "rsds",
            "Surface Downwelling Shortwave Radiation",
            "W m-2",
            1,
            ["swd"],
            copy_of("swd"),
            "`swd`",
        ),
        TRENDYVariable(
            "mrso",
            "Total Soil Moisture Content",
            "kg m-2",
            1,
            ["swc", "si"],
            vars -> dropdims(sum(layer_mass(vars); dims = 3); dims = 3),
            "sum over all soil layers (the full soil column, see soil_column_depth) of (`swc` ρ_liquid + `si` ρ_ice) Δz",
        ),
        TRENDYVariable(
            "mrsos",
            "Moisture in Upper Portion of Soil Column",
            "kg m-2",
            2,
            ["swc", "si"],
            vars -> dropdims(
                sum(
                    layer_mass(vars) .* reshape(
                        top_layer_weights(
                            vars["swc"].z,
                            MRSOS_DEPTH;
                            soil_depth,
                        ),
                        1,
                        1,
                        :,
                        1,
                    );
                    dims = 3,
                );
                dims = 3,
            ),
            "(`swc` ρ_liquid + `si` ρ_ice) Δz summed over the top 0.1 m of soil (as CMIP6 mrsos)",
        ),
        TRENDYVariable(
            "mrros",
            "Surface Runoff",
            "kg m-2 s-1",
            1,
            ["sr"],
            copy_of("sr"; factor = ρ_LIQUID_WATER),
            "`sr` × 1000 kg m-3",
        ),
        TRENDYVariable(
            "mrro",
            "Total Runoff",
            "kg m-2 s-1",
            1,
            ["tr"],
            copy_of("tr"; factor = ρ_LIQUID_WATER),
            "`tr` (surface + subsurface) × 1000 kg m-3",
        ),
        TRENDYVariable(
            "evapotrans",
            "Total Evapo-Transpiration",
            "kg m-2 s-1",
            1,
            ["et"],
            copy_of("et"),
            "`et`, positive upward",
        ),
        TRENDYVariable(
            "cVeg",
            "Carbon in Vegetation",
            "kg m-2",
            1,
            ["cveg"],
            copy_of("cveg"),
            "`cveg`",
        ),
        TRENDYVariable(
            "cSoil",
            "Carbon in Soil (including below-ground litter)",
            "kg m-2",
            1,
            ["soc"],
            vars -> column_integral(vars["soc"]; soil_depth),
            "`soc` integrated over the whole soil column",
        ),
        TRENDYVariable(
            "gpp",
            "Gross Primary Production",
            "kg m-2 s-1",
            1,
            ["gpp"],
            carbon_flux("gpp"),
            "`gpp` × 0.0120107 kg C mol-1",
        ),
        TRENDYVariable(
            "ra",
            "Autotrophic (Plant) Respiration",
            "kg m-2 s-1",
            1,
            ["ra"],
            carbon_flux("ra"),
            "`ra` × 0.0120107 kg C mol-1",
        ),
        TRENDYVariable(
            "npp",
            "Net Primary Production",
            "kg m-2 s-1",
            1,
            ["gpp", "ra"],
            vars -> (vars["gpp"].data .- vars["ra"].data) .* KG_C_PER_MOL_CO2,
            "(`gpp` - `ra`) × 0.0120107 kg C mol-1",
        ),
        TRENDYVariable(
            "rh",
            "Heterotrophic Respiration",
            "kg m-2 s-1",
            1,
            ["hr"],
            carbon_flux("hr"),
            "`hr` × 0.0120107 kg C mol-1",
        ),
        TRENDYVariable(
            "fLuc",
            "CO2 Flux to Atmosphere from Land Use Change",
            "kg m-2 s-1",
            1,
            ["gpp"],
            vars -> vars["gpp"].data .* 0,
            "0 on land: ClimaLand has no land-use change",
        ),
        TRENDYVariable(
            "nbp",
            "Net Biospheric Production",
            "kg m-2 s-1",
            1,
            ["nee"],
            vars -> -vars["nee"].data .* KG_C_PER_MOL_CO2,
            "-`nee` × 0.0120107 kg C mol-1 (positive into the land; no fire, land use or harvest fluxes)",
        ),
        TRENDYVariable(
            "lai",
            "Leaf Area Index",
            "1",
            1,
            ["lai"],
            copy_of("lai"),
            "`lai`",
        ),
        TRENDYVariable(
            "tsl",
            "Temperature of Soil",
            "K",
            2,
            ["tsoil"],
            copy_of("tsoil"),
            "stlayer",
            "`tsoil` at the soil layer centres",
        ),
        TRENDYVariable(
            "msl",
            "Moisture of Soil",
            "kg m-2",
            2,
            ["swc", "si"],
            layer_mass,
            "smlayer",
            "(`swc` ρ_liquid + `si` ρ_ice) Δz of each soil layer",
        ),
        TRENDYVariable(
            "tran",
            "Transpiration",
            "kg m-2 s-1",
            2,
            ["trans"],
            copy_of("trans"),
            "`trans`",
        ),
    ]
end

"""
First-priority TRENDY variables that ClimaLand cannot provide, with the reason.
"""
const TRENDY_VARIABLES_NOT_PROVIDED = [
    "cLitter" => "no separate above-ground litter pool (litter goes to the soil organic carbon)",
    "cProduct" => "no land use, no product pools",
    "fFire" => "no fire model",
    "burntArea" => "no fire model",
    "fFireCveg, fFireLitter, fFireCsoil" => "no fire model",
    "soilr" => "no diagnostic of root respiration in the soil CO2 efflux",
    "cLeaf, cWood, cRoot" => "pools exist in PrognosticCarbonModel but have no diagnostics yet",
    "fVegLitter, fLeafLitter, fWoodLitter, fRootLitter" => "fluxes exist in PrognosticCarbonModel but have no diagnostics yet",
    "fLitterSoil, fVegSoil, cCwd, cSoilPools, rhPools" => "no litter or coarse woody debris pools; one soil carbon pool",
    "fAllocLeaf, fAllocWood, fAllocRoot" => "allocation exists in PrognosticCarbonModel but has no diagnostics yet",
    "landCoverFrac, PFT-level variables (*pft)" => "no PFT tiles; tree, C3 and C4 grass fractions (`ftr`, `fc3g`, `fc4g`) only with ZhouOptimalLAIModel",
]

"""
    trendy_grid(resolution = 1.0)

Longitudes and latitudes of the centres of the regular TRENDY grid of
`resolution` degrees, from -180 to 180 and from -90 to 90.
"""
function trendy_grid(resolution = 1.0)
    isinteger(180 / resolution) ||
        error("180 must be a multiple of the resolution $resolution")
    lon = collect(
        range(-180 + resolution / 2, 180 - resolution / 2; step = resolution),
    )
    lat = collect(
        range(-90 + resolution / 2, 90 - resolution / 2; step = resolution),
    )
    return lon, lat
end

"""
    find_short_name(filepath)

Short name of the diagnostic in a ClimaDiagnostics file `<short_name>_1M_average.nc`.
"""
find_short_name(filepath) = replace(basename(filepath), "_1M_average.nc" => "")

"""
    record_months(ds)

Year and month of each monthly record of the ClimaDiagnostics dataset `ds`,
from the start of its averaging window: `time_bnds` if present, otherwise the
time stamp, which ClimaDiagnostics sets to the start of the window.
"""
function record_months(ds)
    starts = haskey(ds, "time_bnds") ? ds["time_bnds"][1, :] : ds["time"][:]
    starts = Dates.DateTime.(starts)
    all(d -> Dates.day(d) == 1 && Dates.Time(d) == Dates.Time(0), starts) ||
        error(
            "Expected monthly records starting at 00:00 on the first of the month",
        )
    return [(Dates.year(d), Dates.month(d)) for d in starts]
end

"""
    read_clima_variable(filepath)

Read the monthly ClimaDiagnostics file `filepath` as a `ClimaVariable`.
"""
function read_clima_variable(filepath)
    short_name = find_short_name(filepath)
    NCDatasets.NCDataset(filepath) do ds
        v = ds[short_name]
        names = NCDatasets.dimnames(v)
        order = [
            findfirst(==(d), names) for
            d in ("lon", "lat", "z", "time") if d in names
        ]
        length(order) == length(names) ||
            error("Unexpected dimensions $names in $filepath")
        data = permutedims(Float64.(NCDatasets.nomissing(Array(v), NaN)), order)
        z = "z" in names ? Float64.(ds["z"][:]) : nothing
        if !isnothing(z) && z[1] < z[end]
            # Order depths from the surface downward
            z = reverse(z)
            data = reverse(data; dims = 3)
        end
        lon = Float64.(ds["lon"][:])
        lat = Float64.(ds["lat"][:])
        return ClimaVariable(data, lon, lat, z, record_months(ds))
    end
end

"""
    regrid(v::ClimaVariable, lon, lat)

Return `v` on the grid of centres `lon`, `lat`, by bilinear interpolation,
periodic in longitude, ignoring NaN neighbours (NaN only where all four
neighbours are NaN). Returns `v` unchanged if it is already on that grid.
"""
function regrid(v::ClimaVariable, lon, lat)
    src_lon = mod.(v.lon .+ 180, 360) .- 180
    perm = sortperm(src_lon)
    src_lon = src_lon[perm]
    data = selectdim(v.data, 1, perm)
    if length(src_lon) == length(lon) &&
       length(v.lat) == length(lat) &&
       all(isapprox.(src_lon, lon; atol = 1e-6)) &&
       all(isapprox.(v.lat, lat; atol = 1e-6))
        return ClimaVariable(collect(data), lon, lat, v.z, v.months)
    end
    issorted(v.lat) || error("Expected increasing latitudes")
    nlon = length(src_lon)
    out = fill(NaN, length(lon), length(lat), size(v.data)[3:end]...)
    rest = CartesianIndices(size(v.data)[3:end])
    # Periodic longitudes: extend by one point on each side
    ext_lon = [src_lon[end] - 360; src_lon; src_lon[1] + 360]
    ext_index = [nlon; 1:nlon; 1]
    for (j, φ) in enumerate(lat), (i, λ) in enumerate(lon)
        i1 = clamp(searchsortedlast(ext_lon, λ), 1, nlon + 1)
        wλ = (λ - ext_lon[i1]) / (ext_lon[i1 + 1] - ext_lon[i1])
        if φ <= v.lat[1]
            j1, j2, wφ = 1, 1, 0.0
        elseif φ >= v.lat[end]
            j1, j2, wφ = length(v.lat), length(v.lat), 0.0
        else
            j1 = searchsortedlast(v.lat, φ)
            j2 = j1 + 1
            wφ = (φ - v.lat[j1]) / (v.lat[j2] - v.lat[j1])
        end
        stencil = (
            (ext_index[i1], j1, (1 - wλ) * (1 - wφ)),
            (ext_index[i1 + 1], j1, wλ * (1 - wφ)),
            (ext_index[i1], j2, (1 - wλ) * wφ),
            (ext_index[i1 + 1], j2, wλ * wφ),
        )
        for r in rest
            total = 0.0
            weight = 0.0
            for (si, sj, w) in stencil
                x = data[si, sj, r]
                if !isnan(x) && w > 0
                    total += w * x
                    weight += w
                end
            end
            weight > 0 && (out[i, j, r] = total / weight)
        end
    end
    return ClimaVariable(out, lon, lat, v.z, v.months)
end

"""
    cell_edges(centres; periodic = false)

Edges of the cells of the regular grid `centres` (midpoints, with half a spacing
added at both ends; latitudes are clipped to ±90).
"""
function cell_edges(centres; periodic = false)
    Δ = centres[2] - centres[1]
    edges = [centres .- Δ / 2; centres[end] + Δ / 2]
    return periodic ? edges : clamp.(edges, -90, 90)
end

"""
    overlap_matrix(src_edges, dst_edges; period = nothing)

Matrix (destination, source) of the lengths of overlap between source and
destination intervals, periodic with `period` if given.
"""
function overlap_matrix(src_edges, dst_edges; period = nothing)
    shifts = isnothing(period) ? (0.0,) : (-period, 0.0, period)
    O = zeros(length(dst_edges) - 1, length(src_edges) - 1)
    for a in 1:(length(src_edges) - 1),
        s in shifts,
        i in 1:(length(dst_edges) - 1)

        lo = max(src_edges[a] + s, dst_edges[i])
        hi = min(src_edges[a + 1] + s, dst_edges[i + 1])
        hi > lo && (O[i, a] += hi - lo)
    end
    return O
end

"""
    land_fraction_on_grid(src_land, src_lon, src_lat, lon, lat)

Area fraction of each cell of the regular grid `lon`, `lat` covered by the source
cells where `src_land` (longitude, latitude; 1 land, 0 not) is 1. Source and
destination grids are regular in longitude and latitude.
"""
function land_fraction_on_grid(src_land, src_lon, src_lat, lon, lat)
    src_lon = mod.(src_lon .+ 180, 360) .- 180
    perm = sortperm(src_lon)
    src_lon = src_lon[perm]
    src_land = src_land[perm, :]
    issorted(src_lat) || error("Expected increasing latitudes")
    Olon = overlap_matrix(
        cell_edges(src_lon; periodic = true),
        cell_edges(lon; periodic = true);
        period = 360.0,
    )
    # Overlap in sin(latitude), proportional to area
    Olat = overlap_matrix(sind.(cell_edges(src_lat)), sind.(cell_edges(lat)))
    covered = Olon * src_land * transpose(Olat)
    dst_area =
        diff(cell_edges(lon; periodic = true)) .*
        transpose(diff(sind.(cell_edges(lat))))
    return clamp.(covered ./ dst_area, 0, 1)
end

"""
    month_midpoint_days(year, month)

Days since 1700-01-01 of the middle of `month` of `year` in the 365-day
calendar, and the bounds of the month.
"""
function month_midpoint_days(year, month)
    start =
        365 * (year - 1700) +
        sum(DAYS_PER_MONTH_NOLEAP[1:(month - 1)]; init = 0)
    stop = start + DAYS_PER_MONTH_NOLEAP[month]
    return (start + stop) / 2, start, stop
end

"""
    write_coordinates!(ds, lon, lat)

Define the `longitude` and `latitude` dimensions and variables in `ds`.
"""
function write_coordinates!(ds, lon, lat)
    NCDatasets.defVar(
        ds,
        "longitude",
        lon,
        ("longitude",);
        attrib = [
            "standard_name" => "longitude",
            "long_name" => "longitude",
            "units" => "degrees_east",
            "axis" => "X",
        ],
    )
    NCDatasets.defVar(
        ds,
        "latitude",
        lat,
        ("latitude",);
        attrib = [
            "standard_name" => "latitude",
            "long_name" => "latitude",
            "units" => "degrees_north",
            "axis" => "Y",
        ],
    )
end

"""
    write_time!(ds, months)

Define the `time` dimension (365-day calendar, mid-month) and its bounds.
"""
function write_time!(ds, months)
    days = [month_midpoint_days(y, m) for (y, m) in months]
    NCDatasets.defDim(ds, "bnds", 2)
    NCDatasets.defVar(
        ds,
        "time",
        first.(days),
        ("time",);
        attrib = [
            "standard_name" => "time",
            "units" => TRENDY_TIME_UNITS,
            "calendar" => "365_day",
            "axis" => "T",
            "bounds" => "time_bnds",
        ],
    )
    bounds = Float64.(hcat([[d[2], d[3]] for d in days]...))
    NCDatasets.defVar(ds, "time_bnds", bounds, ("bnds", "time"))
end

"""
    to_fill(data)

`data` as `Float32` with NaN replaced by the TRENDY fill value.
"""
to_fill(data) = Float32.(map(x -> isnan(x) ? TRENDY_FILL_VALUE : x, data))

variable_attributes(var::TRENDYVariable) = [
    "long_name" => var.long_name,
    "units" => var.units,
    "_FillValue" => TRENDY_FILL_VALUE,
    "missing_value" => TRENDY_FILL_VALUE,
    "comment" => "ClimaLand: " * var.description,
]

"""
    write_trendy_variable(path, var, data, lon, lat, months, attributes; z = nothing)

Write the TRENDY variable `var` with values `data` (dimensions longitude,
latitude, [layer,] time) to the file `path`.
"""
function write_trendy_variable(
    path,
    var::TRENDYVariable,
    data,
    lon,
    lat,
    months,
    attributes;
    z = nothing,
    soil_depth = nothing,
)
    var_attributes = variable_attributes(var)
    if !isnothing(z) && "swc" in var.inputs
        faces = layer_faces(z; soil_depth)
        push!(var_attributes, "soil_column_depth" => Float32(-faces[end]))
        var.name == "mrsos" &&
            push!(var_attributes, "soil_depth" => Float32(MRSOS_DEPTH))
    end
    NCDatasets.NCDataset(path, "c"; attrib = attributes) do ds
        write_coordinates!(ds, lon, lat)
        write_time!(ds, months)
        dims = ("longitude", "latitude", "time")
        if !isnothing(var.layer_dim)
            NCDatasets.defVar(
                ds,
                var.layer_dim,
                -z,
                (var.layer_dim,);
                attrib = [
                    "long_name" => "depth of the soil layer centre",
                    "units" => "m",
                    "positive" => "down",
                    "axis" => "Z",
                    "bounds" => var.layer_dim * "_bnds",
                ],
            )
            faces = layer_faces(z; soil_depth)
            NCDatasets.defVar(
                ds,
                var.layer_dim * "_bnds",
                permutedims(hcat(-faces[1:(end - 1)], -faces[2:end])),
                ("bnds", var.layer_dim);
                attrib = [
                    "long_name" => "depths of the top and bottom of the soil layer",
                    "units" => "m",
                ],
            )
            dims = ("longitude", "latitude", var.layer_dim, "time")
        end
        NCDatasets.defVar(
            ds,
            var.name,
            to_fill(data),
            dims;
            attrib = var_attributes,
        )
    end
end

"""
    grid_cell_area(lon, lat)

Area (m²) of the cells of the regular grid of centres `lon`, `lat`, as a matrix
(longitude, latitude).
"""
function grid_cell_area(lon, lat)
    Δλ = deg2rad(lon[2] - lon[1])
    Δφ = lat[2] - lat[1]
    band =
        EARTH_RADIUS^2 * Δλ .*
        (sind.(min.(lat .+ Δφ / 2, 90)) .- sind.(max.(lat .- Δφ / 2, -90)))
    return repeat(reshape(band, 1, :), length(lon), 1)
end

"""
    zonal_nbp(nbp, land_fraction, lon, lat, months; allow_partial_years = false)

Annual NBP (PgC yr⁻¹) of the globe, the north (> 30°N), the tropics (30°S–30°N)
and the south (< 30°S), from the monthly `nbp` (kg m⁻² s⁻¹ of land, dimensions
longitude, latitude, time), weighted by the cell area and `land_fraction`. Only
complete years are returned unless `allow_partial_years`.
"""
function zonal_nbp(
    nbp,
    land_fraction,
    lon,
    lat,
    months;
    allow_partial_years = false,
)
    area = grid_cell_area(lon, lat) .* land_fraction
    bands = (trues(length(lat)), lat .> 30, -30 .<= lat .<= 30, lat .< -30)
    rows = Tuple{Int, NTuple{4, Float64}}[]
    for year in unique(first.(months))
        indices = findall(m -> m[1] == year, months)
        length(indices) == 12 || allow_partial_years || continue
        totals = map(bands) do band
            sum(indices) do t
                seconds = DAYS_PER_MONTH_NOLEAP[months[t][2]] * 86400
                values = nbp[:, band, t] .* area[:, band] .* seconds
                sum(x -> isnan(x) ? 0.0 : x, values)
            end * 1e-12
        end
        push!(rows, (year, totals))
    end
    return rows
end

"""
    write_zonal_nbp(path, rows)

Write the rows of `zonal_nbp` as the TRENDY zonal NBP text file.
"""
function write_zonal_nbp(path, rows)
    open(path, "w") do io
        println(io, "Year\tGlobal\tNorth\tTropics\tSouth")
        for (year, totals) in rows
            println(
                io,
                year,
                "\t",
                join((string(round(x; digits = 5)) for x in totals), "\t"),
            )
        end
    end
end

"""
    write_readme_template(path, model; simulations, institution, contact, source_code, commit)

Write a README template for the model submission, to be completed by hand.
"""
function write_readme_template(
    path,
    model;
    simulations,
    institution,
    contact,
    source_code,
    commit,
    resolution,
)
    open(path, "w") do io
        print(
            io,
            """
            $model TRENDY v15 (Global Carbon Budget 2026) submission
            ========================================================

            Model: $model (ClimaLand, CliMA land model)
            Institution: $institution
            Contact: $contact
            Source code: $source_code
            Commit: $commit
            Simulations: $(join(simulations, ", "))
            Grid: regular $(resolution)° longitude-latitude grid, -180 to 180, -90 to 90
            Time: monthly means, 365-day calendar ($TRENDY_TIME_UNITS)
            Units: pools kg m-2 (carbon), fluxes kg m-2 s-1; fill value -99999

            Configuration
            -------------
            [TO COMPLETE: prescribed (MODIS) or prognostic LAI, carbon pools, soil
            carbon model, spin-up length and state, forcing (CRU-JRA v4.10 with the
            protocol diffuse fraction), CO2 (global_co2_ann_1700_2025.txt, annual).]

            Fluxes
            ------
            nbp = gpp - ra - rh - fLuc, positive into the land. ClimaLand has no fire,
            harvest, grazing or land-use change: fFire is not provided and fLuc = 0.
            Global annual NBP: sum over months of nbp × seconds of the month × cell
            area × (1 - oceanCoverFrac); see <MODEL>_zonalNBP_<SIM>.txt (PgC yr-1,
            north > 30°N, tropics 30°S-30°N, south < 30°S).

            Variables not provided
            ----------------------
            """,
        )
        for (names, reason) in TRENDY_VARIABLES_NOT_PROVIDED
            println(io, names, ": ", reason)
        end
    end
end

"""
    write_trendy_output(diagnostics_dir, output_root;
                        model = "ClimaLand",
                        simulation = "S2",
                        resolution = 1.0,
                        variables = nothing,
                        institution = "Climate Modeling Alliance (CliMA), California Institute of Technology",
                        contact = "TO COMPLETE",
                        source_code = "https://github.com/CliMA/ClimaLand.jl",
                        commit = "unknown",
                        soil_depth = nothing,
                        allow_partial_years = false,
                       )

Convert the monthly ClimaDiagnostics files in `diagnostics_dir` into TRENDY files
`<output_root>/<model>/<simulation>/<model>_<simulation>_<var>.nc` on the
regular grid of `resolution` degrees (1 or 0.5), and write
`<output_root>/<model>/<model>_zonalNBP_<simulation>.txt` and, if missing, a
README template `<output_root>/<model>/README_<model>.txt`.

`model` is the TRENDY model name, without underscores: `ClimaLand` for the runs with
prognostic vegetation and `ClimaLand-OBS` for those with prescribed vegetation.
`variables` restricts the TRENDY variables written (default: all that can be
computed from the files present). `oceanCoverFrac` is 1 where the model has no
land and 0 elsewhere. Returns the paths of the files written.
"""
function write_trendy_output(
    diagnostics_dir,
    output_root;
    model = "ClimaLand",
    simulation = "S2",
    resolution = 1.0,
    variables = nothing,
    institution = "Climate Modeling Alliance (CliMA), California Institute of Technology",
    contact = "TO COMPLETE",
    source_code = "https://github.com/CliMA/ClimaLand.jl",
    commit = "unknown",
    soil_depth = nothing,
    allow_partial_years = false,
)
    occursin('_', model) &&
        error("TRENDY model names cannot contain underscores: $model")
    simulation in ("S0", "S1", "S2", "S3", "S4", "S5", "S6") ||
        error("Unknown TRENDY simulation $simulation")
    lon, lat = trendy_grid(resolution)
    files = filter(
        f -> endswith(f, "_1M_average.nc"),
        readdir(diagnostics_dir; join = true),
    )
    available = Dict(find_short_name(f) => f for f in files)
    selected = filter(trendy_variables(; soil_depth)) do var
        (isnothing(variables) || var.name in variables) &&
            all(in(keys(available)), var.inputs)
    end
    isempty(selected) &&
        error("No TRENDY variable can be computed from $diagnostics_dir")
    cache = Dict{String, ClimaVariable}()
    function input(name)
        haskey(cache, name) || (
            cache[name] =
                regrid(read_clima_variable(available[name]), lon, lat)
        )
        return cache[name]
    end
    attributes = [
        "Conventions" => "CF-1.8",
        "title" => "$model $simulation output for TRENDY v15 (Global Carbon Budget 2026)",
        "institution" => institution,
        "contact" => contact,
        "source" => "ClimaLand.jl, the CliMA land model",
        "source_code" => source_code,
        "commit" => commit,
        "model" => model,
        "simulation" => simulation,
        "history" => "$(Dates.now()): converted from ClimaLand diagnostics with experiments/trendy/trendy_output.jl",
    ]
    sim_dir = joinpath(output_root, model, simulation)
    mkpath(sim_dir)
    written = String[]
    # Some diagnostics (e.g. `lai`, from data filled over the ocean) are not masked,
    # so all variables take the land of the first masked diagnostic available
    mask_name = something(
        findfirst(in(keys(available)), ["tair", "gpp", "swc", "et"]),
        0,
    )
    mask_name == 0 &&
        error("No diagnostic to find the land in $diagnostics_dir")
    # Land fraction of each TRENDY cell = area fraction covered by model land cells
    # (exactly 0 or 1 when the diagnostics are on the TRENDY grid). Values are kept
    # wherever the fraction is > 0; global totals need area × (1 - oceanCoverFrac).
    mask_source =
        read_clima_variable(available[["tair", "gpp", "swc", "et"][mask_name]])
    source_land =
        .!isnan.(selectdim(mask_source.data, ndims(mask_source.data), 1))
    ndims(source_land) == 3 && (source_land = source_land[:, :, 1])
    land_fraction = land_fraction_on_grid(
        Float64.(source_land),
        mask_source.lon,
        mask_source.lat,
        lon,
        lat,
    )
    land = land_fraction .> 1e-6
    months = nothing
    nbp = nothing
    for var in selected
        vars = Dict(name => input(name) for name in var.inputs)
        first_input = vars[var.inputs[1]]
        months = first_input.months
        all(v -> v.months == months, values(vars)) ||
            error("Inputs of $(var.name) cover different months")
        data = var.compute(vars)
        data[.!land, ntuple(_ -> :, ndims(data) - 2)...] .= NaN
        path = joinpath(sim_dir, "$(model)_$(simulation)_$(var.name).nc")
        isfile(path) && rm(path)
        write_trendy_variable(
            path,
            var,
            data,
            lon,
            lat,
            months,
            attributes;
            z = first_input.z,
            soil_depth,
        )
        push!(written, path)
        var.name == "nbp" && (nbp = data)
    end

    ocean_path = joinpath(sim_dir, "$(model)_$(simulation)_oceanCoverFrac.nc")
    isfile(ocean_path) && rm(ocean_path)
    NCDatasets.NCDataset(ocean_path, "c"; attrib = attributes) do ds
        write_coordinates!(ds, lon, lat)
        NCDatasets.defVar(
            ds,
            "oceanCoverFrac",
            Float32.(1 .- land_fraction),
            ("longitude", "latitude");
            attrib = [
                "long_name" => "Fractional Ocean Cover",
                "units" => "1",
                "_FillValue" => TRENDY_FILL_VALUE,
                "missing_value" => TRENDY_FILL_VALUE,
                "comment" => "ClimaLand: 1 where the model has no land (land-sea mask of the simulation), 0 elsewhere",
            ],
        )
    end
    push!(written, ocean_path)

    if !isnothing(nbp)
        zonal_path =
            joinpath(output_root, model, "$(model)_zonalNBP_$(simulation).txt")
        rows =
            zonal_nbp(nbp, land_fraction, lon, lat, months; allow_partial_years)
        write_zonal_nbp(zonal_path, rows)
        push!(written, zonal_path)
    end
    readme_path = joinpath(output_root, model, "README_$(model).txt")
    if !isfile(readme_path)
        write_readme_template(
            readme_path,
            model;
            simulations = [simulation],
            institution,
            contact,
            source_code,
            commit,
            resolution,
        )
        push!(written, readme_path)
    end
    return written
end

"""
    trendy_mapping_table()

Markdown table of the TRENDY variables and how they are computed from the
ClimaLand diagnostics, followed by the first-priority variables not provided.
"""
function trendy_mapping_table()
    lines = [
        "| TRENDY | priority | units | ClimaLand diagnostics | conversion |",
        "|---|---|---|---|---|",
    ]
    for var in trendy_variables()
        push!(
            lines,
            "| `$(var.name)` | $(var.priority) | $(var.units) | $(join(var.inputs, ", ")) | $(var.description) |",
        )
    end
    push!(
        lines,
        "| `oceanCoverFrac` | 1 | 1 | land-sea mask | 1 where no model land |",
    )
    push!(lines, "", "Not provided (first priority):", "")
    for (names, reason) in TRENDY_VARIABLES_NOT_PROVIDED
        push!(lines, "- $names: $reason")
    end
    return join(lines, "\n")
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) >= 4 || error(
        "Usage: trendy_output.jl <diagnostics_dir> <output_root> <MODEL> <SIM> [resolution]",
    )
    resolution = length(ARGS) >= 5 ? parse(Float64, ARGS[5]) : 1.0
    written = write_trendy_output(
        ARGS[1],
        ARGS[2];
        model = ARGS[3],
        simulation = ARGS[4],
        resolution,
    )
    foreach(println, written)
end
