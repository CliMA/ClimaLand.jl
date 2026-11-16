# Observation-based targets of phase 1a, on CPU: global and banded land totals (north of
# 60°S, each dataset on its own grid, cos(lat) areas) of
#   - GPP, ER, NEE and Rh from the `inversion_nee` artifact, 2002–2020 means (Pg C/yr);
#   - SOC over 0–1 m and 0–2 m from `soilgrids_ocd_lowres_interp` (Pg C);
#   - live woody carbon from `xusaatchi_woody_carbon` (Pg C).
# Writes <out_dir>/targets.csv (quantity, band, value, unit, source).
#
# Usage: julia --project=.buildkite targets.jl <out_dir>
using NCDatasets
using Dates
using ClimaLand

const R = 6.371e6
const BANDS = (("global", -60, 90), ("north_50N", 50, 90), ("30N_50N", 30, 50), ("tropics_30S_30N", -30, 30), ("south_60S_30S", -60, -30))
out_dir = ARGS[1]
rows = String["quantity,band,value,unit,source"]

function coords(ds)
    lon = Float64.(ds[haskey(ds, "lon") ? "lon" : "longitude"][:])
    lat = Float64.(ds[haskey(ds, "lat") ? "lat" : "latitude"][:])
    return lon, lat
end
cell_area(lat, dlat, dlon) = R^2 * cosd(lat) * deg2rad(dlat) * deg2rad(dlon)

"""Banded integrals of a lon × lat map `x` (per m²), NaN or missing left out."""
function banded(x, lon, lat)
    dlon = abs(lon[2] - lon[1])
    dlat = abs(lat[2] - lat[1])
    out = Dict{String, Float64}()
    for (name, lo, hi) in BANDS
        total = 0.0
        for j in eachindex(lat), i in eachindex(lon)
            (lo <= lat[j] < hi) || continue
            v = x[i, j]
            (ismissing(v) || !isfinite(v)) && continue
            total += v * cell_area(lat[j], dlat, dlon)
        end
        out[name] = total
    end
    return out
end

# inversion fluxes, g C m⁻² month⁻¹ (rh: g C m⁻² day⁻¹)
path = ClimaLand.Artifacts.inversion_nee_dataset_path()
NCDataset(path) do ds
    lon, lat = coords(ds)
    times = ds["time"][:]
    for v in ("gpp", "er", "nee", "rh")
        x = ds[v]
        dims = dimnames(x)
        annual = zeros(length(lon), length(lat))
        counts = zeros(length(lon), length(lat))
        for (m, t) in enumerate(times)
            a = dims[1] in ("lon", "longitude") ? x[:, :, m] : permutedims(x[m, :, :], (1, 2))
            scale = v == "rh" ? Dates.daysinmonth(t) : 1.0
            for I in eachindex(a)
                ismissing(a[I]) && continue
                isfinite(a[I]) || continue
                annual[I] += a[I] * scale
                counts[I] += 1
            end
        end
        nyears = length(times) / 12
        mean_rate = [counts[I] > 0 ? annual[I] / nyears * 1e-3 : NaN for I in eachindex(annual)]  # kg C m⁻² yr⁻¹
        for (band, total) in banded(reshape(mean_rate, size(annual)), lon, lat)
            push!(rows, "$(v)_inversion,$band,$(total * 1e-12),PgC/yr,inversion_nee 2002-2020 mean")
        end
    end
end

# SoilGrids OCD, kg m⁻³ at depth levels
path = ClimaLand.Artifacts.soil_grids_ocd_artifact_path()
NCDataset(path) do ds
    println("SoilGrids variables: ", keys(ds))
    lon, lat = coords(ds)
    zname = first(k for k in keys(ds) if k in ("z", "depth", "lev", "level"))
    z = Float64.(ds[zname][:])
    println("SoilGrids depths ($zname): ", z, " ", get(ds[zname].attrib, "units", ""))
    ocd = ds["ocd"]
    println("ocd dims: ", dimnames(ocd), " size ", size(ocd))
    x = Array{Union{Missing, Float64}}(ocd[:, :, :])
    # depths of the level centres, positive down (m)
    d = abs.(z)
    maximum(d) > 10 && (d = d ./ 100)
    order = sortperm(d)
    d = d[order]
    edges = vcat(0.0, [(d[k] + d[k + 1]) / 2 for k in 1:(length(d) - 1)], 2 * d[end] - (d[end] + d[end - 1]) / 2)
    for zmax in (1.0, 2.0)
        column = zeros(length(lon), length(lat))
        for (kk, k) in enumerate(order)
            top, bottom = edges[kk], min(edges[kk + 1], zmax)
            bottom > top || continue
            column .+= coalesce.(x[:, :, k], NaN) .* (bottom - top)
        end
        for (band, total) in banded(column, lon, lat)
            push!(rows, "SOC_0_$(Int(zmax))m_soilgrids,$band,$(total * 1e-12),PgC,soilgrids_ocd_lowres_interp (level edges at midpoints)")
        end
    end
end

# Xu-Saatchi live woody carbon, kg C m⁻²
path = ClimaLand.Artifacts.xusaatchi_woody_carbon_path()
NCDataset(path) do ds
    lon, lat = coords(ds)
    name = first(k for k in keys(ds) if !(k in ("lon", "lat", "longitude", "latitude", "time")))
    x = ds[name][:, :]
    for (band, total) in banded(x, lon, lat)
        push!(rows, "cwood_xusaatchi,$band,$(total * 1e-12),PgC,xusaatchi_woody_carbon ($name) 2000-2019")
    end
end

mkpath(out_dir)
write(joinpath(out_dir, "targets.csv"), join(rows, "\n") * "\n")
println(join(rows, "\n"))
