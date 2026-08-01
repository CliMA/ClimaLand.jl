## Extract SoilGrids soil properties at each battery site.
##
## The climate-only allocation failed because mean annual precipitation
## anti-correlates with the target for the two hardest cases: the wettest of the
## four critical sites is a grassland and the two driest are boreal forests.
##
## Soil properties are the next candidate - non-PFT, SoilGrids-derived, already
## in the model. This script exists to run the same ordering check FIRST, before
## any of it is wired into the allocation: does any soil property separate
## grassland from forest at similar climate, without inverting the boreal end?
##
## No mechanism is implemented here on purpose. If the ordering fails the same
## way precipitation did, that is the finding.
##
## Usage:  julia --project=.buildkite soil_properties.jl

import ClimaComms
ClimaComms.@import_required_backends
using ClimaCore
using ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.Soil

const FT = Float64

read_sites(csv) = [
    (name = f[4], biome = f[3], lon = parse(FT, f[1]), lat = parse(FT, f[2])) for f in (split(strip(l), ',') for l in readlines(csv)) if
    length(f) == 4 && f[1] != "lon" && !startswith(f[1], "#")
]

"""
    top_mean(field, space)

Mean of a subsurface field over the top ~30 cm - the depth SoilGrids texture is
most meaningful over, and the depth that matters for whether a column can
establish woody vegetation.
"""
function top_mean(field, z)
    vals = Array(parent(field))
    zs = Array(parent(z))
    sel = vec(zs) .> -0.3
    any(sel) || return sum(vec(vals)) / length(vec(vals))
    return sum(vec(vals)[sel]) / count(sel)
end

function main()
    sites = read_sites(joinpath(@__DIR__, "test_sites.csv"))
    println(
        rpad("site", 21),
        rpad("biome", 27),
        rpad("org", 8),
        rpad("sand", 8),
        rpad("gravel", 8),
        rpad("porosity", 10),
        "Ksat (m/s)",
    )
    for s in sites
        lon = s.lon == 0 ? FT(1e-3) : s.lon
        lat = s.lat == 0 ? FT(1e-3) : s.lat
        domain = ClimaLand.Domains.Column(;
            zlim = FT.((-15, 0)),
            longlat = (lon, lat),
            nelements = 15,
            dz_tuple = FT.((3, 0.05)),
        )
        sub = domain.space.subsurface
        z = domain.fields.z
        comp = Soil.soil_composition_parameters(sub, FT)
        vg = Soil.soil_vangenuchten_parameters(sub, FT)
        println(
            rpad(s.name, 21),
            rpad(s.biome, 27),
            rpad(round(top_mean(comp.ν_ss_om, z), digits = 4), 8),
            rpad(round(top_mean(comp.ν_ss_quartz, z), digits = 4), 8),
            rpad(round(top_mean(comp.ν_ss_gravel, z), digits = 4), 8),
            rpad(round(top_mean(vg.ν, z), digits = 4), 10),
            round(top_mean(vg.K_sat, z), sigdigits = 3),
        )
    end
end

main()
