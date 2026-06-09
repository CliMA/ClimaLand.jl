"""
extract_porosity_profiles.jl

Extract the *exact* soil porosity (ν) profile that the NEON calibration runs use,
for every site referenced in the config/ directory.

It reproduces, verbatim, the code path the runs take:
  - site lat/long via the same `_get_neon_site_metadata` helper (ForwardRun.jl)
  - the same `Column` domain (zlim, nelements, dz_tuple, longlat)  (ForwardRun.jl:133-140)
  - ClimaLand's own `Soil.soil_vangenuchten_parameters(subsurface, FT)`, which
    reads + regrids the Gupta et al. (2020) porosity map (the default
    `retention_parameters` inside `Soil.EnergyHydrology` / `LandModel`).

So the ν values printed/saved here are the same ones `land.soil.parameters.ν`
holds in a run (before the `porosity_scale = 1` no-op multiply). No approximation.

Run:
    julia --project=.buildkite \
        experiments/calibrate_neon_pipeline/extract_porosity_profiles.jl

Output:
    experiments/calibrate_neon_pipeline/porosity_profiles.csv
    (columns: site_id, layer, z_m, nu)
"""

using ClimaLand
import ClimaLand.Soil as Soil
import ClimaLand.Domains: Column
import ClimaCore
import ClimaComms
ClimaComms.@import_required_backends
using DelimitedFiles
using Printf

const FT = Float64

# Site metadata helper (lat/long/atmos_h from CSV) — the same file ForwardRun uses.
include(joinpath(pkgdir(ClimaLand), "experiments/calibrate_neon/site_metadata.jl"))

# ── Collect the site IDs actually used, from the config/*.toml files ──────────
function sites_from_configs(config_dir)
    ids = String[]
    for f in readdir(config_dir; join = true)
        endswith(f, ".toml") || continue
        for line in eachline(f)
            m = match(r"NEON-([A-Za-z]+)", line)
            m === nothing && continue
            push!(ids, "NEON-" * uppercase(m.captures[1]))
        end
    end
    return sort(unique(ids))
end

# ── Build the exact run domain (mirrors ForwardRun.jl:133-140) ────────────────
function run_column(long, lat)
    dz_bottom = FT(2)
    dz_top = FT(0.038)
    dz_tuple = (dz_bottom, dz_top)
    nelements = 24
    zmin = FT(-6.2)
    zmax = FT(0)
    return Column(; zlim = (zmin, zmax), nelements = nelements,
        dz_tuple = dz_tuple, longlat = (long, lat))
end

function porosity_profile(site_id)
    md = _get_neon_site_metadata(site_id)
    lat = FT(md.lat)
    long = FT(md.long)
    domain = run_column(long, lat)

    # Identical to the run's default retention_parameters source.
    rp = Soil.soil_vangenuchten_parameters(domain.space.subsurface, FT)
    ν = rp.ν

    z = ClimaCore.Fields.coordinate_field(domain.space.subsurface).z
    z_vals = parent(z)[:, 1]
    ν_vals = parent(ν)[:, 1]
    return (; lat, long, z_vals, ν_vals)
end

# ── Main ──────────────────────────────────────────────────────────────────────
config_dir = joinpath(@__DIR__, "config")
site_ids = sites_from_configs(config_dir)
println("Sites found in configs ($(length(site_ids))): ", join(site_ids, ", "))
println()

rows = Vector{Any}[]
push!(rows, Any["site_id", "lat", "long", "layer", "z_m", "nu"])

for site_id in site_ids
    local prof
    try
        prof = porosity_profile(site_id)
    catch err
        @warn "Failed for $site_id" exception = err
        continue
    end
    nlayer = length(prof.z_vals)
    println("── $site_id  (lat=$(prof.lat), long=$(prof.long)) ──")
    for i in 1:nlayer
        # print top-down (surface first) for readability
        j = nlayer - i + 1
        @printf("    layer %2d  z = %8.4f m   ν = %.4f\n", j, prof.z_vals[j], prof.ν_vals[j])
        push!(rows, Any[site_id, prof.lat, prof.long, j,
            round(prof.z_vals[j]; digits = 5), round(prof.ν_vals[j]; digits = 5)])
    end
    println()
end

out_csv = joinpath(@__DIR__, "porosity_profiles.csv")
open(out_csv, "w") do io
    writedlm(io, rows, ',')
end
println("Wrote: $out_csv")
