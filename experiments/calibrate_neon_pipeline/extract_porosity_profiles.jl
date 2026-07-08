"""
extract_porosity_profiles.jl

Extract the *exact* soil retention parameter profiles that the NEON calibration
runs use, for every site referenced in the config/ directory.

It reproduces, verbatim, the code path the runs take:
  - site lat/long via the same `_get_neon_site_metadata` helper (ForwardRun.jl)
  - the same `Column` domain (zlim, nelements, dz_tuple, longlat)  (ForwardRun.jl:133-140)
  - the same retention-parameter source as the run (see USE_ROSETTA below).

Two retention-parameter sources are supported, selectable via the USE_ROSETTA
switch (or the EXTRACT_USE_ROSETTA env var):
  - USE_ROSETTA = false : `Soil.soil_vangenuchten_parameters(subsurface, FT)`
        reads + regrids the Gupta et al. (2020) maps. (legacy default)
  - USE_ROSETTA = true  : `Soil.rosetta_soil_vangenuchten_parameters(subsurface, FT)`
        the Rosetta-based parameters, matching the current ForwardRun.jl:193 run
        path. Output filenames get a "_Rosetta" suffix.

So the values printed/saved here are the same ones `land.soil.parameters` holds
in a run (porosity before the `porosity_scale = 1` no-op multiply). No approximation.

Parameters extracted (one CSV each):
  - ν     : porosity            -> porosity_profiles.csv
  - θ_r   : residual water content -> theta_r_profiles.csv
  - vg_α  : van Genuchten α (1/m)  -> vg_alpha_profiles.csv
  - vg_n  : van Genuchten n (-)    -> vg_n_profiles.csv

Each CSV: site_id, lat, long, layer, z_m, <value>   (layers top-down).

Run:
    julia --project=.buildkite \
        experiments/calibrate_neon_pipeline/extract_porosity_profiles.jl
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

# ── Retention-parameter source ────────────────────────────────────────────────
# true  -> Rosetta (matches current ForwardRun.jl:193), outputs get "_Rosetta".
# false -> Gupta et al. (2020) (legacy default), outputs have no suffix.
# Override at runtime with:  EXTRACT_USE_ROSETTA=true julia ... extract_porosity_profiles.jl
const USE_ROSETTA = lowercase(get(ENV, "EXTRACT_USE_ROSETTA", "true")) in ("true", "1", "yes")

retention_params(subsurface) =
    USE_ROSETTA ?
    Soil.rosetta_soil_vangenuchten_parameters(subsurface, FT) :
    Soil.soil_vangenuchten_parameters(subsurface, FT)

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

# Returns lat, long, layer depths, and a Dict param_name => per-layer values.
function retention_profiles(site_id)
    md = _get_neon_site_metadata(site_id)
    lat = FT(md.lat)
    long = FT(md.long)
    domain = run_column(long, lat)

    # Identical to the run's retention_parameters source (Gupta or Rosetta).
    rp = retention_params(domain.space.subsurface)

    z = ClimaCore.Fields.coordinate_field(domain.space.subsurface).z
    z_vals = parent(z)[:, 1]

    # hydrology_cm is a Field of vanGenuchten structs; pull α and n per layer.
    hcm = rp.hydrology_cm
    vg_alpha = parent(map(c -> c.α, hcm))[:, 1]
    vg_n = parent(map(c -> c.n, hcm))[:, 1]

    vals = Dict(
        "nu" => parent(rp.ν)[:, 1],
        "theta_r" => parent(rp.θ_r)[:, 1],
        "vg_alpha" => vg_alpha,
        "vg_n" => vg_n,
    )
    return (; lat, long, z_vals, vals)
end

# ── Main ──────────────────────────────────────────────────────────────────────
config_dir = joinpath(@__DIR__, "config")
site_ids = sites_from_configs(config_dir)
println("Sites found in configs ($(length(site_ids))): ", join(site_ids, ", "))
println()

# param key => (output filename, header column name).
# Filenames get a "_Rosetta" suffix when USE_ROSETTA is set, so the two sources
# never overwrite each other.
suffix = USE_ROSETTA ? "_Rosetta" : "Gupta"
outputs = [
    ("nu", "porosity_profiles$(suffix).csv", "nu"),
    ("theta_r", "theta_r_profiles$(suffix).csv", "theta_r"),
    ("vg_alpha", "vg_alpha_profiles$(suffix).csv", "vg_alpha_1_per_m"),
    ("vg_n", "vg_n_profiles$(suffix).csv", "vg_n"),
]
println("Retention source: ", USE_ROSETTA ? "Rosetta" : "Gupta et al. (2020)")

# Initialize each output with a header row.
rows = Dict(key => Vector{Any}[Any["site_id", "lat", "long", "layer", "z_m", col]]
            for (key, _, col) in outputs)

for site_id in site_ids
    local prof
    try
        prof = retention_profiles(site_id)
    catch err
        @warn "Failed for $site_id" exception = err
        continue
    end
    nlayer = length(prof.z_vals)
    println("── $site_id  (lat=$(prof.lat), long=$(prof.long)) ──")
    @printf("    %-6s %10s %10s %10s %10s\n", "layer", "z_m", "nu", "theta_r", "vg_n")
    for i in 1:nlayer
        j = nlayer - i + 1  # surface (top) first
        @printf("    %-6d %10.4f %10.4f %10.4f %10.4f\n",
            j, prof.z_vals[j], prof.vals["nu"][j], prof.vals["theta_r"][j],
            prof.vals["vg_n"][j])
        for (key, _, _) in outputs
            push!(rows[key], Any[site_id, prof.lat, prof.long, j,
                round(prof.z_vals[j]; digits = 5),
                round(prof.vals[key][j]; digits = 6)])
        end
    end
    println()
end

for (key, fname, _) in outputs
    out_csv = joinpath("/kiwi-data/Data/groupMembers/evametz/ClimaLand_Output/", fname)
    open(out_csv, "w") do io
        writedlm(io, rows[key], ',')
    end
    println("Wrote: $out_csv")
end
