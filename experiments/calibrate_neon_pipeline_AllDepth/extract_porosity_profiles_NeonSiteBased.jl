"""
extract_porosity_profiles_NeonSiteBased.jl

Extract the *exact* soil retention parameter profiles that the CURRENT NEON
calibration runs use (run_pipeline.jl → pipeline/ForwardRun.jl), for every site
referenced in the config/ directory.

This is the companion to `extract_porosity_profiles.jl`. That older script pulled
the retention parameters straight from a global map
(`Soil.soil_vangenuchten_parameters` / `rosetta_soil_vangenuchten_parameters`).
The current pipeline no longer uses those values directly: it starts from that
global map but then OVERRIDES ν, θ_r, K_sat, and the van Genuchten (α, n) with
NEON per-site depth profiles read from CSVs in `Neon/Neon_data/` (see
ForwardRun.jl:210-246). This script reproduces that override path verbatim, so
the values it prints/saves are exactly the ones `land.soil.parameters` holds in a
current run (porosity before the `porosity_scale = 1` no-op multiply).

It reproduces, verbatim, the code path the runs take:
  - site lat/long via the same `_get_neon_site_metadata` helper (src/site_metadata.jl)
  - the same `Column` domain (zlim, nelements, dz_tuple, longlat)  (ForwardRun.jl:154-162)
  - the same base map `Soil.soil_vangenuchten_parameters`          (ForwardRun.jl:214)
  - the same NEON per-site CSV override + `read_neon_profile` interp (ForwardRun.jl:219-246)

Parameters extracted (one CSV each), with the "_NeonSiteBased" suffix so they
never overwrite the older script's outputs:
  - ν     : porosity            -> porosity_profiles_NeonSiteBased.csv
  - θ_r   : residual water content -> theta_r_profiles_NeonSiteBased.csv
  - vg_α  : van Genuchten α (1/m)  -> vg_alpha_profiles_NeonSiteBased.csv
  - vg_n  : van Genuchten n (-)    -> vg_n_profiles_NeonSiteBased.csv
  - K_sat : saturated hydraulic conductivity (m/s) -> ksat_profiles_NeonSiteBased.csv

Each CSV: site_id, lat, long, layer, z_m, <value>   (layers top-down).

COMPARISON: if a matching older-pipeline CSV exists (set COMPARE_SUFFIX below;
default "_Rosetta", the source the current pipeline replaced), this script also
writes `<param>_profiles_NeonSiteBased_vs<COMPARE_SUFFIX>.csv` with, per site &
layer, both values, their difference, and % difference — and prints a per-param,
per-site RMSE / mean-abs-difference summary.

Run:
    julia --project=.buildkite \
        experiments/calibrate_neon_pipeline_AllDepth/extract_porosity_profiles_NeonSiteBased.jl
"""

using ClimaLand
import ClimaLand.Soil as Soil
import ClimaLand.Domains: Column
import ClimaCore
import ClimaComms
ClimaComms.@import_required_backends
using DelimitedFiles
using Printf
using CSV
using DataFrames
# Use the SAME linear_interpolation the pipeline uses (ForwardRun.jl:34): the
# ClimaUtilities.Utils one takes (x, y, xq) positionally and returns the value —
# NOT Interpolations.linear_interpolation, which builds an interpolant object.
import ClimaUtilities.Utils: linear_interpolation
using Statistics: mean

const FT = Float64

# ── NEON per-site retention CSV directory (matches ForwardRun.jl:224) ─────────
const NEON_DIR = "/kiwi-data/Data/groupMembers/evametz/Neon/Neon_data"

# ── Output directory (matches the older script) ───────────────────────────────
const OUTPUT_DIR = "/kiwi-data/Data/groupMembers/evametz/ClimaLand_Output/"

# ── Which older-pipeline output to compare against ────────────────────────────
# The current pipeline's base call is `soil_vangenuchten_parameters` (Gupta),
# but its NEON override is what replaced the old direct-Rosetta values, so the
# most meaningful comparison is usually "_Rosetta". Use "" to compare vs the
# legacy no-suffix (Gupta) files, or "Gupta" for the *Gupta.csv files.
# Set to nothing to skip the comparison entirely.
const COMPARE_SUFFIX = get(ENV, "EXTRACT_COMPARE_SUFFIX", "_Rosetta")

# ── Site metadata helper (lat/long/atmos_h from CSV) — same file ForwardRun uses.
include(joinpath(@__DIR__, "src", "site_metadata.jl"))

# ── read_neon_profile: VERBATIM copy of ForwardRun.jl:68-83 ───────────────────
# Read one NEON per-depth retention parameter into a subsurface field.
# Linear-interpolate within the measured range; hold the deepest measured
# value (flat) below it — matching Rosetta's Interpolations.Flat().
function read_neon_profile(csv_path, colname, space, FT)
    df = CSV.read(csv_path, DataFrame)
    valid = .!ismissing.(df[!, colname])
    z_raw = Float64.(df.depth[valid])
    v_raw = Float64.(df[valid, colname])
    si = sortperm(z_raw)                 # ascending z: most-negative → 0
    z = z_raw[si]
    v = v_raw[si]
    z_bot = z[1]                         # deepest (most negative) measurement
    v_bot = v[1]
    zvals = ClimaCore.Fields.coordinate_field(space).z
    return map(zvals) do zc
        val = zc > z_bot ? linear_interpolation(z, v, zc) : v_bot
        FT(val)
    end
end

# ── Collect the site IDs actually used, from the config/*.toml files ──────────
# NOTE on case: the pipeline uses the config's original-case id (e.g. "NEON-cper")
# both as run.site and — critically — as the NEON CSV column prefix, which is
# LOWERCASE ("NEON-cper_nu_[m3/m3]"). The older extract script and its output
# CSVs use the UPPERCASE id ("NEON-CPER") for the site_id column. So we keep the
# original-case id for CSV *column lookups* (read_neon_profile) and uppercase it
# only for the written site_id / comparison join. Config ids are lowercase, so
# this matches both the live pipeline and the old outputs.
function sites_from_configs(config_dir)
    ids = String[]
    for f in readdir(config_dir; join = true)
        endswith(f, ".toml") || continue
        for line in eachline(f)
            m = match(r"NEON-([A-Za-z]+)", line)
            m === nothing && continue
            push!(ids, "NEON-" * m.captures[1])   # preserve original case
        end
    end
    return sort(unique(ids))
end

# ── Build the exact run domain (mirrors ForwardRun.jl:154-162) ────────────────
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

# Returns lat, long, layer depths, and a Dict param_name => per-layer values,
# built by the CURRENT pipeline path: base van Genuchten map overridden with the
# NEON per-site CSV profiles (ForwardRun.jl:214-246).
function retention_profiles(site_id)
    md = _get_neon_site_metadata(site_id)
    lat = FT(md.lat)
    long = FT(md.long)
    domain = run_column(long, lat)
    sp = domain.space.subsurface

    # Base map (same call the pipeline makes before overriding). We build it to
    # match the pipeline exactly even though every field below is overwritten.
    rp = Soil.soil_vangenuchten_parameters(sp, FT)

    # ── NEON per-site overrides (verbatim source/columns from ForwardRun.jl) ──
    α_field = read_neon_profile(
        joinpath(NEON_DIR, "NEON_all_sites_alpha_1_m_2cm_mean.csv"),
        "$(site_id)_alpha_[1/m]", sp, FT)
    n_field = read_neon_profile(
        joinpath(NEON_DIR, "NEON_all_sites_n_-_2cm_mean.csv"),
        "$(site_id)_n_[-]", sp, FT)

    rp.hydrology_cm .=
        ((α, n) -> ClimaLand.Soil.vanGenuchten{FT}(; α, n)).(α_field, n_field)

    rp.K_sat .= read_neon_profile(
        joinpath(NEON_DIR, "NEON_all_sites_Ksat_m_s_2cm_mean.csv"),
        "$(site_id)_Ksat_[m/s]", sp, FT)
    rp.ν .= read_neon_profile(
        joinpath(NEON_DIR, "NEON_all_sites_nu_m3_m3_2cm_mean.csv"),
        "$(site_id)_nu_[m3/m3]", sp, FT)
    rp.θ_r .= read_neon_profile(
        joinpath(NEON_DIR, "NEON_all_sites_theta_r_m3_m3_2cm_mean.csv"),
        "$(site_id)_theta_r_[m3/m3]", sp, FT)

    z = ClimaCore.Fields.coordinate_field(sp).z
    z_vals = parent(z)[:, 1]

    hcm = rp.hydrology_cm
    vg_alpha = parent(map(c -> c.α, hcm))[:, 1]
    vg_n = parent(map(c -> c.n, hcm))[:, 1]

    vals = Dict(
        "nu" => parent(rp.ν)[:, 1],
        "theta_r" => parent(rp.θ_r)[:, 1],
        "vg_alpha" => vg_alpha,
        "vg_n" => vg_n,
        "ksat" => parent(rp.K_sat)[:, 1],
    )
    return (; lat, long, z_vals, vals)
end

# ── Main ──────────────────────────────────────────────────────────────────────
config_dir = joinpath(@__DIR__, "config")
site_ids = sites_from_configs(config_dir)
println("Sites found in configs ($(length(site_ids))): ", join(site_ids, ", "))
println()
println("Retention source: NEON per-site profiles (current pipeline path)")
println("NEON data dir:     $NEON_DIR")
println()

# param key => (output filename, header column name).
outputs = [
    ("nu", "porosity_profiles_NeonSiteBased.csv", "nu"),
    ("theta_r", "theta_r_profiles_NeonSiteBased.csv", "theta_r"),
    ("vg_alpha", "vg_alpha_profiles_NeonSiteBased.csv", "vg_alpha_1_per_m"),
    ("vg_n", "vg_n_profiles_NeonSiteBased.csv", "vg_n"),
    ("ksat", "ksat_profiles_NeonSiteBased.csv", "ksat_m_per_s"),
]

# Initialize each output with a header row.
rows = Dict(key => Vector{Any}[Any["site_id", "lat", "long", "layer", "z_m", col]]
            for (key, _, col) in outputs)

# Keep per-site per-layer values in memory for the comparison step below.
# extracted[key][site_id] => Dict(layer => value)
extracted = Dict(key => Dict{String, Dict{Int, Float64}}() for (key, _, _) in outputs)

for site_id in site_ids
    # site_id (original case, e.g. "NEON-cper") is used for the NEON CSV column
    # lookup inside retention_profiles; display_id (uppercase) is what we write
    # to site_id and join on for the comparison, matching the old outputs.
    display_id = uppercase(site_id)
    local prof
    try
        prof = retention_profiles(site_id)
    catch err
        @warn "Failed for $site_id" exception = err
        continue
    end
    nlayer = length(prof.z_vals)
    println("── $display_id  (lat=$(prof.lat), long=$(prof.long)) ──")
    @printf("    %-6s %10s %10s %10s %10s\n", "layer", "z_m", "nu", "theta_r", "vg_n")
    for i in 1:nlayer
        j = nlayer - i + 1  # surface (top) first
        @printf("    %-6d %10.4f %10.4f %10.4f %10.4f\n",
            j, prof.z_vals[j], prof.vals["nu"][j], prof.vals["theta_r"][j],
            prof.vals["vg_n"][j])
        for (key, _, _) in outputs
            v = round(prof.vals[key][j]; digits = 8)
            push!(rows[key], Any[display_id, prof.lat, prof.long, j,
                round(prof.z_vals[j]; digits = 5), v])
            get!(extracted[key], display_id, Dict{Int, Float64}())[j] = Float64(v)
        end
    end
    println()
end

for (key, fname, _) in outputs
    out_csv = joinpath(OUTPUT_DIR, fname)
    open(out_csv, "w") do io
        writedlm(io, rows[key], ',')
    end
    println("Wrote: $out_csv")
end

# ── Comparison against the older-pipeline outputs ──────────────────────────────
# The older script wrote <param>_profiles<COMPARE_SUFFIX>.csv with the same
# schema (site_id, lat, long, layer, z_m, <value>). We join on (site_id, layer)
# and report the difference. Layers match because both use the identical domain.
if COMPARE_SUFFIX !== nothing
    println()
    println("="^72)
    println("Comparison vs older pipeline: <param>_profiles$(COMPARE_SUFFIX).csv")
    println("="^72)

    # Map each output key to the OLD filename that holds it.
    old_files = Dict(
        "nu"       => "porosity_profiles$(COMPARE_SUFFIX).csv",
        "theta_r"  => "theta_r_profiles$(COMPARE_SUFFIX).csv",
        "vg_alpha" => "vg_alpha_profiles$(COMPARE_SUFFIX).csv",
        "vg_n"     => "vg_n_profiles$(COMPARE_SUFFIX).csv",
        # ksat has no older-script counterpart; skipped in comparison.
    )

    for (key, _, col) in outputs
        haskey(old_files, key) || continue
        old_path = joinpath(OUTPUT_DIR, old_files[key])
        if !isfile(old_path)
            @warn "No older file to compare for '$key' at $old_path — skipping"
            continue
        end
        old_df = CSV.read(old_path, DataFrame)
        # Old value column is the last column; grab by position for robustness.
        old_valcol = names(old_df)[end]

        # Build old lookup: (site_id, layer) => value.
        old_lookup = Dict{Tuple{String, Int}, Float64}()
        for r in eachrow(old_df)
            old_lookup[(String(r.site_id), Int(r.layer))] = Float64(r[old_valcol])
        end

        cmp_rows = Vector{Any}[Any["site_id", "layer", "z_m",
            "new_$(col)", "old_$(col)", "diff_new_minus_old", "pct_diff"]]

        # Per-site accumulators for the summary.
        site_sqerr = Dict{String, Vector{Float64}}()

        # Reconstruct z per (site, layer) from the freshly written rows.
        z_lookup = Dict{Tuple{String, Int}, Float64}()
        for row in rows[key][2:end]   # skip header
            z_lookup[(String(row[1]), Int(row[4]))] = Float64(row[5])
        end

        for (site_id, layers) in sort(collect(extracted[key]); by = first)
            for (layer, newv) in sort(collect(layers); by = first)
                haskey(old_lookup, (site_id, layer)) || continue
                oldv = old_lookup[(site_id, layer)]
                diff = newv - oldv
                pct = oldv == 0 ? NaN : 100 * diff / oldv
                z = get(z_lookup, (site_id, layer), NaN)
                push!(cmp_rows, Any[site_id, layer, round(z; digits = 5),
                    round(newv; digits = 8), round(oldv; digits = 8),
                    round(diff; digits = 8),
                    isnan(pct) ? "NaN" : round(pct; digits = 4)])
                push!(get!(site_sqerr, site_id, Float64[]), diff^2)
            end
        end

        cmp_name = replace(old_files[key], COMPARE_SUFFIX * ".csv" => "") *
                   "_NeonSiteBased_vs$(COMPARE_SUFFIX).csv"
        cmp_path = joinpath(OUTPUT_DIR, cmp_name)
        open(cmp_path, "w") do io
            writedlm(io, cmp_rows, ',')
        end

        # Summary: overall + per-site RMSE / mean-abs-diff.
        all_sq = reduce(vcat, values(site_sqerr); init = Float64[])
        if isempty(all_sq)
            println("\n[$col] no overlapping (site, layer) rows with $(old_files[key])")
        else
            overall_rmse = sqrt(mean(all_sq))
            println("\n[$col]  overall RMSE(new-old) = $(round(overall_rmse; digits = 6))" *
                    "   (over $(length(all_sq)) layer-values)")
            @printf("    %-12s %12s %12s\n", "site", "rmse", "n_layers")
            for site_id in sort(collect(keys(site_sqerr)))
                sq = site_sqerr[site_id]
                @printf("    %-12s %12.6f %12d\n",
                    site_id, sqrt(mean(sq)), length(sq))
            end
        end
        println("  Wrote comparison: $cmp_path")
    end
end

println()
println("Done.")
