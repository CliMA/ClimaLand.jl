"""
NEON soil-CO₂ depth lookup helper.

Single source of truth for "which measurement depths does a site have, and which
ClimaLand model layer does each map to". The messy parsing (positive-zOffset sign
fix, duplicate-row dedup, across-plot median, nearest-layer matching) lives in the
validated Python notebook

    /home/evametz/Scripts/CliMA/neon_depth_layer_mapping.ipynb

which exports the CSV read here. Regenerate that CSV (rerun the notebook) if the
depth source or the model grid changes.

CSV columns: site, depth_code, z_obs_m, model_layer, z_layer_m, mismatch_m
  - site         : bare NEON code (e.g. "CPER"), matched via _neon_site_key(SITE_ID)
  - depth_code   : "501" / "502" / "503"
  - z_obs_m      : across-plot median sensor depth (m, negative below surface)
  - model_layer  : 1-based ClimaLand subsurface layer index (argmin |z_layer - z_obs|)
  - z_layer_m    : center depth of that model layer (m)
  - mismatch_m   : z_layer_m - z_obs_m (diagnostic)

Requires `_neon_site_key` from site_metadata.jl — include that first.
"""

using DelimitedFiles

const NEON_DEPTH_LOOKUP_CSV =
    get(ENV, "NEON_DEPTH_LOOKUP_CSV",
        "/home/evametz/Scripts/CliMA/neon_depth_layer_mapping.csv")

"""
    neon_depths_for_site(SITE_ID; codes = ["501", "502", "503"])

Return the per-depth lookup for `SITE_ID`, restricted to `codes` and ordered to
match `codes`. Each entry is a named tuple
`(; code::String, z_obs_m::Float64, model_layer::Int)`.

`codes` is the explicit depth-code list the run calibrates on; the returned order
defines the stacking order of the observation vector and the G matrix, so callers
on the obs side and the model side MUST pass the same `codes` in the same order.

Errors if the CSV is missing, the site is absent, or a requested code has no row.
"""
function neon_depths_for_site(SITE_ID; codes = ["501", "502", "503"])
    isfile(NEON_DEPTH_LOOKUP_CSV) ||
        error("NEON depth lookup CSV not found: $NEON_DEPTH_LOOKUP_CSV (rerun the notebook).")

    (data, columns) = DelimitedFiles.readdlm(NEON_DEPTH_LOOKUP_CSV, ','; header = true)
    header = vec(String.(columns))
    i_site = findfirst(==("site"), header)
    i_code = findfirst(==("depth_code"), header)
    i_zobs = findfirst(==("z_obs_m"), header)
    i_layer = findfirst(==("model_layer"), header)
    (i_site === nothing || i_code === nothing || i_zobs === nothing ||
     i_layer === nothing) &&
        error("Unexpected columns in $NEON_DEPTH_LOOKUP_CSV: $(header)")

    key = _neon_site_key(SITE_ID)

    out = NamedTuple[]
    for code in string.(codes)
        row_idx = findfirst(axes(data, 1)) do i
            uppercase(strip(string(data[i, i_site]))) == key &&
                strip(string(data[i, i_code])) == code
        end
        row_idx === nothing && error(
            "No depth-lookup row for site $key (from $SITE_ID), code $code " *
            "in $NEON_DEPTH_LOOKUP_CSV",
        )

        z_obs_m = Float64(data[row_idx, i_zobs])
        model_layer = Int(data[row_idx, i_layer])
        push!(out, (; code = code, z_obs_m = z_obs_m, model_layer = model_layer))
    end
    return out
end
