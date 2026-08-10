# Print (or write) every FLUXNET2015 site ID known to the `fluxnet2015`
# artifact, one per line, sorted. Sized for use as a Slurm/PBS array index
# source (see `calibrate_all_fluxnet.sbatch`, `calibrate_all_fluxnet.pbs`).
#
# The `fluxnet2015` artifact is HPC-only (no download URL), so this only works
# on a machine where it has been pre-staged (central HPC, Derecho).
#
# Usage:
#   # Print to stdout:
#   julia --project=.buildkite experiments/integrated/generic_site/list_fluxnet_sites.jl
#   # Write to a file (avoids shell redirection — useful when copy-pasting):
#   julia --project=.buildkite experiments/integrated/generic_site/list_fluxnet_sites.jl path/to/site_ids.txt
#
# Or programmatically:
#   include("list_fluxnet_sites.jl"); ids = list_fluxnet_sites()

import ClimaLand

"""
    list_fluxnet_sites() -> Vector{String}

Return the sorted list of FLUXNET2015 site IDs found in the `fluxnet2015`
artifact, parsed from the per-site directory names of the form
`FLX_<SITE_ID>_FLUXNET2015_FULLSET_<years>_<version>`.
"""
function list_fluxnet_sites()
    root = ClimaLand.Artifacts.fluxnet2015_data_path()
    site_ids = String[]
    for d in readdir(root)
        isdir(joinpath(root, d)) || continue
        startswith(d, "FLX_") || continue
        parts = split(d, '_')
        length(parts) >= 2 || continue
        push!(site_ids, parts[2])
    end
    return sort!(unique(site_ids))
end

if abspath(PROGRAM_FILE) == @__FILE__
    sites = list_fluxnet_sites()
    if length(ARGS) >= 1
        out = ARGS[1]
        mkpath(dirname(abspath(out)))
        open(out, "w") do io
            for s in sites
                println(io, s)
            end
        end
        @info "Wrote $(length(sites)) site IDs to $out"
    else
        for s in sites
            println(s)
        end
    end
end
