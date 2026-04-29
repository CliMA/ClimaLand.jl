# Run the AU-Cum-style single-site calibration (calibrate_site.jl) for every
# FLUXNET2015 site, sequentially, in this Julia process. Each site that errors
# is logged and skipped — the loop does not abort.
#
# This is the simple path: useful for smoke tests, small subsets, or running
# overnight on a single workstation/login node. For the central HPC, prefer
# `calibrate_all_fluxnet.sbatch` (Slurm job array) so sites run in parallel.
#
# Usage:
#   julia --project=.buildkite experiments/integrated/generic_site/calibrate_all_fluxnet.jl
#
# Optional ENV vars (forwarded to calibrate_site.jl):
#   CAL_DURATION_DAYS   default "365"
#   N_ITERS             default "5"
#   SITES               comma-separated subset, e.g. "AU-Cum,US-MOz" — if unset,
#                       runs all sites returned by list_fluxnet_sites().

include(joinpath(@__DIR__, "list_fluxnet_sites.jl"))

const CALIBRATE_SCRIPT = joinpath(@__DIR__, "calibrate_site.jl")

function _selected_sites()
    raw = get(ENV, "SITES", "")
    isempty(raw) && return list_fluxnet_sites()
    return strip.(split(raw, ','))
end

function calibrate_all()
    sites = _selected_sites()
    @info "Calibrating $(length(sites)) sites" first = first(sites) last = last(sites)

    failures = Tuple{String,String}[]
    for (i, site) in enumerate(sites)
        @info "[$i/$(length(sites))] === $site ==="
        try
            withenv("SITE_ID" => site, "LOCAL" => "false") do
                Base.include(Main, CALIBRATE_SCRIPT)
            end
        catch err
            msg = sprint(showerror, err)
            @error "Site $site failed" exception = (err, catch_backtrace())
            push!(failures, (site, first(split(msg, '\n'))))
        end
    end

    @info "Done. successes=$(length(sites) - length(failures)) failures=$(length(failures))"
    if !isempty(failures)
        @info "Failures:"
        for (site, msg) in failures
            println("  $site\t$msg")
        end
    end
    return failures
end

if abspath(PROGRAM_FILE) == @__FILE__
    calibrate_all()
end
