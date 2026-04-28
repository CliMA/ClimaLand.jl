# Smoke-test FluxnetSimulations.generic_site_simulation against the four
# bundled (`fluxnet_sites` artifact) sites: US-MOz, US-NR1, US-Ha1, US-Var.
#
# Each site is run for 7 days using local-mode coordinates (from the per-site
# Val{} dispatchers). Pass = `solve!` returns without throwing.
#
# Usage:
#   julia --project=experiments experiments/integrated/generic_site/test_local_artifact_sites.jl

import ClimaComms
ClimaComms.@import_required_backends
using Dates
using ClimaLand
import ClimaLand.FluxnetSimulations as FluxnetSimulations
using ClimaLand.Simulations: solve!

const FT = Float64
const SITES = ["US-MOz", "US-NR1", "US-Ha1", "US-Var"]
const DURATION = Day(parse(Int, get(ENV, "SIM_DURATION_DAYS", "7")))

results = Dict{String, Symbol}()
errors = Dict{String, Any}()

for site_ID in SITES
    @info "==> $site_ID"
    try
        site_val = FluxnetSimulations.replace_hyphen(site_ID)
        (; time_offset, lat, long) =
            FluxnetSimulations.get_location(FT, Val(site_val))
        (; atmos_h) = FluxnetSimulations.get_fluxtower_height(FT, Val(site_val))

        result = FluxnetSimulations.generic_site_simulation(
            site_ID;
            lat,
            long,
            time_offset,
            atmos_h,
            duration = DURATION,
        )
        @time solve!(result.simulation)
        results[site_ID] = :ok
    catch e
        results[site_ID] = :failed
        errors[site_ID] = e
        @error "Failure at $site_ID" exception = (e, catch_backtrace())
    end
end

println()
println("=== Summary ===")
for site_ID in SITES
    println("  $site_ID: $(results[site_ID])")
end

n_failed = count(==(:failed), values(results))
if n_failed > 0
    println("\n$n_failed site(s) failed.")
    exit(1)
else
    println("\nAll $(length(SITES)) sites passed.")
end
