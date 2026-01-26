using Dates
using ClimaAnalysis

include("data_sources.jl")

# Get ERA5 observation data
era5_obs = get_era5_obs_var_dict()
lhf_var = era5_obs["lhf"]

# Check available dates
dates = lhf_var.dims["time"]
println("ERA5 LHF Observation Data:")
println("  Total dates available: $(length(dates))")
println("  First date: $(dates[1])")
println("  Last date: $(dates[end])")
println("  Date interval: $(dates[2] - dates[1])")

# Check if specific dates exist
test_dates = [
    DateTime(2015, 6, 1),
    DateTime(2015, 6, 15),
    DateTime(2015, 6, 30),
    DateTime(2015, 7, 1),
]

println("\nTesting specific dates:")
for d in test_dates
    exists = d in dates
    println("  $d: $(exists ? "✅ EXISTS" : "❌ NOT FOUND")")
end

# Find nearest dates to 2015-06-30
idx = searchsortedfirst(dates, DateTime(2015, 6, 30))
println("\nNearest dates to 2015-06-30:")
if idx > 1
    println("  Before: $(dates[idx-1])")
end
if idx <= length(dates)
    println("  At/After: $(dates[idx])")
end