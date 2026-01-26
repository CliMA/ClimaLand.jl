using JLD2
using Statistics

println("Loading trait data from calculated_traits.jld2...")
data = load("calculated_traits.jld2")

# Extract data
lons = data["lons"]
lats = data["lats"]
ψx50 = data["ψx50"]
kx = data["kx"]
ΠR = data["ΠR"]

println("Total points: $(length(lons))")

# Filter out NaN values for statistics
valid_mask = .!isnan.(ψx50)
n_land = count(valid_mask)
n_ocean = count(.!valid_mask)

println("Land points: $n_land")
println("Ocean points: $n_ocean")

println("\nψx50 (MPa) statistics (land only):")
valid_ψx50 = ψx50[valid_mask]
println("  min: $(minimum(valid_ψx50)), max: $(maximum(valid_ψx50))")
println("  mean: $(mean(valid_ψx50)), median: $(median(valid_ψx50))")

println("\nkx statistics (land only):")
valid_kx = kx[valid_mask]
println("  min: $(minimum(valid_kx)), max: $(maximum(valid_kx))")
println("  mean: $(mean(valid_kx)), median: $(median(valid_kx))")

println("\nΠR statistics (land only):")
valid_ΠR = ΠR[valid_mask]
println("  min: $(minimum(valid_ΠR)), max: $(maximum(valid_ΠR))")
println("  mean: $(mean(valid_ΠR)), median: $(median(valid_ΠR))")

# Export for Python visualization
using DelimitedFiles

# Export ψx50
println("\nExporting ψx50 data...")
open("traits_psi50_data.csv", "w") do io
    println(io, "lon,lat,psi50")
    for i in 1:length(lons)
        if valid_mask[i]
            println(io, "$(lons[i]),$(lats[i]),$(ψx50[i])")
        end
    end
end

# Export kx
println("Exporting kx data...")
open("traits_kx_data.csv", "w") do io
    println(io, "lon,lat,kx")
    for i in 1:length(lons)
        if valid_mask[i]
            println(io, "$(lons[i]),$(lats[i]),$(kx[i])")
        end
    end
end

# Export ΠR
println("Exporting ΠR data...")
open("traits_pr_data.csv", "w") do io
    println(io, "lon,lat,pr")
    for i in 1:length(lons)
        if valid_mask[i]
            println(io, "$(lons[i]),$(lats[i]),$(ΠR[i])")
        end
    end
end

println("\n✓ Exported CSV files for Python visualization:")
println("  - traits_psi50_data.csv")
println("  - traits_kx_data.csv")
println("  - traits_pr_data.csv")
