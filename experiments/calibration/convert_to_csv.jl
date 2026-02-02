using JLD2

# Load data
data = JLD2.load("calculated_traits_with_uncertainty.jld2")

# Write to CSV manually
open("traits_with_uncertainty.csv", "w") do io
    # Header
    println(io, "lat,lon,psi50_mean,psi50_std,kx_mean,kx_std,pr_mean,pr_std,aridity")
    
    # Data
    n = length(data["lats"])
    for i in 1:n
        println(io, "$(data["lats"][i]),$(data["lons"][i]),$(data["ψx50_mean"][i]),$(data["ψx50_std"][i]),$(data["kx_mean"][i]),$(data["kx_std"][i]),$(data["ΠR_mean"][i]),$(data["ΠR_std"][i]),$(data["aridity"][i])")
    end
end

println("✓ Converted to CSV: traits_with_uncertainty.csv")
