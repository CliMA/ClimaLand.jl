using JLD2

"""
Inspect the CWD field data to understand its structure.
"""
function inspect_cwd_data(cwd_file::String)
    @info "Loading CWD field from" cwd_file
    
    # Load the data
    data = JLD2.load(cwd_file)
    lons = data["lons"]
    lats = data["lats"]
    CWD_spatial = data["CWD_spatial"]
    
    println("\n=== Data Structure ===")
    println("Longitude range: $(minimum(lons)) to $(maximum(lons))")
    println("Latitude range: $(minimum(lats)) to $(maximum(lats))")
    println("Number of longitudes: $(length(lons))")
    println("Number of latitudes: $(length(lats))")
    println("CWD array size: $(size(CWD_spatial))")
    println("\n=== CWD Statistics ===")
    println("Minimum: $(minimum(CWD_spatial))")
    println("Maximum: $(maximum(CWD_spatial))")
    println("Mean: $(sum(CWD_spatial)/length(CWD_spatial))")
    println("NaN count: $(sum(isnan.(CWD_spatial)))")
    
    println("\n=== First few coordinates ===")
    println("First 5 lons: $(lons[1:min(5, length(lons))])")
    println("First 5 lats: $(lats[1:min(5, length(lats))])")
end

if abspath(PROGRAM_FILE) == @__FILE__
    cwd_file = joinpath(@__DIR__, "cwd_era5_2015_2023.jld2")
    inspect_cwd_data(cwd_file)
end