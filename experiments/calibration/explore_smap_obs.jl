using JLD2

# Load the SMAP observation file
filepath = "land_observation_vector.jld2"
data = JLD2.load(filepath)

println("=" ^ 80)
println("SMAP Observation File Structure")
println("=" ^ 80)
println()

# Show all top-level keys
println("Top-level keys:")
for key in sort(collect(keys(data)))
    println("  • $key")
end
println()

# Examine each key in detail
for key in sort(collect(keys(data)))
    println("-" ^ 80)
    println("Key: $key")
    println("-" ^ 80)
    
    value = data[key]
    println("Type: $(typeof(value))")
    
    if value isa Array
        println("Array size: $(size(value))")
        println("Element type: $(eltype(value))")
        if !isempty(value)
            println("First element type: $(typeof(value[1]))")
            if length(value) <= 5
                println("Contents:")
                for (i, elem) in enumerate(value)
                    println("  [$i]: $elem")
                end
            else
                println("First 3 elements:")
                for i in 1:3
                    println("  [$i]: $(value[i])")
                end
            end
        end
    elseif value isa String || value isa Number
        println("Value: $value")
    else
        println("Value: $value")
    end
    println()
end

# Deep dive into observation_vector
println("=" ^ 80)
println("Deep Dive: observation_vector")
println("=" ^ 80)

if haskey(data, "observation_vector")
    obs_vec = data["observation_vector"]
    println("Number of periods: $(length(obs_vec))")
    println()
    
    for (i, obs) in enumerate(obs_vec)
        println("Period $i:")
        println("  Type: $(typeof(obs))")
        
        # Check if it's an EKP.Observation object
        if hasproperty(obs, :samples)
            println("  samples: $(typeof(obs.samples))")
            if obs.samples isa Array && !isempty(obs.samples)
                println("    Number of sample sets: $(length(obs.samples))")
                for (j, sample) in enumerate(obs.samples)
                    println("    Sample set $j:")
                    println("      Type: $(typeof(sample))")
                    println("      Length: $(length(sample))")
                    if length(sample) > 0
                        println("      First few values: $(sample[1:min(5, length(sample))])")
                        println("      Value range: [$(minimum(sample)), $(maximum(sample))]")
                        n_nan = sum(isnan.(sample))
                        println("      NaN count: $n_nan / $(length(sample))")
                    end
                end
            end
        end
        
        if hasproperty(obs, :covariances)
            println("  covariances: $(typeof(obs.covariances))")
            if obs.covariances isa Array && !isempty(obs.covariances)
                println("    Number of covariance matrices: $(length(obs.covariances))")
                for (j, cov) in enumerate(obs.covariances)
                    println("    Covariance $j:")
                    println("      Type: $(typeof(cov))")
                    println("      Size: $(size(cov))")
                end
            end
        end
        
        if hasproperty(obs, :names)
            println("  names: $(obs.names)")
        end
        println()
    end
end

println("=" ^ 80)
println("Summary Statistics")
println("=" ^ 80)

if haskey(data, "observation_vector")
    obs_vec = data["observation_vector"]
    total_obs = 0
    for obs in obs_vec
        if hasproperty(obs, :samples) && obs.samples isa Array
            for sample in obs.samples
                total_obs += length(sample)
            end
        end
    end
    println("Total observation count: $total_obs")
end

if haskey(data, "sample_date_ranges")
    println("Date ranges:")
    for (i, range) in enumerate(data["sample_date_ranges"])
        println("  Period $i: $(range[1]) to $(range[2])")
    end
end

println()
println("File created: $(data["created"])")
