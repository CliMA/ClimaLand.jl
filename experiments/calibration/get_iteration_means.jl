#!/usr/bin/env julia
using TOML
using Statistics

iteration = 2
param_names = ["βkx_base", "βkx_coord", "βψx50_base", "βψx50_slope", "βΠR_base", "βΠR_slope"]

values = Dict(name => Float64[] for name in param_names)

for member in 1:13
    member_dir = "land_model/iteration_$(lpad(iteration, 3, '0'))/member_$(lpad(member, 3, '0'))"
    param_file = joinpath(member_dir, "parameters.toml")
    
    if isfile(param_file)
        params = TOML.parsefile(param_file)
        for name in param_names
            if haskey(params, name) && haskey(params[name], "value")
                push!(values[name], params[name]["value"])
            end
        end
    end
end

println("Parameter means from iteration $(lpad(iteration, 3, '0')):")
println("="^60)
for name in param_names
    if !isempty(values[name])
        μ = mean(values[name])
        σ = std(values[name])
        println("  $name:")
        println("    mean = $μ")
        println("    std  = $σ")
    end
end
