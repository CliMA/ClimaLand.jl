"""
Usage:
    Manually specify `sol_name` to be the folder where the implicit
        solution information is stored (output of imex_testing.jl).
    Manually specify  filenames in `files` array, where the first
        file is the explicit solution you want to compare error to.
    Manually specify `dts` to be the timesteps used in each simulation
        in `files`.
    This file will output a .csv file containing the mean and median
        times to solution, and RMSE compared to the explicit solution
        for each implicit solution.
"""

using JLD2
using Statistics
using DelimitedFiles

rmse(v1, v2) = sqrt(mean((v1 .- v2) .^ 2))

sol_name = "ssp333-mod-imp_res-iters_10"
output_file = sol_name * ".csv"
files = [
    "ssp333-mod-exp_res/rre_sol-exp-mod_picard-iters_1-dt_1.jld2",
    sol_name * "/rre_sol-imp-mod_picard-iters_10-dt_1.jld2",
    sol_name * "/rre_sol-imp-mod_picard-iters_10-dt_4.jld2",
    sol_name * "/rre_sol-imp-mod_picard-iters_10-dt_16.jld2",
]

dts = [1; 1; 4; 16; 64]

header = Array{String}(undef, 1, length(files))
errors = Array{Float64}(undef, 1, length(files))
mean_times = Array{Float64}(undef, 1, length(files))
median_times = Array{Float64}(undef, 1, length(files))
labels = ["dt", "RMSE (vs exp soln)", "mean time (ns)", "median time (ns)"]

for i in eachindex(files)
    jldopen(files[i], "r") do file
        if i == 1
            global exp_sol = file["ϑ_l"]
        end
        # compute RMSE compared to explicit solution
        err = rmse(exp_sol, file["ϑ_l"])
        dt = dts[i]

        global header[i] = "$dt"
        global errors[i] = err
        global mean_times[i] = mean(file["benchmark"].times)
        global median_times[i] = median(file["benchmark"].times)
    end
end

out = vcat.([header, errors, mean_times, median_times])
out = hcat.(labels, out)

open(output_file, "w") do io
    writedlm(io, out, ",")
end
