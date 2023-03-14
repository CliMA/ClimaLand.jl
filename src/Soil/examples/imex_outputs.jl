"""
Usage: pass reference explicit solution as first file to compare
second file to it.
"""

using JLD2
using Statistics
using DelimitedFiles

rmse(v1, v2) = sqrt(mean((v1 .- v2) .^ 2))

files = [
    "exp_res/rre_sol-exp-true_picard-iters_2-dt_0p25.jld2",
    "imp_res-iters_2/rre_sol-imp-true_picard-iters_2-dt_0p25.jld2",
    "imp_res-iters_2/rre_sol-imp-true_picard-iters_2-dt_1.jld2",
    "imp_res-iters_2/rre_sol-imp-true_picard-iters_2-dt_4.jld2",
    "imp_res-iters_2/rre_sol-imp-true_picard-iters_2-dt_16.jld2",
    "imp_res-iters_2/rre_sol-imp-true_picard-iters_2-dt_64.jld2",
    "imp_res-iters_2/rre_sol-imp-true_picard-iters_2-dt_82.jld2",
    "imp_res-iters_2/rre_sol-imp-true_picard-iters_2-dt_96.jld2",
    "imp_res-iters_2/rre_sol-imp-true_picard-iters_2-dt_112.jld2",
    "imp_res-iters_2/rre_sol-imp-true_picard-iters_2-dt_120.jld2",
]

# dts for iters_2, iters_10:
dts = [0.25; 0.25; 1; 4; 16; 64; 82; 96; 112; 120]
# dts for iters_1:
# dts = [0.25; 0.25; 1; 4; 16; 64; 128; 136; 144; 152; 160; 168; 176]

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


output_file = "imp_res-iters_2.csv"
open(output_file, "w") do io
    writedlm(io, out, ",")
end

# jldopen(files[2], "r") do file
#     global imp_res1 = file["ϑ_l"]
#     global imp_z1 = file["z"]
#     global imp_b1 = file["benchmark"]
# end

# check error
# @assert maximum(abs.(imp_res .- exp_res)) < 1e-3
# @show mean(exp_b.times), median(exp_b.times)
# @show mean(imp_b.times), median(imp_b.times)
