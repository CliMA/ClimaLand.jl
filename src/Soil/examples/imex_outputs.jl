"""
Usage: pass reference explicit solution as first file to compare
second file to it.
"""

using JLD2
using Statistics


files = ARGS

jldopen(files[1], "r") do file
    global exp_res = file["ϑ_l"]
    global exp_z = file["z"]
    global exp_b = file["benchmark"]
end

jldopen(files[2], "r") do file
    global imp_res = file["ϑ_l"]
    global imp_z = file["z"]
    global imp_b = file["benchmark"]
end

# check error
@assert maximum(abs.(imp_res .- exp_res)) < 1e-3
@show mean(exp_b.times), median(exp_b.times)
@show mean(imp_b.times), median(imp_b.times)


# jldopen(f, "r") do file
#     global exp_res = file["ϑ_l"]
#     global exp_z = file["z"]
#     global exp_b = file["benchmark"]
# end
