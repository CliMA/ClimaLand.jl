#= run these lines for test / development

ARGS = ["US-MOz"]
include("integrated/fluxnet/setup.jl")
include("integrated/fluxnet/run_fluxnet.jl")
include("integrated/fluxnet/inputs_dataframe.jl")

=#

"""
    getoutput(variable::Symbol, variables::Symbol...; result = sv.saveval, depth = 1)

Return a vector of FT corresponding to the variable of interest at all times.
By default, get output from sv.saveval, but user can specify e.g., result = sol.u
By default, get surface value, but user can specify depth for e.g., soil temperature
"""
function getoutput(
    variable::Symbol,
    variables::Symbol...;
    result = sv.saveval,
    depth = 1,
)
    for v in (variable, variables...)
        result = getproperty.(result, v)
    end
    return [parent(r)[depth] for r in result]
end
# example: getoutput(:SW_out)
# example 2: getoutput(:canopy, :conductance, :gs)
# example 3: getoutput(:soil, :T; result = sol.u, depth = 2)

# List of output that we want
output_list = (
    (:SW_out,),
    (:LW_out,),
    (:canopy, :conductance, :gs),
    (:canopy, :conductance, :transpiration),
    (:canopy, :autotrophic_respiration, :Ra),
    (:canopy, :photosynthesis, :GPP),
    (:canopy, :hydraulics, :β),
    (:canopy, :hydraulics, :area_index, :leaf),
    (:canopy, :energy, :lhf),
    (:soil, :turbulent_fluxes, :shf),
    (:soil, :turbulent_fluxes, :lhf),
    (:soil, :T),
    (:soil, :θ_l),
    (:soil, :turbulent_fluxes, :vapor_flux),
)

using DataFrames
output_vectors = [getoutput(args...) for args in output_list]
# Extract the last symbol from each tuple for column names
column_names = [names[end] for names in output_list]
# Create a dictionary with simplified column names and corresponding vectors
data_dict = Dict(zip(column_names, output_vectors))
# Create a DataFrame from the dictionary
climalsm = DataFrame(data_dict)
# now if I want for example GPP, I can just do df.GPP

index_t_start = 120 * 48 # we shouldn't hardcode that 120 in ozark_simulation.jl
index_t_end = 120 * 48 + (N_days - N_spinup_days) * 48
model_dt = inputs.DateTime[index_t_start:index_t_end]

insertcols!(climalsm, 1, :DateTime => model_dt)
