export getoutput, make_output_df

"""
    getoutput(sv, variable::Symbol, variables::Symbol...; result = sv.saveval, depth = 1)

Return a vector of FT corresponding to the variable of interest at all times.
By default, get output from sv.saveval, but user can specify e.g., result = sol.u
By default, get surface value, but user can specify depth for e.g., soil temperature

example 1: 
julia> getoutput(sv, 1, :SW_out)

example 2: 
julia> getoutput(sv, 1, :canopy, :conductance, :gs)

example 3: 
julia> getoutput(sol, 1, :soil, :ϑ_l; result = sol.u)
"""
function getoutput(
    sv, # or sol for prognostic variables
    depth, # 1 is surface
    variable::Symbol,
    variables::Symbol...;
    result = sv.saveval, # or sol.u for prognostic variables
)
    for v in (variable, variables...)
        result = getproperty.(result, v)
    end
    return [parent(r)[depth] for r in result]
end

"""
    make_output_df(sv, inputs)

Return a dataframe containing climaland outputs
"""
function make_output_df(
    site_ID,
    sv,
    inputs;
    setup = make_setup(site_ID),
    timestepper = make_timestepper(setup),
)
    # List of output that we want
    output_list = vcat(
        (1, :SW_out),
        (1, :LW_out),
        (1, :canopy, :conductance, :gs),
        (1, :canopy, :conductance, :transpiration),
        (1, :canopy, :autotrophic_respiration, :Ra),
        (1, :canopy, :photosynthesis, :GPP),
        (1, :canopy, :hydraulics, :β),
        (1, :canopy, :hydraulics, :area_index, :leaf),
        (1, :canopy, :energy, :lhf),
        (1, :soil, :turbulent_fluxes, :shf),
        (1, :soil, :turbulent_fluxes, :lhf),
        collect(map(i -> (i, :soil, :T), 1:20)), # 20 shouldn't be hard-coded, but an arg, equal to n layers
        collect(map(i -> (i, :soil, :θ_l), 1:20)),
        (1, :soil, :turbulent_fluxes, :vapor_flux_liq),
        (1, :canopy, :sif, :SIF),
    )

    output_vectors = [getoutput(sv, args...) for args in output_list]
    # Extract the last symbol from each tuple for column names
    column_names =
        [Symbol(names[end], "_depth_", names[1]) for names in output_list]
    # remove "_1"
    for i in 1:length(column_names)
        name = string(column_names[i])
        if endswith(name, "_1")
            column_names[i] = Symbol(chop(name, tail = 8))
        end
    end
    # Create a dictionary with simplified column names and corresponding vectors
    data_dict = Dict(zip(column_names, output_vectors))
    # Create a DataFrame from the dictionary
    climaland = DataFrame(data_dict)
    # now if I want for example GPP, I can just do df.GPP

    index_t_start = Int(setup.t0 / (3600 * 24) * 48)
    index_t_end = Int(
        index_t_start + (timestepper.N_days - timestepper.N_spinup_days) * 48,
    )
    model_dt = inputs.DateTime[index_t_start:index_t_end]

    insertcols!(climaland, 1, :DateTime => model_dt)
    return climaland
end
