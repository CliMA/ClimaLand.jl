"""
    dalec_811_parmin()

Return minimum parameter bounds for the dalec 811 model.
"""
function dalec_811_parmin()
    parmin=[0.00001, 0.2, 0.01, 0.01, 1.001, 0.000025, 0.0001, 0.0001,
    0.0000001, 0.018, 5, 365.25, 0.01, 365.25/12, 365.25, 365.25/12,
     5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 10, 1, 1, 1, 0.01, 0.01, 0.01, 0.01, 1.001, 0.01]
     return parmin
end

"""
    dalec_811_parmax()

Return maximum parameter bounds for the dalec 811 model.
"""
function dalec_811_parmax()
   parmax = [0.01, 0.8, 0.5, 1, 8, 0.001, 0.01, 0.01, 0.001, 0.08, 50, 365.25*4,
    0.5, 100, 365.25 * 4, 150, 200, 2000.0, 2000.0, 2000.0, 100000.0, 2000.0,
     200000.0, 50, 100000, 10000, 10000, 1, 1, 1, 1, 8, 1]
   return parmax
end

"""
    dalec_811_parnames()

Return trainable parameter names for the dalec 811 model.
"""
function dalec_811_parnames()
   parnames = fieldnames(DALEC811Parameters)
   return parnames
end

"""
    check_dalec_811_parameter_bounds(params::Vector{FT}) where {FT <: AbstractFloat}

Validate if a parameter vector is within the parameter bounds of the dalec 811 model.
"""
function check_dalec_811_parameter_bounds(params::Vector{FT}) where {FT <: AbstractFloat}
    parnames = dalec_811_parnames()
    parmax = dalec_811_parmax()
    parmin = dalec_811_parmin()
    
    npar = length(params)
    
    if npar != length(parnames)
        error("The length of the param vector should be 33 for dalec 811 model!")
    end
    
    for i in 1:npar
        if (params[i] < parmin[i] || params[i] > parmax[i])
            error("Parameter no. " * string(i) * ": " * string(parnames[i]) * " is out of bound with min: " * string(parmin[i]) * " and max: " * string(parmax[i]) * ".") 

        end
    end
end

"""
    load_initial_condition!(model::DALECModel, Y::AbstractVector{FT}) where {FT<: AbstractFloat}

Load initia parameters into the Y vector.
"""
function load_initial_condition!(model::DALECModel, Y::AbstractVector{FT}) where {FT<: AbstractFloat}
    @. Y.dalec811.next_labile_pool = model.parameters.Clab
    @. Y.dalec811.next_foliar_pool = model.parameters.Cfol
    @. Y.dalec811.next_root_pool = model.parameters.Croot
    @. Y.dalec811.next_wood_pool = model.parameters.Cwood
    @. Y.dalec811.next_litter_pool = model.parameters.Clitter
    @. Y.dalec811.next_som_pool = model.parameters.Csom
    @. Y.dalec811.next_water_pool = model.parameters.initial_water
    return nothing
end
 