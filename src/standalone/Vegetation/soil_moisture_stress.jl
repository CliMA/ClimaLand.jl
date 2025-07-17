"""

"""

abstract type AbstractSoilMoistureStressModel{FT} <: AbstractCanopyComponent{FT} end

Base.@kwdef struct TuzetMoistureStressParameters{
    FT <: AbstractFloat,
}
    "Sensitivity to low water pressure, in the moisture stress factor, (Pa^{-1})"
    sc::FT
    "Reference water pressure for the moisture stress factor (Pa)"
    pc::FT
end

Base.eltype(::TuzetMoistureStressParameters{FT}) where {FT} = FT

"""
    TuzetMoistureStressModel{FT, TMSP <: TuzetMoistureStressParameters{FT}}
    <: AbstractSoilMoistureStressModel{FT}

"""

struct TuzetMoistureStressModel{FT, TMSP <: TuzetMoistureStressParameters{FT}} <:
       AbstractSoilMoistureStressModel{FT}
    parameters::TMSP
end

function TuzetMoistureStressModel{FT}(
    parameters::TuzetMoistureStressParameters{FT},
) where {FT <: AbstractFloat}
    return TuzetMoistureStressModel{eltype(parameters), typeof(parameters)}(parameters)
end

ClimaLand.name(model::AbstractSoilMoistureStressModel) = :soil_moisture_stress
ClimaLand.auxiliary_vars(model::TuzetMoistureStressModel) = (:βm)
ClimaLand.auxiliary_types(model::TuzetMoistureStressModel{FT}) where {FT} = (FT,)
ClimaLand.auxiliary_domain_names(::TuzetMoistureStressModel) = (:surface,) 

function compute_tuzet_moisture_stress(
    parameters::TuzetMoistureStressParameters{FT},
    p_leaf::FT
) where {FT}
    (; sc, pc) = parameters
    β = min(FT(1), (1 + exp(sc * pc)) / (1 + exp(sc * (pc - p_leaf))))
    return β
end

function update_soil_moisture_stress!(p, Y, model::TuzetMoistureStressModel, canopy)
    # unpack parameters and constants
    (; sc, pc) = model.parameters
    grav = LP.grav(earth_param_set)
    ρ_water = LP.ρ_cloud_liq(earth_param_set)
    n_stem = canopy.hydraulics.n_stem
    n_leaf = canopy.hydraulics.n_leaf
    i_end = n_stem + n_leaf

    ψ = p.canopy.hydraulics.ψ
    p_leaf = @. lazy(ψ.:($$i_end) * ρ_water * grav) # converts height to hydrostatic pressure 

    # Compute the moisture stress factor
    @. p.canopy.soil_moisture_stress.βm = compute_tuzet_moisture_stress(model.parameters, p_leaf)
end

