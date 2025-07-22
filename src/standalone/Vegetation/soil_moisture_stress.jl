export TuzetMoistureStressParameters, TuzetMoistureStressModel,
    NoMoistureStressModel, 
    PiecewiseMoistureStressParameters, PiecewiseMoistureStressModel

# define an abstract type for all soil moisture stress models
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
Base.broadcastable(x::TuzetMoistureStressParameters) = tuple(x)

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
ClimaLand.auxiliary_vars(model::TuzetMoistureStressModel) = (:βm,)
ClimaLand.auxiliary_types(model::TuzetMoistureStressModel{FT}) where {FT} = (FT,)
ClimaLand.auxiliary_domain_names(::TuzetMoistureStressModel) = (:surface,) 


"""
compute_tuzet_moisture_stress(
    parameters::TuzetMoistureStressParameters{FT},
    p_leaf::FT
) where {FT}

This function computes the soil moisture stress factor using the leaf water potential (Pa) and two
parameters `sc` (sensitivity to low water potential, Pa^-1) and `pc` (reference water potential, Pa).
"""
function compute_tuzet_moisture_stress(
    parameters::TuzetMoistureStressParameters{FT},
    p_leaf::FT
) where {FT}
    (; sc, pc) = parameters
    β = min(FT(1), (1 + exp(sc * pc)) / (1 + exp(sc * (pc - p_leaf))))
    return β
end

function update_soil_moisture_stress!(p, Y, model::TuzetMoistureStressModel, canopy)
    # unpack constants
    earth_param_set = canopy.parameters.earth_param_set
    grav = LP.grav(earth_param_set)
    ρ_water = LP.ρ_cloud_liq(earth_param_set)
    n_stem = canopy.hydraulics.n_stem
    n_leaf = canopy.hydraulics.n_leaf
    i_end = n_stem + n_leaf

    ψ = p.canopy.hydraulics.ψ
    p_leaf = @. lazy(ψ.:($$i_end) * ρ_water * grav) # converts head to hydrostatic pressure 

    # Compute the moisture stress factor
    @. p.canopy.soil_moisture_stress.βm = compute_tuzet_moisture_stress(model.parameters, p_leaf)
end


struct NoMoistureStressModel{FT} <: AbstractSoilMoistureStressModel{FT} end

ClimaLand.auxiliary_vars(model::NoMoistureStressModel) = (:βm,)
ClimaLand.auxiliary_types(model::NoMoistureStressModel{FT}) where {FT} = (FT,)
ClimaLand.auxiliary_domain_names(::NoMoistureStressModel) = (:surface,) 

function update_soil_moisture_stress!(p, Y, model::NoMoistureStressModel, canopy)
    @. p.canopy.soil_moisture_stress.βm = 1.0
end


Base.@kwdef struct PiecewiseMoistureStressParameters{
    FT <: AbstractFloat,
    TC <: Union{FT, ClimaCore.Fields.Field},
    TW <: Union{FT, ClimaCore.Fields.Field},
    CP <: Union{FT, ClimaCore.Fields.Field}
}
    """Field capacity volumetric water content (m^3 / m^3)"""
    θ_c::TC
    """Wilting point volumetric water content (m^3 / m^3)"""
    θ_w::TW
    """Curvature parameter (unitless)"""
    c::CP
end

Base.eltype(::PiecewiseMoistureStressParameters{FT}) where {FT} = FT
Base.broadcastable(x::PiecewiseMoistureStressParameters) = tuple(x) 

"""
    PiecewiseMoistureStressModel{FT, TMSP <: PiecewiseMoistureStressParameters{FT}}
    <: AbstractSoilMoistureStressModel{FT}

"""
struct PiecewiseMoistureStressModel{FT, TMSP <: PiecewiseMoistureStressParameters{FT}} <:
       AbstractSoilMoistureStressModel{FT}
    parameters::TMSP
end

function PiecewiseMoistureStressModel{FT}(
    parameters::PiecewiseMoistureStressParameters{FT},
) where {FT <: AbstractFloat}
    return PiecewiseMoistureStressModel{eltype(parameters), typeof(parameters)}(parameters)
end

ClimaLand.auxiliary_vars(model::PiecewiseMoistureStressModel) = (:βm, :ϑ_root,)
ClimaLand.auxiliary_types(model::PiecewiseMoistureStressModel{FT}) where {FT} = (FT, FT,)
ClimaLand.auxiliary_domain_names(::PiecewiseMoistureStressModel) = (:surface, :surface,)


"""
compute_piecewise_moisture_stress(
    parameters::PiecewiseMoistureStressParameters{FT},
    θ::FT,
) where {FT}

This function computes the soil moisture stress factor using the volumetric water content θ (m^3/m^3) 
and three parameters: `θ_c` (field capacity, m^3/m^3), `θ_w` (wilting point, m^3/m^3), and `c` (curvature 
parameter, unitless). We follow the piecewise formulation defined in Egea et al. (2011). 

Citation: https://doi.org/10.1016/j.agrformet.2011.05.019
"""
function compute_piecewise_moisture_stress(
    parameters::PiecewiseMoistureStressParameters{FT},
    θ::FT,
) where {FT}
    (; θ_c, θ_w, c) = parameters

    return max(
        FT(0),
        min(
            FT(1),
            ((θ - θ_w) / (θ_c - θ_w))^c,
        )
    )
end

function update_soil_moisture_stress!(p, Y, model::PiecewiseMoistureStressModel, canopy)
    if :soil in canopy.boundary_conditions.prognostic_land_components 
        ϑ_l = Y.soil.ϑ_l
        z = ClimaCore.Fields.coordinate_field(axes(Y.soil.ϑ_l)).z

        # normalized distribution for root density 
        root_distribution = @. lazy(Canopy.PlantHydraulics.root_distribution(
            z, canopy.hydraulics.parameters.rooting_depth))
        
        # compute the root zone-averaged volumetric water content
        ClimaCore.Operators.column_integral_definite!(p.canopy.soil_moisture_stress.ϑ_root,
            ϑ_l .* root_distribution)

        # Note: right now, if we do a column integral over root_distribution, we get something that isn't 
        # quite normalized 
        ϑ_root = p.canopy.soil_moisture_stress.ϑ_root
    else
        ϑ_root = canopy.boundary_conditions.ground.ϑ_root
    end

    @. p.canopy.soil_moisture_stress.βm = compute_piecewise_moisture_stress(model.parameters, ϑ_root)
end
