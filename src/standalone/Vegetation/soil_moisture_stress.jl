export TuzetMoistureStressParameters, 
    TuzetMoistureStressModel,
    NoMoistureStressModel, 
    PiecewiseMoistureStressParameters, 
    PiecewiseMoistureStressModel, 
    PiecewiseMoistureStressParametersFromHydrology

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
    """Intrinsic moisture stress factor (unitless)"""
    β0::FT
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
    piecewise_moisture_stress_params_from_hydrology(
        hydrology_cm,
        ν,
        θ_r;
        c = 1.0,
        β0 = 1.0,
        ψ_fc = -3.36,  # Field capacity matric potential (m H2O, equivalent to -33 kPa)
        ψ_wp = -152.9  # Wilting point matric potential (m H2O, equivalent to -1500 kPa)
    )

Helper function to create `PiecewiseMoistureStressParameters` by calculating field capacity (θ_c) 
and wilting point (θ_w) from the soil hydraulic model `hydrology_cm`.

# Arguments
- `hydrology_cm`: Soil hydraulic model (e.g., vanGenuchten or BrooksCorey)
- `ν`: Soil porosity (m³/m³)
- `θ_r`: Residual water content (m³/m³)
- `c`: Curvature parameter (unitless), default 1.0
- `β0`: Intrinsic moisture stress factor (unitless), default 1.0
- `ψ_fc`: Field capacity matric potential (m H2O), default -3.36 m (-33 kPa)
- `ψ_wp`: Wilting point matric potential (m H2O), default -152.9 m (-1500 kPa)

The function uses the `inverse_matric_potential` function to calculate the effective 
saturation at field capacity and wilting point, then converts to volumetric water content.
"""
function PiecewiseMoistureStressParametersFromHydrology(
    hydrology_cm,
    ν::FT,
    θ_r::FT;
    c::FT = FT(1.0),
    β0::FT = FT(1.0),
    ψ_fc::FT = FT(-3.36734694),    # -33 kPa, equivalent in m H2O
    ψ_wp::FT = FT(-153.06122449),  # -1500 kPa, equivalent in m H2O
    verbose::Bool = false
) where {FT <: AbstractFloat}
    # Calculate effective saturation at field capacity and wilting point
    S_fc = ClimaLand.Soil.inverse_matric_potential(hydrology_cm, ψ_fc)
    S_wp = ClimaLand.Soil.inverse_matric_potential(hydrology_cm, ψ_wp)

    θ_c = S_fc * (ν - θ_r) + θ_r
    θ_w = S_wp * (ν - θ_r) + θ_r

    verbose && println("Constructing PiecewiseMoistureStressParameters: field capacity θ_c = $θ_c, wilting point θ_w = $θ_w")
    return PiecewiseMoistureStressParameters(
        θ_c = θ_c,
        θ_w = θ_w,
        c = c,
        β0 = β0,
    )
end

"""
compute_piecewise_moisture_stress(
    parameters::PiecewiseMoistureStressParameters{FT},
    θ::FT,
) where {FT}

This function computes the soil moisture stress factor using the volumetric water content θ (m^3/m^3) 
and four parameters: `θ_c` (field capacity, m^3/m^3), `θ_w` (wilting point, m^3/m^3), `c` (curvature 
parameter, unitless), and `β0` (intrinsic moisture stress factor, unitless). This is a modification to
the piecewise formulation defined in Egea et al. (2011) where the entire functional form is multiplied by
`β0` to allow for some "intrinsic moisture stress" even in well-watered conditions.

Citation: https://doi.org/10.1016/j.agrformet.2011.05.019
"""
function compute_piecewise_moisture_stress(
    parameters::PiecewiseMoistureStressParameters{FT},
    θ::FT,
) where {FT}
    (; θ_c, θ_w, c, β0) = parameters

    return β0 * max(
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
        # norm = @. lazy(Canopy.PlantHydraulics.root_distribution_CDF(depth, canopy.hydraulics.parameters.rooting_depth))
        ϑ_l_root_distribution = @. lazy(ϑ_l * Canopy.PlantHydraulics.root_distribution(
            z, canopy.hydraulics.parameters.rooting_depth))
        
        # compute the root zone-averaged volumetric water content
        ClimaCore.Operators.column_integral_definite!(p.canopy.soil_moisture_stress.ϑ_root, ϑ_l_root_distribution)

        # Note: right now, if we do a column integral over root_distribution, we get something that isn't 
        # quite normalized 
        ϑ_root = p.canopy.soil_moisture_stress.ϑ_root
    else
        error("Soil moisture stress for prescribed ground conditions is not implemented yet.")
    end

    @. p.canopy.soil_moisture_stress.βm = compute_piecewise_moisture_stress(model.parameters, ϑ_root)
end
