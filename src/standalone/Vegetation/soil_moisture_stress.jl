export TuzetMoistureStressParameters, TuzetMoistureStressModel
export NoMoistureStressModel
export PiecewiseMoistureStressParameters,
    PiecewiseMoistureStressModel, PiecewiseMoistureStressParametersFromHydrology

export compute_piecewise_moisture_stress
export compute_tuzet_moisture_stress

# define an abstract type for all soil moisture stress models
abstract type AbstractSoilMoistureStressModel{FT} <: AbstractCanopyComponent{FT} end


"""
    TuzetMoistureStressParameters{FT <: AbstractFloat}

The required parameters for the Tuzet plant moisture stress model.
"""
Base.@kwdef struct TuzetMoistureStressParameters{FT <: AbstractFloat}
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

An implementation of the Tuzet moisture stress function. 

TUZET, A., PERRIER, A. and LEUNING, R. (2003), A coupled model of stomatal conductance, 
    photosynthesis and transpiration. Plant, Cell & Environment, 26: 1097-1116. 
    https://doi.org/10.1046/j.1365-3040.2003.01035.x
"""
struct TuzetMoistureStressModel{
    FT,
    TMSP <: TuzetMoistureStressParameters{FT},
} <: AbstractSoilMoistureStressModel{FT}
    parameters::TMSP
end

"""
    TuzetMoistureStressModel{FT}(
        parameters::TuzetMoistureStressParameters{FT},
    ) where {FT <: AbstractFloat}

Outer constructor for TuzetMoistureStressModel.
"""
function TuzetMoistureStressModel{FT}(
    parameters::TuzetMoistureStressParameters{FT},
) where {FT <: AbstractFloat}
    return TuzetMoistureStressModel{eltype(parameters), typeof(parameters)}(
        parameters,
    )
end

ClimaLand.name(model::AbstractSoilMoistureStressModel) = :soil_moisture_stress
ClimaLand.auxiliary_vars(model::TuzetMoistureStressModel) = (:βm,)
ClimaLand.auxiliary_types(model::TuzetMoistureStressModel{FT}) where {FT} =
    (FT,)
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
    p_leaf::FT,
) where {FT}
    (; sc, pc) = parameters
    β = min(FT(1), (1 + exp(sc * pc)) / (1 + exp(sc * (pc - p_leaf))))
    return β
end


"""
    update_soil_moisture_stress!(
        p,
        Y,
        model::TuzetMoistureStressModel,
        canopy,
    )    

Computes and updates in place the soil moisture stress for the Tuzet formulation. 
"""
function update_soil_moisture_stress!(
    p,
    Y,
    model::TuzetMoistureStressModel,
    canopy,
)
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
    @. p.canopy.soil_moisture_stress.βm =
        compute_tuzet_moisture_stress(model.parameters, p_leaf)
end

"""
    struct NoMoistureStressModel{FT} <: AbstractSoilMoistureStressModel{FT} end

A constructor for NoMoistureStressModel, which sets βm = 1.0 always. 
"""
struct NoMoistureStressModel{FT} <: AbstractSoilMoistureStressModel{FT} end

Base.eltype(::NoMoistureStressModel{FT}) where {FT} = FT

ClimaLand.auxiliary_vars(model::NoMoistureStressModel) = (:βm,)
ClimaLand.auxiliary_types(model::NoMoistureStressModel{FT}) where {FT} = (FT,)
ClimaLand.auxiliary_domain_names(::NoMoistureStressModel) = (:surface,)

function update_soil_moisture_stress!(
    p,
    Y,
    model::NoMoistureStressModel,
    canopy,
)
    FT = eltype(model)
    @. p.canopy.soil_moisture_stress.βm = FT(1.0)
end


"""
    struct PiecewiseMoistureStressParameters{
        F <: Union{AbstractFloat, ClimaCore.Fields.Field},
    }

The required parameters for the piecewise moisture stress model.
"""
Base.@kwdef struct PiecewiseMoistureStressParameters{
    F <: Union{AbstractFloat, ClimaCore.Fields.Field},
}
    """Field capacity volumetric water content (m^3 / m^3)"""
    θ_c::F
    """Wilting point volumetric water content (m^3 / m^3)"""
    θ_w::F
    """Curvature parameter (unitless)"""
    c::F
    """Intrinsic moisture stress factor (unitless)"""
    β0::F
    """Depth of lowest soil cell face (m), negative"""
    soil_zmin::F
end


"""
    PiecewiseMoistureStressParameters(
        ::Type{FT};
        θ_c, 
        θ_w,
        c, 
        β0, 
        soil_zmin
    ) where {FT <: AbstractFloat}
    
An outer constructor for PiecewiseMoistureStressParameters. The parameters should fall within
the following ranges:
- θ_c, θ_w should be in [0,ν] where ν is the porosity of the soil 
- c should be positive
- β should be in (0,1]
- soil_zmin should be negative
"""
function PiecewiseMoistureStressParameters(
    ::Type{FT};
    θ_c,
    θ_w,
    c,
    β0,
    soil_zmin = FT(-1.0),
) where {FT <: AbstractFloat}
    return PiecewiseMoistureStressParameters{FT}(θ_c, θ_w, c, β0, soil_zmin)
end

Base.eltype(::PiecewiseMoistureStressParameters{FT}) where {FT} = FT
Base.broadcastable(x::PiecewiseMoistureStressParameters) = tuple(x)

"""
    PiecewiseMoistureStressModel{FT, TMSP <: PiecewiseMoistureStressParameters{FT}}
    <: AbstractSoilMoistureStressModel{FT}

"""
struct PiecewiseMoistureStressModel{
    FT,
    TMSP <: PiecewiseMoistureStressParameters{FT},
} <: AbstractSoilMoistureStressModel{FT}
    parameters::TMSP
end

function PiecewiseMoistureStressModel{FT}(
    parameters::PiecewiseMoistureStressParameters{FT},
) where {FT <: AbstractFloat}
    return PiecewiseMoistureStressModel{eltype(parameters), typeof(parameters)}(
        parameters,
    )
end

ClimaLand.auxiliary_vars(model::PiecewiseMoistureStressModel) = (:βm, :ϑ_root)
ClimaLand.auxiliary_types(model::PiecewiseMoistureStressModel{FT}) where {FT} =
    (FT, FT)
ClimaLand.auxiliary_domain_names(::PiecewiseMoistureStressModel) =
    (:surface, :surface)

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
and wilting point (θ_w) from the soil hydraulic closure model `hydrology_cm`.

# Arguments
- `hydrology_cm`: Soil hydraulic closure model (e.g., vanGenuchten or BrooksCorey)
- `ν`: Soil porosity (m^3/m^3)
- `θ_r`: Residual water content (m^3/m^3)
- `soil_zmin`: Depth of the lowest soil cell face (m), must be negative
- `c`: Curvature parameter (unitless), default 1.0
- `β0`: Intrinsic moisture stress factor (unitless), default 1.0
- `ψ_fc`: Field capacity matric potential (m H2O), default -3.36 m (-33 kPa)
- `ψ_wp`: Wilting point matric potential (m H2O), default -152.9 m (-1500 kPa)

The function uses the `inverse_matric_potential` function to calculate the effective 
saturation at field capacity and wilting point, then converts to volumetric water content.
"""
function PiecewiseMoistureStressParametersFromHydrology(
    ::Type{FT},
    hydrology_cm,
    ν::Union{FT, ClimaCore.Fields.Field},
    θ_r::Union{FT, ClimaCore.Fields.Field},
    soil_zmin::Union{FT, ClimaCore.Fields.Field};
    c::FT = FT(1.0),
    β0::FT = FT(1.0),
    ψ_fc::FT = FT(-3.36734694),    # -33 kPa, equivalent in m H2O
    ψ_wp::FT = FT(-153.06122449),  # -1500 kPa, equivalent in m H2O
    verbose::Bool = false,
) where {FT <: AbstractFloat}
    if c <= 0
        throw(
            ArgumentError("Curvature parameter `c` must be greater than zero"),
        )
    end
    if β0 <= 0
        throw(
            ArgumentError(
                "Intrinsic moisture stress factor `β0` must be greater than zero",
            ),
        )
    end

    c = FT(c)
    β0 = FT(β0)
    ψ_fc = FT(ψ_fc)
    ψ_wp = FT(ψ_wp)

    # Calculate effective saturation at field capacity and wilting point
    if hydrology_cm isa ClimaCore.Fields.Field
        S_fc = @. ClimaLand.Soil.inverse_matric_potential(hydrology_cm, ψ_fc)
        S_wp = @. ClimaLand.Soil.inverse_matric_potential(hydrology_cm, ψ_wp)
        θ_c = @. S_fc * (ν - θ_r) + θ_r
        θ_w = @. S_wp * (ν - θ_r) + θ_r
    else
        S_fc = ClimaLand.Soil.inverse_matric_potential(hydrology_cm, ψ_fc)
        S_wp = ClimaLand.Soil.inverse_matric_potential(hydrology_cm, ψ_wp)
        θ_c = S_fc * (ν - θ_r) + θ_r
        θ_w = S_wp * (ν - θ_r) + θ_r
    end

    verbose && println(
        "Constructing PiecewiseMoistureStressParameters: field capacity θ_c = $θ_c, wilting point θ_w = $θ_w",
    )
    return PiecewiseMoistureStressParameters(
        θ_c = θ_c,
        θ_w = θ_w,
        c = c,
        β0 = β0,
        soil_zmin = soil_zmin,
    )
end

"""
compute_piecewise_moisture_stress(
    θ_c::FT, 
    θ_w::FT,
    c::FT,
    β0::FT,
    θ::FT,
) where {FT}

This function computes at a point the soil moisture stress factor using the volumetric water content θ (m^3/m^3) 
and four parameters: `θ_c` (field capacity, m^3/m^3), `θ_w` (wilting point, m^3/m^3), `c` (curvature parameter, 
unitless), and `β0` (intrinsic moisture stress factor, unitless). This is a modification to
the piecewise formulation defined in Egea et al. (2011) where the entire functional form is multiplied by
`β0` to allow for some "intrinsic moisture stress" even in well-watered conditions.

Citation: https://doi.org/10.1016/j.agrformet.2011.05.019
"""
function compute_piecewise_moisture_stress(
    θ_c::FT,
    θ_w::FT,
    c::FT,
    β0::FT,
    θ::FT,
) where {FT}
    # avoid taking e.g. sqrt of negative numbers for rational c
    arg = max((θ - θ_w) / (θ_c - θ_w), FT(0))

    return β0 * max(FT(0), min(FT(1), arg^c))
end

"""
    update_soil_moisture_stress!(
        p,
        Y,
        model::PiecewiseMoistureStressModel,
        canopy,
    )

This updates the soil moisture stress factor according to the piecewise soil moisture stress model. 
"""
function update_soil_moisture_stress!(
    p,
    Y,
    model::PiecewiseMoistureStressModel,
    canopy,
)
    # compute root-integrated soil water content
    if :soil in canopy.boundary_conditions.prognostic_land_components
        # unpack soil moisture field, parameters
        ϑ_l = Y.soil.ϑ_l
        (; θ_c, θ_w, c, β0, soil_zmin) = model.parameters
        z = ClimaCore.Fields.coordinate_field(axes(ϑ_l)).z

        # normalized distribution for root density 
        norm = @. lazy(
            Canopy.PlantHydraulics.root_distribution_CDF(
                soil_zmin,
                canopy.hydraulics.parameters.rooting_depth,
            ),
        )

        # per soil element
        βm = @. lazy(compute_piecewise_moisture_stress(θ_c, θ_w, c, β0, ϑ_l))
        βm_root_distribution = @. lazy(
            βm * Canopy.PlantHydraulics.root_distribution(
                z,
                canopy.hydraulics.parameters.rooting_depth,
            ) / norm,
        )

        # compute the root zone-averaged βm
        ClimaCore.Operators.column_integral_definite!(
            p.canopy.soil_moisture_stress.βm,
            βm_root_distribution,
        )
    else
        (; θ_c, θ_w, c, β0, soil_zmin) = model.parameters
        @. p.canopy.soil_moisture_stress.βm =
            compute_piecewise_moisture_stress(θ_c, θ_w, c, β0, p.drivers.θ)
    end
end
