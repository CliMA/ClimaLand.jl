export TuzetMoistureStressParameters,
    TuzetMoistureStressModel,
    NoMoistureStressModel,
    PiecewiseMoistureStressParameters,
    PiecewiseMoistureStressModel,
    PiecewiseMoistureStressParametersFromHydrology,
    compute_piecewise_moisture_stress,
    compute_tuzet_moisture_stress,
    update_piecewise_soil_moisture_stress!

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

βm = min(1, (1 + exp(sc * pc)) / (1 + exp(sc * (pc - p_leaf)))), where sc and pc
are stored in the parameter struct.

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

This pointwise function computes the soil moisture stress factor using the leaf water potential (Pa) and two
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

The parameters should fall within
the following ranges:
- θ_high>θ_low should be in (θ_r, ν] where ν is the porosity of the soil
- c should be positive
"""
Base.@kwdef struct PiecewiseMoistureStressParameters{
    FT,
    F <: Union{AbstractFloat, ClimaCore.Fields.Field},
} where {FT}
    """Field capacity volumetric water content or porosity (m^3 / m^3)"""
    θ_high::F
    """Wilting point volumetric water content or residual water fraction(m^3 / m^3)"""
    θ_low::F
    """Curvature parameter (unitless)"""
    c::FT
end

Base.eltype(::PiecewiseMoistureStressParameters{FT}) where {FT} = FT
Base.broadcastable(x::PiecewiseMoistureStressParameters) = tuple(x)

"""
    PiecewiseMoistureStressModel{FT, TMSP <: PiecewiseMoistureStressParameters{FT}}
    <: AbstractSoilMoistureStressModel{FT}

An implementation of a piecewise moisture stress model, taking the form

βm(z) = min(1, [(θ(z)-θ_low)/(θ_high - θ_low)]^c), where θ_high,
θ_low, and c are parameters, and θ(z) is the soil moisture at z.

See  Egea et al. (2011) for details.

Citation: https://doi.org/10.1016/j.agrformet.2011.05.019
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

ClimaLand.auxiliary_vars(model::PiecewiseMoistureStressModel) = (:βm,)
ClimaLand.auxiliary_types(model::PiecewiseMoistureStressModel{FT}) where {FT} =
    (FT,)
ClimaLand.auxiliary_domain_names(::PiecewiseMoistureStressModel) = (:surface,)

"""
    compute_piecewise_moisture_stress(
        θ_high::FT,
        θ_low::FT,
        c::FT,
        θ::FT,
    ) where {FT}

This function computes at a point the soil moisture stress factor using the volumetric water content θ (m^3/m^3)
and four parameters: `θ_high` (field capacity, m^3/m^3), `θ_low` (wilting point, m^3/m^3), `c` (curvature parameter,
unitless). See  Egea et al. (2011).

Citation: https://doi.org/10.1016/j.agrformet.2011.05.019
"""
function compute_piecewise_moisture_stress(
    θ_high::FT,
    θ_low::FT,
    c::FT,
    θ::FT,
) where {FT}
    # avoid taking e.g. sqrt of negative numbers for rational c
    arg = max((θ - θ_low) / (θ_high - θ_low), FT(0))
    return min(FT(1), arg^c)
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
    ground = canopy.boundary_conditions.ground
    update_piecewise_soil_moisture_stress!(ground, p, Y, model, canopy)
end

"""
    update_piecewise_soil_moisture_stress!(ground::PrescribedGroundConditions, p, Y, model, canopy)

Updates the soil moisture stress using the piecewise model for Prescribed
GroundConditions (p.drivers.θ prescribed).
"""
function update_piecewise_soil_moisture_stress!(
    ground::PrescribedGroundConditions,
    p,
    Y,
    model,
    canopy,
)
    @error(
        "You cannot use the PiecewiseSoilMoistureStress model with a prescribed soil yet."
    )
    #(; θ_high, θ_low, c,) = model.parameters
    # Interpret p.drivers.θ as the root zone value.
    #@. p.canopy.soil_moisture_stress.βm =
    #    compute_piecewise_moisture_stress(θ_high, θ_low, c, p.drivers.θ)
end
