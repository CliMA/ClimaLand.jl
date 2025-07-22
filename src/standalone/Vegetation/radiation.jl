export BeerLambertParameters,
    BeerLambertModel,
    TwoStreamParameters,
    TwoStreamModel,
    canopy_radiant_energy_fluxes!,
    ConstantGFunction,
    CLMGFunction

abstract type AbstractRadiationModel{FT} <: AbstractCanopyComponent{FT} end

abstract type AbstractGFunction{FT <: AbstractFloat} end

"""
    ConstantGFunction

A type for a constant G function, which is used to represent the leaf angle
distribution function in the radiative transfer models.
"""
struct ConstantGFunction{FT} <: AbstractGFunction{FT}
    "Leaf angle distribution value (unitless)"
    ld::FT
end

# Make the ConstantGFunction broadcastable
Base.broadcastable(G::ConstantGFunction) = tuple(G)

"""
    CLMGFunction

A type for a G function that is parameterized by the cosine of the
solar zenith angle,
following the CLM approach to parameterizing the leaf angle distribution function.
"""
struct CLMGFunction{FT} <: AbstractGFunction{FT}
    "Leaf orientation index (unitless)"
    χl::FT
end

# Make the CLMGFunction broadcastable
Base.broadcastable(G::CLMGFunction) = tuple(G)

"""
    BeerLambertParameters{
        FT <: AbstractFloat,
        G <: Union{AbstractGFunction, ClimaCore.Fields.Field},
        F <: Union{FT, ClimaCore.Fields.Field},
    }

The required parameters for the Beer-Lambert radiative transfer model.
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct BeerLambertParameters{
    FT <: AbstractFloat,
    G <: Union{AbstractGFunction, ClimaCore.Fields.Field},
    F <: Union{FT, ClimaCore.Fields.Field},
}
    "PAR leaf reflectance (unitless)"
    α_PAR_leaf::F
    "NIR leaf reflectance"
    α_NIR_leaf::F
    "Emissivity of the canopy"
    ϵ_canopy::FT
    "Clumping index following Braghiere (2021) (unitless)"
    Ω::FT
    "Typical wavelength per PAR photon (m)"
    λ_γ_PAR::FT
    "Leaf angle distribution function"
    G_Function::G
end

Base.eltype(::BeerLambertParameters{FT}) where {FT} = FT

struct BeerLambertModel{FT, BLP <: BeerLambertParameters{FT}} <:
       AbstractRadiationModel{FT}
    parameters::BLP
end

function BeerLambertModel{FT}(
    parameters::BeerLambertParameters{FT},
) where {FT <: AbstractFloat}
    return BeerLambertModel{eltype(parameters), typeof(parameters)}(parameters)
end

"""
    TwoStreamParameters{FT <: AbstractFloat}

The required parameters for the two-stream radiative transfer model.
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct TwoStreamParameters{
    FT <: AbstractFloat,
    G <: Union{AbstractGFunction, ClimaCore.Fields.Field},
    F <: Union{FT, ClimaCore.Fields.Field},
}
    "PAR leaf reflectance (unitless)"
    α_PAR_leaf::F
    "PAR leaf element transmittance"
    τ_PAR_leaf::F
    "NIR leaf reflectance"
    α_NIR_leaf::F
    "NIR leaf element transmittance"
    τ_NIR_leaf::F
    "Emissivity of the canopy"
    ϵ_canopy::FT
    "Clumping index following Braghiere 2021 (unitless)"
    Ω::F
    "Typical wavelength per PAR photon (m)"
    λ_γ_PAR::FT
    "Number of layers to partition the canopy into when integrating the
    absorption over the canopy vertically. Unrelated to the number of layers in
    the vertical discretization of the canopy for the plant hydraulics model.
    (Constant, and should eventually move to ClimaParams)"
    n_layers::UInt64
    "Leaf angle distribution function"
    G_Function::G
end

Base.eltype(::TwoStreamParameters{FT}) where {FT} = FT

struct TwoStreamModel{FT, TSP <: TwoStreamParameters{FT}} <:
       AbstractRadiationModel{FT}
    parameters::TSP
end

function TwoStreamModel{FT}(
    parameters::TwoStreamParameters{FT},
) where {FT <: AbstractFloat}
    return TwoStreamModel{eltype(parameters), typeof(parameters)}(parameters)
end

"""
    compute_PAR!(par,
        model::AbstractRadiationModel,
        solar_radiation::ClimaLand.PrescribedRadiativeFluxes,
        p,
        t,
    )

Updates `par` with the estimated PAR (W/,m^2) given the input solar radiation
for a radiative transfer model.

The estimated PAR is half of the incident shortwave radiation.
"""
function compute_PAR!(
    par,
    model::AbstractRadiationModel,
    solar_radiation::ClimaLand.AbstractRadiativeDrivers,
    p,
    t,
)
    @. par = p.drivers.SW_d / 2
end

"""
    compute_NIR!(nir,
        model::AbstractRadiationModel,
        solar_radiation::ClimaLand.PrescribedRadiativeFluxes,
        p,
        t,
    )

Update `nir` with the estimated NIR (W/m^2) given the input solar radiation
for a radiative transfer model.

The estimated PNIR is half of the incident shortwave radiation.
"""
function compute_NIR!(
    nir,
    model::AbstractRadiationModel,
    solar_radiation::ClimaLand.AbstractRadiativeDrivers,
    p,
    t,
)
    @. nir = p.drivers.SW_d / 2
end

# Make radiation models broadcastable
Base.broadcastable(RT::AbstractRadiationModel) = tuple(RT)

ClimaLand.name(model::AbstractRadiationModel) = :radiative_transfer
ClimaLand.auxiliary_vars(model::Union{BeerLambertModel, TwoStreamModel}) =
    (:nir_d, :par_d, :nir, :par, :LW_n, :SW_n, :ϵ)
ClimaLand.auxiliary_types(
    model::Union{BeerLambertModel{FT}, TwoStreamModel{FT}},
) where {FT} = (
    FT,
    FT,
    NamedTuple{(:abs, :refl, :trans), Tuple{FT, FT, FT}},
    NamedTuple{(:abs, :refl, :trans), Tuple{FT, FT, FT}},
    FT,
    FT,
    FT,
)
ClimaLand.auxiliary_domain_names(::Union{BeerLambertModel, TwoStreamModel}) =
    (:surface, :surface, :surface, :surface, :surface, :surface, :surface)

"""
    canopy_radiant_energy_fluxes!(p::NamedTuple,
                                  ground::PrescribedGroundConditions
                                  canopy,
                                  radiation::PrescribedRadiativeFluxes,
                                  earth_param_set::PSE,
                                  Y::ClimaCore.Fields.FieldVector,
                                  t,
                                 ) where {PSE}


Computes and stores the net long and short wave radiation, in W/m^2, over all bands,
absorbed by the canopy when the canopy is run in standalone mode, with only
a :canopy model as a prognostic component,
with PrescribedGroundConditions.

LW and SW net radiation are stored in `p.canopy.radiative_transfer.LW_n`
and `p.canopy.radiative_transfer.SW_n`.
"""
function canopy_radiant_energy_fluxes!(
    p::NamedTuple,
    ground::PrescribedGroundConditions,
    canopy,
    radiation::PrescribedRadiativeFluxes,
    earth_param_set::PSE,
    Y::ClimaCore.Fields.FieldVector,
    t,
) where {PSE}
    FT = eltype(earth_param_set)
    par_d = p.canopy.radiative_transfer.par_d
    nir_d = p.canopy.radiative_transfer.nir_d
    f_abs_par = p.canopy.radiative_transfer.par.abs
    f_abs_nir = p.canopy.radiative_transfer.nir.abs
    @. p.canopy.radiative_transfer.SW_n = f_abs_par * par_d + f_abs_nir * nir_d
    ϵ_canopy = p.canopy.radiative_transfer.ϵ # this takes into account LAI/SAI
    # Long wave: use ground conditions from the ground driver
    T_ground::FT = ground.T(t)
    ϵ_ground = ground.ϵ
    _σ = FT(LP.Stefan(earth_param_set))
    LW_d = p.drivers.LW_d
    T_canopy = canopy_temperature(canopy.energy, canopy, Y, p)
    LW_d_canopy = @. (1 - ϵ_canopy) * LW_d + ϵ_canopy * _σ * T_canopy^4
    LW_u_ground = @. ϵ_ground * _σ * T_ground^4 + (1 - ϵ_ground) * LW_d_canopy
    @. p.canopy.radiative_transfer.LW_n =
        ϵ_canopy * LW_d - 2 * ϵ_canopy * _σ * T_canopy^4 +
        ϵ_canopy * LW_u_ground
end

## For interfacing with ClimaParams

"""
    function TwoStreamParameters(FT::AbstractFloat;
        ld = (_) -> 0.5,
        α_PAR_leaf = 0.3,
        τ_PAR_leaf = 0.2,
        α_NIR_leaf = 0.4,
        τ_NIR_leaf = 0.25,
        Ω = 1,
        n_layers = UInt64(20),
        kwargs...
    )
    function TwoStreamParameters(toml_dict;
        ld = (_) -> 0.5,
        α_PAR_leaf = 0.3,
        τ_PAR_leaf = 0.2,
        α_NIR_leaf = 0.4,
        τ_NIR_leaf = 0.25,
        Ω = 1,
        n_layers = UInt64(20),
        kwargs...
    )

Floating-point and toml dict based constructor supplying default values
for the TwoStreamParameters struct. Additional parameter values can be directly set via kwargs.
"""
TwoStreamParameters(::Type{FT}; kwargs...) where {FT <: AbstractFloat} =
    TwoStreamParameters(CP.create_toml_dict(FT); kwargs...)

function TwoStreamParameters(
    toml_dict::CP.AbstractTOMLDict;
    G_Function = ConstantGFunction(CP.float_type(toml_dict)(0.5)),
    α_PAR_leaf::F = 0.3,
    τ_PAR_leaf::F = 0.2,
    α_NIR_leaf::F = 0.4,
    τ_NIR_leaf::F = 0.25,
    Ω = 1,
    n_layers = UInt64(20),
    kwargs...,
) where {F}
    name_map = (;
        :wavelength_per_PAR_photon => :λ_γ_PAR,
        :canopy_emissivity => :ϵ_canopy,
    )

    parameters = CP.get_parameter_values(toml_dict, name_map, "Land")
    FT = CP.float_type(toml_dict)
    # default value for keyword args must be converted manually
    # automatic conversion not possible to Union types
    α_PAR_leaf = FT.(α_PAR_leaf)
    τ_PAR_leaf = FT.(τ_PAR_leaf)
    α_NIR_leaf = FT.(α_NIR_leaf)
    τ_NIR_leaf = FT.(τ_NIR_leaf)
    return TwoStreamParameters{FT, typeof(G_Function), typeof(α_PAR_leaf)}(;
        G_Function,
        α_PAR_leaf,
        τ_PAR_leaf,
        α_NIR_leaf,
        τ_NIR_leaf,
        Ω,
        n_layers,
        parameters...,
        kwargs...,
    )
end

"""
    function BeerLambertParameters(FT::AbstractFloat;
        ld = (_) -> 0.5,
        α_PAR_leaf = 0.1,
        α_NIR_leaf = 0.4,
        Ω = 1,
        kwargs...
    )
    function BeerLambertParameters(toml_dict;
        ld = (_) -> 0.5,
        α_PAR_leaf = 0.1,
        α_NIR_leaf = 0.4,
        Ω = 1,
        kwargs...
    )

Floating-point and toml dict based constructor supplying default values
for the BeerLambertParameters struct. Additional parameter values can be directly set via kwargs.
"""
BeerLambertParameters(::Type{FT}; kwargs...) where {FT <: AbstractFloat} =
    BeerLambertParameters(CP.create_toml_dict(FT); kwargs...)

function BeerLambertParameters(
    toml_dict::CP.AbstractTOMLDict;
    G_Function = ConstantGFunction(CP.float_type(toml_dict)(0.5)),
    α_PAR_leaf::F = 0.1,
    α_NIR_leaf::F = 0.4,
    Ω = 1,
    kwargs...,
) where {F}
    name_map = (;
        :wavelength_per_PAR_photon => :λ_γ_PAR,
        :canopy_emissivity => :ϵ_canopy,
    )

    parameters = CP.get_parameter_values(toml_dict, name_map, "Land")
    FT = CP.float_type(toml_dict)
    # default value for keyword args must be converted manually
    # automatic conversion not possible to Union types
    α_PAR_leaf = FT.(α_PAR_leaf)
    α_NIR_leaf = FT.(α_NIR_leaf)
    return BeerLambertParameters{FT, typeof(G_Function), typeof(α_PAR_leaf)}(;
        G_Function,
        α_PAR_leaf,
        α_NIR_leaf,
        Ω,
        parameters...,
        kwargs...,
    )
end
