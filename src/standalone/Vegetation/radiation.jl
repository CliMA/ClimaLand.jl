export BeerLambertParameters,
    BeerLambertModel,
    TwoStreamParameters,
    TwoStreamModel,
    canopy_radiant_energy_fluxes!,
    ConstantGFunction,
    CLMGFunction,
    ground_albedo_PAR,
    ground_albedo_NIR,
    canopy_sw_rt_beer_lambert,
    canopy_sw_rt_two_stream,
    extinction_coeff,
    compute_G

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
    FF <: Union{FT, ClimaCore.Fields.Field},
}
    "PAR leaf reflectance (unitless)"
    α_PAR_leaf::F
    "NIR leaf reflectance"
    α_NIR_leaf::F
    "Emissivity of the canopy"
    ϵ_canopy::FT
    "Clumping index following Braghiere (2021) (unitless)"
    Ω::FF
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
    T_ground = p.drivers.T_ground
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


"""
    ground_albedo_PAR(prognostic_land_components::Val{(:canopy,)}, ground::PrescribedGroundConditions, _...)

Returns the ground albedo in the PAR for a `PrescribedGroundConditions` driver. In this case,
the `prognostic_land_components` only contain `:canopy`, because the canopy is being run in standalone
mode.
"""
function ground_albedo_PAR(
    prognostic_land_components::Val{(:canopy,)},
    ground::PrescribedGroundConditions,
    _...,
)
    return ground.α_PAR
end

"""
    ground_albedo_NIR(prognostic_land_components::Val{(:canopy,)}, ground::PrescribedGroundConditions, _...)

Returns the ground albedo in the NIR for a `PrescribedGroundConditions` driver. In this case,
the `prognostic_land_components` only contain `:canopy`, because the canopy is being run in standalone
mode.
"""
function ground_albedo_NIR(
    prognostic_land_components::Val{(:canopy,)},
    ground::PrescribedGroundConditions,
    _...,
)
    return ground.α_NIR
end


## For interfacing with ClimaParams

"""
    function TwoStreamParameters(
        toml_dict::CP.ParamDict;
        G_Function,
        α_PAR_leaf,
        τ_PAR_leaf,
        α_NIR_leaf,
        τ_NIR_leaf,
        Ω,
        n_layers = UInt64(20),
        ϵ_canopy = toml_dict["canopy_emissivity"],
    )

TOML dict based constructor supplying default values for the
`TwoStreamParameters` struct.
"""
function TwoStreamParameters(
    toml_dict::CP.ParamDict;
    G_Function,
    α_PAR_leaf,
    τ_PAR_leaf,
    α_NIR_leaf,
    τ_NIR_leaf,
    Ω,
    n_layers = UInt64(20),
    ϵ_canopy = toml_dict["canopy_emissivity"],
)
    FT = CP.float_type(toml_dict)
    λ_γ_PAR = toml_dict["wavelength_per_PAR_photon"]
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
        ϵ_canopy,
        λ_γ_PAR,
    )
end

"""
    function BeerLambertParameters(
        toml_dict::CP.ParamDict;
        G_Function,
        α_PAR_leaf,
        α_NIR_leaf,
        Ω,
        ϵ_canopy = toml_dict["canopy_emissivity"],
    )

TOML dict based constructor supplying default values for the
`BeerLambertParameters` struct. Additional parameter values can be directly set
via kwargs.
"""
function BeerLambertParameters(
    toml_dict::CP.ParamDict;
    G_Function,
    α_PAR_leaf,
    α_NIR_leaf,
    Ω,
    ϵ_canopy = toml_dict["canopy_emissivity"],
)
    FT = CP.float_type(toml_dict)
    # default value for keyword args must be converted manually
    # automatic conversion not possible to Union types
    α_PAR_leaf = FT.(α_PAR_leaf)
    α_NIR_leaf = FT.(α_NIR_leaf)
    Ω = FT.(Ω)
    λ_γ_PAR = toml_dict["wavelength_per_PAR_photon"]
    return BeerLambertParameters{
        FT,
        typeof(G_Function),
        typeof(α_PAR_leaf),
        typeof(Ω),
    }(;
        G_Function,
        α_PAR_leaf,
        α_NIR_leaf,
        Ω,
        ϵ_canopy,
        λ_γ_PAR,
    )
end

# 1. Radiative transfer

"""
    compute_G(
        G::ConstantGFunction,
        _,
    )

Returns the constant leaf angle distribution value for the given G function.
Takes in an arbitrary value for the cosine of the solar zenith angle, which is not used.
"""
function compute_G(G::ConstantGFunction, _)
    return G.ld
end

"""
    compute_G(
        G::CLMGFunction,
        cosθs,
    )

Returns the leaf angle distribution value for CLM G function as a function of the
cosine of the solar zenith angle and the leaf orientation index.
See section 3.1 of https://www2.cesm.ucar.edu/models/cesm2/land/CLM50_Tech_Note.pdf.

Note that the zenith angle is defined ∈ [0,2π), so to prevent a negative value
of G when the sun is below the horizon, we clip cosθs >= 0.
"""
function compute_G(G::CLMGFunction, cosθs::FT) where {FT}
    χl = G.χl
    ϕ1 = 0.5 - 0.633 * χl - 0.33 * χl^2
    ϕ2 = 0.877 * (1 - 2 * ϕ1)
    return FT(ϕ1 + ϕ2 * max(cosθs, 0))
end

"""
    compute_fractional_absorbances!(
        p,
        RT::BeerLambertModel{FT},
        LAI,
        α_soil_PAR,
        α_soil_NIR,
    )

Computes the PAR and NIR fractional absorbances, reflectances, and tranmittances
for a canopy in the case of the
Beer-Lambert model. The absorbances are a function of the radiative transfer
model, as well as the leaf area index, the clumping index,
the cosine of the zenith angle, the leaf angle distribution,
the extinction coefficient, and the
soil albedo in the PAR and NIR bands. Returns a
NamedTuple of NamedTuple, of the form:
(; par = (; refl = , trans = , abs = ),  nir = (; refl = , trans = , abs = ))
"""
function compute_fractional_absorbances!(
    p,
    RT::BeerLambertModel{FT},
    LAI,
    α_soil_PAR,
    α_soil_NIR,
) where {FT}
    RTP = RT.parameters
    cosθs = p.drivers.cosθs
    @. p.canopy.radiative_transfer.par = canopy_sw_rt_beer_lambert(
        RTP.G_Function,
        cosθs,
        RTP.Ω,
        RTP.α_PAR_leaf,
        LAI,
        α_soil_PAR,
    )
    @. p.canopy.radiative_transfer.nir = canopy_sw_rt_beer_lambert(
        RTP.G_Function,
        cosθs,
        RTP.Ω,
        RTP.α_NIR_leaf,
        LAI,
        α_soil_NIR,
    )
end

"""
    compute_fractional_absorbances!(p,
        RT::TwoStreamModel{FT},
        LAI,
        α_soil_PAR,
        α_soil_NIR,
    )

Computes the PAR and NIR fractional absorbances, reflectances, and tranmittances
for a canopy in the case of the
Two-stream model. The absorbances are a function of the radiative transfer
model, as well as the leaf area index, the clumping index,
the cosine of the zenith angle, the leaf angle distribution,
the extinction coefficient, and the
soil albedo in the PAR and NIR bands.

This model also depends on the diffuse fraction.
Returns a
NamedTuple of NamedTuple, of the form:
(; par = (; refl = , trans = , abs = ),  nir = (; refl = , trans = , abs = ))
"""
function compute_fractional_absorbances!(
    p,
    RT::TwoStreamModel{FT},
    LAI,
    α_soil_PAR,
    α_soil_NIR,
) where {FT}
    RTP = RT.parameters
    cosθs = p.drivers.cosθs
    frac_diff = p.drivers.frac_diff
    @. p.canopy.radiative_transfer.par = canopy_sw_rt_two_stream(
        RTP.G_Function,
        RTP.Ω,
        RTP.n_layers,
        RTP.α_PAR_leaf,
        RTP.τ_PAR_leaf,
        LAI,
        cosθs,
        α_soil_PAR,
        frac_diff,
    )
    @. p.canopy.radiative_transfer.nir = canopy_sw_rt_two_stream(
        RTP.G_Function,
        RTP.Ω,
        RTP.n_layers,
        RTP.α_NIR_leaf,
        RTP.τ_NIR_leaf,
        LAI,
        cosθs,
        α_soil_NIR,
        frac_diff,
    )
end

"""
    canopy_sw_rt_beer_lambert(
        G_Function,
        cosθs::FT,
        Ω::FT,
        α_leaf::FT,
        LAI::FT,
        α_soil::FT,
    )

Computes the absorbed, reflected, and transmitted flux fractions by radiation band.

This applies the Beer-Lambert law, which is a function of leaf reflectance
(`α_leaf`), the leaf angle distribution and zenith angle (defined via `G_Function`, and `cosθs`), leaf area index (`LAI`),
and the albedo of the soil (`α_soil`).

Returns a tuple of reflected, absorbed, and transmitted radiation fractions.
"""
function canopy_sw_rt_beer_lambert(
    G_Function,
    cosθs::FT,
    Ω::FT,
    α_leaf::FT,
    LAI::FT,
    α_soil::FT,
) where {FT}
    K = extinction_coeff(G_Function, cosθs)
    AR = (1 - α_leaf) * (1 - exp(-K * LAI * Ω)) * (1 - α_soil)
    TR = exp(-K * LAI * Ω)
    RR = FT(1) - AR - TR * (1 - α_soil)
    return (; abs = AR, refl = RR, trans = TR)
end

"""
    canopy_sw_rt_two_stream(
        G_Function,
        Ω::FT,
        n_layers::UInt64,
        SW_d::FT,
        α_leaf::FT,
        τ_leaf::FT,
        LAI::FT,
        cosθs::FT,
        α_soil::FT,
        frac_diff::FT,
    )

Computes the absorbed, reflected, and transmitted flux fractions by radiation band.

This applies the two-stream radiative transfer solution which takes into account
the impacts of scattering within the canopy. The function takes in all
parameters from the parameter struct of a TwoStreamModel, along with the
incident radiation, LAI, extinction coefficient K, soil albedo from the
canopy soil_driver, the cosine of the solar zenith angle, and τ.

Returns a tuple of reflected, absorbed, and transmitted radiation fractions.
"""
function canopy_sw_rt_two_stream(
    G_Function,
    Ω::FT,
    n_layers::UInt64,
    α_leaf::FT,
    τ_leaf::FT,
    LAI::FT,
    cosθs::FT,
    α_soil::FT,
    frac_diff::FT,
) where {FT}
    α_soil = max(eps(FT), α_soil) # this prevents division by zero, below.
    cosθs = max(eps(FT), cosθs) # The insolations package returns θs > π/2 (nighttime), but this assumes cosθs >0
    G = compute_G(G_Function, cosθs)
    K = extinction_coeff(G_Function, cosθs)
    # Compute μ̄, the average inverse diffuse optical length per LAI
    μ̄ = 1 / max(2G, eps(FT))

    # Clip this to prevent dividing by zero; the sum must also be < 1. Note that using eps(FT)
    # as the threshold leads to numerical errors, so we use 1e-4
    ω = min(max(α_leaf + τ_leaf, FT(1e-4)), 1 - FT(1e-4))

    # Compute aₛ, the single scattering albedo
    aₛ =
        0.5 * ω * (1 - cosθs * log((abs(cosθs) + 1) / max(abs(cosθs), eps(FT))))

    # Compute β₀, the direct upscattering parameter
    β₀ = (1 / ω) * aₛ * (1 + μ̄ * K) / (μ̄ * K)

    # Compute β, the diffuse upscattering parameter
    diff = α_leaf - τ_leaf
    # With uniform distribution, Dickinson integral becomes following:
    c²θ̄ = pi * G / 4
    β = 0.5 * (ω + diff * c²θ̄) / ω

    # Compute coefficients for two-stream solution
    b = 1 - ω + ω * β
    c = ω * β
    d = ω * β₀ * μ̄ * K
    f = ω * μ̄ * K * (1 - β₀)
    h = √(b^2 - c^2) / μ̄
    σ = (μ̄ * K)^2 + c^2 - b^2

    u₁ = b - c / α_soil
    u₂ = b - c * α_soil
    u₃ = f + c * α_soil

    s₁ = exp(-h * LAI * Ω)
    s₂ = exp(-K * LAI * Ω)

    p₁ = b + μ̄ * h
    p₂ = b - μ̄ * h
    p₃ = b + μ̄ * K
    p₄ = b - μ̄ * K

    d₁ = p₁ * (u₁ - μ̄ * h) / s₁ - p₂ * (u₁ + μ̄ * h) * s₁
    d₂ = (u₂ + μ̄ * h) / s₁ - (u₂ - μ̄ * h) * s₁

    # h coefficients for direct upward flux
    h₁ = -d * p₄ - c * f
    h₂ =
        1 / d₁ * (
            (d - h₁ / σ * p₃) * (u₁ - μ̄ * h) / s₁ -
            p₂ * s₂ * (d - c - h₁ / σ * (u₁ + μ̄ * K))
        )
    h₃ =
        -1 / d₁ * (
            (d - h₁ / σ * p₃) * (u₁ + μ̄ * h) * s₁ -
            p₁ * s₂ * (d - c - h₁ / σ * (u₁ + μ̄ * K))
        )

    # h coefficients for direct downward flux
    h₄ = -f * p₃ - c * d
    h₅ =
        -1 / d₂ *
        (h₄ * (u₂ + μ̄ * h) / (σ * s₁) + (u₃ - h₄ / σ * (u₂ - μ̄ * K)) * s₂)
    h₆ =
        1 / d₂ *
        (h₄ / σ * (u₂ - μ̄ * h) * s₁ + (u₃ - h₄ / σ * (u₂ - μ̄ * K)) * s₂)

    # h coefficients for diffuse upward flux
    h₇ = c * (u₁ - μ̄ * h) / (d₁ * s₁)
    h₈ = -c * s₁ * (u₁ + μ̄ * h) / d₁

    # h coefficients for diffuse downward flux
    h₉ = (u₂ + μ̄ * h) / (d₂ * s₁)
    h₁₀ = -s₁ * (u₂ - μ̄ * h) / d₂

    # Compute the LAI per layer for this canopy
    Lₗ = LAI / n_layers

    # Initialize the fraction absorbed value and layer counter
    F_abs = 0
    i = 0

    # Total light reflected from top of canopy
    F_refl = 0

    # Total light transmitted through the canopy on the downward pass
    F_trans = 0

    # Intialize vars to save computed fluxes from each layer for the next layer
    I_dir_up_prev = 0
    I_dir_dn_prev = 0
    I_dif_up_prev = 0
    I_dif_dn_prev = 0


    # Compute F_abs in each canopy layer
    while i <= n_layers

        # Compute cumulative LAI at this layer
        L = i * Lₗ

        # Compute the direct fluxes into/out of the layer
        I_dir_up =
            h₁ * exp(-K * L * Ω) / σ +
            h₂ * exp(-h * L * Ω) +
            h₃ * exp(h * L * Ω)
        I_dir_dn =
            h₄ * exp(-K * L * Ω) / σ +
            h₅ * exp(-h * L * Ω) +
            h₆ * exp(h * L * Ω)

        # Add collimated radiation to downward flux
        I_dir_dn += exp(-K * L * Ω)

        # Compute the diffuse fluxes into/out of the layer
        I_dif_up = h₇ * exp(-h * L * Ω) + h₈ * exp(h * L * Ω)
        I_dif_dn = h₉ * exp(-h * L * Ω) + h₁₀ * exp(h * L * Ω)

        # Energy balance giving radiation absorbed in the layer
        if i == 0
            I_dir_abs = 0
            I_dif_abs = 0
        else
            I_dir_abs = I_dir_up - I_dir_up_prev - I_dir_dn + I_dir_dn_prev
            I_dif_abs = I_dif_up - I_dif_up_prev - I_dif_dn + I_dif_dn_prev
        end

        if i == 1 # prev = top of layer
            F_refl =
                (1 - frac_diff) * I_dir_up_prev + (frac_diff) * I_dif_up_prev
        end
        if i == n_layers # not prev = bottom of layer
            F_trans = (1 - frac_diff) * I_dir_dn + (frac_diff) * I_dif_dn
        end


        # Add radiation absorbed in the layer to total absorbed radiation
        F_abs += (1 - frac_diff) * I_dir_abs + (frac_diff) * I_dif_abs

        # Save input/output values to compute energy balance of next layer
        I_dir_up_prev = I_dir_up
        I_dir_dn_prev = I_dir_dn
        I_dif_up_prev = I_dif_up
        I_dif_dn_prev = I_dif_dn

        # Move on to the next layer
        i += 1
    end

    # Reflected is reflected from the land surface. F_abs
    # refers to radiation absorbed by the canopy on either pass.
    # F_trans refers to the downward pass only. Note that
    # (1-α_soil)*F_trans + F_abs = total absorbed fraction by land, so the following
    # must hold
    # @assert (1 - α_soil) * FT(F_trans) + FT(F_abs) + FT(F_refl) ≈ 1
    # This is tested in test/standalone/Vegetation/test_two_stream.jl
    return (; abs = FT(F_abs), refl = FT(F_refl), trans = FT(F_trans))
end

"""
    extinction_coeff(G_Function,
                     cosθs::FT) where {FT}

Computes the vegetation extinction coefficient (`K`), as a function
of the cosine of the sun zenith angle (`cosθs`),
and the leaf angle distribution function(`G_Function`).

In the two-stream scheme, values of K ~ 1/epsilon can lead to numerical issues.
Here we clip it to 1e6.
"""
function extinction_coeff(G_Function, cosθs::FT) where {FT}
    G = compute_G(G_Function, cosθs)
    K = min(G / max(cosθs, eps(FT)), FT(1e6))
    return K
end

"""
   update_radiative_transfer!(p, Y, t, radiative_transfer::AbstractRadiationModel, canopy)

Updates the following cache variables in place:
- downwelling PAR and NIR in W/m^2: p.canopy.radiative_transfer.par_d, .nir_d
- absorbed, reflected, and transmitted fractions of PAR and NIR: p.canopy.radiative_transfer.par, .nir
- canopy emissivity: p.canopy.radiative_transfer.ϵ

This implies that all concrete types of AbstractRadiationModel must
set these variables.
"""
function update_radiative_transfer!(
    p,
    Y,
    t,
    radiative_transfer::AbstractRadiationModel,
    canopy,
)
    par_d = p.canopy.radiative_transfer.par_d
    nir_d = p.canopy.radiative_transfer.nir_d
    area_index = p.canopy.hydraulics.area_index
    LAI = area_index.leaf
    SAI = area_index.stem
    bc = canopy.boundary_conditions

    # update radiative transfer
    (; G_Function, Ω, λ_γ_PAR) = radiative_transfer.parameters
    @. p.canopy.radiative_transfer.ϵ =
        radiative_transfer.parameters.ϵ_canopy * (1 - exp(-(LAI + SAI))) #from CLM 5.0, Tech note 4.20
    compute_PAR!(par_d, radiative_transfer, bc.radiation, p, t)
    compute_NIR!(nir_d, radiative_transfer, bc.radiation, p, t)

    compute_fractional_absorbances!(
        p,
        radiative_transfer,
        LAI,
        ground_albedo_PAR(
            Val(bc.prognostic_land_components),
            bc.ground,
            Y,
            p,
            t,
        ),
        ground_albedo_NIR(
            Val(bc.prognostic_land_components),
            bc.ground,
            Y,
            p,
            t,
        ),
    )
end
