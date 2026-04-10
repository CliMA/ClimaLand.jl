export AbstractCanopyInterceptionModel,
    NoCanopyInterception,
    CanopyInterceptionModel,
    CanopyInterceptionParameters

"""
    AbstractCanopyInterceptionModel{FT}

An abstract type for canopy interception models. Both `NoCanopyInterception` and
`CanopyInterceptionModel` are subtypes of this abstract type.

Canopy interception models compute the interception, storage, drainage, and
wet-canopy evaporation of precipitation by the canopy.
"""
abstract type AbstractCanopyInterceptionModel{FT} <:
              AbstractCanopyComponent{FT} end

ClimaLand.name(::AbstractCanopyInterceptionModel) = :interception

ClimaLand.auxiliary_vars(::AbstractCanopyInterceptionModel) =
    (:W_int, :f_wet, :throughfall)
ClimaLand.auxiliary_types(::AbstractCanopyInterceptionModel{FT}) where {FT} =
    (FT, FT, FT)
ClimaLand.auxiliary_domain_names(::AbstractCanopyInterceptionModel) =
    (:surface, :surface, :surface)

"""
    NoCanopyInterception{FT} <: AbstractCanopyInterceptionModel{FT}

A no-op interception model: all precipitation passes through the canopy
unmodified. This is the default for backward compatibility.
"""
struct NoCanopyInterception{FT} <: AbstractCanopyInterceptionModel{FT} end

Base.eltype(::NoCanopyInterception{FT}) where {FT} = FT

"""
    update_interception!(p, Y, t, model::NoCanopyInterception, canopy)

No-op update: sets `W_int = 0`, `f_wet = 0`, and `throughfall = P_liq`.
"""
function update_interception!(
    p,
    Y,
    t,
    model::NoCanopyInterception{FT},
    canopy,
) where {FT}
    @. p.canopy.interception.W_int = FT(0)
    @. p.canopy.interception.f_wet = FT(0)
    @. p.canopy.interception.throughfall = p.drivers.P_liq
end

## Active interception model

"""
    CanopyInterceptionParameters{FT <: AbstractFloat}

Parameters for the Rutter et al. (1971, 1975) canopy interception model,
with the Deardorff (1978) wet fraction formulation.

$(DocStringExtensions.FIELDS)
"""
struct CanopyInterceptionParameters{FT <: AbstractFloat}
    "Extinction coefficient for free throughfall (Beer-Lambert analogy), unitless"
    k_int::FT
    "Storage capacity per unit leaf area index [m]"
    S_leaf::FT
    "Storage capacity per unit stem area index [m]"
    S_stem::FT
    "Rutter drainage rate coefficient at canopy saturation [m/s]"
    D_s::FT
    "Rutter drainage exponent [1/m]"
    b::FT
end

Base.eltype(::CanopyInterceptionParameters{FT}) where {FT} = FT

"""
    CanopyInterceptionParameters{FT}(;
        k_int = FT(0.5),
        S_leaf = FT(0.0002),
        S_stem = FT(0.00025),
        D_s = FT(2.4e-5),
        b = FT(3700),
    ) where {FT}

Construct `CanopyInterceptionParameters` with default values from
Rutter et al. (1971, 1975).
"""
function CanopyInterceptionParameters{FT}(;
    k_int = FT(0.5),
    S_leaf = FT(0.0002),
    S_stem = FT(0.00025),
    D_s = FT(2.4e-5),
    b = FT(3700),
) where {FT}
    return CanopyInterceptionParameters{FT}(k_int, S_leaf, S_stem, D_s, b)
end

"""
    CanopyInterceptionParameters(toml_dict::CP.ParamDict)

Construct `CanopyInterceptionParameters` from a TOML parameter dictionary.
"""
function CanopyInterceptionParameters(toml_dict::CP.ParamDict)
    FT = CP.float_type(toml_dict)
    return CanopyInterceptionParameters{FT}(
        toml_dict["interception_k_int"],
        toml_dict["interception_S_leaf"],
        toml_dict["interception_S_stem"],
        toml_dict["interception_D_s"],
        toml_dict["interception_b"],
    )
end

"""
    CanopyInterceptionModel{FT, CIP} <: AbstractCanopyInterceptionModel{FT}

A Rutter-type running water balance interception model (Rutter et al. 1971, 1975)
with Deardorff (1978) wet fraction formulation.

The canopy intercepts a fraction of incoming precipitation proportional to
`1 - exp(-k_int * LAI)` (Beer-Lambert analogy). Intercepted water is stored
on the canopy up to a maximum capacity `S_max = S_leaf * LAI + S_stem * SAI`.
Drainage follows the Rutter exponential: `D = D_s * exp(b * (W_int - S_max))`
when `W_int > 0`. The wet fraction of the canopy `f_wet = min(1, (W_int / S_max)^(2/3))`
evaporates at the potential rate (zero stomatal resistance), while the dry
fraction transpires through stomata as usual.

$(DocStringExtensions.FIELDS)
"""
struct CanopyInterceptionModel{
    FT,
    CIP <: CanopyInterceptionParameters{FT},
} <: AbstractCanopyInterceptionModel{FT}
    "Interception parameters"
    parameters::CIP
end

function CanopyInterceptionModel{FT}(
    parameters::CanopyInterceptionParameters{FT},
) where {FT <: AbstractFloat}
    return CanopyInterceptionModel{eltype(parameters), typeof(parameters)}(
        parameters,
    )
end

ClimaLand.prognostic_vars(::CanopyInterceptionModel) = (:W_int,)
ClimaLand.prognostic_types(::CanopyInterceptionModel{FT}) where {FT} = (FT,)
ClimaLand.prognostic_domain_names(::CanopyInterceptionModel) = (:surface,)

# W_int is prognostic, so only f_wet and throughfall are auxiliary
ClimaLand.auxiliary_vars(::CanopyInterceptionModel) = (:f_wet, :throughfall)
ClimaLand.auxiliary_types(::CanopyInterceptionModel{FT}) where {FT} = (FT, FT)
ClimaLand.auxiliary_domain_names(::CanopyInterceptionModel) =
    (:surface, :surface)

"""
    update_interception!(p, Y, t, model::CanopyInterceptionModel, canopy)

Update the auxiliary variables for canopy interception:
- `f_wet`: wet canopy fraction (Deardorff 1978)
- `throughfall`: free throughfall + drainage (m/s)
"""
function update_interception!(
    p,
    Y,
    t,
    model::CanopyInterceptionModel{FT},
    canopy,
) where {FT}
    params = model.parameters
    LAI = p.canopy.biomass.area_index.leaf
    SAI = p.canopy.biomass.area_index.stem
    W_int = Y.canopy.interception.W_int
    P_liq = p.drivers.P_liq

    # Maximum canopy storage capacity
    # Drainage (Rutter exponential): D = D_s * exp(b * max(W_int - S_max, 0))
    # Throughfall = free throughfall + drainage
    @. p.canopy.interception.f_wet = min(
        FT(1),
        (max(W_int, FT(0)) / max(params.S_leaf * LAI + params.S_stem * SAI, eps(FT)))^FT(2 / 3),
    )
    @. p.canopy.interception.throughfall =
        exp(-params.k_int * LAI) * P_liq + ifelse(
            W_int > FT(0),
            params.D_s *
            exp(params.b * max(W_int - (params.S_leaf * LAI + params.S_stem * SAI), FT(0))),
            FT(0),
        )
end

"""
    total_liq_water_vol_per_area_interception!(surface_field, ::NoCanopyInterception, Y)

No-op: no intercepted water to add.
"""
function total_liq_water_vol_per_area_interception!(
    surface_field,
    ::NoCanopyInterception,
    Y,
)
    return nothing
end

"""
    total_liq_water_vol_per_area_interception!(surface_field, ::CanopyInterceptionModel, Y)

Adds the intercepted water volume to `surface_field`.
"""
function total_liq_water_vol_per_area_interception!(
    surface_field,
    ::CanopyInterceptionModel,
    Y,
)
    @. surface_field += Y.canopy.interception.W_int
    return nothing
end

"""
    ClimaLand.make_compute_exp_tendency(
        model::CanopyInterceptionModel{FT},
        canopy,
    ) where {FT}

Creates the explicit tendency function for canopy interception.

The canopy water balance is:
    dW_int/dt = interception - drainage - wet_canopy_evaporation

where:
- interception = (1 - exp(-k_int * LAI)) * P_liq
- drainage = D_s * exp(b * max(W_int - S_max, 0)) when W_int > 0
- wet_canopy_evaporation = f_wet * vapor_flux (from turbulent fluxes)
"""
function ClimaLand.make_compute_exp_tendency(
    model::CanopyInterceptionModel{FT},
    canopy,
) where {FT}
    function compute_exp_tendency!(dY, Y, p, t)
        params = model.parameters
        LAI = p.canopy.biomass.area_index.leaf
        SAI = p.canopy.biomass.area_index.stem
        P_liq = p.drivers.P_liq
        W_int = Y.canopy.interception.W_int

        # dW_int/dt = interception - drainage - E_c
        @. dY.canopy.interception.W_int =
            (1 - exp(-params.k_int * LAI)) * P_liq - ifelse(
                W_int > FT(0),
                params.D_s *
                exp(params.b * max(W_int - (params.S_leaf * LAI + params.S_stem * SAI), FT(0))),
                FT(0),
            ) - p.canopy.interception.f_wet * p.canopy.turbulent_fluxes.vapor_flux
    end
    return compute_exp_tendency!
end
