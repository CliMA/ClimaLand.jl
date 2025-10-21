export OptimalLAIParameters,
    OptimalLAIModel,
    AbstractLAIModel,
    update_optimal_LAI!,
    make_OptimalLAI_callback

"""
    AbstractLAIModel{FT} <: AbstractCanopyComponent{FT}

An abstract type for LAI models that predict leaf area index dynamically.
"""
abstract type AbstractLAIModel{FT} <: AbstractCanopyComponent{FT} end

ClimaLand.name(model::AbstractLAIModel) = :lai_model

"""
    OptimalLAIParameters{FT<:AbstractFloat}

The required parameters for the optimal LAI model based on Zhou et al. (2025).
This model predicts LAI by optimizing the balance between carbon assimilation
and the cost of maintaining leaf area.

Reference:
Zhou, B., Cai, W., Zhu, Z., Wang, H., Harrison, S. P., & Prentice, C. (2025).
A General Model for the Seasonal to Decadal Dynamics of Leaf Area.
Global Change Biology, 31(1), e70125. https://doi.org/10.1111/gcb.70125

$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct OptimalLAIParameters{FT <: AbstractFloat}
    """Maintenance cost parameter (unitless). Represents the marginal cost of
    maintaining additional leaf area, typically ~ 0.3"""
    m::FT
    """Minimum assimilation threshold (μmol CO2 m^-2 s^-1). Below this threshold,
    LAI is limited to prevent unrealistic growth in low-light conditions, typically ~ 12.227"""
    z::FT
    """Timescale parameter used in EMA for acclimation of optimal LAI (unitless).
    Setting this to 0 represents no incorporation of past values. Since we update
    the EMA equation once per day, α = 1 - 1 day/τ where τ is the acclimation
    timescale in days. Typically ~ 0.933 (corresponding to ~15 day timescale)"""
    α::FT
end

Base.eltype(::OptimalLAIParameters{FT}) where {FT} = FT
Base.broadcastable(x::OptimalLAIParameters) = tuple(x)

"""
    OptimalLAIModel{FT, OLAP <: OptimalLAIParameters{FT}} <: AbstractLAIModel{FT}

An implementation of the optimal LAI model from Zhou et al. (2025).
This model dynamically predicts LAI based on optimizing the balance between
photosynthetic carbon assimilation and the marginal cost of maintaining leaf area.

The model updates LAI at local noon using an exponential moving average (EMA) to
account for acclimation.

Reference:
Zhou, B., Cai, W., Zhu, Z., Wang, H., Harrison, S. P., & Prentice, C. (2025).
A General Model for the Seasonal to Decadal Dynamics of Leaf Area.
Global Change Biology, 31(1), e70125. https://doi.org/10.1111/gcb.70125
"""
struct OptimalLAIModel{FT, OLAP <: OptimalLAIParameters{FT}} <:
       AbstractLAIModel{FT}
    parameters::OLAP
end

Base.eltype(::OptimalLAIModel{FT, OLAP}) where {FT, OLAP} = FT
Base.broadcastable(x::OptimalLAIModel) = tuple(x)

"""
    OptimalLAIModel{FT}(
        parameters::OptimalLAIParameters{FT}
    ) where {FT <: AbstractFloat}

Outer constructor for OptimalLAIModel.
"""
function OptimalLAIModel{FT}(
    parameters::OptimalLAIParameters{FT},
) where {FT <: AbstractFloat}
    return OptimalLAIModel{eltype(parameters), typeof(parameters)}(parameters)
end

"""
    OptimalLAIParameters(
        toml_dict::CP.ParamDict;
        m = toml_dict["optimal_lai_m"],
        z = toml_dict["optimal_lai_z"],
        α = toml_dict["optimal_lai_alpha"],
    )

Creates an `OptimalLAIParameters` object from a TOML dictionary with default values.
"""
function OptimalLAIParameters(
    toml_dict::CP.ParamDict;
    m = toml_dict["optimal_lai_m"],
    z = toml_dict["optimal_lai_z"],
    α = toml_dict["optimal_lai_alpha"],
)
    FT = CP.float_type(toml_dict)

    return OptimalLAIParameters{FT}(
        m = FT(m),
        z = FT(z),
        α = FT(α),
    )
end

"""
    ClimaLand.auxiliary_vars(model::OptimalLAIModel)

Defines the auxiliary variables for the OptimalLAIModel:
- `L`: The optimal LAI (m²/m²) updated at local noon using EMA
"""
ClimaLand.auxiliary_vars(model::OptimalLAIModel) = (:L,)

"""
    ClimaLand.auxiliary_types(model::OptimalLAIModel{FT}) where {FT}

Defines the types of auxiliary variables for OptimalLAIModel.
"""
ClimaLand.auxiliary_types(model::OptimalLAIModel{FT}) where {FT} = (FT,)

"""
    ClimaLand.auxiliary_domain_names(::OptimalLAIModel)

Defines the domain names for auxiliary variables (surface domain).
"""
ClimaLand.auxiliary_domain_names(::OptimalLAIModel) = (:surface,)

# Note: The core LAI computation functions (compute_Ao, compute_L_max, g, dgdL, 
# compute_L_opt, compute_L) are defined in pmodel.jl and imported here for use.
# This avoids code duplication while maintaining a clean module structure.

"""
    update_optimal_LAI!(
        L,
        A,
        cosθs,
        canopy,
        local_noon_mask,
        parameters::OptimalLAIParameters,
    )

Updates the optimal LAI in place at local noon following Zhou et al. (2025).

This is the main update function that:
1. Computes the extinction coefficient from the radiative transfer model
2. Extracts parameters (m, z, α) from the parameters struct
3. Inverts canopy assimilation to get top-of-canopy assimilation (Ao)
4. Computes and applies the optimal LAI update

Args:
- `L`: Current LAI field (m²/m²) - updated in place
- `A`: Canopy-integrated net assimilation rate field (mol CO2 m^-2 s^-1)
- `cosθs`: Cosine of solar zenith angle field (unitless)
- `canopy`: Parent canopy model
- `local_noon_mask`: Field indicating local noon (0 or 1 at each point)
- `parameters`: OptimalLAIParameters containing m, z, and α
"""
function update_optimal_LAI!(L, A, cosθs, canopy, local_noon_mask, parameters)
    radiation = canopy.radiative_transfer
    k = @. lazy(
        radiation.parameters.Ω *
        extinction_coeff(radiation.parameters.G_Function, cosθs),
    )
    FT = eltype(cosθs)
    (; m, z, α) = parameters
    Ao = @. lazy(compute_Ao(A, k, L)) # with current L, compute Ao from A
    # Now update L
    @. L = compute_L(L, Ao, m, k, α, z, local_noon_mask)
end

"""
    set_initial_LAI!(p, Y0, model::OptimalLAIModel, canopy)

Sets the initial LAI to a physically reasonable value. This is called before
the simulation starts to initialize the LAI field.

By default, initializes LAI to a moderate value (3.0 m²/m²). In the future,
this could be made more sophisticated by computing an initial optimal LAI based
on initial conditions, following Zhou et al. (2025).

Args:
- `p`: Auxiliary state cache
- `Y0`: Initial prognostic state
- `model`: OptimalLAIModel instance
- `canopy`: Parent canopy model
"""
function set_initial_LAI!(p, Y0, model::OptimalLAIModel, canopy)
    FT = eltype(model.parameters)
    # Initialize with a reasonable LAI value
    # Could be made more sophisticated in the future
    initial_LAI = FT(3.0)
    p.canopy.lai_model.L .= initial_LAI
end

"""
    call_update_optimal_LAI(
        p,
        Y,
        t;
        canopy,
        dt,
        local_noon,
        model::OptimalLAIModel,
    )

Wrapper function that computes the local noon mask and calls update_optimal_LAI!.
This is the function that gets called by the callback.

Args:
- `p`: Auxiliary state cache
- `Y`: Prognostic state
- `t`: Current simulation time (seconds UTC)
- `canopy`: Parent canopy model
- `dt`: Timestep for local noon window (seconds)
- `local_noon`: Time of local noon in seconds UTC (field)
- `model`: OptimalLAIModel instance
"""
function call_update_optimal_LAI(
    p,
    Y,
    t;
    canopy,
    dt,
    local_noon,
    model::OptimalLAIModel,
)
    # Get local noon mask (1 if within window, 0 otherwise)
    local_noon_mask = @. lazy(get_local_noon_mask(t, dt, local_noon))

    # Get current LAI
    L = p.canopy.lai_model.L

    # Get canopy-integrated net assimilation
    A = p.canopy.photosynthesis.An  # mol CO2 m^-2 s^-1

    # Get cosine of solar zenith angle
    cosθs = p.drivers.cosθs

    # Update LAI using the optimal LAI model
    update_optimal_LAI!(L, A, cosθs, canopy, local_noon_mask, model.parameters)
end

"""
    make_OptimalLAI_callback(
        ::Type{FT},
        t0::ITime,
        dt,
        canopy,
        model::OptimalLAIModel,
        longitude = nothing,
    ) where {FT <: AbstractFloat}

Constructs an IntervalBasedCallback for the optimal LAI model that updates LAI
at local noon using an exponential moving average, following Zhou et al. (2025).

The callback checks for local noon every dt using the provided longitude. The time
of local noon is expressed in seconds UTC and neglects the effects of obliquity
and eccentricity (maximum error ~20 minutes).

Args:
- `FT`: Floating-point type (e.g., Float32, Float64)
- `t0`: Initial time (ITime with epoch in UTC)
- `dt`: Timestep (seconds)
- `canopy`: Parent canopy model
- `model`: OptimalLAIModel instance
- `longitude`: Optional longitude in degrees (default: inferred from domain)

Returns:
- IntervalBasedCallback that updates LAI at local noon
"""
function make_OptimalLAI_callback(
    ::Type{FT},
    t0::ITime,
    dt,
    canopy,
    model::OptimalLAIModel,
    longitude = nothing,
) where {FT <: AbstractFloat}
    # Helper function to convert date to seconds after midnight
    function seconds_after_midnight(date)
        return FT(
            Hour(date).value * 3600 +
            Minute(date).value * 60 +
            Second(date).value,
        )
    end

    # Get longitude from domain if not provided
    if isnothing(longitude)
        try
            longitude = get_long(canopy.domain.space.surface)
        catch e
            error(
                "Longitude must be provided explicitly if the domain does not " *
                "have axes that specify longitude: $e",
            )
        end
    end

    # Compute local noon time in seconds UTC (constant throughout year)
    seconds_in_a_day = IP.day(IP.InsolationParameters(FT))
    start_t = seconds_after_midnight(date(t0))
    local_noon = @. seconds_in_a_day * (FT(1 / 2) - longitude / 360)

    # Define the callback function
    affect! = (integrator) -> call_update_optimal_LAI(
        integrator.p,
        integrator.u,
        (float(integrator.t) + start_t) % seconds_in_a_day,
        canopy = canopy,
        dt = dt,
        local_noon = local_noon,
        model = model,
    )

    return IntervalBasedCallback(
        dt,      # period of this callback
        t0,      # simulation start
        dt,      # integration timestep
        affect!,
    )
end

"""
    get_model_callbacks(
        component::OptimalLAIModel{FT},
        canopy;
        t0,
        Δt,
    ) where {FT}

Creates the callback tuple for the OptimalLAIModel. This is called by the parent
canopy model to set up callbacks for all components.

Args:
- `component`: OptimalLAIModel instance
- `canopy`: Parent canopy model
- `t0`: Initial time
- `Δt`: Timestep

Returns:
- Tuple containing the OptimalLAI callback
"""
function get_model_callbacks(
    component::OptimalLAIModel{FT},
    canopy;
    t0,
    Δt,
) where {FT}
    optimal_lai_cb = make_OptimalLAI_callback(FT, t0, Δt, canopy, component)
    return (optimal_lai_cb,)
end
