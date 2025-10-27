export AbstractLAIModel,
    OptimalLAIModel,
    OptimalLAIParameters,
    update_LAI,
    initialize_LAI!

"""
    AbstractLAIModel{FT <: AbstractFloat}

An abstract type for LAI (Leaf Area Index) model parameterizations.
"""
abstract type AbstractLAIModel{FT <: AbstractFloat} <: AbstractCanopyComponent{FT} end

"""
    OptimalLAIParameters{FT<:AbstractFloat}

The required parameters for the optimal LAI model based on Zhou et al. (2025).

# References
Zhou et al. (2025) "A General Model for the Seasonal to Decadal Dynamics of Leaf Area"
Global Change Biology. https://onlinelibrary.wiley.com/doi/pdf/10.1111/gcb.70125

$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct OptimalLAIParameters{FT <: AbstractFloat}
    """Light extinction coefficient (dimensionless), typically 0.5"""
    k::FT
    """Unit cost of constructing and maintaining leaves (mol m⁻² yr⁻¹), globally fitted as 12.227 mol m⁻² yr⁻¹"""
    z::FT
    """Ratio of leaf-internal to ambient CO₂ partial pressure (dimensionless), typically from the P-model"""
    chi::FT
    """Fraction of annual precipitation used by plants (dimensionless), varies with aridity index"""
    f0::FT
    """Dimensionless parameter representing departure from square-wave LAI dynamics, globally fitted as 0.771"""
    sigma::FT
    """Smoothing factor for exponential moving average (dimensionless, 0-1). Set to 0.067 for ~15 days of memory"""
    alpha::FT
end

Base.eltype(::OptimalLAIParameters{FT}) where {FT} = FT

# make these custom structs broadcastable as tuples
Base.broadcastable(x::OptimalLAIParameters) = tuple(x)

"""
    OptimalLAIModel{FT, OLPT <: OptimalLAIParameters{FT}} <: AbstractLAIModel{FT}

An implementation of the optimal LAI model from Zhou et al. (2025).

This model predicts seasonal to decadal dynamics of leaf area index based on
optimality principles, balancing energy and water constraints.

# References
Zhou et al. (2025) "A General Model for the Seasonal to Decadal Dynamics of Leaf Area"
Global Change Biology. https://onlinelibrary.wiley.com/doi/pdf/10.1111/gcb.70125
"""
struct OptimalLAIModel{FT, OLPT <: OptimalLAIParameters{FT}} <:
       AbstractLAIModel{FT}
    "Required parameters for the optimal LAI model"
    parameters::OLPT
end

Base.eltype(::OptimalLAIModel{FT, OLPT}) where {FT, OLPT} = FT

"""
    OptimalLAIModel{FT}(parameters::OptimalLAIParameters{FT})

Outer constructor for the OptimalLAIModel struct.
"""
function OptimalLAIModel{FT}(
    parameters::OptimalLAIParameters{FT},
) where {FT <: AbstractFloat}
    return OptimalLAIModel{FT, typeof(parameters)}(parameters)
end

"""
    ClimaLand.auxiliary_vars(model::OptimalLAIModel)
    ClimaLand.auxiliary_types(model::OptimalLAIModel)
    ClimaLand.auxiliary_domain_names(model::OptimalLAIModel)

Defines the auxiliary variable for the OptimalLAIModel: the leaf area index (LAI).
"""
ClimaLand.auxiliary_vars(model::OptimalLAIModel) = (:LAI,)
ClimaLand.auxiliary_types(model::OptimalLAIModel{FT}) where {FT} = (FT,)
ClimaLand.auxiliary_domain_names(::OptimalLAIModel) = (:surface,)

ClimaLand.name(::AbstractLAIModel) = :lai_model

"""
    initialize_LAI!(p, model::OptimalLAIModel, initial_value::FT = FT(1.0)) where {FT}

Initialize the LAI field to a given initial value. Default is 1.0 m² m⁻².
"""
function initialize_LAI!(p, model::OptimalLAIModel, initial_value = eltype(model)(1.0))
    p.canopy.lai_model.LAI .= initial_value
end

"""
    compute_L_max(Ao_annual, P_annual, D_growing, k, z, ca, chi, f0)

Compute seasonal maximum leaf area index (LAI_max) based on energy and water limitations.

This implements Equations 11-12 from Zhou et al. (2025), which predicts LAI_max as the
minimum of an energy-limited rate (maximizing GPP) and a water-limited rate (maximizing
the use of available precipitation).

# Arguments
- `Ao_annual::FT`: Annual total potential GPP (mol m⁻² yr⁻¹). This is the integral of daily
  A₀ over the year. A₀ is the GPP that would be achieved if fAPAR = 1.
- `P_annual::FT`: Annual total precipitation (mol m⁻² yr⁻¹). Note: 1 mol H₂O ≈ 18 g
- `D_growing::FT`: Mean vapor pressure deficit during the growing season when T > 0°C (Pa)
- `k::FT`: Light extinction coefficient (dimensionless), typically 0.5
- `z::FT`: Unit cost of constructing and maintaining leaves (mol m⁻² yr⁻¹), globally fitted
  as 12.227 mol m⁻² yr⁻¹
- `ca::FT`: Ambient CO₂ partial pressure (Pa)
- `chi::FT`: Ratio of leaf-internal to ambient CO₂ partial pressure (dimensionless),
  typically from the P-model
- `f0::FT`: Fraction of annual precipitation used by plants (dimensionless), varies with
  aridity index (Equation 6)

# Returns
- `LAI_max::FT`: Seasonal maximum leaf area index (m² m⁻², dimensionless)

# Example values
```julia
# IMPORTANT NOTE: This equation as published produces LAI values that are typically
# much smaller than observed forest LAI when calculated for individual sites.
# The paper (Zhou et al. 2025) reports good global-scale agreement with satellite data,
# but there may be discrepancies at local scales. The equation is implemented exactly
# as specified in Equations 11-12 of the paper.

# Example 1: Dry temperate conditions
Ao_annual = 100.0   # mol m⁻² yr⁻¹ (typical temperate forest)
P_annual = 30000.0    # mol m⁻² yr⁻¹ (~540 mm precipitation)
D_growing = 2500.0     # Pa (mean VPD during growing season)
k = 0.5               # dimensionless
z = 12.227            # mol m⁻² yr⁻¹
ca = 40.0             # Pa (~400 ppm)
chi = 0.7             # dimensionless
f0 = 0.62             # dimensionless

LAI_max = compute_L_max(Ao_annual, P_annual, D_growing, k, z, ca, chi, f0)
# Returns: 1.6 m² m⁻² (strongly water-limited)

# Example 2: Wet temperate conditions
P_annual = 150000.0   # mol m⁻² yr⁻¹ (~2700 mm precipitation)
D_growing = 600.0     # Pa (low VPD, humid)
chi = 0.78            # dimensionless
f0 = 0.75             # dimensionless

LAI_max = compute_L_max(Ao_annual, P_annual, D_growing, k, z, ca, chi, f0)
# Returns: 2.8 m² m⁻² (water-limited)
```

# References
Zhou et al. (2025) Global Change Biology, Equations 11-12
"""
function compute_L_max(
    Ao_annual::FT,  # mol m⁻² yr⁻¹
    P_annual::FT,   # mol m⁻² yr⁻¹
    D_growing::FT,  # Pa
    k::FT,          # dimensionless
    z::FT,          # mol m⁻² yr⁻¹
    ca::FT,         # Pa
    chi::FT,        # dimensionless
    f0::FT,         # dimensionless
) where {FT}
    # Energy-limited fAPAR (Equation 11, first part)
    fAPAR_energy = 1 - z / (k * Ao_annual)

    # Water-limited fAPAR (Equation 11, second part)
    fAPAR_water = (ca * (1 - chi) / (1.6 * D_growing)) * (f0 * P_annual / Ao_annual)

    # Take minimum of energy and water limitations (Equation 11)
    fAPAR_max = min(fAPAR_energy, fAPAR_water)

    # Ensure fAPAR is in valid range [0, 1]
    fAPAR_max = max(0, min(1, fAPAR_max))

    # Convert fAPAR to LAI using Beer's law (Equation 12)
    # LAI_max = -(1/k) * ln(1 - fAPAR_max)
    LAI_max = -(1 / k) * log(1 - fAPAR_max)

    return LAI_max
end

fAPAR_max_fun(k, LAI_max) = 1 - exp(-k * LAI_max)

"""
    compute_m(GSL, LAI_max, Ao_annual, sigma)

Compute the parameter m, which represents the ratio of steady-state LAI to steady-state GPP.

This implements Equation 20 from Zhou et al. (2025). The parameter m quantifies the
relationship between LAI and GPP dynamics, representing the extent to which seasonal LAI
dynamics depart from a "square wave" (where maximum LAI would be maintained throughout
the growing season).

# Arguments
- `GSL::FT`: Growing season length (days). Defined as the length of continuous period
  above 0°C longer than 5 days.
- `LAI_max::FT`: Seasonal maximum leaf area index (m² m⁻², dimensionless)
- `Ao_annual::FT`: Annual total potential GPP (mol m⁻² yr⁻¹). This is the integral of daily
  A₀ over the year.
- `sigma::FT`: Dimensionless parameter representing departure from square-wave LAI dynamics.
  Globally fitted as σ = 0.771
- `k::FT`: Light extinction coefficient (dimensionless)

# Returns
- `m::FT`: Parameter relating steady-state LAI to steady-state GPP (dimensionless, units
  work out as: days × m²m⁻² / (mol m⁻² yr⁻¹ × dimensionless) with implicit conversion)

# Example values
```julia
GSL = 180.0          # days (6-month growing season)
LAI_max = 1.6        # m² m⁻²
Ao_annual = 100.0     # mol m⁻² yr⁻¹
sigma = 0.771        # dimensionless
k = 0.5

m = compute_m(GSL, LAI_max, Ao_annual, sigma, k)
# Returns: ~4 (dimensionless, after unit conversion)
```

# References
Zhou et al. (2025) Global Change Biology, Equation 20
"""
function compute_m(
    GSL::FT,        # days
    LAI_max::FT,    # m² m⁻² (dimensionless)
    Ao_annual::FT,     # mol m⁻² yr⁻¹
    sigma::FT,      # dimensionless
    k::FT,
) where {FT}
    # Equation 20: m = (σ × GSL × LAI_max) / (A₀_sum × fAPAR_max)
    fAPAR_max = fAPAR_max_fun(k, LAI_max)
    m = (sigma * GSL * LAI_max) / (Ao_annual * fAPAR_max)
    return m
end

"""
    lambertw0(x)

Compute the principal branch (W₀) of the Lambert W function for x ∈ [-1/e, ∞).

The Lambert W function satisfies W(x)·exp(W(x)) = x. This implementation uses
Halley's method for fast convergence, typically requiring only 2-3 iterations.

# Arguments
- `x::Real`: Input value, must be ≥ -1/e ≈ -0.36788

# Returns
- `W::Float64`: Lambert W₀(x), the principal branch value

# Algorithm
Uses Halley's method with an appropriate initial guess:
- For x near -1/e: use series expansion
- For x ∈ [-1/e, -0.1]: use fitted approximation
- For x ∈ [-0.1, 10]: use log-based approximation
- For x > 10: use asymptotic expansion

# References
Corless et al. (1996) "On the Lambert W function"
"""
function lambertw0(x::T) where {T<:Real}
    # Check domain
    min_x = -1 / ℯ
    if x < min_x
        throw(DomainError(x, "Lambert W₀ is only defined for x ≥ -1/e ≈ -0.36788"))
    end

    # Special cases
    if x == 0
        return zero(T)
    elseif abs(x - min_x) < 1e-10
        return -one(T)
    end

    # Choose initial guess based on the region
    if x < -0.32  # Near the branch point -1/e
        # Series expansion near -1/e
        p = sqrt(2 * (ℯ * x + 1))
        w = -1 + p - p^2/3 + p^3*11/72
    elseif x < 0
        # For x ∈ [-0.32, 0], use a rational approximation
        w = x / (1 + x)  # Simple approximation that's good enough for Halley
    elseif x < 2.5
        # For small positive x, start with x as the guess (works well up to ~2.5)
        w = x * (1 - x/3)  # Slightly better than just x
    elseif x < 10
        # Log-based approximation (safe since x >= 2.5)
        l1 = log(x)
        l2 = log(l1)
        w = l1 - l2 + l2 / l1
    else
        # Asymptotic expansion for large x
        l1 = log(x)
        l2 = log(l1)
        w = l1 - l2 + l2 / l1 + l2 * (l2 - 2) / (2 * l1 * l1)
    end

    # Halley's method refinement (typically converges in 2-3 iterations)
    for _ in 1:10  # Maximum iterations
        ew = exp(w)
        wew = w * ew
        f = wew - x

        # Check convergence
        if abs(f) < 1e-14 * (1 + abs(x))
            break
        end

        # Halley's method update
        # w_new = w - f / (f' - f * f'' / (2 * f'))
        # where f = w*exp(w) - x
        # f' = exp(w) * (w + 1)
        # f'' = exp(w) * (w + 2)
        w1 = w + 1
        denom = ew * w1 - f * (w + 2) / (2 * w1)

        if abs(denom) < 1e-20
            break
        end

        w = w - f / denom
    end

    return w
end



"""
    compute_steady_state_LAI(Ao_daily, m, k, LAI_max)

Compute steady-state LAI from daily potential GPP using the Lambert W function solution.

This implements Equations 13-15 from Zhou et al. (2025). The steady-state LAI (L_s) is
the LAI that would be in equilibrium with GPP if weather conditions were held constant.
Given daily meteorological conditions, this is computed on a daily basis.

# Arguments
- `Ao_daily::FT`: Daily potential GPP (mol m⁻² day⁻¹). This is the GPP that would be
  achieved if fAPAR = 1, calculated from LUE × PPFD.
- `m::FT`: Parameter relating steady-state LAI to steady-state GPP (dimensionless), from
  `compute_m()`
- `k::FT`: Light extinction coefficient (dimensionless), typically 0.5
- `LAI_max::FT`: Seasonal maximum LAI constraint (m² m⁻², dimensionless)

# Returns
- `L_steady::FT`: Steady-state leaf area index (m² m⁻², dimensionless). Always ≥ 0.

# Example values
```julia
Ao_daily = 0.4      # mol m⁻² day⁻¹ (summer day in temperate forest)
m = 7.0           # dimensionless (from compute_m)
k = 0.5             # dimensionless
LAI_max = 3.0       # m² m⁻²

L_steady = compute_steady_state_LAI(Ao_daily, m, k, LAI_max)
# Returns: ~1.4 m² m⁻² (close to LAI_max during peak growing season)
```

# Notes
The solution uses the Lambert W₀ function: L_s = min{μ + (1/k)W₀[-kμ exp(-kμ)], LAI_max}
where μ = m × A₀. The result is constrained to be non-negative and below LAI_max.

# References
Zhou et al. (2025) Global Change Biology, Equations 13-15
"""
function compute_steady_state_LAI(
    Ao_daily::FT,  # mol m⁻² day⁻¹
    m::FT,         # dimensionless
    k::FT,         # dimensionless
    LAI_max::FT,   # m² m⁻² (dimensionless)
) where {FT}
    # μ = m × A₀ (Equation 15)
    mu = m * Ao_daily

    # Compute argument for Lambert W function
    arg = -k * mu * exp(-k * mu)

    # Check if argument is in valid range for W₀ branch: [-1/e, 0]
    # If outside this range, use boundary solution
    if arg < -1 / exp(1)
        # Beyond valid range; use maximum possible LAI
        L_s = LAI_max
    else
        # Equation 15: L_s = μ + (1/k) × W₀[-k μ exp(-k μ)]
        # Using our custom lambertw0 function (W₀ is the principal branch)
        w_val = lambertw0(arg)
        L_s = mu + (1 / k) * w_val

        # Take minimum with LAI_max (Equation 15)
        L_s = min(L_s, LAI_max)
    end

    # Ensure non-negative (should be guaranteed mathematically, but enforce for numerical stability)
    L_s = max(zero(FT), L_s)

    return L_s
end

"""
    update_LAI!(LAI_prev, L_steady, alpha, local_noon_mask)

Update LAI using exponential weighted moving average to represent lag between carbon
allocation and steady-state LAI.

This implements Equation 16 from Zhou et al. (2025). The exponential moving average
represents the time lag (days to months) for photosynthate allocation to leaves and
leaf development.

# Arguments
- `LAI_prev::FT`: LAI from previous time step (m² m⁻², dimensionless)
- `L_steady::FT`: Current steady-state LAI (m² m⁻², dimensionless), from
  `compute_steady_state_LAI()`
- `alpha::FT`: Smoothing factor (dimensionless, 0-1). Set to 0.067 for ~15 days of memory,
  meaning LAI[t] = 0.067 × L_steady[t] + 0.933 × LAI[t-1]
- `local_noon_mask::FT`: A mask (0 or 1) indicating whether the current time is within the local noon window.

# Returns
- `LAI_new::FT`: Updated leaf area index (m² m⁻², dimensionless). Always ≥ 0.

# Example values
```julia
LAI_prev = 3.5      # m² m⁻² (LAI from yesterday)
L_steady = 4.0      # m² m⁻² (steady-state LAI for today)
alpha = 0.067       # dimensionless (15-day memory)
local_noon_mask = 1.0  # update at local noon

LAI_new = update_LAI!(LAI_prev, L_steady, alpha, local_noon_mask)
# Returns: ~3.53 m² m⁻² (slowly approaching steady-state)
```

# Notes
The parameter α controls the response time:
- α = 0.067 ≈ 15 days of memory (paper default)
- α = 0.1 ≈ 10 days of memory (faster response)
- α = 0.033 ≈ 30 days of memory (slower response)

The time scale τ ≈ 1/α days.

# References
Zhou et al. (2025) Global Change Biology, Equation 16
"""
function update_LAI!(
        LAI_prev::FT,   # m² m⁻² (dimensionless)
        L_steady::FT,   # m² m⁻² (dimensionless)
        alpha::FT,      # dimensionless (0-1)
        local_noon_mask::FT,
    ) where {FT}
    if local_noon_mask == FT(1.0)
        # Equation 16: LAI_sim = α × L_s[t] + (1-α) × LAI_sim[t-1]
        LAI_new = alpha * L_steady + (1 - alpha) * LAI_prev

        # Ensure non-negative (important for numerical stability, especially at initialization)
        LAI_new = max(zero(FT), LAI_new)
        return LAI_new
    else
        return LAI_prev
    end
end

function compute_Ao_daily(A::FT, k::FT, L::FT) where {FT}
    Ao_daily_inst = A/max(1-exp(-k*L), eps(FT))
    Ao_daily = Ao_daily_inst * 3600 * 8 # 3600 sec in an hour, 8 hour of sunlight
    # Note: would be better to integrate daily A instead of assuming noon * 8 hours...
    return Ao_daily
end

function update_optimal_LAI(
        local_noon_mask::FT,
        A::FT,
        L::FT; # m2 m-2
        # Ao_daily = 100.0, # this one should be computed from A?
        k = FT(0.5),
        Ao_annual = FT(100.0), # mol CO2 m-2 y-1, for an average forest
        P_annual = FT(60000.0), # ~ 1000 mm precipitation per year,
        D_growing = FT(1000.0), # mean VPD - also some average value,
        z = FT(12.227), # mol m-2 yr-1, param in the paper,
        ca = FT(40.0), # co2 concentration (400ppm),
        chi = FT(0.7), # dimensionless,
        f0 = FT(0.62), # dimensionless,
        GSL = FT(180.0), # days, growing season length
        sigma = FT(0.771), # dimensionless
        alpha = FT(0.067), # dimensionless (15-day memory)
    ) where {FT}
    Ao_daily = compute_Ao_daily(A, k, L)
    LAI_max = compute_L_max(Ao_annual, P_annual, D_growing, k, z, ca, chi, f0)
    m = compute_m(GSL, LAI_max, Ao_annual, sigma, k)
    L_steady = compute_steady_state_LAI(Ao_daily, m, k, LAI_max)
    L = update_LAI!(L, L_steady, alpha, local_noon_mask)
    return L
end
