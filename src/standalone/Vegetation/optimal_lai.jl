export AbstractLAIModel,
    OptimalLAIModel,
    OptimalLAIParameters,
    update_LAI!,
    initialize_LAI!,
    set_historical_cache!,
    make_OptimalLAI_callback,
    call_update_optimal_LAI

"""
    AbstractLAIModel{FT <: AbstractFloat}

An abstract type for LAI (Leaf Area Index) model parameterizations.
"""
abstract type AbstractLAIModel{FT <: AbstractFloat} <:
              AbstractCanopyComponent{FT} end

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
    OptimalLAIParameters{FT}(toml_dict::CP.ParamDict) where {FT}

Creates an `OptimalLAIParameters` object from a TOML parameter dictionary.
"""
function OptimalLAIParameters{FT}(toml_dict::CP.ParamDict) where {FT}
    return OptimalLAIParameters{FT}(
        k = FT(toml_dict["optimal_lai_k"]),
        z = FT(toml_dict["optimal_lai_z"]),
        chi = FT(toml_dict["optimal_lai_chi"]),
        f0 = FT(toml_dict["optimal_lai_f0"]),
        sigma = FT(toml_dict["optimal_lai_sigma"]),
        alpha = FT(toml_dict["optimal_lai_alpha"]),
    )
end

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

Defines the auxiliary variables for the OptimalLAIModel:
- `LAI`: leaf area index (m² m⁻²)
- `A0_daily`: daily potential GPP from previous day (mol CO₂ m⁻² day⁻¹)
- `A0_annual`: annual potential GPP from previous year (mol CO₂ m⁻² yr⁻¹)
- `A0_daily_acc`: accumulator for current day's potential GPP (mol CO₂ m⁻² day⁻¹)
- `A0_annual_acc`: accumulator for current year's potential GPP (mol CO₂ m⁻² yr⁻¹)
- `last_day_of_year`: day of year of last update (for detecting year change)
"""
ClimaLand.auxiliary_vars(model::OptimalLAIModel) =
    (:LAI, :A0_daily, :A0_annual, :A0_daily_acc, :A0_annual_acc, :last_day_of_year)
ClimaLand.auxiliary_types(model::OptimalLAIModel{FT}) where {FT} =
    (FT, FT, FT, FT, FT, FT)
ClimaLand.auxiliary_domain_names(::OptimalLAIModel) =
    (:surface, :surface, :surface, :surface, :surface, :surface)

ClimaLand.name(::AbstractLAIModel) = :lai_model

"""
    initialize_LAI!(p, model::OptimalLAIModel; initial_LAI, initial_A0_daily, initial_A0_annual)

Initialize the LAI and A0 fields to given initial values.

# Arguments
- `p`: auxiliary state
- `model::OptimalLAIModel`: the optimal LAI model
- `initial_LAI`: initial LAI value (m² m⁻², default 1.0)
- `initial_A0_daily`: initial daily potential GPP (mol CO₂ m⁻² day⁻¹, default 0.5)
- `initial_A0_annual`: initial annual potential GPP (mol CO₂ m⁻² yr⁻¹, default 258.0)
"""
function initialize_LAI!(
    p,
    model::OptimalLAIModel;
    initial_LAI = eltype(model)(1.0),
    initial_A0_daily = eltype(model)(0.5),
    initial_A0_annual = eltype(model)(258.0),
)
    FT = eltype(model)
    p.canopy.lai_model.LAI .= initial_LAI
    p.canopy.lai_model.A0_daily .= initial_A0_daily
    p.canopy.lai_model.A0_annual .= initial_A0_annual
    p.canopy.lai_model.A0_daily_acc .= FT(0)
    p.canopy.lai_model.A0_annual_acc .= FT(0)
    p.canopy.lai_model.last_day_of_year .= FT(0)
end

"""
    set_historical_cache!(p, Y0, model::OptimalLAIModel, canopy; A0_annual, A0_daily, P_annual, D_growing, ca, GSL)

The optimal LAI model requires initialization of LAI and A0 values before the simulation.
This method assumes that the LAI is initially in equilibrium with the initial conditions
of the simulation.

# Arguments
- `A0_annual`: Initial annual potential GPP (mol CO₂ m⁻² yr⁻¹), default 258.0
- `A0_daily`: Initial daily potential GPP (mol CO₂ m⁻² day⁻¹), default 0.5
- `P_annual`: Annual precipitation (mol H₂O m⁻² yr⁻¹), default 60000.0 (~1080 mm)
- `D_growing`: Mean VPD during growing season (Pa), default 1000.0
- `ca`: Ambient CO₂ partial pressure (Pa), default 40.0 (~400 ppm)
- `GSL`: Growing season length (days), default 180.0

An alternative to this approach is to initialize the initial values to some reasonable values
based on a spun-up simulation or observational data.
"""
function set_historical_cache!(
    p,
    Y0,
    model::OptimalLAIModel,
    canopy;
    A0_annual = eltype(model.parameters)(258.0),
    A0_daily = eltype(model.parameters)(0.5),
    P_annual = eltype(model.parameters)(60000.0),
    D_growing = eltype(model.parameters)(1000.0),
    ca = eltype(model.parameters)(40.0),
    GSL = eltype(model.parameters)(240.0),
)
    parameters = model.parameters
    FT = eltype(parameters)

    L = p.canopy.lai_model.LAI

    # Initialize LAI with dummy value which will be overwritten
    fill!(L, FT(1.0))

    # Initialize A0 variables
    fill!(p.canopy.lai_model.A0_daily, A0_daily)
    fill!(p.canopy.lai_model.A0_annual, A0_annual)
    fill!(p.canopy.lai_model.A0_daily_acc, FT(0))
    fill!(p.canopy.lai_model.A0_annual_acc, FT(0))
    fill!(p.canopy.lai_model.last_day_of_year, FT(0))

    # Set local_noon_mask to 1 to force update for initialization
    local_noon_mask = FT(1)

    # Compute initial LAI based on initial conditions
    # Note: ca here is already in Pa (as expected by the optimal LAI model)
    @. L = update_optimal_LAI(
        local_noon_mask,
        A0_daily,
        L,
        parameters.k,
        A0_annual,
        P_annual,
        D_growing,
        parameters.z,
        ca,
        parameters.chi,
        parameters.f0,
        GSL,
        parameters.sigma,
        parameters.alpha,
    )
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
    fAPAR_energy = FT(1) - z / (k * Ao_annual)

    # Water-limited fAPAR (Equation 11, second part)
    fAPAR_water =
        (ca * (FT(1) - chi) / (FT(1.6) * D_growing)) *
        (f0 * P_annual / Ao_annual)

    # Take minimum of energy and water limitations (Equation 11)
    fAPAR_max = min(fAPAR_energy, fAPAR_water)

    # Ensure fAPAR is in valid range [0, 1]
    fAPAR_max = max(FT(0), min(FT(1), fAPAR_max))

    # Convert fAPAR to LAI using Beer's law (Equation 12)
    # LAI_max = -(1/k) * ln(1 - fAPAR_max)
    LAI_max = -(FT(1) / k) * log(FT(1) - fAPAR_max)

    return LAI_max
end

fAPAR_max_fun(k::FT, LAI_max::FT) where {FT} = FT(1) - exp(-k * LAI_max)

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

const MINARG = -inv(Base.MathConstants.e)

"""
    _lambertw0_initial_guess(x::T) where {T<:AbstractFloat}

Provide a robust initial guess for the Lambert W₀ function for use in iterative solvers.

# Arguments
- `x::T`: Input value, should be ≥ -1/e

# Returns
- Initial guess for W₀(x)

# Algorithm
- For x > 1: uses log(x) - log(log(x)) approximation
- For x < -0.32 (near -1/e): uses series expansion for accurate convergence near branch point
- For -0.32 ≤ x ≤ 1: uses max(x, -0.3) as a simple starting point
"""
@inline function _lambertw0_initial_guess(x::T) where {T <: AbstractFloat}
    if x > one(T)
        return log(x) - log(max(log(x), T(1e-6)))
    elseif x < T(-0.32)
        # Near the branch point -1/e, use series expansion
        # This handles the singular behavior at x = -1/e where W(x) = -1
        p = sqrt(T(2) * (T(ℯ) * x + one(T)))
        return -one(T) + p - p^2 / T(3) + p^3 * T(11) / T(72)
    else
        return max(x, T(-0.3))
    end
end

"""
    lambertw0(x::T; maxiter::Int = 16) where {T<:AbstractFloat}

Compute the principal branch (W₀) of the Lambert W function for x ∈ [-1/e, ∞).

This is a GPU-device-friendly implementation using a fixed number of Halley iterations.
The Lambert W function satisfies W(x)·exp(W(x)) = x.

# Arguments
- `x::T`: Input value, must be ≥ -1/e ≈ -0.36788
- `maxiter::Int`: Maximum number of Halley iterations (default: 16)

# Returns
- `W::T`: Lambert W₀(x), the principal branch value, or NaN for invalid inputs

# Algorithm
Uses Halley's method with a fixed number of iterations for GPU compatibility:
- No dynamic memory allocation
- No conditional breaks (runs all iterations)
- Broadcastable for use with CuArrays: lambertw0.(cuarray)

# Device Compatibility
This implementation is designed to work on both CPU and GPU:
- All operations are scalar and supported on CUDA.jl
- No array allocations or dynamic loops
- Type-generic over AbstractFloat (Float32, Float64)

# References
Corless et al. (1996) "On the Lambert W function"
"""
@inline function lambertw0(x::T; maxiter::Int = 16) where {T <: AbstractFloat}
    if !(isfinite(x)) || x < T(MINARG)
        return T(NaN)
    end
    w = _lambertw0_initial_guess(x)
    for i in 1:maxiter
        ew = exp(w)
        f = w * ew - x
        # Halley denominator
        # Special case: when w ≈ -1, both numerator and denominator approach 0
        # This happens at the branch point x = -1/e, where W(-1/e) = -1
        w_plus_1 = w + one(T)
        if abs(w_plus_1) < eps(T)
            # Already at or very near the solution w = -1, no update needed
            Δ = zero(T)
        else
            two_w_plus_2 = T(2) * w_plus_1
            if abs(two_w_plus_2) < eps(T)
                # Near w = -1, use Newton's method instead of Halley
                Δ = f / (ew * w_plus_1)
            else
                denom = ew * w_plus_1 - (w + T(2)) * f / two_w_plus_2
                if abs(denom) < eps(T)
                    Δ = f / (ew * w_plus_1)
                else
                    Δ = f / denom
                end
            end
        end
        w -= Δ
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
    if arg < -FT(1) / FT(exp(1))
        # Beyond valid range; use maximum possible LAI
        L_s = LAI_max
    else
        # Equation 15: L_s = μ + (1/k) × W₀[-k μ exp(-k μ)]
        # Using our custom lambertw0 function (W₀ is the principal branch)
        w_val = lambertw0(arg)
        L_s = mu + (FT(1) / k) * w_val

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

"""
    update_optimal_LAI(local_noon_mask, A0_daily, L, k, A0_annual, P_annual, D_growing, z, ca, chi, f0, GSL, sigma, alpha)

Update LAI using the optimal LAI model with precomputed daily and annual potential GPP.

# Arguments
- `local_noon_mask::FT`: Mask (0 or 1) indicating if it's local noon
- `A0_daily::FT`: Daily potential GPP (mol CO₂ m⁻² day⁻¹)
- `L::FT`: Current LAI (m² m⁻²)
- `k::FT`: Light extinction coefficient
- `A0_annual::FT`: Annual potential GPP (mol CO₂ m⁻² yr⁻¹)
- Other parameters as in compute_L_max and compute_m

# Returns
Updated LAI value.
"""
function update_optimal_LAI(
    local_noon_mask::FT,
    A0_daily::FT,
    L::FT, # m2 m-2
    k::FT,
    A0_annual::FT, # mol CO2 m-2 y-1
    P_annual::FT, # mol H2O m-2 y-1
    D_growing::FT, # Pa, mean VPD during growing season
    z::FT, # mol m-2 yr-1, leaf construction cost
    ca::FT, # Pa, CO2 partial pressure
    chi::FT, # dimensionless, ci/ca ratio
    f0::FT, # dimensionless, precip fraction
    GSL::FT, # days, growing season length
    sigma::FT, # dimensionless
    alpha::FT, # dimensionless (15-day memory)
) where {FT}
    LAI_max = compute_L_max(A0_annual, P_annual, D_growing, k, z, ca, chi, f0)
    m = compute_m(GSL, LAI_max, A0_annual, sigma, k)
    L_steady = compute_steady_state_LAI(A0_daily, m, k, LAI_max)
    L = update_LAI!(L, L_steady, alpha, local_noon_mask)
    return L
end

"""
    compute_PPFD(par_d::FT, λ_γ_PAR::FT, lightspeed::FT, planck_h::FT, N_a::FT) where {FT}

Compute photosynthetic photon flux density (PPFD) from PAR downwelling flux.
This is the total incoming PAR (not absorbed), in units of mol photons m⁻² s⁻¹.
"""
function compute_PPFD(
    par_d::FT,
    λ_γ_PAR::FT,
    lightspeed::FT,
    planck_h::FT,
    N_a::FT,
) where {FT}
    energy_per_mole_photon_par = planck_h * lightspeed * N_a / λ_γ_PAR
    return par_d / energy_per_mole_photon_par
end

"""
    call_update_optimal_LAI(p, Y, t, current_date; canopy, dt, local_noon, P_annual, D_growing, GSL)

Updates LAI and accumulates A0 at each timestep. At local noon, finalizes daily A0
and updates LAI. On year change (Jan 1), finalizes annual A0.

The function:
1. Computes instantaneous A0 using the P-model with fAPAR=1, β=1
2. Accumulates A0 to daily accumulator (converted to mol CO₂ m⁻² per timestep)
3. At local noon: finalizes daily A0, adds to annual accumulator, updates LAI
4. On Jan 1: finalizes annual A0
"""
function call_update_optimal_LAI(
    p,
    Y,
    t,
    current_date;
    canopy,
    dt,
    local_noon,
    P_annual,
    D_growing,
    GSL,
)
    FT = eltype(canopy.lai_model.parameters)

    # Compute local noon mask
    local_noon_mask = @. lazy(get_local_noon_mask(t, dt, local_noon))

    # Get P-model parameters and constants for computing A0
    pmodel_parameters = canopy.photosynthesis.parameters
    pmodel_constants = canopy.photosynthesis.constants
    is_c3 = canopy.photosynthesis.is_c3

    # Get drivers for A0 computation
    T_canopy = canopy_temperature(canopy.energy, canopy, Y, p)
    P_air = p.drivers.P
    ca = p.drivers.c_co2  # mol/mol
    earth_param_set = canopy.earth_param_set

    # Compute VPD (clipped to avoid numerical issues)
    VPD = @. lazy(
        max(
            ClimaLand.vapor_pressure_deficit(
                p.drivers.T,
                p.drivers.P,
                p.drivers.q,
                LP.thermodynamic_parameters(earth_param_set),
            ),
            sqrt(eps(FT)),
        ),
    )

    # Compute PPFD from PAR downwelling (total incoming, not absorbed)
    par_d = p.canopy.radiative_transfer.par_d
    λ_γ_PAR = canopy.radiative_transfer.parameters.λ_γ_PAR
    PPFD = @. lazy(
        compute_PPFD(
            par_d,
            λ_γ_PAR,
            pmodel_constants.lightspeed,
            pmodel_constants.planck_h,
            pmodel_constants.N_a,
        ),
    )

    # Compute instantaneous A0 (kg C m⁻² s⁻¹)
    A0_inst = @. lazy(
        compute_A0_potential(
            is_c3,
            pmodel_parameters,
            pmodel_constants,
            T_canopy,
            P_air,
            VPD,
            ca,
            PPFD,
        ),
    )

    # Convert A0 from kg C m⁻² s⁻¹ to mol CO₂ m⁻² per timestep
    # Mc = 0.0120107 kg/mol C
    dt_seconds = FT(float(dt))
    Mc = pmodel_constants.Mc
    A0_per_timestep = @. lazy(A0_inst * dt_seconds / Mc)

    # Accumulate to daily A0
    @. p.canopy.lai_model.A0_daily_acc += A0_per_timestep

    # Get parameters from the LAI model
    parameters = canopy.lai_model.parameters

    # Current day of year (scalar)
    current_doy = FT(Dates.dayofyear(current_date))

    # Convert ca from mol/mol to Pa for LAI computation
    ca_Pa = @. lazy(ca * P_air)

    # At local noon: finalize daily A0, update annual accumulator, update LAI
    # Inline the update logic to avoid tuple broadcasting issues.
    # Compute finalized A0 values before any resets.
    A0_daily_final = @. p.canopy.lai_model.A0_daily_acc
    A0_annual_final = @. ifelse(
        local_noon_mask == FT(1) &&
        current_doy < p.canopy.lai_model.last_day_of_year,
        p.canopy.lai_model.A0_annual_acc,
        p.canopy.lai_model.A0_annual,
    )

    # Update LAI at noon using the finalized values
    @. p.canopy.lai_model.LAI = ifelse(
        local_noon_mask == FT(1),
        update_optimal_LAI(
            FT(1),  # At noon, we always update
            A0_daily_final,
            p.canopy.lai_model.LAI,
            parameters.k,
            A0_annual_final, # already accounts for year change
            P_annual,
            D_growing,
            parameters.z,
            ca_Pa,
            parameters.chi,
            parameters.f0,
            GSL,
            parameters.sigma,
            parameters.alpha,
        ),
        p.canopy.lai_model.LAI,
    )

    # Update A0_annual if year changed
    @. p.canopy.lai_model.A0_annual = ifelse(
        local_noon_mask == FT(1) &&
        current_doy < p.canopy.lai_model.last_day_of_year,
        A0_annual_final,
        p.canopy.lai_model.A0_annual,
    )

    # Update A0_daily at noon
    @. p.canopy.lai_model.A0_daily = ifelse(
        local_noon_mask == FT(1),
        A0_daily_final,
        p.canopy.lai_model.A0_daily,
    )

    # Reset A0_annual_acc if year changed, otherwise add new daily value at noon
    @. p.canopy.lai_model.A0_annual_acc = ifelse(
        local_noon_mask == FT(1) &&
        current_doy < p.canopy.lai_model.last_day_of_year,
        FT(0),
        ifelse(
            local_noon_mask == FT(1),
            p.canopy.lai_model.A0_annual_acc + A0_daily_final,
            p.canopy.lai_model.A0_annual_acc,
        ),
    )

    # Reset A0_daily_acc at noon
    @. p.canopy.lai_model.A0_daily_acc = ifelse(
        local_noon_mask == FT(1),
        FT(0),
        p.canopy.lai_model.A0_daily_acc,
    )

    # Update last_day_of_year at noon
    @. p.canopy.lai_model.last_day_of_year = ifelse(
        local_noon_mask == FT(1),
        current_doy,
        p.canopy.lai_model.last_day_of_year,
    )
end

"""
    update_A0_and_LAI_at_noon(local_noon_mask, A0_daily, A0_annual, A0_daily_acc, A0_annual_acc, last_doy, current_doy, L, k, P_annual, D_growing, z, ca, chi, f0, GSL, sigma, alpha)

At local noon: finalize daily A0, add to annual, reset daily accumulator, and update LAI.
On year change (current_doy < last_doy, indicating Jan 1): finalize annual A0.

Returns tuple: (A0_daily, A0_annual, A0_daily_acc, A0_annual_acc, last_doy, L)
"""
function update_A0_and_LAI_at_noon(
    local_noon_mask::FT,
    A0_daily::FT,
    A0_annual::FT,
    A0_daily_acc::FT,
    A0_annual_acc::FT,
    last_doy::FT,
    current_doy::FT,
    L::FT,
    k::FT,
    P_annual::FT,
    D_growing::FT,
    z::FT,
    ca::FT,
    chi::FT,
    f0::FT,
    GSL::FT,
    sigma::FT,
    alpha::FT,
) where {FT}
    if local_noon_mask == FT(1)
        # Check for year change (day of year decreased = new year started)
        year_changed = current_doy < last_doy

        if year_changed
            # Finalize annual A0 before processing the new year
            A0_annual = A0_annual_acc
            A0_annual_acc = FT(0)
        end

        # Finalize daily A0
        A0_daily = A0_daily_acc
        A0_daily_acc = FT(0)

        # Add today's A0 to annual accumulator
        A0_annual_acc += A0_daily

        # Update last day of year
        last_doy = current_doy

        # Update LAI using the optimal LAI model
        L = update_optimal_LAI(
            local_noon_mask,
            A0_daily,
            L,
            k,
            A0_annual,
            P_annual,
            D_growing,
            z,
            ca,
            chi,
            f0,
            GSL,
            sigma,
            alpha,
        )

        return A0_daily, A0_annual, A0_daily_acc, A0_annual_acc, last_doy, L
    else
        return A0_daily, A0_annual, A0_daily_acc, A0_annual_acc, last_doy, L
    end
end

"""
    make_OptimalLAI_callback(::Type{FT}, t0::ITime, dt, canopy; P_annual, D_growing, GSL, longitude) where {FT <: AbstractFloat}

This constructs an IntervalBasedCallback for the optimal LAI model that:
1. Computes and accumulates potential GPP (A₀) at each timestep
2. Updates LAI using an exponential moving average at local noon
3. Tracks daily and annual A₀ sums

We check for local noon using the provided `longitude` every dt.
The time of local noon is expressed in seconds UTC and neglects the effects of obliquity and eccentricity, so
it is constant throughout the year.

# Arguments
- `FT`: The floating-point type used in the model (e.g., `Float32`, `Float64`).
- `t0`: ITime, with epoch in UTC.
- `dt`: timestep
- `canopy`: the canopy object containing the optimal LAI model parameters.
- `P_annual`: Annual total precipitation (mol H₂O m⁻² yr⁻¹), default FT(60000.0) (~1080 mm)
- `D_growing`: Mean vapor pressure deficit during the growing season (Pa), default FT(1000.0)
- `GSL`: Growing season length (days), default FT(180.0)
- `longitude`: optional longitude in degrees for local noon calculation (default is `nothing`, which means
    that it will be inferred from the canopy domain).

# Notes
- A₀ (potential GPP) is computed dynamically using the P-model with fAPAR=1 and β=1
- Daily A₀ is accumulated over each day and finalized at local noon
- Annual A₀ is accumulated and reset on January 1
"""
function make_OptimalLAI_callback(
    ::Type{FT},
    t0,
    dt,
    canopy;
    P_annual = FT(60000.0),
    D_growing = FT(1000.0),
    GSL = FT(240.0),
    longitude = nothing,
) where {FT <: AbstractFloat}
    function seconds_after_midnight(d)
        return FT(
            Hour(d).value * 3600 +
            Minute(d).value * 60 +
            Second(d).value,
        )
    end

    if isnothing(longitude)
        try
            longitude = get_long(canopy.domain.space.surface)
        catch e
            error(
                "Longitude must be provided explicitly if the domain you are working on does not \
                  have axes that specify longitude $e",
            )
        end
    end

    # this computes the time of local noon in seconds UTC without considering the
    # effects of obliquity and orbital eccentricity, so it is constant throughout the year
    # the max error is on the order of 20 minutes
    seconds_in_a_day = IP.day(IP.InsolationParameters(FT))
    start_t = seconds_after_midnight(date(t0))
    start_date = date(t0)
    local_noon = @. seconds_in_a_day * (FT(1 / 2) - longitude / 360) # allocates, but only on init

    affect! =
        (integrator) -> begin
            # Compute current date from t0 and elapsed time
            elapsed_days = floor(Int, float(integrator.t) / seconds_in_a_day)
            current_date = start_date + Dates.Day(elapsed_days)

            call_update_optimal_LAI(
                integrator.p,
                integrator.u,
                (float(integrator.t) + start_t) % (seconds_in_a_day), # current time in seconds UTC
                current_date;
                canopy = canopy,
                dt = dt,
                local_noon = local_noon,
                P_annual = P_annual,
                D_growing = D_growing,
                GSL = GSL,
            )
        end

    return IntervalBasedCallback(
        dt,         # period of this callback
        t0,         # simulation start
        dt,         # integration timestep
        affect!;
    )
end

"""
    get_model_callbacks(component::OptimalLAIModel, canopy; t0, Δt, P_annual, D_growing, GSL)

Creates the optimal LAI callback and returns it as a single element tuple of model callbacks.

The additional parameters (P_annual, D_growing, GSL) can be provided
as keyword arguments. If not provided, default values will be used.

Note: A₀ (potential GPP) is now computed dynamically from the P-model with fAPAR=1 and β=1,
so no Ao_annual parameter is needed.
"""
function get_model_callbacks(
    component::OptimalLAIModel{FT},
    canopy;
    t0,
    Δt,
    P_annual = FT(60000.0),
    D_growing = FT(1000.0),
    GSL = FT(240.0),
) where {FT}
    lai_cb = make_OptimalLAI_callback(
        FT,
        t0,
        Δt,
        canopy;
        P_annual = P_annual,
        D_growing = D_growing,
        GSL = GSL,
    )
    return (lai_cb,)
end
