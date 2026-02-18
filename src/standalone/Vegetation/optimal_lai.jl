export OptimalLAIParameters,
    compute_L_max,
    compute_m,
    lambertw0,
    compute_steady_state_LAI,
    compute_LAI,
    update_optimal_LAI,
    call_update_optimal_LAI,
    make_OptimalLAI_callback

"""
    OptimalLAIParameters{FT<:AbstractFloat}

The required parameters for the optimal LAI model based on Zhou et al. (2025).

Water limitation is handled through the f0*P/A0 term following Zhou et al. (2025) Equation 11,
where P is annual precipitation and A0 is annual potential GPP.

# References
Zhou et al. (2025) "A General Model for the Seasonal to Decadal Dynamics of Leaf Area"
Global Change Biology. https://onlinelibrary.wiley.com/doi/pdf/10.1111/gcb.70125

$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct OptimalLAIParameters{FT <: AbstractFloat}
    """Light extinction coefficient (dimensionless), typically 0.5"""
    k::FT
    """Unit cost of constructing and maintaining leaves (mol m^-2 yr^-1), globally fitted as 12.227 mol m^-2 yr^-1"""
    z::FT
    """Dimensionless parameter representing departure from square-wave LAI dynamics, globally fitted as 0.771"""
    sigma::FT
    """Smoothing factor for exponential moving average (dimensionless, 0-1). Set to 0.067 for ~15 days of memory"""
    alpha::FT
    """Fraction of annual precipitation available for transpiration (dimensionless, 0-1).
    Following Zhou et al. (2025), f0 = 0.65 at the energy-water limitation transition.
    In arid regions, f0 can be lower: f0 = 0.65 * exp(-0.604 * ln^2(AI/1.9)) where AI is aridity index.
    Default value 0.65 assumes optimal water use efficiency."""
    f0::FT
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
        sigma = FT(toml_dict["optimal_lai_sigma"]),
        alpha = FT(toml_dict["optimal_lai_alpha"]),
        f0 = FT(toml_dict["optimal_lai_f0"]),
    )
end

"""
    compute_L_max(Ao_annual, k, z, precip_annual, f0, ca_pa, chi, vpd_gs)

Compute seasonal maximum leaf area index (LAI_max) based on annual potential GPP
and water availability, following Zhou et al. (2025) Equation 11.

LAI_max is determined by the minimum of energy-limited and water-limited fAPAR:
- Energy-limited: fAPAR_energy = 1 - z/(k*A0)
- Water-limited: fAPAR_water = f0*P/A0 * (ca(1-chi))/(1.6*D)

# Arguments
- `Ao_annual::FT`: Annual total potential GPP (mol CO2 m^-2 yr^-1).
- `k::FT`: Light extinction coefficient (dimensionless), typically 0.5
- `z::FT`: Unit cost of constructing and maintaining leaves (mol m^-2 yr^-1), 12.227
- `precip_annual::FT`: Mean annual precipitation (mol H2O m^-2 yr^-1)
- `f0::FT`: Fraction of precipitation available for transpiration (dimensionless), 0.65
- `ca_pa::FT`: Ambient CO2 partial pressure (Pa), typically ~40 Pa at 400 ppm
- `chi::FT`: Optimal ratio of intercellular to ambient CO2 (dimensionless), typically 0.7-0.8
- `vpd_gs::FT`: Mean vapor pressure deficit during growing season (Pa)

# Returns
- `LAI_max::FT`: Seasonal maximum leaf area index (m^2 m^-2)

# Notes
Following Zhou et al. (2025) Equation 11:
```
fAPAR_max = min{1 - z/(k*A0), f0*P/A0 * (ca(1-chi))/(1.6*D)}
```
The first term is energy-limited (carbon gain vs leaf cost trade-off).
The second term is water-limited (precipitation constrains transpiration, scaled by
intrinsic water use efficiency iWUE = ca(1-chi)/(1.6*D)).

The iWUE factor converts water flux to carbon flux:
- ca(1-chi): CO2 drawdown from ambient to intercellular (Pa)
- 1.6*D: VPD adjusted for CO2/H2O diffusivity ratio (Pa)

# References
Zhou et al. (2025) Global Change Biology, Equation 11
"""
function compute_L_max(
    Ao_annual::FT,      # mol CO2 m^-2 yr^-1
    k::FT,              # dimensionless
    z::FT,              # mol m^-2 yr^-1
    precip_annual::FT,  # mol H2O m^-2 yr^-1
    f0::FT,             # dimensionless
    ca_pa::FT,          # Pa
    chi::FT,            # dimensionless
    vpd_gs::FT,         # Pa
) where {FT}
    # Handle edge case: very small or zero Ao_annual (e.g., polar regions)
    # When Ao_annual ~ 0, z / (k * Ao_annual) -> Inf, causing numerical issues.
    # Use ifelse for GPU compatibility.
    Ao_annual_safe = max(Ao_annual, eps(FT))

    # Energy-limited fAPAR (Equation 11, first term)
    # Plants optimize leaf area to maximize carbon gain minus construction cost
    fAPAR_energy = FT(1) - z / (k * Ao_annual_safe)

    # Water-limited fAPAR (Equation 11, second term)
    # fAPAR_water = f0 * P / A0 * (ca(1-chi)) / (1.6*D)
    # The iWUE factor (ca(1-chi))/(1.6*D) converts water flux to carbon flux
    # Guard against zero VPD
    vpd_safe = max(vpd_gs, eps(FT))
    iWUE_factor = (ca_pa * (FT(1) - chi)) / (FT(1.6) * vpd_safe)
    fAPAR_water = f0 * precip_annual / Ao_annual_safe * iWUE_factor

    # fAPAR_max is the minimum of energy and water constraints (Equation 11)
    fAPAR_max = min(fAPAR_energy, fAPAR_water)

    # Ensure fAPAR is in valid range [0, 1]
    fAPAR_max = max(FT(0), min(FT(1), fAPAR_max))

    # Convert fAPAR to LAI using Beer's law (Equation 12)
    # fAPAR = 1 - exp(-k * LAI)  ->  LAI = -(1/k) * ln(1 - fAPAR)
    # Guard against fAPAR_max = 1 which would give -log(0) = Inf
    fAPAR_max_safe = min(fAPAR_max, FT(1) - eps(FT))
    LAI_max = -(FT(1) / k) * log(FT(1) - fAPAR_max_safe)

    return LAI_max
end

fAPAR_max_fun(k::FT, LAI_max::FT) where {FT} = FT(1) - exp(-k * LAI_max)

"""
    compute_m(GSL, LAI_max, Ao_annual, sigma, k)

Compute the parameter m, which represents the ratio of steady-state LAI to steady-state GPP.

This implements Equation 20 from Zhou et al. (2025). The parameter m quantifies the
relationship between LAI and GPP dynamics, representing the extent to which seasonal LAI
dynamics depart from a "square wave" (where maximum LAI would be maintained throughout
the growing season).

# Arguments
- `GSL::FT`: Growing season length (days). Defined as the length of continuous period
  above 0C longer than 5 days.
- `LAI_max::FT`: Seasonal maximum leaf area index (m^2 m^-2, dimensionless)
- `Ao_annual::FT`: Annual total potential GPP (mol m^-2 yr^-1). This is the integral of daily
  A0 over the year.
- `sigma::FT`: Dimensionless parameter representing departure from square-wave LAI dynamics.
  Globally fitted as sigma = 0.771
- `k::FT`: Light extinction coefficient (dimensionless)

# Returns
- `m::FT`: Parameter relating steady-state LAI to steady-state GPP (dimensionless, units
  work out as: days * m^2 m^-2 / (mol m^-2 yr^-1 * dimensionless) with implicit conversion)

# References
Zhou et al. (2025) Global Change Biology, Equation 20
"""
function compute_m(
    GSL::FT,        # days
    LAI_max::FT,    # m^2 m^-2 (dimensionless)
    Ao_annual::FT,  # mol m^-2 yr^-1
    sigma::FT,      # dimensionless
    k::FT,
) where {FT}
    # Equation 20: m = (sigma * GSL * LAI_max) / (A0_sum * fAPAR_max)
    fAPAR_max = fAPAR_max_fun(k, LAI_max)

    # Guard against division by zero: when LAI_max ~ 0, fAPAR_max ~ 0,
    # but the numerator (sigma * GSL * LAI_max) is also ~ 0, so m ~ 0 naturally.
    fAPAR_max_safe = max(fAPAR_max, eps(FT))
    Ao_annual_safe = max(Ao_annual, eps(FT))
    m = (sigma * GSL * LAI_max) / (Ao_annual_safe * fAPAR_max_safe)
    return m
end

const MINARG = -inv(Base.MathConstants.e)

"""
    _lambertw0_initial_guess(x::T) where {T<:AbstractFloat}

Provide a robust initial guess for the Lambert W0 function for use in iterative solvers.

# Arguments
- `x::T`: Input value, should be >= -1/e

# Returns
- Initial guess for W0(x)

# Algorithm
- For x > 1: uses log(x) - log(log(x)) approximation
- For x < -0.32 (near -1/e): uses series expansion for accurate convergence near branch point
- For -0.32 <= x <= 1: uses max(x, -0.3) as a simple starting point
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
    lambertw0(x::T; maxiter::Int = 8) where {T<:AbstractFloat}

Compute the principal branch (W0) of the Lambert W function for x in [-1/e, Inf).

This is a GPU-device-friendly implementation using a fixed number of Halley iterations.
The Lambert W function satisfies W(x)*exp(W(x)) = x.

# Arguments
- `x::T`: Input value, must be >= -1/e ~ -0.36788
- `maxiter::Int`: Maximum number of Halley iterations (default: 8; Halley's method has cubic convergence, so 8 is generous)

# Returns
- `W::T`: Lambert W0(x), the principal branch value, or NaN for invalid inputs

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
@inline function lambertw0(x::T; maxiter::Int = 8) where {T <: AbstractFloat}
    # In our usage, arg = -k*mu*exp(-k*mu) with k > 0, mu >= 0,
    # so x is always in [-1/e, 0]. This check is a safety net.
    if !(isfinite(x)) || x < T(MINARG)
        return T(NaN)
    end
    w = _lambertw0_initial_guess(x)
    for i in 1:maxiter
        ew = exp(w)
        f = w * ew - x
        # Halley denominator
        # Special case: when w ~ -1, both numerator and denominator approach 0
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
- `Ao_daily::FT`: Daily potential GPP (mol m^-2 day^-1). This is the GPP that would be
  achieved if fAPAR = 1, calculated from LUE * PPFD.
- `m::FT`: Parameter relating steady-state LAI to steady-state GPP (dimensionless), from
  `compute_m()`
- `k::FT`: Light extinction coefficient (dimensionless), typically 0.5
- `LAI_max::FT`: Seasonal maximum LAI constraint (m^2 m^-2, dimensionless)

# Returns
- `L_steady::FT`: Steady-state leaf area index (m^2 m^-2, dimensionless). Always >= 0.

# Notes
The solution uses the Lambert W0 function: L_s = min{mu + (1/k)W0[-k*mu*exp(-k*mu)], LAI_max}
where mu = m * A0. The result is constrained to be non-negative and below LAI_max.

# References
Zhou et al. (2025) Global Change Biology, Equations 13-15
"""
function compute_steady_state_LAI(
    Ao_daily::FT,  # mol m^-2 day^-1
    m::FT,         # dimensionless
    k::FT,         # dimensionless
    LAI_max::FT,   # m^2 m^-2 (dimensionless)
) where {FT}
    # mu = m * A0 (Equation 15)
    mu = m * Ao_daily

    # Compute argument for Lambert W function
    arg = -k * mu * exp(-k * mu)

    # Check if argument is in valid range for W0 branch: [-1/e, 0]
    # If outside this range, use boundary solution
    if arg < -FT(1) / FT(exp(1))
        # Beyond valid range; use maximum possible LAI
        L_s = LAI_max
    else
        # Equation 15: L_s = mu + (1/k) * W0[-k mu exp(-k mu)]
        # Using our custom lambertw0 function (W0 is the principal branch)
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
    compute_LAI(LAI_prev, L_steady, alpha, local_noon_mask)

Compute updated LAI using exponential weighted moving average to represent lag between carbon
allocation and steady-state LAI.

This implements Equation 16 from Zhou et al. (2025). The exponential moving average
represents the time lag (days to months) for photosynthate allocation to leaves and
leaf development.

# Arguments
- `LAI_prev::FT`: LAI from previous time step (m^2 m^-2, dimensionless)
- `L_steady::FT`: Current steady-state LAI (m^2 m^-2, dimensionless), from
  `compute_steady_state_LAI()`
- `alpha::FT`: Smoothing factor (dimensionless, 0-1). Set to 0.067 for ~15 days of memory,
  meaning LAI[t] = 0.067 * L_steady[t] + 0.933 * LAI[t-1]
- `local_noon_mask::FT`: A mask (0 or 1) indicating whether the current time is within the local noon window.

# Returns
- `LAI_new::FT`: Updated leaf area index (m^2 m^-2, dimensionless). Always >= 0.

# Notes
The parameter alpha controls the response time:
- alpha = 0.067 ~ 15 days of memory (paper default)
- alpha = 0.1 ~ 10 days of memory (faster response)
- alpha = 0.033 ~ 30 days of memory (slower response)

The time scale tau ~ 1/alpha days.

# References
Zhou et al. (2025) Global Change Biology, Equation 16
"""
function compute_LAI(
    LAI_prev::FT,   # m^2 m^-2 (dimensionless)
    L_steady::FT,   # m^2 m^-2 (dimensionless)
    alpha::FT,      # dimensionless (0-1)
    local_noon_mask::FT,
) where {FT}
    if local_noon_mask == FT(1.0)
        # Equation 16: LAI_sim = alpha * L_s[t] + (1-alpha) * LAI_sim[t-1]
        LAI_new = alpha * L_steady + (1 - alpha) * LAI_prev

        # L_steady is already clipped to >= 0 in compute_steady_state_LAI,
        # so this is a safety net for initialization edge cases only.
        LAI_new = max(zero(FT), LAI_new)
        return LAI_new
    else
        return LAI_prev
    end
end

"""
    update_optimal_LAI(local_noon_mask, A0_daily, L, k, A0_annual, z, GSL, sigma, alpha, precip_annual, f0, ca_pa, chi, vpd_gs)

Update LAI using the optimal LAI model with precomputed daily and annual potential GPP.

# Arguments
- `local_noon_mask::FT`: Mask (0 or 1) indicating if it's local noon
- `A0_daily::FT`: Daily potential GPP (mol CO2 m^-2 day^-1), with moisture stress factor beta
- `L::FT`: Current LAI (m^2 m^-2)
- `k::FT`: Light extinction coefficient
- `A0_annual::FT`: Annual potential GPP (mol CO2 m^-2 yr^-1)
- `z::FT`: Unit cost of constructing and maintaining leaves (mol m^-2 yr^-1)
- `GSL::FT`: Growing season length (days)
- `sigma::FT`: Dimensionless parameter for LAI dynamics
- `alpha::FT`: Smoothing factor for exponential moving average (~15-day memory)
- `precip_annual::FT`: Mean annual precipitation (mol H2O m^-2 yr^-1)
- `f0::FT`: Fraction of precipitation available for transpiration (dimensionless)
- `ca_pa::FT`: Ambient CO2 partial pressure (Pa)
- `chi::FT`: Optimal ratio of intercellular to ambient CO2 (dimensionless)
- `vpd_gs::FT`: Mean vapor pressure deficit during growing season (Pa)

# Returns
Updated LAI value.

# Notes
Following Zhou et al. (2025):
- A0_daily uses moisture stress factor beta to drive daily LAI dynamics
- A0_annual uses actual moisture stress factor beta for LAI_max computation
- Water limitation enters LAI_max through the f0*P/A0 * (ca(1-chi))/(1.6*D) term (Equation 11)
"""
function update_optimal_LAI(
    local_noon_mask::FT,
    A0_daily::FT,
    L::FT, # m2 m-2
    k::FT,
    A0_annual::FT, # mol CO2 m-2 y-1
    z::FT, # mol m-2 yr-1, leaf construction cost
    GSL::FT, # days, growing season length
    sigma::FT, # dimensionless
    alpha::FT, # dimensionless (~15-day memory)
    precip_annual::FT, # mol H2O m-2 yr-1, mean annual precipitation
    f0::FT, # dimensionless, fraction of precip for transpiration
    ca_pa::FT, # Pa, ambient CO2 partial pressure
    chi::FT, # dimensionless, optimal ci/ca ratio
    vpd_gs::FT, # Pa, mean VPD during growing season
) where {FT}
    LAI_max =
        compute_L_max(A0_annual, k, z, precip_annual, f0, ca_pa, chi, vpd_gs)
    m = compute_m(GSL, LAI_max, A0_annual, sigma, k)
    L_steady = compute_steady_state_LAI(A0_daily, m, k, LAI_max)
    L = compute_LAI(L, L_steady, alpha, local_noon_mask)
    return L
end


"""
    call_update_optimal_LAI(p, Y, t, current_date; canopy, dt, local_noon)

Updates LAI and accumulates potential GPP (A0) at each timestep.

Every timestep: accumulates instantaneous potential GPP into the daily accumulator.
At local noon: finalizes daily A0, adds it to annual accumulator, updates LAI.
Every 365 days: finalizes annual A0 from the accumulator.

Uses air temperature (not canopy temperature) for A0 computation, since canopy
temperature includes energy balance feedbacks that should not affect potential GPP.

GSL (Growing Season Length) is read from p.canopy.biomass.GSL, which supports spatially
varying values initialized via set_historical_cache!.
"""
function call_update_optimal_LAI(p, Y, t, current_date; canopy, dt, local_noon)
    FT = eltype(canopy.biomass.parameters)

    # Compute local noon mask
    local_noon_mask = @. lazy(get_local_noon_mask(t, dt, local_noon))

    # Get P-model parameters and constants for computing A0
    pmodel_parameters = canopy.photosynthesis.parameters
    pmodel_constants = canopy.photosynthesis.constants
    is_c3 = canopy.photosynthesis.is_c3

    # Get drivers for A0 computation
    # Use air temperature (not canopy temperature) for potential GPP computation.
    # Canopy temperature includes energy balance feedbacks that should not affect A0.
    T_air = p.drivers.T
    P_air = p.drivers.P
    ca = p.drivers.c_co2  # mol/mol
    earth_param_set = canopy.earth_param_set

    # Get soil moisture stress factor (beta)
    βm = p.canopy.soil_moisture_stress.βm

    # Compute VPD (clipped to avoid numerical issues)
    VPD = @. lazy(
        max(
            Thermodynamics.vapor_pressure_deficit(
                LP.thermodynamic_parameters(earth_param_set),
                p.drivers.T,
                p.drivers.P,
                p.drivers.q,
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

    # Compute instantaneous potential GPP (kg C m^-2 s^-1)
    # A single computation serves both daily and annual accumulation
    dt_seconds = FT(float(dt))
    Mc = pmodel_constants.Mc
    @. p.canopy.biomass.A0_daily_acc +=
        compute_A0_daily(
            is_c3,
            pmodel_parameters,
            pmodel_constants,
            T_air,
            P_air,
            VPD,
            ca,
            PPFD,
            βm,
        ) * dt_seconds / Mc

    # Get parameters from the biomass model
    parameters = canopy.biomass.parameters

    # At local noon: finalize daily A0, update annual accumulator, update LAI
    # Snapshot the daily accumulator before any resets
    A0_daily_final = @. p.canopy.biomass.A0_daily_acc

    # Use days_since_reset to track accumulation period (365-day rolling window).
    # This ensures the first annual update happens exactly 365 days after
    # simulation start, and subsequent updates every 365 days thereafter.
    days_since_reset = p.canopy.biomass.days_since_reset

    # At noon: finalize annual A0 when we have accumulated 365 days
    A0_annual_final = @. ifelse(
        local_noon_mask == FT(1) && days_since_reset >= FT(365),
        p.canopy.biomass.A0_annual_acc,
        p.canopy.biomass.A0_annual,
    )

    # Compute chi for water limitation term (using growing season mean VPD)
    chi = @. lazy(
        compute_chi(
            pmodel_parameters,
            pmodel_constants,
            T_air,
            P_air,
            p.canopy.biomass.vpd_gs,
            ca,
        ),
    )

    # Update LAI at noon using the finalized values
    @. p.canopy.biomass.area_index.leaf = ifelse(
        local_noon_mask == FT(1),
        update_optimal_LAI(
            FT(1),  # At noon, we always update
            A0_daily_final,
            p.canopy.biomass.area_index.leaf,
            parameters.k,
            A0_annual_final,
            parameters.z,
            p.canopy.biomass.GSL,
            parameters.sigma,
            parameters.alpha,
            p.canopy.biomass.precip_annual,
            p.canopy.biomass.f0,
            ca * P_air,  # ca_pa: CO2 partial pressure (Pa)
            chi,
            p.canopy.biomass.vpd_gs,
        ),
        p.canopy.biomass.area_index.leaf,
    )

    # Update A0_annual when 365 days have accumulated
    @. p.canopy.biomass.A0_annual = ifelse(
        local_noon_mask == FT(1) && days_since_reset >= FT(365),
        A0_annual_final,
        p.canopy.biomass.A0_annual,
    )

    # Update A0_daily at noon (stores the daily sum for diagnostics)
    @. p.canopy.biomass.A0_daily = ifelse(
        local_noon_mask == FT(1),
        A0_daily_final,
        p.canopy.biomass.A0_daily,
    )

    # At noon: add daily sum to annual accumulator, then reset daily accumulator
    # When 365 days reached: start fresh with today's value
    @. p.canopy.biomass.A0_annual_acc = ifelse(
        local_noon_mask == FT(1) && days_since_reset >= FT(365),
        A0_daily_final,  # Reset and start with today's value
        ifelse(
            local_noon_mask == FT(1),
            p.canopy.biomass.A0_annual_acc + A0_daily_final,
            p.canopy.biomass.A0_annual_acc,
        ),
    )

    # Update days_since_reset counter at noon
    @. p.canopy.biomass.days_since_reset = ifelse(
        local_noon_mask == FT(1) && days_since_reset >= FT(365),
        FT(1),  # Reset to 1 (today's value already added)
        ifelse(
            local_noon_mask == FT(1),
            days_since_reset + FT(1),
            days_since_reset,
        ),
    )

    # Reset daily accumulator at noon
    @. p.canopy.biomass.A0_daily_acc =
        ifelse(local_noon_mask == FT(1), FT(0), p.canopy.biomass.A0_daily_acc)
end

"""
    make_OptimalLAI_callback(::Type{FT}, t0::ITime, dt, canopy; longitude) where {FT <: AbstractFloat}

This constructs an IntervalBasedCallback for the optimal LAI model that:
1. Computes and accumulates potential GPP (A0) at each timestep
2. Updates LAI using an exponential moving average at local noon
3. Tracks daily and annual A0 sums

We check for local noon using the provided `longitude` every dt.
The time of local noon is expressed in seconds UTC and neglects the effects of obliquity and eccentricity, so
it is constant throughout the year.

# Arguments
- `FT`: The floating-point type used in the model (e.g., `Float32`, `Float64`).
- `t0`: ITime, with epoch in UTC.
- `dt`: timestep
- `canopy`: the canopy object containing the optimal LAI model parameters.
- `longitude`: optional longitude in degrees for local noon calculation (default is `nothing`, which means
    that it will be inferred from the canopy domain).

# Notes
- Daily A0 is computed with fAPAR=1 and moisture stress factor beta - drives L_steady
- Annual A0 is computed with fAPAR=1 and actual moisture stress factor beta - used for LAI_max
- Water limitation enters LAI_max through the f0*P/A0 * (ca(1-chi))/(1.6*D) term (Equation 11)
- Daily A0 is accumulated over each day and finalized at local noon
- Annual A0 is accumulated and reset on January 1
- GSL (Growing Season Length) is read from p.canopy.biomass.GSL, which supports spatially
  varying values initialized via set_historical_cache!.
"""
function make_OptimalLAI_callback(
    ::Type{FT},
    t0,
    dt,
    canopy;
    longitude = nothing,
) where {FT <: AbstractFloat}
    function seconds_after_midnight(d)
        return FT(Hour(d).value * 3600 + Minute(d).value * 60 + Second(d).value)
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
            )
        end

    return IntervalBasedCallback(
        dt,         # period of this callback
        t0,         # simulation start
        dt,         # integration timestep
        affect!;
    )
end
